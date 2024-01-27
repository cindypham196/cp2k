!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Environment independent subroutines for TIDFT
!> \par History
!>       2021.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_methods
   USE almo_scf_methods,                ONLY: almo_scf_t_to_proj,&
                                              orthogonalize_mos
   USE almo_scf_types,                  ONLY: almo_scf_env_type,&
                                              almo_scf_history_type,&
                                              almo_mat_dim_aobasis,&
                                              almo_mat_dim_occ
   USE bibliography,                    ONLY: Kolafa2004,&
                                              Kuhne2007,&
                                              cite_reference
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_lu_decompose,&
                                              cp_fm_adjugate
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, dbcsr_frobenius_norm, &
        dbcsr_get_block_p, dbcsr_get_diag, dbcsr_get_info, dbcsr_get_stored_coordinates, &
        dbcsr_init_random, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_print, dbcsr_release, &
        dbcsr_reserve_block2d, dbcsr_scale_by_vector, dbcsr_set, dbcsr_set_diag, dbcsr_transposed, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_work_create, &
        dbcsr_complete_redistribute, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm, dbcsr_norm, dbcsr_dot
   USE domain_submatrix_methods,        ONLY: &
        add_submatrices, construct_dbcsr_from_submatrices, construct_submatrices, &
        copy_submatrices, copy_submatrix_data, init_submatrices, multiply_submatrices, &
        print_submatrices, release_submatrices
   USE domain_submatrix_types,          ONLY: domain_map_type,&
                                              domain_submatrix_type,&
                                              select_row,&
                                              select_row_col
   USE fermi_utils,                     ONLY: FermiFixed
   USE input_constants,                 ONLY: almo_domain_layout_molecular,&
                                              almo_mat_distr_atomic,&
                                              almo_scf_diag,&
                                              spd_inversion_dense_cholesky,&
                                              spd_inversion_ls_hotelling,&
                                              spd_inversion_ls_taylor
   USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                              invert_Taylor,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: binomial
   USE tidft_qs,                        ONLY: matrix_layout
   USE tidft_types,                     ONLY: tidft_env_type
   USE util,                            ONLY: sort
   USE machine,                         ONLY: m_flush
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_methods'

   PUBLIC stitch_columns, determinant_local, &
          split_columns, matrix_invert_local, &
          gradVAR_scale_split_increment, &
          precMO_to_precVAR_increment, &
          tidft_optimizer_mainvar_multiply, &
          tidft_optimizer_mainvar_multiply_diff, &
          tidft_optimizer_mainvar_dot, &
          tidft_optimizer_mainvar_add, &
          tidft_optimizer_mainvar_copy, &
          tidft_optimizer_mainvar_set, &
          tidft_optimizer_mainvar_print, &
          dbcsr_set_element, &
          project_out, &
          tidft_optimizer_gradient_norm, &
          tidft_mainvar_orthonormalize, &
          swap_columns

CONTAINS

! **************************************************************************************************
!> \brief Mutliply all mainvar matrices (index i) by different NxN matrices
!>        C_i = beta C_i + alpha A_i**tranA . B_i**tranB
!> \param ncols_var ...
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_multiply_diff(transA, transB, alpha, matA, matB, &
      beta, matC, filter_eps, ncols_var)

      CHARACTER, INTENT(IN)                              :: transA, transB
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)        :: matA
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: matB
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: matC
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta, &
                                                            filter_eps
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var

      INTEGER                                            :: istate, ivar, &
                                                            nstates, nvars, &
                                                            nstates2, nvars2

      nstates = SIZE(matB, 1)
      nvars = SIZE(matB, 2)
      CPASSERT(nvars.EQ.3)

      nstates2 = SIZE(matC, 1)
      CPASSERT(nstates.EQ.nstates2)
      nvars2 = SIZE(matC, 2)
      CPASSERT(nvars.EQ.nvars2)

      nstates2 = SIZE(matB, 1)
      CPASSERT(nstates.EQ.nstates2)
      nvars2 = SIZE(matB, 2)
      CPASSERT(nvars.EQ.nvars2)

      DO ivar = 1, nvars

         IF (ncols_var(ivar).GT.0) THEN

            DO istate = 1, nstates

                  CALL dbcsr_multiply(transA, transB, &
                          alpha, &
                          matA(istate,ivar), &
                          matB(istate,ivar), &
                          beta, &
                          matC(istate,ivar), &
                          filter_eps=filter_eps)

            ENDDO ! stateI

         ENDIF

      ENDDO ! ivar

   END SUBROUTINE tidft_optimizer_mainvar_multiply_diff

! **************************************************************************************************
!> \brief Mutliply all mainvar matrices (index i) by the same NxN matrix
!>        C_i = beta C_i + alpha A**tranA . B_i**tranB
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_multiply(transA, transB, alpha, matA, matB, &
      beta, matC, filter_eps, ncols_var)

      CHARACTER, INTENT(IN)                              :: transA, transB
      TYPE(dbcsr_type), INTENT(IN)                       :: matA
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: matB
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: matC
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta, &
                                                            filter_eps
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var

      INTEGER                                            :: istate, ivar, &
                                                            nstates, nvars, &
                                                            nstates2, nvars2

      nstates = SIZE(matB, 1)
      nvars = SIZE(matB, 2)
      CPASSERT(nvars.EQ.3)

      nstates2 = SIZE(matC, 1)
      CPASSERT(nstates.EQ.nstates2)
      nvars2 = SIZE(matC, 2)
      CPASSERT(nvars.EQ.nvars2)

      DO ivar = 1, nvars

         IF (ncols_var(ivar).GT.0) THEN

            DO istate = 1, nstates

                  CALL dbcsr_multiply(transA, transB, &
                          alpha, &
                          matA, &
                          matB(istate,ivar), &
                          beta, &
                          matC(istate,ivar), &
                          filter_eps=filter_eps)

            ENDDO ! stateI

         ENDIF

      ENDDO ! ivar

   END SUBROUTINE tidft_optimizer_mainvar_multiply

! **************************************************************************************************
!> \brief Dot product of matrices describing main variables: 
!>        dot_product = sum_states ( 2 Tr(L1**T . R1) + Tr(L2**T.R2) + Tr (L3**T.R3) )
!>        The subroutine is made tensorially consistent
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_dot(mat_left, mat_right, &
      mat_left_covariant, mat_right_covariant, metric, metric_inv, &
      dot_product, ncols_var)

      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: mat_left, mat_right
      LOGICAL, INTENT(IN)                                :: mat_left_covariant, &
                                                            mat_right_covariant
      TYPE(dbcsr_type), OPTIONAL, INTENT(IN)             :: metric
      TYPE(dbcsr_type), OPTIONAL, INTENT(IN)             :: metric_inv
      REAL(KIND=dp), INTENT(OUT)                         :: dot_product
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var
      !REAL(KIND=dp), OPTIONAL, INTENT(IN)                :: eps_filter

      INTEGER                                            :: istate, ivar, &
                                                            nstates, nvars, &
                                                            nstates2, nvars2
      REAL(KIND=dp)                                      :: prod, factor
      LOGICAL                                            :: mat_left_contravariant, &
                                                            mat_right_contravariant
      TYPE(dbcsr_type)                                   :: mat_temp

      mat_left_contravariant = .NOT.mat_left_covariant
      mat_right_contravariant = .NOT.mat_right_covariant

      IF (mat_left_contravariant .AND. mat_right_contravariant) THEN
         CPASSERT(PRESENT(metric))
         !CPASSERT(PRESENT(eps_filter))
      ENDIF

      IF (mat_left_covariant .AND. mat_right_covariant) THEN
         CPASSERT(PRESENT(metric_inv))
         !CPASSERT(PRESENT(eps_filter))
      ENDIF

      nstates = SIZE(mat_left, 1)
      nvars = SIZE(mat_left, 2)
      CPASSERT(nvars.EQ.3)

      nstates2 = SIZE(mat_right, 1)
      CPASSERT(nstates.EQ.nstates2)
      nvars2 = SIZE(mat_right, 2)
      CPASSERT(nvars.EQ.nvars2)

      dot_product = 0.0_dp
      DO ivar = 1, nvars

         IF (ncols_var(ivar).GT.0) THEN

            !RZK-WARNING: this factor of 2 worries me. Should it be 1?
            IF (ivar.EQ.1) THEN
               factor = 1.0_dp
            ELSE
               factor = 1.0_dp
            ENDIF

            DO istate = 1, nstates

               IF (mat_left_contravariant .AND. mat_right_contravariant) THEN

                  CALL dbcsr_create(mat_temp, template=mat_right(istate,ivar))
                  CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          metric, &
                          mat_right(istate,ivar), &
                          0.0_dp, mat_temp) !, &
                          !filter_eps=eps_filter)
                  CALL dbcsr_dot(mat_left(istate,ivar), mat_temp, prod)
                  CALL dbcsr_release(mat_temp)

               ELSE IF (mat_left_covariant .AND. mat_right_covariant) THEN

                  CALL dbcsr_create(mat_temp, template=mat_right(istate,ivar))
                  CALL dbcsr_multiply("N", "N", 1.0_dp, & 
                          metric_inv, &
                          mat_right(istate,ivar), &
                          0.0_dp, mat_temp) !, &
                          !filter_eps=eps_filter)
                  CALL dbcsr_dot(mat_left(istate,ivar), mat_temp, prod)
                  CALL dbcsr_release(mat_temp)

               ELSE

                  CALL dbcsr_dot(mat_left(istate,ivar), mat_right(istate,ivar), prod)

               ENDIF

               dot_product = dot_product + factor * prod

            ENDDO ! stateI

         ENDIF

      ENDDO ! ivar

   END SUBROUTINE tidft_optimizer_mainvar_dot

! **************************************************************************************************
!> \brief Add matrices describing main variables: 
!>        mat_dest = alpha * mat_dest + beta * mat_source
!> \param mat_dest ...
!> \param mat_source ...
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_add(mat_dest, mat_source, alpha, beta, ncols_var)

      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: mat_dest
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: mat_source
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var

      INTEGER                                            :: istate, ivar, &
                                                            nstates, nvars, &
                                                            nstates2, nvars2

      nstates = SIZE(mat_source, 1)
      nvars = SIZE(mat_source, 2)
      CPASSERT(nvars.EQ.3)

      nstates2 = SIZE(mat_dest, 1)
      CPASSERT(nstates.EQ.nstates2)
      nvars2 = SIZE(mat_dest, 2)
      CPASSERT(nvars.EQ.nvars2)

      DO istate = 1, nstates
         DO ivar = 1, nvars
            IF (ncols_var(ivar).GT.0) THEN
               CALL dbcsr_add(mat_dest(istate,ivar), mat_source(istate,ivar), &
                       alpha, beta)
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

   END SUBROUTINE tidft_optimizer_mainvar_add

! **************************************************************************************************
!> \brief 
!> \param mat_dest ...
!> \param mat_source ...
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_print(mat, ncols_var)

      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: mat
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var

      INTEGER                                            :: istate, ivar, &
                                                            nstates, nvars

      nstates = SIZE(mat, 1)
      nvars = SIZE(mat, 2)
      CPASSERT(nvars.EQ.3)

      DO istate = 1, nstates
         DO ivar = 1, nvars
            IF (ncols_var(ivar).GT.0) THEN
               CALL dbcsr_print(mat(istate,ivar))
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

   END SUBROUTINE tidft_optimizer_mainvar_print

! **************************************************************************************************
!> \brief Copy matrices describing main variables
!> \param mat_dest ...
!> \param mat_source ...
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_set(mat, val, ncols_var)

      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: mat
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var
      REAL(KIND=dp), INTENT(IN)                          :: val

      INTEGER                                            :: istate, ivar, &
                                                            nstates, nvars

      nstates = SIZE(mat, 1)
      nvars = SIZE(mat, 2)
      CPASSERT(nvars.EQ.3)

      DO istate = 1, nstates
         DO ivar = 1, nvars
            IF (ncols_var(ivar).GT.0) THEN
               CALL dbcsr_set(mat(istate,ivar), val)
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

   END SUBROUTINE tidft_optimizer_mainvar_set

! **************************************************************************************************
!> \brief Copy matrices describing main variables
!> \param mat_dest ...
!> \param mat_source ...
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_copy(mat_dest, mat_source, ncols_var)

      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: mat_dest
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: mat_source
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var

      INTEGER                                            :: istate, ivar, &
                                                            nstates, nvars, &
                                                            nstates2, nvars2

      nstates = SIZE(mat_source, 1)
      nvars = SIZE(mat_source, 2)
      CPASSERT(nvars.EQ.3)

      nstates2 = SIZE(mat_dest, 1)
      CPASSERT(nstates.EQ.nstates2)
      nvars2 = SIZE(mat_dest, 2)
      CPASSERT(nvars.EQ.nvars2)

      DO istate = 1, nstates
         DO ivar = 1, nvars
            IF (ncols_var(ivar).GT.0) THEN
               CALL dbcsr_copy(mat_dest(istate,ivar), mat_source(istate,ivar))
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

   END SUBROUTINE tidft_optimizer_mainvar_copy

! **************************************************************************************************
!> \brief Compute the gradient norm
!> \param mat_grad ...
!> \param ncols_var ...
!> \param grad_norm ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_gradient_norm(mat_grad, which_norm, grad_norm, &
      ncols_var, metric_inv)

      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: mat_grad
      INTEGER, INTENT(IN)                                :: which_norm
      REAL(KIND=dp), INTENT(OUT)                         :: grad_norm
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var
      TYPE(dbcsr_type), OPTIONAL, INTENT(IN)             :: metric_inv

      INTEGER                                            :: istate, ivar, &
                                                            nstates, nvars, &
                                                            ncols, nbasis, &
                                                            itmp
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: grad_norm_array
      REAL(KIND=dp)                                      :: denom

      nstates = SIZE(mat_grad, 1)
      nvars = SIZE(mat_grad, 2)
      CPASSERT(nvars.EQ.3)

      SELECT CASE (which_norm)
      CASE (dbcsr_norm_maxabsnorm)

         ALLOCATE(grad_norm_array(nstates,nvars))
         grad_norm_array=0.0_dp

         DO istate = 1, nstates
            DO ivar = 1, nvars
               IF (ncols_var(ivar).GT.0) THEN
                  CALL dbcsr_norm(matrix=mat_grad(istate,ivar), &
                          which_norm=dbcsr_norm_maxabsnorm, &
                          norm_scalar=grad_norm_array(istate,ivar))
               ENDIF
            ENDDO ! ivar
         ENDDO ! stateI
         grad_norm = MAXVAL(grad_norm_array)

         DEALLOCATE(grad_norm_array)

      CASE (dbcsr_norm_frobenius)

         CPASSERT(PRESENT(metric_inv))
         CALL tidft_optimizer_mainvar_dot( &
                 mat_left = mat_grad, &
                 mat_right = mat_grad, &
                 mat_left_covariant = .TRUE., &
                 mat_right_covariant = .TRUE., &
                 metric_inv = metric_inv, &
                 dot_product = grad_norm, &
                 ncols_var = ncols_var) !, &

         ncols = 0
         DO istate = 1, nstates
            DO ivar = 1, nvars
               IF (ncols_var(ivar).GT.0) THEN
                  CALL dbcsr_get_info(mat_grad(istate,ivar), nfullrows_total=nbasis)
                  CALL dbcsr_get_info(mat_grad(istate,ivar), nfullcols_total=itmp)
                  ncols = ncols + itmp
               ENDIF
            ENDDO ! ivar
         ENDDO ! stateI

         denom = SQRT( 1.0_dp*INT(nbasis, KIND=int_8) ) * SQRT( 1.0_dp*INT(ncols, KIND=int_8) )
         IF ( denom .NE. 0.0_dp) THEN
            grad_norm = SQRT(ABS(grad_norm)) / denom
         ELSE
            CPABORT("Zero electronic variables")
         ENDIF

      END SELECT

   END SUBROUTINE tidft_optimizer_gradient_norm

! **************************************************************************************************
!> \brief Converts MO-matrices to the corresponding matrices of independent variables
!> \param mat_g_out ...
!> \param ncols_var ...
!> \param mat_g_in ...
!> \param idet ...
!> \param ispin ...
!> \param nspins ...
!> \param triplet_correction ...
!> \param energy_term ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE gradVAR_scale_split_increment(mat_g_out, ncols_var, mat_g_in, mat_s, &
      idet, ispin, nspins, triplet_correction, energy_term)

      TYPE(dbcsr_type), DIMENSION(3), INTENT(INOUT)      :: mat_g_out
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_g_in
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_s
      INTEGER, INTENT(IN)                                :: idet, ispin, nspins
      LOGICAL                                            :: triplet_correction, &
                                                            energy_term

      REAL(KIND=dp)                                      :: triplet_correction_energy_det1_multiplier 
      INTEGER                                            :: varA, varB, &
                                                            ncolsA, ncolsB
      TYPE(dbcsr_type)                                   :: tempNVarA, tempNVarB, &
                                                            mat_var23


      triplet_correction_energy_det1_multiplier = 1.0_dp
      IF (triplet_correction .AND. energy_term) THEN
         IF (idet.EQ.1) triplet_correction_energy_det1_multiplier = 2.0_dp
         IF (idet.EQ.2) triplet_correction_energy_det1_multiplier = -1.0_dp
      ENDIF

      IF (nspins .EQ. 1) THEN

         ! idet is 1 since closed-shell calculations cannot have more determinants
         CPASSERT(idet.EQ.1)

         ! scale the term using spin factor 2.0
         CALL dbcsr_add(mat_g_out(1), &
                 mat_g_in, &
                 1.0_dp, 2.0_dp*triplet_correction_energy_det1_multiplier)

      ELSE

         IF (idet .EQ. 2) THEN ! this is the "spin-flipped" determinant

            varA = 1
            ncolsA = ncols_var(varA)

            IF (ispin .EQ. 1) THEN
            ! Alpha T-matrix: contains Var1, Var2, Var3. Only ncols(Var1) can be zero

               IF (ncolsA .NE. 0) THEN
               ! separate into v1, v2 and v3 using two splits
               ! the convention: the rightmost colum is the "spin-flipped" electron

                  ! first split
                  CALL dbcsr_create(tempNVarA, template=mat_g_out(varA))
                  CALL matrix_layout( &
                          matrix_ao_template = mat_s, &
                          matrix_new = mat_var23 , &
                          name_new = 'm_var23', &
                          orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                          suggested_mo_blk_size = 2, &
                          symmetry_new = dbcsr_type_no_symmetry, &
                          nmos = (/-1, 2/) ) ! 2 singly-occupied orbitals - alpha and beta
                  CALL split_columns( &
                          mat_g_in, &
                          tempNVarA, &
                          mat_var23 )
                  CALL dbcsr_add(mat_g_out(varA), &
                          tempNVarA, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)
                  CALL dbcsr_release(tempNVarA)

                  ! second split
                  CALL dbcsr_create(tempNVarA, template=mat_g_out(2))
                  CALL dbcsr_create(tempNVarB, template=mat_g_out(3))
                  CALL split_columns( &
                          mat_var23, &
                          tempNVarA, &
                          tempNVarB )
                  CALL dbcsr_add(mat_g_out(2), &
                          tempNVarA, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)
                  CALL dbcsr_add(mat_g_out(3), &
                          tempNVarB, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)
                  CALL dbcsr_release(tempNVarA)
                  CALL dbcsr_release(tempNVarB)
                  CALL dbcsr_release(mat_var23)

               ELSE ! ncolsA==0
               ! no doubly-occupied MOs: separate into v2 and v3
               ! the convention: the rightmost colum is the "spin-flipped" electron

                  CALL dbcsr_create(tempNVarA, template=mat_g_out(2))
                  CALL dbcsr_create(tempNVarB, template=mat_g_out(3))
                  CALL split_columns( &
                          mat_g_in, &
                          tempNVarA, &
                          tempNVarB )
                  CALL dbcsr_add(mat_g_out(2), &
                          tempNVarA, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)
                  CALL dbcsr_add(mat_g_out(3), &
                          tempNVarB, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)
                  CALL dbcsr_release(tempNVarA)
                  CALL dbcsr_release(tempNVarB)

               ENDIF ! ncolsA decision

            ELSE ! beta electrons
            ! Beta T-matrix: contains only Var1. ncols(Var1) can be zero

               IF (ncolsA .NE. 0) THEN
               ! no singly-occupied MOs

                  CALL dbcsr_add(mat_g_out(varA), &
                          mat_g_in, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)

               ELSE ! ncolsA==0
               ! zero beta electrons
               ! this is allowed but such situations should be caught on the upper levels

                  CPABORT("Cannot end up here: zero T matrices must be filtered out")

               ENDIF

            ENDIF

         ELSE ! idet check
         ! idet.NE.2 - regular determinants

            varA = 1
            IF (ispin .EQ. 1) THEN
               varB = 2
            ELSE ! beta electrons
               varB = 3
            ENDIF
            ncolsA = ncols_var(varA)
            ncolsB = ncols_var(varB)

            IF (ncolsA .NE. 0) THEN ! There are doubly-occupied states

               IF (ncolsB .NE. 0) THEN
               ! Split the matrix since both parts are non-zero

                  CALL dbcsr_create(tempNVarA, template=mat_g_out(varA))
                  CALL dbcsr_create(tempNVarB, template=mat_g_out(varB))
                  CALL split_columns( &
                          mat_g_in, &
                          tempNVarA, &
                          tempNVarB )
                  CALL dbcsr_add(mat_g_out(varA), &
                          tempNVarA, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)
                  CALL dbcsr_add(mat_g_out(varB), &
                          tempNVarB, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)
                  CALL dbcsr_release(tempNVarA)
                  CALL dbcsr_release(tempNVarB)

               ELSE ! ncolsB==0
               ! restricted closed-shell cases with no unpairings can be treated here
               ! but with the 2.0 multiplier

                  ! all MOs are doubly-occupied
                  CALL dbcsr_add(mat_g_out(varA), &
                          mat_g_in, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)

               ENDIF

            ELSE ! ncolsA==0
            ! no paired electrons

               IF (ncolsB .NE. 0) THEN
               ! there are only unpaired electrons

                  CALL dbcsr_add(mat_g_out(varB), &
                          mat_g_in, &
                          1.0_dp, 1.0_dp*triplet_correction_energy_det1_multiplier)

               ELSE ! ncolsB==0
               ! no electrons
               ! this is possible but such situations should be caught on the upper levels

                  CPABORT("Warning: empty matrices should be eliminated above")

               ENDIF

            ENDIF

         ENDIF ! idet check

      ENDIF ! nspins check

   END SUBROUTINE gradVAR_scale_split_increment

! **************************************************************************************************
!> \brief Converts MO-matrices to the corresponding matrices of independent variables
!> \param mat_g_out ...
!> \param ncols_var ...
!> \param mat_g_in ...
!> \param idet ...
!> \param ispin ...
!> \param nspins ...
!> \param triplet_correction ...
!> \param energy_term ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE precMO_to_precVAR_increment(mat_prec_out, ncols_var, mat_prec_in, &
      idet, ispin)

      TYPE(dbcsr_type), DIMENSION(3), INTENT(INOUT)      :: mat_prec_out
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_prec_in
      INTEGER, INTENT(IN)                                :: idet, ispin

      INTEGER                                            :: varA, varB, &
                                                            ncolsA, ncolsB

      IF (idet .EQ. 2) THEN ! this is the "spin-flipped" determinant

         varA = 1
         ncolsA = ncols_var(varA)

         IF (ncolsA .NE. 0) THEN
            CALL dbcsr_add(mat_prec_out(varA), mat_prec_in, 1.0_dp, 1.0_dp)
         ENDIF ! ncolsA != 0

         IF (ispin .EQ. 1) THEN
         ! Alpha matrix: contains Var1, Var2, Var3. 

            CALL dbcsr_add(mat_prec_out(2), mat_prec_in, 1.0_dp, 1.0_dp)
            CALL dbcsr_add(mat_prec_out(3), mat_prec_in, 1.0_dp, 1.0_dp)

         ENDIF
         ! Beta T-matrix: contains only Var1

      ELSE ! idet check
      ! idet.NE.2 - regular determinants

         varA = 1
         IF (ispin .EQ. 1) THEN
            varB = 2
         ELSE ! beta electrons
            varB = 3
         ENDIF
         ncolsA = ncols_var(varA)
         ncolsB = ncols_var(varB)

         IF (ncolsA .NE. 0) THEN
            CALL dbcsr_add(mat_prec_out(varA), mat_prec_in, 1.0_dp, 1.0_dp)
         ENDIF

         IF (ncolsB .NE. 0) THEN
            CALL dbcsr_add(mat_prec_out(varB), mat_prec_in, 1.0_dp, 1.0_dp)
         ENDIF

      ENDIF ! idet check

   END SUBROUTINE precMO_to_precVAR_increment

! **************************************************************************************************
!> \brief Computes the determinant of a square matrix via LU-decomposition
!> \param matrix_a ...
!> \param determinant ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE matrix_invert_local(matrix_a, n0)

      REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT)       :: matrix_a
      INTEGER, OPTIONAL, INTENT(IN)                      :: n0

      CHARACTER(len=*), PARAMETER :: routineN = 'matrix_invert_local'

      INTEGER                                            :: handle, i, info, n, &
                                                            m, lda, j
      !REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: a

      CALL timeset(routineN, handle)

      IF (.NOT.PRESENT(n0)) THEN
         n = SIZE(matrix_a,1)
         m = SIZE(matrix_a,2)
         IF (m .NE. n) CPABORT("Matrix is not square")
      ELSE
         n = n0
      ENDIF

      lda = SIZE(matrix_a, 1)
      CALL dpotrf('U', n, matrix_a, lda, info)
      IF (info .NE. 0) CPABORT("Cholesky decomposition failed")

      CALL dpotri('U', n, matrix_a, lda, info)
      IF (info .NE. 0) CPABORT("Cholesky inversion failed")

      DO i = 1, n
         DO j = i+1, n
            matrix_a(j,i) = matrix_a(i,j)
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Computes the determinant of a square matrix via LU-decomposition
!> \param matrix_a ...
!> \param determinant ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE determinant_local(matrix_a, determinant, n0)

      REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)          :: matrix_a
      REAL(KIND=dp), INTENT(OUT)                         :: determinant
      INTEGER, OPTIONAL, INTENT(IN)                      :: n0

      CHARACTER(len=*), PARAMETER :: routineN = 'determinant_local'

      INTEGER                                            :: handle, i, info, n, &
                                                            m, lda
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipivot
      REAL(KIND=dp)                                      :: determinant_temp
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: a

      CALL timeset(routineN, handle)

      IF (.NOT.PRESENT(n0)) THEN
         n = SIZE(matrix_a,1)
         m = SIZE(matrix_a,2)
         IF (m .NE. n) CPABORT("Matrix is not square")
      ELSE
         n = n0
      ENDIF

      ALLOCATE(a(n,n))
      a(1:n,1:n)=matrix_a(1:n,1:n)
      
      ALLOCATE (ipivot(n))
      ipivot(:) = 0

      lda = SIZE(a, 1)
      CALL dgetrf(n, n, a, lda, ipivot, info)
      IF (info .NE. 0) CPABORT("DGETRF failed to compute LU decomposition")

      determinant_temp = 1.0_dp
      DO i = 1, n
         IF (ipivot(i) .NE. i) THEN
            determinant_temp = -determinant_temp*a(i, i)
         ELSE
            determinant_temp = determinant_temp*a(i, i)
         END IF
      END DO

      DEALLOCATE (ipivot)
      DEALLOCATE(a)

      determinant = determinant_temp

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief 
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE dbcsr_set_element(mat, val, iblkrow, irow, iblkcol, icol)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat
      REAL(KIND=dp), INTENT(IN)                          :: val
      INTEGER, INTENT(IN)                                :: iblkrow, irow, &
                                                            iblkcol, icol

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_set_element'

      INTEGER                                            :: handle, mynode, &
                                                            hold
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: p_block, p_new_block
      LOGICAL                                            :: found
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat, distribution=dist)
      CALL dbcsr_distribution_get(dist, mynode=mynode)

      CALL dbcsr_get_stored_coordinates(mat, iblkrow, iblkcol, hold)

      IF ( hold .EQ. mynode ) THEN

         CALL dbcsr_get_block_p(mat, iblkrow, iblkcol, p_block, found)

         IF (found) THEN
            p_block(irow,icol) = val
         ELSE
            ! block does not exist
            CALL dbcsr_work_create(mat, work_mutable=.TRUE.)
            NULLIFY (p_new_block)
            CALL dbcsr_reserve_block2d(mat, iblkrow, iblkcol, p_new_block)
            CPASSERT(ASSOCIATED(p_new_block))
            p_new_block(:,:) = 0.0_dp
            p_new_block(irow,icol) = val
            CALL dbcsr_finalize(mat)
         ENDIF
   
      ENDIF ! mynode

      CALL timestop(handle)

   END SUBROUTINE dbcsr_set_element

! **************************************************************************************************
!> \brief Swap columns of two dbcsr matrices
!>        The row blocks and distribution of the the matrices has to be the same.
!>        The columns to swap can be specified using either:
!>        1) global indices of the columns (1,...,n_total_columns) OR
!>        2) block numbers and indices within blocks,
!>        but not both.
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE swap_columns(matrix_A, matrix_B, global_indx_A, blk_cols_A, blk_indx_A, &
      global_indx_B, blk_cols_B, blk_indx_B, AtoB, BtoA)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: blk_cols_A, blk_indx_A, &
                                                            blk_cols_B, blk_indx_B, &
                                                            global_indx_A, &
                                                            global_indx_B
      LOGICAL, INTENT(IN), OPTIONAL                      :: AtoB, BtoA

      CHARACTER(LEN=*), PARAMETER :: routineN = 'swap_columns'

      INTEGER                                            :: handle, n_A, n_B, nmos, &
                                                            imo, iblock_row, iblock_col, &
                                                            unit_nr, &
                                                            current_sum, &
                                                            iblock, &
                                                            nblocks
      TYPE(cp_logger_type), POINTER                      :: logger
      LOGICAL                                            :: AtoB_local, BtoA_local, &
                                                            found, &
                                                            global_indx_A_convert, &
                                                            blk_indx_A_complete, &
                                                            global_indx_B_convert, &
                                                            blk_indx_B_complete
      TYPE(dbcsr_distribution_type)                      :: dist_new, dist_qs
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: blk_cols_A_local, &
                                                            blk_indx_A_local, &
                                                            blk_cols_B_local, &
                                                            blk_indx_B_local
      INTEGER, DIMENSION(:), POINTER                     :: blk_distr, blk_sizes, &
                                                            col_distr_new, col_sizes_new, &
                                                            row_distr_new, row_sizes_new
      TYPE(dbcsr_type)                                   :: matrix_A_selected_cols
      TYPE(dbcsr_type)                                   :: matrix_B_selected_cols
      TYPE(dbcsr_type)                                   :: matrix_B_selected_cols_copy
      TYPE(dbcsr_iterator_type)                          :: iter
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block, p_block

      CALL timeset(routineN, handle)

      ! RZK-WARNING: the code of this subroutine is redundant
      ! It would be better to split it into several subroutines:
      !  A. selected_columns_matrix_create
      !  B. copy_data_to_selected_columns_matrix
      !  C. copy_data_from_selected_columns_matrix

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      AtoB_local = .TRUE.
      IF (PRESENT(AtoB)) AtoB_local = AtoB

      BtoA_local = .TRUE.
      IF (PRESENT(BtoA)) BtoA_local = BtoA

      ! Check that only one type of indices is specified
      global_indx_A_convert = .FALSE.
      blk_indx_A_complete = .FALSE.
      IF (PRESENT(global_indx_A)) global_indx_A_convert = .TRUE.
      IF (PRESENT(blk_cols_A) .AND. PRESENT(blk_indx_A)) blk_indx_A_complete = .TRUE.
      IF (PRESENT(global_indx_A) .AND. ( PRESENT(blk_cols_A) .OR. PRESENT(blk_indx_A) )) THEN
         CPABORT("Cannot specify global and block indices")
      ENDIF
      IF ( .NOT. PRESENT(global_indx_A) ) THEN
         IF ( .NOT. blk_indx_A_complete ) THEN
            CPABORT("Specify global or 2 block indices for matrixA")
         ENDIF
      ENDIF

      ! Check that only one type of indices is specified
      global_indx_B_convert = .FALSE.
      blk_indx_B_complete = .FALSE.
      IF (PRESENT(global_indx_B)) global_indx_B_convert = .TRUE.
      IF (PRESENT(blk_cols_B) .AND. PRESENT(blk_indx_B)) blk_indx_B_complete = .TRUE.
      IF (PRESENT(global_indx_B) .AND. ( PRESENT(blk_cols_B) .OR. PRESENT(blk_indx_B) )) THEN
         CPABORT("Cannot specify global and block indices")
      ENDIF
      IF ( .NOT. PRESENT(global_indx_B) ) THEN
         IF ( .NOT. blk_indx_B_complete ) THEN
            CPABORT("Specify global or 2 block indices for matrixB")
         ENDIF
      ENDIF

      CALL dbcsr_get_info(matrix_A, nblkrows_total=n_A)
      CALL dbcsr_get_info(matrix_B, nblkrows_total=n_B)
      CPASSERT(n_B.EQ.n_A)

      CALL dbcsr_get_info(matrix_A, nfullrows_total=n_A)
      CALL dbcsr_get_info(matrix_B, nfullrows_total=n_B)
      CPASSERT(n_B.EQ.n_A)

      ! Create block column indices for A
      IF (global_indx_A_convert) THEN

         nmos = SIZE(global_indx_A)

         ALLOCATE(blk_cols_A_local(nmos))
         ALLOCATE(blk_indx_A_local(nmos))

         CALL dbcsr_get_info(matrix_A, row_blk_size=blk_sizes)

         nblocks = SIZE(blk_sizes)

         ! Convert global indices to block indices
         MOloopA: DO imo = 1, nmos
            current_sum = 0
            DO iblock = 1, nblocks
               current_sum = current_sum + blk_sizes(iblock)
               IF (global_indx_A(imo) .LE. current_sum) THEN
                  blk_cols_A_local(imo) = iblock
                  blk_indx_A_local(imo) = global_indx_A(imo) - (current_sum - blk_sizes(iblock))
                  CYCLE MOloopA
               ENDIF
            ENDDO
            CPABORT("Global index not converted for A")
         ENDDO MOloopA

      ELSE

         CPASSERT(SIZE(blk_cols_A).EQ.SIZE(blk_indx_A))
         nmos = SIZE(blk_cols_A)

         ALLOCATE(blk_cols_A_local(nmos))
         ALLOCATE(blk_indx_A_local(nmos))

         blk_cols_A_local(:) = blk_cols_A(:)
         blk_indx_A_local(:) = blk_indx_A(:)

      ENDIF

      ! Create block column indices for B
      IF (global_indx_B_convert) THEN

         nmos = SIZE(global_indx_B)

         ALLOCATE(blk_cols_B_local(nmos))
         ALLOCATE(blk_indx_B_local(nmos))

         CALL dbcsr_get_info(matrix_B, row_blk_size=blk_sizes)
         
         nblocks = SIZE(blk_sizes)

         ! Convert global indices to block indices
         MOloopB: DO imo = 1, nmos
            current_sum = 0
            DO iblock = 1, nblocks
               current_sum = current_sum + blk_sizes(iblock)
               IF (global_indx_B(imo) .LE. current_sum) THEN
                  blk_cols_B_local(imo) = iblock
                  blk_indx_B_local(imo) = global_indx_B(imo) - (current_sum - blk_sizes(iblock))
                  CYCLE MOloopB
               ENDIF
            ENDDO
            CPABORT("Global index not converted for B")
         ENDDO MOloopB

      ELSE

         CPASSERT(SIZE(blk_cols_B).EQ.SIZE(blk_indx_B))
         nmos = SIZE(blk_cols_B)

         ALLOCATE(blk_cols_B_local(nmos))
         ALLOCATE(blk_indx_B_local(nmos))

         blk_cols_B_local(:) = blk_cols_B(:)
         blk_indx_B_local(:) = blk_indx_B(:)

      ENDIF

      CPASSERT(SIZE(blk_cols_A_local).EQ.SIZE(blk_cols_B_local))
      CPASSERT(SIZE(blk_indx_A_local).EQ.SIZE(blk_indx_B_local))

      CALL dbcsr_get_info(matrix_A, distribution=dist_qs)

      ! structure of the AO dimension can be copied from one of the matrices
      CALL dbcsr_distribution_get(dist_qs, row_dist=blk_distr)
      CALL dbcsr_get_info(matrix_A, row_blk_size=blk_sizes)
      ALLOCATE (row_sizes_new(SIZE(blk_sizes)))
      ALLOCATE (row_distr_new(SIZE(blk_distr)))
      row_sizes_new(:) = blk_sizes(:)
      row_distr_new(:) = blk_distr(:)

      ! structure for the MO dimension is created from scratch for each matrix
      ALLOCATE (col_sizes_new(nmos))
      ALLOCATE (col_distr_new(nmos))
      col_sizes_new = (/ (1, imo = 1, nmos) /) ! array of 1s, length nmos

      ! now, data specific to matrix_A
      !col_distr_new = (cores1D of all swapped columns) - use get_stored_coordinates for each MO block
      DO imo = 1, nmos
         CALL dbcsr_get_stored_coordinates(matrix_A, 1, blk_cols_A_local(imo), col_distr_new(imo))
      ENDDO

      ! Create the distribution and then the matrix of selected columns
      CALL dbcsr_distribution_new(dist_new, template=dist_qs, &
                                  row_dist=row_distr_new, col_dist=col_distr_new, &
                                  reuse_arrays=.TRUE.)
      CALL dbcsr_create(matrix_A_selected_cols, "A_cols", &
                        dist_new, dbcsr_type_no_symmetry, &
                        row_sizes_new, col_sizes_new, reuse_arrays=.TRUE.)
      CALL dbcsr_finalize(matrix_A_selected_cols)
      CALL dbcsr_distribution_release(dist_new)

!WRITE(*,*)"DONE CREATE MATRIX OF THE SLECTED COLUMNS MAT_A (mat_T)"
!CALL dbcsr_print(matrix_A_selected_cols)

      !CINDY: add the AO copy from mat_B
      ! structure of the AO dimension can be copied from one of the matrices
      CALL dbcsr_distribution_get(dist_qs, row_dist=blk_distr)
      CALL dbcsr_get_info(matrix_B, row_blk_size=blk_sizes)
      ALLOCATE(row_sizes_new(SIZE(blk_sizes)))
      ALLOCATE(row_distr_new(SIZE(blk_distr)))
      row_sizes_new(:) = blk_sizes(:)
      row_distr_new(:) = blk_distr(:)

      !CINDY: add the MO for mat_B 
      ! structure for the MO dimension is created from scratch for each matrix
      ALLOCATE(col_sizes_new(nmos))
      ALLOCATE(col_distr_new(nmos))
      col_sizes_new = (/ (1, imo = 1, nmos) /) ! array of 1s, length nmos

      ! now, data specific to matrix_B
      DO imo = 1, nmos
         CALL dbcsr_get_stored_coordinates(matrix_B, 1, blk_cols_B_local(imo), col_distr_new(imo))
      ENDDO

      ! Create the distribution and then the matrix of selected columns
      CALL dbcsr_distribution_new(dist_new, template=dist_qs, &
                                  row_dist=row_distr_new, col_dist=col_distr_new, &
                                  reuse_arrays=.TRUE.)
      CALL dbcsr_create(matrix_B_selected_cols, "B_cols", &
                        dist_new, dbcsr_type_no_symmetry, &
                        row_sizes_new, col_sizes_new, reuse_arrays=.TRUE.)
      CALL dbcsr_finalize(matrix_B_selected_cols)
      CALL dbcsr_distribution_release(dist_new)

!WRITE(*,*)"DONE CREATE MATRIX OF THE SLECTED COLUMNS MAT_B (mat_V)"
!CALL dbcsr_print(matrix_B_selected_cols)

! CINDY: Here if I deallocate these arrays there will be segmentation error?!
!      DEALLOCATE(row_sizes_new)
!      DEALLOCATE(row_distr_new)
!      DEALLOCATE(col_sizes_new)
!      DEALLOCATE(col_distr_new)

      ! copy data into the newly created matrix of selected columns (same dist as the source)
      ! then redistribute into the other matrix of selected columns (same dist as the destination)
      ! copy data from the other matrix into the destination
      IF (AtoB_local) THEN

         CALL dbcsr_work_create(matrix_A_selected_cols, work_mutable=.TRUE.)
         CALL dbcsr_iterator_start(iter, matrix_A)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p)

            DO imo = 1, nmos

               IF (blk_cols_A_local(imo) .EQ. iblock_col) THEN
               ! this block contains needed info

                  NULLIFY (p_new_block)
                  CALL dbcsr_reserve_block2d(matrix_A_selected_cols, iblock_row, imo, p_new_block)
                  CPASSERT(ASSOCIATED(p_new_block))
                  p_new_block(:,1) = data_p(:,blk_indx_A_local(imo))

               ENDIF ! mynode

            ENDDO

         ENDDO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_finalize(matrix_A_selected_cols)

      ENDIF ! AtoB

      IF (BtoA_local) THEN

         CALL dbcsr_work_create(matrix_B_selected_cols, work_mutable=.TRUE.)
         CALL dbcsr_iterator_start(iter, matrix_B)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p)

            DO imo = 1, nmos

               IF (blk_cols_B_local(imo) .EQ. iblock_col) THEN
               ! this block contains needed info

                  NULLIFY (p_new_block)
                  CALL dbcsr_reserve_block2d(matrix_B_selected_cols, iblock_row, imo, p_new_block)
                  CPASSERT(ASSOCIATED(p_new_block))
                  p_new_block(:,1) = data_p(:,blk_indx_B_local(imo))

               ENDIF ! mynode

            ENDDO

         ENDDO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_finalize(matrix_B_selected_cols)

         CALL dbcsr_create(matrix_B_selected_cols_copy, &
                 template = matrix_B_selected_cols)
         CALL dbcsr_copy(matrix_B_selected_cols_copy, matrix_B_selected_cols)

!WRITE(*,*) "DONE COPY DATA INTO THE NEW MAT_B (mat_V)"   
!CALL dbcsr_print(matrix_B_selected_cols_copy)
      ENDIF ! BtoA

      ! REturn to the next step of the AtoB process
      IF (AtoB_local) THEN

         ! redistribute data so it can be copied into matrix_B
         CALL dbcsr_set(matrix_B_selected_cols, 0.0_dp)
         CALL dbcsr_filter(matrix_B_selected_cols,HUGE(1.0_dp))
         CALL dbcsr_complete_redistribute( &
                 matrix_A_selected_cols, &
                 matrix_B_selected_cols, &
                 keep_sparsity=.FALSE.)

         ! now the data can be copied to matrix_B
         CALL dbcsr_work_create(matrix_B, work_mutable=.TRUE.)
         CALL dbcsr_iterator_start(iter, matrix_B_selected_cols)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p)

            CALL dbcsr_get_block_p(matrix_B, iblock_row, blk_cols_B_local(iblock_col), p_block, found)

            IF (found) THEN
               p_block(:, blk_indx_B_local(iblock_col)) = data_p(:, 1)
            ELSE
               ! block does not exist
               !CALL dbcsr_work_create(mat, work_mutable=.TRUE.)
               NULLIFY (p_new_block)
               CALL dbcsr_reserve_block2d(matrix_B, iblock_row, blk_cols_B_local(iblock_col), p_new_block)
               CPASSERT(ASSOCIATED(p_new_block))
               p_new_block(:,:) = 0.0_dp
               p_new_block(:, blk_indx_B_local(iblock_col)) = data_p(:, 1)
               !CALL dbcsr_finalize(mat)
            ENDIF

         ENDDO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_finalize(matrix_B)

      ENDIF ! AtoB

      ! Copy B to A if necessary
      IF (BtoA_local) THEN

         ! redistribute data so it can be copied into matrix_A
         CALL dbcsr_set(matrix_A_selected_cols, 0.0_dp)
         CALL dbcsr_filter(matrix_A_selected_cols,HUGE(1.0_dp))
         CALL dbcsr_complete_redistribute( &
                 matrix_B_selected_cols_copy, &
                 matrix_A_selected_cols, &
                 keep_sparsity=.FALSE.)

         ! now the data can be copied to matrix_A
         CALL dbcsr_work_create(matrix_A, work_mutable=.TRUE.)
         CALL dbcsr_iterator_start(iter, matrix_A_selected_cols)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p)

            CALL dbcsr_get_block_p(matrix_A, iblock_row, blk_cols_A_local(iblock_col), p_block, found)

            IF (found) THEN
               p_block(:, blk_indx_A_local(iblock_col)) = data_p(:, 1)
            ELSE
               ! block does not exist
               !CALL dbcsr_work_create(mat, work_mutable=.TRUE.)
               NULLIFY (p_new_block)
               CALL dbcsr_reserve_block2d(matrix_A, iblock_row, blk_cols_A_local(iblock_col), p_new_block)
               CPASSERT(ASSOCIATED(p_new_block))
               p_new_block(:,:) = 0.0_dp
               p_new_block(:, blk_indx_A_local(iblock_col)) = data_p(:, 1)
               !CALL dbcsr_finalize(mat)
            ENDIF

         ENDDO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_finalize(matrix_A)
!WRITE(*,*)"DONE SWAPPING VIRTUAL COLS MAT_V TO OCCUPIED COLS MAT_T"

      ENDIF ! BtoA

      DEALLOCATE(blk_cols_A_local)
      DEALLOCATE(blk_indx_A_local)
      DEALLOCATE(blk_cols_B_local)
      DEALLOCATE(blk_indx_B_local)

      CALL dbcsr_release(matrix_A_selected_cols)
      CALL dbcsr_release(matrix_B_selected_cols)

      CALL timestop(handle)

   END SUBROUTINE swap_columns

! **************************************************************************************************
!> \brief Join columns of two matrices to produce a new matrix.
!>        The row blocks and distribution of the the three matrices has
!>        to be the same.
!> \par History
!>       2021.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE stitch_columns(matrix_left, matrix_right, matrix_stitched)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_left, matrix_right
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_stitched

      CHARACTER(LEN=*), PARAMETER :: routineN = 'stitch_columns'

      INTEGER                                            :: handle, iblock_col, iblock_row, &
                                                            nblocks_offset, n_left, &
                                                            n_right, n_combined, hold, &
                                                            mynode
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block
      TYPE(dbcsr_iterator_type)                          :: iter, iter_right
      TYPE(dbcsr_type)                                   :: matrix_intermediate
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_left, nblkrows_total=n_left)
      CALL dbcsr_get_info(matrix_right, nblkrows_total=n_right)
      CPASSERT(n_left.EQ.n_right)
      CALL dbcsr_get_info(matrix_stitched, nblkrows_total=n_combined)
      CPASSERT(n_left.EQ.n_combined)

      CALL dbcsr_get_info(matrix_left, nfullrows_total=n_left)
      CALL dbcsr_get_info(matrix_right, nfullrows_total=n_right)
      CPASSERT(n_left.EQ.n_right)
      CALL dbcsr_get_info(matrix_stitched, nfullrows_total=n_combined)
      CPASSERT(n_left.EQ.n_combined)

      CALL dbcsr_get_info(matrix_left, nfullcols_total=n_left)
      CALL dbcsr_get_info(matrix_right, nfullcols_total=n_right)
      CALL dbcsr_get_info(matrix_stitched, nfullcols_total=n_combined)
      CPASSERT((n_left + n_right) .EQ. n_combined)

      CALL dbcsr_get_info(matrix_left, nblkcols_total=nblocks_offset)

      CALL create_matrix_column_stitched(matrix_left, matrix_right, matrix_intermediate)
      CALL dbcsr_get_info(matrix_intermediate, distribution=dist)
      CALL dbcsr_distribution_get(dist, mynode=mynode)
      CALL dbcsr_work_create(matrix_intermediate, work_mutable=.TRUE.)

      CALL dbcsr_iterator_start(iter, matrix_left)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p)

         CALL dbcsr_get_stored_coordinates(matrix_intermediate, iblock_row, iblock_col, hold)

         IF (hold .EQ. mynode) THEN
            NULLIFY (p_new_block)
            CALL dbcsr_reserve_block2d(matrix_intermediate, iblock_row, iblock_col, p_new_block)
            CPASSERT(ASSOCIATED(p_new_block))
            p_new_block = data_p
         ELSE
            CPABORT("Left part of the intermediate matrix is incorrect")
         ENDIF ! mynode

      ENDDO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_iterator_start(iter_right, matrix_right)

      DO WHILE (dbcsr_iterator_blocks_left(iter_right))

         CALL dbcsr_iterator_next_block(iter_right, iblock_row, iblock_col, data_p)

         CALL dbcsr_get_stored_coordinates(matrix_intermediate, iblock_row, iblock_col, hold)

         IF (hold .EQ. mynode) THEN
            NULLIFY (p_new_block)
            CALL dbcsr_reserve_block2d(matrix_intermediate, iblock_row, &
                    iblock_col + nblocks_offset, p_new_block)
            CPASSERT(ASSOCIATED(p_new_block))
            p_new_block = data_p
         ELSE
            CPABORT("Right part of the intermediate matrix is incorrect")
         ENDIF ! mynode

      ENDDO
      CALL dbcsr_iterator_stop(iter_right)

      CALL dbcsr_finalize(matrix_intermediate)
      !CALL dbcsr_filter(matrix_intermediate, eps_filter)

      !CALL dbcsr_print(matrix_intermediate)
      CALL dbcsr_set(matrix_stitched, 0.0_dp)
      CALL dbcsr_filter(matrix_stitched,HUGE(1.0_dp))
      CALL dbcsr_complete_redistribute(matrix_intermediate, matrix_stitched, &
              keep_sparsity=.FALSE.)

      CALL dbcsr_release(matrix_intermediate)

      CALL timestop(handle)

   END SUBROUTINE stitch_columns

! **************************************************************************************************
!> \brief Split columns of a matrix to produce two matrices.
!>        The row blocks and distribution of the the three matrices has
!>        to be the same.
!> \par History
!>       2021.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE split_columns(matrix_stitched, matrix_left, matrix_right)

      TYPE(dbcsr_type), INTENT(IN)                     :: matrix_stitched
      TYPE(dbcsr_type), INTENT(INOUT)                  :: matrix_left, matrix_right

      CHARACTER(LEN=*), PARAMETER :: routineN = 'split_columns'

      INTEGER                                            :: handle, iblock_col, iblock_row, &
                                                            nblocks_offset, n_left, &
                                                            n_right, n_combined
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block
      TYPE(dbcsr_iterator_type)                          :: iter, iter_right
      TYPE(dbcsr_type)                                   :: matrix_intermediate

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_left, nblkrows_total=n_left)
      CALL dbcsr_get_info(matrix_right, nblkrows_total=n_right)
      CPASSERT(n_left.EQ.n_right)
      CALL dbcsr_get_info(matrix_stitched, nblkrows_total=n_combined)
      CPASSERT(n_left.EQ.n_combined)

      CALL dbcsr_get_info(matrix_left, nfullrows_total=n_left)
      CALL dbcsr_get_info(matrix_right, nfullrows_total=n_right)
      CPASSERT(n_left.EQ.n_right)
      CALL dbcsr_get_info(matrix_stitched, nfullrows_total=n_combined)
      CPASSERT(n_left.EQ.n_combined)

      CALL dbcsr_get_info(matrix_left, nfullcols_total=n_left)
      CALL dbcsr_get_info(matrix_right, nfullcols_total=n_right)
      CALL dbcsr_get_info(matrix_stitched, nfullcols_total=n_combined)
      CPASSERT((n_left + n_right) .EQ. n_combined)

      CALL create_matrix_column_stitched(matrix_left, matrix_right, matrix_intermediate)
      CALL dbcsr_set(matrix_intermediate, 0.0_dp)
      CALL dbcsr_complete_redistribute(matrix_stitched, matrix_intermediate, &
              keep_sparsity=.FALSE.)

      CALL dbcsr_get_info(matrix_left, nblkcols_total=nblocks_offset)
      CALL dbcsr_set(matrix_left, 0.0_dp)
      CALL dbcsr_set(matrix_right, 0.0_dp)
      CALL dbcsr_filter(matrix_left,HUGE(1.0_dp))
      CALL dbcsr_filter(matrix_right,HUGE(1.0_dp))
      CALL dbcsr_work_create(matrix_left, work_mutable=.TRUE.)
      CALL dbcsr_work_create(matrix_right, work_mutable=.TRUE.)

      CALL dbcsr_iterator_start(iter, matrix_intermediate)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p)

         NULLIFY (p_new_block)
         IF (iblock_col.LE.nblocks_offset) THEN
            CALL dbcsr_reserve_block2d(matrix_left, iblock_row, iblock_col, p_new_block)
         ELSE
            CALL dbcsr_reserve_block2d(matrix_right, iblock_row, iblock_col-nblocks_offset, p_new_block)
         ENDIF
         CPASSERT(ASSOCIATED(p_new_block))
         p_new_block = data_p

      ENDDO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_finalize(matrix_left)
      CALL dbcsr_finalize(matrix_right)

      CALL dbcsr_release(matrix_intermediate)

      CALL timestop(handle)

   END SUBROUTINE split_columns

   SUBROUTINE create_matrix_column_stitched(matrix_left, matrix_right, matrix_stitched)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_left, matrix_right
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_stitched

      INTEGER                                            :: dimen, handle, hold, iatom, iblock_col, &
                                                            iblock_row, imol, mynode, natoms, &
                                                            nblkrows_tot, nlength, nmols, row, &
                                                            ncores
      INTEGER, DIMENSION(:), POINTER                     :: blk_distr_row, blk_distr_left, &
                                                            blk_distr_right, blk_sizes_row, &
                                                            blk_sizes_left, blk_sizes_right, &
                                                            block_sizes_row_new, &
                                                            block_distr_row_new, &
                                                            block_sizes_col_new, &
                                                            block_distr_col_new
      TYPE(dbcsr_distribution_type)                      :: dist_right, dist_left, dist_new

      ! obtain distribution from the left and right matrices
      CALL dbcsr_get_info(matrix_left, distribution=dist_left)
      CALL dbcsr_get_info(matrix_right, distribution=dist_right)

      CALL dbcsr_distribution_get(dist_left, row_dist=blk_distr_row)
      CALL dbcsr_distribution_get(dist_left, col_dist=blk_distr_left)
      CALL dbcsr_distribution_get(dist_right, col_dist=blk_distr_right)

      CALL dbcsr_get_info(matrix_left, row_blk_size=blk_sizes_row)
      CALL dbcsr_get_info(matrix_left, col_blk_size=blk_sizes_left)
      CALL dbcsr_get_info(matrix_right, col_blk_size=blk_sizes_right)

      ALLOCATE (block_sizes_row_new(SIZE(blk_sizes_row)))
      ALLOCATE (block_distr_row_new(SIZE(blk_distr_row)))
      ALLOCATE (block_sizes_col_new(SIZE(blk_sizes_left)+SIZE(blk_sizes_right)))
      ALLOCATE (block_distr_col_new(SIZE(blk_distr_left)+SIZE(blk_distr_right)))
      block_sizes_row_new(:) = blk_sizes_row(:)
      block_distr_row_new(:) = blk_distr_row(:)
      block_sizes_col_new(1:SIZE(blk_sizes_left)) = blk_sizes_left(:)
      block_sizes_col_new(SIZE(blk_sizes_left)+1:) = blk_sizes_right(:)
      block_distr_col_new(1:SIZE(blk_distr_left)) = blk_distr_left(:)
      block_distr_col_new(SIZE(blk_distr_left)+1:) = blk_distr_right(:)

      ! Create the distribution
      CALL dbcsr_distribution_new(dist_new, template=dist_left, &
              row_dist=block_distr_row_new, col_dist=block_distr_col_new, &
                                  reuse_arrays=.TRUE.)

      ! Create the matrix
      CALL dbcsr_create(matrix_stitched, 'stitched', &
                        dist_new, dbcsr_type_no_symmetry, &
                        block_sizes_row_new, block_sizes_col_new, reuse_arrays=.TRUE.)
      CALL dbcsr_finalize(matrix_stitched)

      CALL dbcsr_distribution_release(dist_new)

      IF (ASSOCIATED(block_sizes_row_new)) DEALLOCATE(block_sizes_row_new)
      IF (ASSOCIATED(block_distr_row_new)) DEALLOCATE(block_distr_row_new)
      IF (ASSOCIATED(block_sizes_col_new)) DEALLOCATE(block_sizes_col_new)
      IF (ASSOCIATED(block_distr_col_new)) DEALLOCATE(block_distr_col_new)

   END SUBROUTINE create_matrix_column_stitched

! **************************************************************************************************
!> \brief Orthonormalize main variables
!> \param qs_env ...
!> \param tidft_env ...
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_mainvar_orthonormalize(tidft_env, only_state)
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env
      INTEGER, INTENT(IN), OPTIONAL                      :: only_state

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_mainvar_orthonormalize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iaspc, ispin, istore, naspc, &
                                                            nspins, unit_nr, nmos_stateI, &
                                                            stateI, ivar, nstates, idet, &
                                                            ndets, slaterI, stateJ, &
                                                            nmos_stateJ, my_state_start, &
                                                            my_state_end
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      LOGICAL                                            :: aspc_guess, has_unit_metric
      REAL(KIND=dp)                                      :: alpha, cs_pos, energy, kTS_sum, &
                                                            random_percentile, &
                                                            energy_singlet, energy_triplet, &
                                                            state_energy, obj_energy, &
                                                            interstate_nonorthogonality_penalty, &
                                                            intrastate_nonorthogonality_penalty, &
                                                            DetDiagonalSigma, penalty_stateI, &
                                                            determinant_spd
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: tempNO1, mat_proj_var1
      TYPE(dbcsr_type)                                   :: temp2, mat_var23, mat_var23_sigma
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: fake_e_per_block
      INTEGER                                            :: group, nsize
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diagonal

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! to re-use some almo_scf subroutines set
      ! nocc_of_domain to the negative number
      ALLOCATE(fake_e_per_block(1))
      fake_e_per_block(:) = -1

      nspins = tidft_env%nspins
      nstates = tidft_env%nstates
      ndets=1
      IF (tidft_env%triplet_correction) ndets = 2

      my_state_start = 1
      my_state_end = nstates
      IF (PRESENT(only_state)) THEN
         my_state_start = only_state
         my_state_end = only_state
      ENDIF

      ! matrices of independent variables
      ! unrestricted case: no. of electron pairs is zero
      DO stateI = my_state_start, my_state_end

         CALL dbcsr_create(mat_proj_var1, template=tidft_env%mat_s)

         DO ivar = 1, 3

            ! Orthogonalize vectors describing variables in several steps
            ! * Orthogonalize doubly occupied orbitals
            ! * Create the projector R1 onto the doubly occupied subspace
            ! * Project out the doubly occupied subspace from both singly-occupied
            !   alpha and singly-occupied beta subspaces
            ! * Orthogonalize singly occupied orbitals

            nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

            IF ( nmos_stateI .GT. 0 ) THEN

               ! Project out the doubly occupied subspace if any
               IF (ivar .NE. 1) THEN

                  IF (tidft_env%e_states(stateI)%slater_det(1)%ncols_var(1) .NE. 0) THEN

                     CALL dbcsr_create(tempNO1, &
                        template=tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))

                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   tidft_env%mat_s, &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                                   0.0_dp, tempNO1, &
                                   filter_eps=tidft_env%eps_filter)

                     CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                   mat_proj_var1, &
                                   tempNO1, &
                                   1.0_dp, tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                                   filter_eps=tidft_env%eps_filter)

                     CALL dbcsr_release(tempNO1)

                  ENDIF

               ENDIF

               IF (tidft_env%triplet_correction .AND. ivar .GT. 1) THEN

                  ! The singly-occupied alpha and beta electrons must be orthogonalized simultaneusly
                  ! Do this for ivar.EQ.3, when R1 is projected out from both singly-occupied states

                  ! Do nothing for ivar.EQ.2
                  IF (ivar .EQ. 3) THEN

                     CALL matrix_layout( &
                             matrix_ao_template = tidft_env%mat_s, &
                             matrix_new = mat_var23 , &
                             name_new = 'm_var23', &
                             orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                             suggested_mo_blk_size = 2, &
                             symmetry_new = dbcsr_type_no_symmetry, &
                             nmos = (/-1, 2/) ) ! 2 singly-occupied orbitals - alpha and beta
                     CALL matrix_layout( &
                             matrix_ao_template = tidft_env%mat_s, &
                             matrix_new = mat_var23_sigma , &
                             name_new = 'm_var23_sigma', &
                             orbital_type_new = (/almo_mat_dim_occ, almo_mat_dim_occ/), &
                             suggested_mo_blk_size = 2, &
                             symmetry_new = dbcsr_type_no_symmetry, &
                             nmos = (/2, 2/) ) ! these are not MOs but independent vars

                     CALL stitch_columns( &
                             tidft_env%e_states(stateI)%slater_det(1)%mat_var(2), &
                             tidft_env%e_states(stateI)%slater_det(1)%mat_var(3), &
                             mat_var23 &
                          )

                     CALL orthogonalize_mos(ket=mat_var23, &
                                            overlap=mat_var23_sigma, &
                                            metric=tidft_env%mat_s, &
                                            retain_locality=.FALSE., &
                                            only_normalize=.FALSE., &
                                            nocc_of_domain=fake_e_per_block, &
                                            eps_filter=tidft_env%eps_filter, &
                                            order_lanczos=tidft_env%order_lanczos, &
                                            eps_lanczos=tidft_env%eps_lanczos, &
                                            max_iter_lanczos=tidft_env%max_iter_lanczos, &
                                            smear=tidft_env%smear)
                     CALL split_columns( &
                             mat_var23, &
                             tidft_env%e_states(stateI)%slater_det(1)%mat_var(2), &
                             tidft_env%e_states(stateI)%slater_det(1)%mat_var(3) &
                          )
                     CALL dbcsr_release(mat_var23)
                     CALL dbcsr_release(mat_var23_sigma)

                  ENDIF

               ELSE ! not triplet_correction

                  CALL orthogonalize_mos(ket=tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                                         overlap=tidft_env%mat_var_sigma(ivar), &
                                         metric=tidft_env%mat_s, &
                                         retain_locality=.FALSE., &
                                         only_normalize=.FALSE., &
                                         nocc_of_domain=fake_e_per_block, &
                                         eps_filter=tidft_env%eps_filter, &
                                         order_lanczos=tidft_env%order_lanczos, &
                                         eps_lanczos=tidft_env%eps_lanczos, &
                                         max_iter_lanczos=tidft_env%max_iter_lanczos, &
                                         smear=tidft_env%smear)

               ENDIF ! triplet correction

               IF (ivar .EQ. 1) THEN
                  CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                                tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                                0.0_dp, mat_proj_var1, &
                                filter_eps=tidft_env%eps_filter)
               ENDIF

            ENDIF ! non-zero columns

         ENDDO ! end ivar

         CALL dbcsr_release(mat_proj_var1)

      ENDDO  ! end stateI

      DEALLOCATE(fake_e_per_block)

      CALL timestop(handle)

   END SUBROUTINE tidft_mainvar_orthonormalize

! **************************************************************************************************
!> \brief IF (post_multiply == TRUE)
!>           temp = in * projector_orbitals * sigma_inv(projector_orbitals)
!>        ELSE
!>           temp = in
!>        ENDIF
!>        IF (covariant_input == TRUE)
!>           out = (I - S * projector(projector_orbitals)) * temp
!>        ELSE
!>           out = (I - projector(projector_orbitals) * S) * temp
!>        ENDIF
!>        Matrices should be compatible in size, blocks, distribution
!>        If projector_orbitals is NxK, then
!>           template_sigma is KxK
!>           mat_S is NxN
!>           mat_in is Nx(N/M) if post_multiply == T/F
!>           mat_out is Nx(K/M) if post_multiply == T/F
!>           M - is any integer
!> \param project_out...
!> \par History
!>       2023.01 created [Cindy Pham]
!> \author Cindy Pham, Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE project_out(mat_in, mat_out, projector_orbitals, template_sigma, &
                 mat_s, covariant_input, eps_filter, para_env, blacs_env, post_multiply)

      TYPE(dbcsr_type), INTENT(IN)                                 :: mat_in
      TYPE(dbcsr_type), INTENT(INOUT)                              :: mat_out
      TYPE(dbcsr_type), INTENT(IN)                                 :: projector_orbitals, &
                                                                      template_sigma, mat_s
      LOGICAL                                                      :: covariant_input
      REAL(KIND=dp), INTENT(IN)                                    :: eps_filter
      TYPE(cp_para_env_type), POINTER                              :: para_env
      TYPE(cp_blacs_env_type), POINTER                             :: blacs_env
      LOGICAL, OPTIONAL                                            :: post_multiply

      CHARACTER(LEN=*), PARAMETER :: routineN = 'project_out'

      TYPE(dbcsr_type)                                             :: tempNN, ST
      INTEGER, DIMENSION(:), ALLOCATABLE                           :: fake_e_per_block
      TYPE(dbcsr_type)                                             :: sigma, sigma_inv, &
                                                                      theta_siginv
      LOGICAL                                                      :: my_post_multiply

      my_post_multiply = .FALSE.
      IF(PRESENT(post_multiply)) THEN
         my_post_multiply = post_multiply
      ENDIF

      ! Re-use some almo_scf subroutines set nocc_of_domain to the negative number
      ALLOCATE(fake_e_per_block(1))
      fake_e_per_block(:) = -1

      ! mat_in = T
      ! mat_out = T~
      ! projector_orbital = T'

      CALL dbcsr_create(tempNN, template=mat_s)
      CALL dbcsr_create(sigma, template = template_sigma)
      CALL dbcsr_create(sigma_inv, template = template_sigma)

      ! R = T'.(sigma_inv).T'^t
      CALL almo_scf_t_to_proj(t=projector_orbitals, &
                              p=tempNN, &
                              eps_filter=eps_filter, &
                              nocc_of_domain= fake_e_per_block, &
                              s=mat_s, &
                              sigma = sigma,&
                              sigma_inv = sigma_inv,&
                              algorithm = 1, & ! use Cholesky
                              orthog_orbs=.FALSE., &
                              use_guess=.FALSE., &
                              smear=.FALSE., &
                              para_env=para_env, &
                              blacs_env=blacs_env)

      IF (my_post_multiply) THEN

         CALL dbcsr_create(theta_siginv, template = projector_orbitals)

         ! NxK = NxK.KxK
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                              projector_orbitals, &
                              sigma_inv, &
                              0.0_dp, theta_siginv, &
                              filter_eps=eps_filter)

         ! NxK = NxN.NxK
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                              mat_in, &
                              theta_siginv, &
                              0.0_dp, mat_out, &
                              filter_eps=eps_filter)

         CALL dbcsr_release(theta_siginv)

      ELSE

         CALL dbcsr_copy(mat_out, mat_in)

      ENDIF

      CALL dbcsr_release(sigma)
      CALL dbcsr_release(sigma_inv)
      CALL dbcsr_create(ST, template = mat_out)

      IF (covariant_input) THEN
         ! T~ = T - S.R.T

         ! R.T
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                              tempNN, &
                              mat_out, &
                              0.0_dp, ST, &
                              filter_eps=eps_filter)

         ! S.RT
         CALL dbcsr_multiply("N", "N", -1.0_dp, &
                              mat_s, &
                              ST, &
                              1.0_dp, mat_out, &
                              filter_eps=eps_filter)

      ELSE
         ! T~ = T - R.S.T

         ! S.T
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                              mat_s, &
                              mat_out, &
                              0.0_dp, ST, &
                              filter_eps=eps_filter)

         ! R.ST
         CALL dbcsr_multiply("N", "N", -1.0_dp, &
                              tempNN, &
                              ST, &
                              1.0_dp, mat_out, &
                              filter_eps=eps_filter)

      ENDIF

      CALL dbcsr_release(ST)

      DEALLOCATE(fake_e_per_block)

      CALL dbcsr_release(tempNN)

   END SUBROUTINE project_out

END MODULE tidft_methods

