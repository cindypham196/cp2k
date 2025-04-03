!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Interface to (sca)lapack for the Cholesky based procedures
!> \author  Rustam Z. Khaliullin
! **************************************************************************************************
MODULE cp_dbcsr_singular

   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_dbcsr_operations,             ONLY: copy_cfm_to_dbcsr,&
                                              copy_dbcsr_to_cfm,&
                                              copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_power,&
                                              cp_fm_syevx
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE cp_dbcsr_api,                    ONLY: &
                                              dbcsr_create, dbcsr_multiply, dbcsr_add_on_diag, &
                                              dbcsr_release, dbcsr_set, dbcsr_type, &
                                              dbcsr_type_no_symmetry, dbcsr_get_info, dbcsr_get_diag, &
                                              dbcsr_print, dbcsr_set_diag, dbcsr_dot, dbcsr_desymmetrize
   USE kinds,                           ONLY: dp
   USE cp_fm_basic_linalg,              ONLY: cp_fm_lu_decompose,&
                                              cp_fm_adjugate
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_singular'

   ! Public subroutines

   PUBLIC :: cp_dbcsr_pseudoinverse, &
             cp_dbcsr_adjugate, &
             cp_dbcsr_determinant_generalized

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param para_env ...
!> \param blacs_env ...
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_pseudoinverse(matrix, zero_eps, min_eigenvalues_to_nullify, &
                 spd_expected, para_env, blacs_env)

      TYPE(dbcsr_type)                                   :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: zero_eps
      LOGICAL, INTENT(IN)                                :: spd_expected
      INTEGER, OPTIONAL, INTENT(IN)                      :: min_eigenvalues_to_nullify
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_pseudoinverse'

      INTEGER                                            :: handle, &
                                                            ncols, ieigen, &
                                                            nrows, &
                                                            my_min_eigenvalues_to_nullify
      TYPE(dbcsr_type)                                   :: eigenvectors, &
                                                            temp_mat_diag, &
                                                            temp_mat
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: eigenvalues

      CALL timeset(routineN, handle)

      my_min_eigenvalues_to_nullify = 0
      IF (PRESENT(min_eigenvalues_to_nullify)) &
         my_min_eigenvalues_to_nullify = min_eigenvalues_to_nullify

      CALL dbcsr_create(eigenvectors, &
              template=matrix, matrix_type = dbcsr_type_no_symmetry)

      CALL dbcsr_get_info(matrix, nfullrows_total=nrows, nfullcols_total=ncols)
      CPASSERT(nrows.EQ.ncols)
      ALLOCATE(eigenvalues(nrows))

      CALL cp_dbcsr_syevd(matrix, eigenvectors, eigenvalues, para_env, blacs_env)

      ! The code assumes that eigenvalues are sorted from lowest(1) to highest(nrows)
      ! If matrix is expected to be spd,
      !    a. check that there are no negative eigenvalues lower than -ABS(zero_eps)
      !    b. issue a warning if the removal of the specified minimum number of eigenvalues
      !          still leaves small eigenvalues
      IF (spd_expected) THEN
         ! a.
         IF (eigenvalues(1) .LT. -ABS(zero_eps)) THEN
            CPWARN("There are negative eigenvalues")
         ENDIF
         ! b.
         IF ( PRESENT(min_eigenvalues_to_nullify) ) THEN
            IF ( eigenvalues(my_min_eigenvalues_to_nullify+1) .LE. ABS(zero_eps)) THEN
               CPWARN("Min number of eigenvalues does not remove all small eigenvalues")
            ENDIF
         ENDIF
      ENDIF

      DO ieigen = 1, nrows
         IF ( eigenvalues(ieigen) .LE. ABS(zero_eps) .OR. &
              ieigen .LE. my_min_eigenvalues_to_nullify ) THEN
            eigenvalues(ieigen) = 0.0_dp
         ELSE
            eigenvalues(ieigen) = 1.0_dp / eigenvalues(ieigen)
         ENDIF
      END DO

      CALL dbcsr_create(temp_mat_diag, &
              template=matrix, matrix_type = dbcsr_type_no_symmetry)
      CALL dbcsr_create(temp_mat, &
              template=matrix, matrix_type = dbcsr_type_no_symmetry)

      CALL dbcsr_set(temp_mat_diag, 0.0_dp)
      CALL dbcsr_add_on_diag(temp_mat_diag, 1.0_dp)
      CALL dbcsr_set_diag(temp_mat_diag, eigenvalues)

      ! temp = diagonal.U^T
      CALL dbcsr_multiply("N", "T", 1.0_dp, &
              temp_mat_diag, &
              eigenvectors, &
              0.0_dp, temp_mat) !, &
              !filter_eps=eps_filter)
      ! matrix = U.temp
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
              eigenvectors, &
              temp_mat, &
              0.0_dp, matrix) !, &
              !filter_eps=eps_filter)

      CALL dbcsr_release(temp_mat_diag)
      CALL dbcsr_release(temp_mat)
      DEALLOCATE(eigenvalues)
      CALL dbcsr_release(eigenvectors)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_pseudoinverse

! *****************************************************************************
!> \brief Computes the adjugate of a square DBCSR matrix
!> \param matrix - input: initial matrix, output: the adjugate
!> \param det - determinant
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_adjugate(matrix, det, para_env, blacs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(INOUT)                       :: det
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_adjugate'

      INTEGER                                            :: handle, &
                                                            nfullcols_total, &
                                                            nfullrows_total
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_matrix

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrix, fm_struct, name="fm_matrix")
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix, fm_matrix)

      ! Do the work - compute adjugate
      CALL cp_fm_adjugate(fm_matrix, det)

      CALL copy_fm_to_dbcsr(fm_matrix, matrix)

      CALL cp_fm_release(fm_matrix)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_adjugate

! *****************************************************************************
!> \brief Computes the determinant of a square matrix
!> \param matrix ...
!> \param det - determinant
!> \par History
!>       2022.01 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_determinant_generalized(matrix, det, para_env, blacs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(INOUT)                       :: det
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_determinant_generalized'

      INTEGER                                            :: handle, &
                                                            nfullcols_total, &
                                                            nfullrows_total
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_matrix

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrix, fm_struct, name="fm_matrix")
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix, fm_matrix)

      CALL cp_fm_lu_decompose(fm_matrix, det, correct_sign = .TRUE.)

      CALL cp_fm_release(fm_matrix)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_determinant_generalized

END MODULE cp_dbcsr_singular

