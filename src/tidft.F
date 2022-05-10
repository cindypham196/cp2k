!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Basic routines for time-independent DFT
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Cindy Pham, Rustam Z. Khaliullin
! **************************************************************************************************
MODULE tidft
   USE almo_scf_methods,                ONLY: almo_scf_p_blk_to_t_blk,&
                                              almo_scf_t_rescaling,&
                                              almo_scf_t_to_proj,&
                                              distribute_domains,&
                                              orthogonalize_mos,&
                                              get_overlap
   USE almo_scf_optimizer,              ONLY: almo_scf_block_diagonal,&
                                              almo_scf_construct_nlmos,&
                                              almo_scf_xalmo_eigensolver,&
                                              almo_scf_xalmo_pcg,&
                                              almo_scf_xalmo_trustr
   USE almo_scf_qs,                     ONLY: almo_dm_to_almo_ks,&
                                              almo_scf_construct_quencher,&
                                              calculate_w_matrix_almo,&
                                              construct_qs_mos,&
                                              init_almo_ks_matrix_via_qs,&
                                              matrix_almo_create,&
                                              matrix_qs_to_almo
   USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                              almo_mat_dim_occ,&
                                              almo_mat_dim_virt,&
                                              almo_mat_dim_virt_disc,&
                                              almo_mat_dim_virt_full,&
                                              optimizer_options_type,&
                                              print_optimizer_options
   USE tidft_types,                     ONLY: tidft_env_type
   USE tidft_methods,                   ONLY: stitch_columns, &
                                              determinant_generalized
   USE tidft_env_methods,               ONLY: tidft_env_deallocate,&
                                              tidft_init
   USE tidft_qs,                        ONLY: tidft_dm_to_tidft_ks
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Khaliullin2013,&
                                              Kolafa2004,&
                                              Scheiber2018,&
                                              Staub2019,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_retain
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_retain
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_local,                     ONLY: dbcsr_init_random
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_binary_read, dbcsr_checksum, dbcsr_copy, dbcsr_create, &
        dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_info, dbcsr_get_stored_coordinates, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_nblkcols_total, &
        dbcsr_nblkrows_total, dbcsr_p_type, dbcsr_release, dbcsr_reserve_block2d, dbcsr_scale, &
        dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_work_create, &
        dbcsr_print, dbcsr_get_diag
   USE domain_submatrix_methods,        ONLY: init_submatrices,&
                                              release_submatrices
   USE input_constants,                 ONLY: &
        almo_deloc_none, almo_deloc_scf, almo_deloc_x, almo_deloc_x_then_scf, &
        almo_deloc_xalmo_1diag, almo_deloc_xalmo_scf, almo_deloc_xalmo_x, almo_deloc_xk, &
        almo_domain_layout_molecular, almo_mat_distr_atomic, almo_mat_distr_molecular, &
        almo_scf_diag, almo_scf_dm_sign, almo_scf_pcg, almo_scf_skip, almo_scf_trustr, &
        atomic_guess, molecular_guess, optimizer_diis, optimizer_lin_eq_pcg, optimizer_pcg, &
        optimizer_trustr, restart_guess, smear_fermi_dirac, virt_full, virt_number, virt_occ_size, &
        xalmo_case_block_diag, xalmo_case_fully_deloc, xalmo_case_normal, xalmo_trial_r0_out
   USE input_section_types,             ONLY: section_vals_type
   USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                              determinant,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathlib,                         ONLY: binomial
   USE message_passing,                 ONLY: mp_sum
   USE molecule_types,                  ONLY: get_molecule_set_info,&
                                              molecule_type
   USE mscfg_types,                     ONLY: get_matrix_from_submatrices,&
                                              molecular_scf_guess_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_initial_guess,                ONLY: calculate_atomic_block_dm,&
                                              calculate_mopac_dm
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              scf_env_release
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft'

   PUBLIC :: tidft_entry

   LOGICAL, PARAMETER :: debug_mode = .FALSE.
   LOGICAL, PARAMETER :: safe_mode = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief The entry point into TIDFT
!> \param qs_env   pointer to the QS environment
!> \param calc_forces   calculate forces?
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Cindy Pham, Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_entry(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_entry', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(tidft_env_type), POINTER                      :: tidft_env

      CALL timeset(routineN, handle)

      ! get a pointer to the environment
      CALL get_qs_env(qs_env, tidft_env=tidft_env)

      tidft_env%calc_forces = calc_forces
      CALL tidft_init(qs_env, tidft_env)

      !!! create the initial guess for the MOs of all states
      CALL tidft_initial_guess(qs_env, tidft_env)

      !!! perform SCF for block diagonal ALMOs
      !!CALL almo_scf_main(qs_env, almo_scf_env)

      !!! allow electron delocalization
      !!CALL almo_scf_delocalization(qs_env, almo_scf_env)

      !!! construct NLMOs
      !!CALL construct_nlmos(qs_env, almo_scf_env)

      !!! do post scf processing
      !!CALL almo_scf_post(qs_env, almo_scf_env)

      CALL tidft_env_deallocate(tidft_env)

      CALL timestop(handle)

   END SUBROUTINE tidft_entry

! **************************************************************************************************
!> \brief create the initial guess for MOs of all states
!> \param qs_env ...
!> \param tidft_env ...
!> \par History
!>       2021.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_initial_guess(qs_env, tidft_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_initial_guess', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: handle, iaspc, ispin, istore, naspc, &
                                                            nspins, unit_nr, nmos_stateI, &
                                                            stateI, ivar, nstates, idet, &
                                                            ndets, slaterI, stateJ, &
                                                            nmos_stateJ
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
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao 
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_type)                                   :: tempNO1
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: fake_e_per_block
      INTEGER                                            :: group, nsize
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diagonal

      CALL timeset(routineN, handle)

      NULLIFY (rho, rho_ao)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! get basic quantities from the qs_env
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      has_unit_metric=has_unit_metric, &
                      para_env=para_env, &
                      nelectron_spin=nelectron_spin, &
                      !mscfg_env=mscfg_env, &
                      rho=rho)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      nspins = tidft_env%nspins
      nstates = tidft_env%nstates
      ndets=1
      IF (tidft_env%triplet_correction) ndets = 2

      ! create an initial guess
      !SELECT CASE (tidft_env%tidft_guess)
      !CASE (atomic_guess, ionic_guess)

         ! use isolated atoms to generate approximate DM for the GS
         IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical .OR. &
             dft_control%qs_control%xtb) THEN
            CALL calculate_mopac_dm(rho_ao, &
                                    matrix_s(1)%matrix, has_unit_metric, &
                                    dft_control, particle_set, atomic_kind_set, qs_kind_set, &
                                    nspins, nelectron_spin, &
                                    para_env)
         ELSE
            CALL calculate_atomic_block_dm(rho_ao, &
                                           matrix_s(1)%matrix, particle_set, atomic_kind_set, qs_kind_set, &
                                           nspins, nelectron_spin, unit_nr, para_env)
         ENDIF

         ! matrices of independent variables
         ! unrestricted case: no. of electron pairs is zero
         DO stateI = 1, nstates

            DO ivar = 1, 3

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

               IF ( nmos_stateI .GT. 0 ) THEN

                  ! Generate random matrix
                  CALL dbcsr_init_random(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                          mini_seed=stateI+ivar*2)

                  IF (stateI .EQ. 1) THEN ! ground state is treated differently

                     !CALL dbcsr_print(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))

                     CALL dbcsr_create(tempNO1, &
                        template=tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))

                     ! Obtain GS orbitals by projecting random MOs onto the generated GS DM
                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   tidft_env%mat_s, &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                                   0.0_dp, tempNO1,filter_eps=tidft_env%eps_filter)

                     IF (ivar .EQ. 3) THEN
                        ispin = 2 !RZK-WARNING: do we always have the second DM when we need it?
                     ELSE
                        ispin = 1
                     ENDIF
                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   rho_ao(ispin)%matrix, &
                                   tempNO1, &
                                   0.0_dp, tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                                   filter_eps=tidft_env%eps_filter)

                     CALL dbcsr_release(tempNO1)

                  ELSE ! excited states
                     ! RZK-WARNING: introduce a new user-defined INPUT-FILE variable that specifies 
                     ! the percient of the random component in the guess
                     random_percentile = 0.05_dp
                     ! Obtain ES by slightly perturbing generated GS orbitals
                     CALL dbcsr_add(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                        tidft_env%e_states(1)%slater_det(1)%mat_var(ivar), &
                        random_percentile, 1.0_dp - random_percentile)
                  ENDIF

               ENDIF ! non-zero columns

            ENDDO ! end ivar

         ENDDO  ! end stateI

!         ! CINDY's code starts here
!         ! Combine X1-X3 to generate MOs mat_t(alpha) and mat_t(beta) for ROKS and UKS
!         DO stateI = 1, nstates
!            DO slaterI = ndets, 1, -1 ! do T before S to avoid overwrite in the X1 U X2 step
!               DO ispin = 1, nspins
!                  ! Do we have paired e-? X1
!                  IF (tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(1) .NE. 0) THEN
!                     WRITE(*,*) "State, determinant, sipn: ", stateI, slaterI, ispin
!                     IF (tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(ispin+1) .NE. 0) THEN
!                        ! combine paired and unpaired electrons
!                        IF (slaterI .EQ. 2)  THEN ! Triplet
!                           IF (ispin .EQ. 1) THEN ! T_alpha
!                                 ! Use the singlet matrix to store intermediate X1 U X2 matrix
!                                 CALL stitch_columns(tidft_env%e_states(stateI)%slater_det(1)%mat_var(1), & ! X1 U X2
!                                                      tidft_env%e_states(stateI)%slater_det(1)%mat_var(2), &
!                                                      tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin))
!                                 CALL stitch_columns(tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin), & ! X1 U X2 U X3
!                                                      tidft_env%e_states(stateI)%slater_det(1)%mat_var(3), &
!                                                      tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin))
!                           ELSE ! T_beta = X1, copy mat_var(1) into mat_t
!                                 CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin), &
!                                                 tidft_env%e_states(stateI)%slater_det(1)%mat_var(1))
!                           ENDIF
!                        ELSE ! Singlet state, here can do for both beta and alpha
!                              CALL stitch_columns(tidft_env%e_states(stateI)%slater_det(1)%mat_var(1), & ! X1 U X2 and X1 U X3
!                                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(ispin+1), &
!                                                   tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin))
!                        ENDIF
!                     ELSE !copy mat_var(1) into mat_t
!                           CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin), &
!                                           tidft_env%e_states(stateI)%slater_det(1)%mat_var(1))
!                     ENDIF
!
!                  ELSE ! No paired e-. Do we have unpaired e-? n_alpha and n_beta are different
!                     IF (tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(ispin+1) .EQ. 0) THEN
!                        IF (ispin .EQ. 1) CPABORT("Number of alpha electrons cannot be zero!")
!                     ELSE ! No paired electrons, only unpaired. This also takes care of UKS
!
!
!                        !RZK-WARNING: Should deal with the triplet state here!
!                        IF (slaterI .EQ. 2)  THEN ! Triplet
!                           IF (ispin .EQ. 1) THEN ! T_alpha
!                           ! Use the singlet matrix to store intermediate X1 U X2 matrix
!                              CALL stitch_columns(tidft_env%e_states(stateI)%slater_det(1)%mat_var(2), & ! X1 U X2
!                                                tidft_env%e_states(stateI)%slater_det(1)%mat_var(3), &
!                                                tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin))
!                           ELSE ! no T_beta for the triplet state, do nothing
!                           ENDIF
!                        ELSE ! Singlet state, here can do for both beta and alpha
!                           CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin), &
!                               tidft_env%e_states(stateI)%slater_det(1)%mat_var(ispin+1))
!                        ENDIF
!
!                        
!                     ENDIF
!                  ENDIF
!               ENDDO ! end ispin
!            ENDDO ! end slaterI
!         ENDDO  ! end stateI
         
         ! CINDY: generate orbitals from mat_var for GS and excited states
         ! deal with cases when some mat_var do not exist

         !! RZK: templates for looping over the key matrices
         !! orbital matrices
         !DO idet = 1, ndets
         !   DO ispin = 1, nspins
         !      DO stateI = 1, nstates
         !         nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)
         !         IF ( nmos_stateI .GT. 0) THEN
         !            DO stateJ = stateI, nstates
         !               nmos_stateJ = tidft_env%e_states(stateJ)%slater_det(idet)%nelectrons_ab(ispin)
         !               IF ( nmos_stateJ .GT. 0) THEN
         !               ENDIF ! zero MOs for state J
         !            ENDDO ! stateJ
         !         ENDIF ! zero MOs
         !      ENDDO ! stateI
         !   ENDDO ! ispin
         !ENDDO ! idet

         !! matrices of independent variables
         !DO stateI = 1, nstates
         !   DO ivar = 1, 3
         !      nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)
         !      IF ( nmos_stateI .GT. 0) THEN
         !         tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar)
         !      ENDIF
         !   ENDDO ! ivar
         !ENDDO ! stateI

         ! This RZK code converts independent variables into MO matrices
         ! Should work for unrestricted case
         ! Reverse order of determinants allows to use empty singlet determinants
         ! to temporarily store the triplet determinant
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) "NStates", nstates
            WRITE (unit_nr, *) "NSpins", nspins
            WRITE (unit_nr, *) "NDets", ndets
         ENDIF

         DO idet = ndets, 1, -1
            DO stateI = 1, nstates
               DO ispin = 1, nspins

                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, *) "idet,istate,ispin", idet, stateI, ispin
                     WRITE (unit_nr, *) "pairs", tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(1)
                     WRITE (unit_nr, *) "single es", tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(ispin+1)
                     WRITE (unit_nr, *) "total", tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)
                  ENDIF

                  IF (tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(1) .EQ. 0) THEN
                     IF (tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(ispin+1) .EQ. 0) THEN
                        ! there are no electrons
                        IF (ispin .EQ. 1) THEN
                           CPABORT("No alpha electrons.")
                        !ELSE
                        !   CPWARN("No beta electrons. MO matrix is not created.")
                        ENDIF
                     ELSE
                        ! there are only unpaired electrons
                        IF (tidft_env%triplet_correction .AND. idet.EQ.2) THEN
                           ! combine alpha and beta into alpha
                           CPASSERT(ispin.EQ.1) ! for triplets, cannot have 0 pairs and non-zero unpaired beta electrons
                           CALL stitch_columns( &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(2), &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(3), &
                                   tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin) &
                                )
                        ELSE
                           CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(ispin+1))
                        ENDIF
                     ENDIF
                  ELSE
                     IF (tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(ispin+1) .EQ. 0) THEN
                        ! there are only paired electrons
                        IF (tidft_env%triplet_correction .AND. idet.EQ.2) THEN
                           CPASSERT(ispin.EQ.2) ! for triplets, cannot have zero unpaired alpha electrons
                           ! for beta spin simply copy pairs
                        ENDIF
                        CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                tidft_env%e_states(stateI)%slater_det(1)%mat_var(1))
                     ELSE
                        IF (tidft_env%triplet_correction .AND. idet.EQ.2) THEN
                           ! for beta spin do nothing - it is already taken care of
                           CPASSERT(ispin.EQ.1) ! for triplets, cannot have non-zero unpaired beta electrons
                           ! for alpha spin combine: pairs + alpha + beta
                           ! use mat_t of singlet as a temporary storage
                           CALL stitch_columns( &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(1), &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(ispin+1), &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin) &
                                )
                           CALL stitch_columns( &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin), &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(ispin+2), &
                                   tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin) &
                                )
                        ELSE
                           ! combine paired and unpared
                           CALL stitch_columns( &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(1), &
                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(ispin+1), &
                                   tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin) &
                                )
                        ENDIF ! triplet correction
                     ENDIF ! zero unpaired?
                  ENDIF ! zero paired?
               ENDDO ! ispin
            ENDDO ! stateI
         ENDDO ! idet

      !CASE (restart_guess)
      ! Read from file

      !DEFAULT

      !   CPABORT("Illegal initial guess")

      !END SELECT

      ! Use MO matrices to compute Kohn-Sham matrices and the energy
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A38)') "Initial guess energy"
      ENDIF

      ALLOCATE(fake_e_per_block(1))
      fake_e_per_block(:) = -1

      obj_energy = 0.0_dp
      DO stateI = 1, nstates

         DO idet = 1, ndets

            DO ispin = 1, nspins

               ! Create projectors only if there are electrons of this spin
               IF (tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin).GT.0) THEN

                  !RZK-warning : orthogonalizing orbitals will destroy their connection to mat_var
                  !CALL orthogonalize_mos(ket=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                  !                       overlap=tidft_env%mat_sigma(stateI,stateI,idet,ispin), &
                  !                       metric=tidft_env%mat_s, &
                  !                       retain_locality=.FALSE., &
                  !                       only_normalize=.FALSE., &
                  !                       nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                  !                       eps_filter=tidft_env%eps_filter, &
                  !                       order_lanczos=tidft_env%order_lanczos, &
                  !                       eps_lanczos=tidft_env%eps_lanczos, &
                  !                       max_iter_lanczos=tidft_env%max_iter_lanczos, &
                  !                       smear=tidft_env%smear)
   
                  ! compute projectors from electronic states
                  CALL almo_scf_t_to_proj(t=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                          p=tidft_env%e_states(stateI)%slater_det(idet)%mat_proj(ispin), &
                                          eps_filter=tidft_env%eps_filter, &
                                          orthog_orbs=.FALSE., &
                                          ! use unallocated array, it will be properly ignored inside
                                          nocc_of_domain=fake_e_per_block, &
                                          s=tidft_env%mat_s, &
                                          ! sigma and sigma_inv will be overwritten
                                          sigma=tidft_env%mat_sigma(stateI,stateI,idet,ispin), &
                                          sigma_inv=tidft_env%mat_sigma_inv(stateI,idet,ispin), &
                                          use_guess=.FALSE., &
                                          smear=tidft_env%smear, &
                                          algorithm=tidft_env%sigma_inv_algorithm, &
                                          !eps_lanczos=tidft_env%eps_lanczos, &
                                          !max_iter_lanczos=tidft_env%max_iter_lanczos, &
                                          !inv_eps_factor=tidft_env%matrix_iter_eps_error_factor, &
                                          para_env=tidft_env%para_env, &
                                          blacs_env=tidft_env%blacs_env)

               ENDIF

            ENDDO ! ispin

            ! compute dm from the projector(s)
            IF (nspins == 1) THEN
               CALL dbcsr_scale(tidft_env%e_states(stateI)%slater_det(idet)%mat_proj(1), 2.0_dp)
            ENDIF

            CALL tidft_dm_to_tidft_ks(qs_env, &
                                      tidft_env%e_states(stateI)%slater_det(idet)%mat_proj, &
                                      tidft_env%e_states(stateI)%slater_det(idet)%mat_ks, &
                                      energy, &
                                      tidft_env%eps_filter, &
                                      smear=tidft_env%smear, &
                                      kTS_sum=kTS_sum)

            ! re-compute projector from the dm
            IF (nspins == 1) THEN
               CALL dbcsr_scale(tidft_env%e_states(stateI)%slater_det(idet)%mat_proj(1), 0.5_dp)
            ENDIF

            IF (unit_nr > 0 .AND. ndets.EQ.2) THEN
               WRITE (unit_nr, '(T2,A28,I10,F40.10)') "Determinant", idet, energy
            ENDIF

            IF (idet.EQ.1) THEN
               energy_singlet = energy
               energy_triplet = energy
            ENDIF
            IF (idet.EQ.2) energy_triplet = energy

         ENDDO ! idet

         state_energy = 2*energy_singlet-energy_triplet
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A28,I10,F40.10)') "State ", stateI, state_energy
         ENDIF

         obj_energy = obj_energy + state_energy

      ENDDO ! stateI

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A38,F40.10)') "Sum of all states", obj_energy
         WRITE (unit_nr, '()')
      ENDIF

      !! compute nonorthogonality penalty
      tidft_env%mat_phi(:,:,:) = 0.0_dp
      intrastate_nonorthogonality_penalty = 0.0_dp
      DO stateI = 1, nstates
         DO idet = 1, ndets
            IF (tidft_env%triplet_correction .AND. idet.EQ.1) cycle ! only triplet determinant is used
            DO ispin = 1, nspins
               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)
               IF ( nmos_stateI .GT. 0) THEN
                  DO stateJ = stateI, nstates
                     nmos_stateJ = tidft_env%e_states(stateJ)%slater_det(idet)%nelectrons_ab(ispin)
                     IF ( nmos_stateJ .GT. 0) THEN
                        ! compute sigma if necessary: (stateI,stateI) is already computed with the energy
                        IF (stateI.NE.stateJ) THEN
                           CALL get_overlap(bra=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                            ket=tidft_env%e_states(stateJ)%slater_det(idet)%mat_t(ispin), &
                                            overlap=tidft_env%mat_sigma(stateI,stateJ,idet,ispin), &
                                            metric=tidft_env%mat_s, &
                                            retain_overlap_sparsity=.FALSE., &
                                            eps_filter=tidft_env%eps_filter)
                           !CALL dbcsr_print(tidft_env%mat_sigma(stateI,stateJ,idet,ispin))
                        ENDIF
                        ! compute the determinant and save it to the Phi matrix
                        CALL determinant_generalized( &
                                matrix = tidft_env%mat_sigma(stateI,stateJ,idet,ispin), &
                                det = determinant_spd, & !tidft_env%mat_phi(stateI,stateJ,ispin), &
                                para_env = tidft_env%para_env, &
                                blacs_env = tidft_env%blacs_env)
                        tidft_env%mat_phi(stateI,stateJ,ispin) = 1.0_dp
                        tidft_env%mat_phi(stateJ,stateI,ispin)=tidft_env%mat_phi(stateI,stateJ,ispin)
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, *) "StateI, StateJ, Slater, Spin: ", stateI, stateJ, idet, ispin, &
                                              "Determinant: ", tidft_env%mat_phi(stateI,stateJ,ispin)
                        ENDIF

                        IF (stateI.EQ.stateJ) THEN

                           IF (debug_mode) THEN
                              CALL determinant(tidft_env%mat_sigma(stateI,stateJ,idet,ispin), &
                                      determinant_spd, &
                                      tidft_env%eps_filter)
                              IF (unit_nr > 0) THEN
                                 WRITE (unit_nr, *) "StateI, StateJ, Slater, Spin: ", stateI, stateJ, idet, ispin, &
                                                    "SPD determi: ", determinant_spd
                              ENDIF
                           ENDIF

                           ! compute the product of the diagonal elements
                           CALL dbcsr_get_info(tidft_env%mat_sigma(stateI,stateJ,idet,ispin), &
                                   nfullrows_total=nsize, group=group)
                           ALLOCATE (diagonal(nsize))
                           CALL dbcsr_get_diag(tidft_env%mat_sigma(stateI,stateJ,idet,ispin), diagonal)
                           CALL mp_sum(diagonal, group)
                           DetDiagonalSigma = PRODUCT(diagonal)
                           DEALLOCATE(diagonal)

                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "StateI, StateJ, Slater, Spin: ", stateI, stateJ, idet, ispin, &
                                                 "Det. diagon: ", DetDiagonalSigma
                              WRITE (unit_nr, *) "StateI, StateJ, Slater, Spin: ", stateI, stateJ, idet, ispin, &
                                                 "Det. s*sd-1: ", &
                                                 tidft_env%mat_phi(stateI,stateJ,ispin) / DetDiagonalSigma
                           ENDIF

                           penalty_stateI = - tidft_env%intrastate_nonorthogonality_penalty_coeff * &
                              ( LOG(tidft_env%mat_phi(stateI,stateI,ispin)) - LOG(DetDiagonalSigma) )

                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, '(T2,A28,I10,F40.10)') "State ", stateI, penalty_stateI
                           ENDIF

                           intrastate_nonorthogonality_penalty = intrastate_nonorthogonality_penalty + &
                              penalty_stateI

                        ENDIF

                     ENDIF ! zero MOs for state J
                  ENDDO ! stateJ
               ENDIF ! zero MOs for state I
            ENDDO ! ispin
         ENDDO ! idet
      ENDDO ! stateI

      ! compute interstate nonorthogonality penalty
      interstate_nonorthogonality_penalty = 0.0_dp

      DEALLOCATE(fake_e_per_block)

      CALL timestop(handle)

   END SUBROUTINE tidft_initial_guess

!tiidft! **************************************************************************************************
!tiidft!> \brief store a history of matrices for later use in almo_scf_initial_guess
!tiidft!> \param almo_scf_env ...
!tiidft!> \par History
!tiidft!>       2016.11 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE almo_scf_store_extrapolation_data(almo_scf_env)
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_store_extrapolation_data', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: handle, ispin, istore, unit_nr
!tiidft      LOGICAL :: delocalization_uses_extrapolation
!tiidft      TYPE(cp_logger_type), POINTER                      :: logger
!tiidft      TYPE(dbcsr_type)                                   :: matrix_no_tmp1, matrix_no_tmp2, &
!tiidft                                                            matrix_no_tmp3, matrix_no_tmp4
!tiidft
!tiidft      CALL timeset(routineN, handle)
!tiidft
!tiidft      ! get a useful output_unit
!tiidft      logger => cp_get_default_logger()
!tiidft      IF (logger%para_env%ionode) THEN
!tiidft         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
!tiidft      ELSE
!tiidft         unit_nr = -1
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%almo_history%nstore > 0) THEN
!tiidft
!tiidft         almo_scf_env%almo_history%istore = almo_scf_env%almo_history%istore + 1
!tiidft
!tiidft         DO ispin = 1, SIZE(almo_scf_env%matrix_t_blk)
!tiidft
!tiidft            istore = MOD(almo_scf_env%almo_history%istore - 1, almo_scf_env%almo_history%nstore) + 1
!tiidft
!tiidft            IF (almo_scf_env%almo_history%istore == 1) THEN
!tiidft               CALL dbcsr_create(almo_scf_env%almo_history%matrix_t(ispin), &
!tiidft                                 template=almo_scf_env%matrix_t_blk(ispin), &
!tiidft                                 matrix_type=dbcsr_type_no_symmetry)
!tiidft            ENDIF
!tiidft            CALL dbcsr_copy(almo_scf_env%almo_history%matrix_t(ispin), &
!tiidft                            almo_scf_env%matrix_t_blk(ispin))
!tiidft
!tiidft            IF (almo_scf_env%almo_history%istore <= almo_scf_env%almo_history%nstore) THEN
!tiidft               CALL dbcsr_create(almo_scf_env%almo_history%matrix_p_up_down(ispin, istore), &
!tiidft                                 template=almo_scf_env%matrix_s(1), &
!tiidft                                 matrix_type=dbcsr_type_no_symmetry)
!tiidft            ENDIF
!tiidft
!tiidft            CALL dbcsr_create(matrix_no_tmp1, template=almo_scf_env%matrix_t_blk(ispin), &
!tiidft                              matrix_type=dbcsr_type_no_symmetry)
!tiidft            CALL dbcsr_create(matrix_no_tmp2, template=almo_scf_env%matrix_t_blk(ispin), &
!tiidft                              matrix_type=dbcsr_type_no_symmetry)
!tiidft
!tiidft            ! compute contra-covariant density matrix
!tiidft            CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1), &
!tiidft                                almo_scf_env%matrix_t_blk(ispin), &
!tiidft                                0.0_dp, matrix_no_tmp1, &
!tiidft                                filter_eps=almo_scf_env%eps_filter)
!tiidft            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_no_tmp1, &
!tiidft                                almo_scf_env%matrix_sigma_inv_0deloc(ispin), &
!tiidft                                0.0_dp, matrix_no_tmp2, &
!tiidft                                filter_eps=almo_scf_env%eps_filter)
!tiidft            CALL dbcsr_multiply("N", "T", 1.0_dp, &
!tiidft                                almo_scf_env%matrix_t_blk(ispin), &
!tiidft                                matrix_no_tmp2, &
!tiidft                                0.0_dp, almo_scf_env%almo_history%matrix_p_up_down(ispin, istore), &
!tiidft                                filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft            CALL dbcsr_release(matrix_no_tmp1)
!tiidft            CALL dbcsr_release(matrix_no_tmp2)
!tiidft
!tiidft         ENDDO
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft      ! exrapolate xalmos?
!tiidft      delocalization_uses_extrapolation = &
!tiidft         .NOT. ((almo_scf_env%deloc_method .EQ. almo_deloc_none) .OR. &
!tiidft                (almo_scf_env%deloc_method .EQ. almo_deloc_xalmo_1diag))
!tiidft      IF (almo_scf_env%xalmo_history%nstore > 0 .AND. &
!tiidft          delocalization_uses_extrapolation) THEN
!tiidft
!tiidft         almo_scf_env%xalmo_history%istore = almo_scf_env%xalmo_history%istore + 1
!tiidft
!tiidft         DO ispin = 1, SIZE(almo_scf_env%matrix_t)
!tiidft
!tiidft            istore = MOD(almo_scf_env%xalmo_history%istore - 1, almo_scf_env%xalmo_history%nstore) + 1
!tiidft
!tiidft            IF (almo_scf_env%xalmo_history%istore == 1) THEN
!tiidft               CALL dbcsr_create(almo_scf_env%xalmo_history%matrix_t(ispin), &
!tiidft                                 template=almo_scf_env%matrix_t(ispin), &
!tiidft                                 matrix_type=dbcsr_type_no_symmetry)
!tiidft            ENDIF
!tiidft            CALL dbcsr_copy(almo_scf_env%xalmo_history%matrix_t(ispin), &
!tiidft                            almo_scf_env%matrix_t(ispin))
!tiidft
!tiidft            IF (almo_scf_env%xalmo_history%istore <= almo_scf_env%xalmo_history%nstore) THEN
!tiidft               !CALL dbcsr_init(almo_scf_env%xalmo_history%matrix_x(ispin, istore))
!tiidft               !CALL dbcsr_create(almo_scf_env%xalmo_history%matrix_x(ispin, istore), &
!tiidft               !        template=almo_scf_env%matrix_t(ispin), &
!tiidft               !        matrix_type=dbcsr_type_no_symmetry)
!tiidft               CALL dbcsr_create(almo_scf_env%xalmo_history%matrix_p_up_down(ispin, istore), &
!tiidft                                 template=almo_scf_env%matrix_s(1), &
!tiidft                                 matrix_type=dbcsr_type_no_symmetry)
!tiidft            ENDIF
!tiidft
!tiidft            CALL dbcsr_create(matrix_no_tmp3, template=almo_scf_env%matrix_t(ispin), &
!tiidft                              matrix_type=dbcsr_type_no_symmetry)
!tiidft            CALL dbcsr_create(matrix_no_tmp4, template=almo_scf_env%matrix_t(ispin), &
!tiidft                              matrix_type=dbcsr_type_no_symmetry)
!tiidft
!tiidft            ! compute contra-covariant density matrix
!tiidft            CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1), &
!tiidft                                almo_scf_env%matrix_t(ispin), &
!tiidft                                0.0_dp, matrix_no_tmp3, &
!tiidft                                filter_eps=almo_scf_env%eps_filter)
!tiidft            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_no_tmp3, &
!tiidft                                almo_scf_env%matrix_sigma_inv(ispin), &
!tiidft                                0.0_dp, matrix_no_tmp4, &
!tiidft                                filter_eps=almo_scf_env%eps_filter)
!tiidft            CALL dbcsr_multiply("N", "T", 1.0_dp, &
!tiidft                                almo_scf_env%matrix_t(ispin), &
!tiidft                                matrix_no_tmp4, &
!tiidft                                0.0_dp, almo_scf_env%xalmo_history%matrix_p_up_down(ispin, istore), &
!tiidft                                filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft            ! store the difference between t and t0
!tiidft            !CALL dbcsr_copy(almo_scf_env%xalmo_history%matrix_x(ispin, istore),&
!tiidft            !        almo_scf_env%matrix_t(ispin))
!tiidft            !CALL dbcsr_add(almo_scf_env%xalmo_history%matrix_x(ispin, istore),&
!tiidft            !        almo_scf_env%matrix_t_blk(ispin),1.0_dp,-1.0_dp)
!tiidft
!tiidft            CALL dbcsr_release(matrix_no_tmp3)
!tiidft            CALL dbcsr_release(matrix_no_tmp4)
!tiidft
!tiidft         ENDDO
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft      CALL timestop(handle)
!tiidft
!tiidft   END SUBROUTINE almo_scf_store_extrapolation_data
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief Initializes the ALMO SCF copy of the AO overlap matrix
!tiidft!>        and all necessary functions (sqrt, inverse...)
!tiidft!> \param matrix_s ...
!tiidft!> \param almo_scf_env ...
!tiidft!> \par History
!tiidft!>       2011.06 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE almo_scf_init_ao_overlap(matrix_s, almo_scf_env)
!tiidft      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_init_ao_overlap', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: handle, unit_nr
!tiidft      TYPE(cp_logger_type), POINTER                      :: logger
!tiidft
!tiidft      CALL timeset(routineN, handle)
!tiidft
!tiidft      ! get a useful output_unit
!tiidft      logger => cp_get_default_logger()
!tiidft      IF (logger%para_env%ionode) THEN
!tiidft         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
!tiidft      ELSE
!tiidft         unit_nr = -1
!tiidft      ENDIF
!tiidft
!tiidft      ! make almo copy of S
!tiidft      ! also copy S to S_blk (i.e. to S with the domain structure imposed)
!tiidft      IF (almo_scf_env%orthogonal_basis) THEN
!tiidft         CALL dbcsr_set(almo_scf_env%matrix_s(1), 0.0_dp)
!tiidft         CALL dbcsr_add_on_diag(almo_scf_env%matrix_s(1), 1.0_dp)
!tiidft         CALL dbcsr_set(almo_scf_env%matrix_s_blk(1), 0.0_dp)
!tiidft         CALL dbcsr_add_on_diag(almo_scf_env%matrix_s_blk(1), 1.0_dp)
!tiidft      ELSE
!tiidft         CALL matrix_qs_to_almo(matrix_s, almo_scf_env%matrix_s(1), &
!tiidft                                almo_scf_env%mat_distr_aos, .FALSE.)
!tiidft         CALL dbcsr_copy(almo_scf_env%matrix_s_blk(1), &
!tiidft                         almo_scf_env%matrix_s(1), keep_sparsity=.TRUE.)
!tiidft      ENDIF
!tiidft
!tiidft      CALL dbcsr_filter(almo_scf_env%matrix_s(1), almo_scf_env%eps_filter)
!tiidft      CALL dbcsr_filter(almo_scf_env%matrix_s_blk(1), almo_scf_env%eps_filter)
!tiidft
!tiidft      IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_diag) THEN
!tiidft         CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_s_blk_sqrt(1), &
!tiidft                                        almo_scf_env%matrix_s_blk_sqrt_inv(1), &
!tiidft                                        almo_scf_env%matrix_s_blk(1), &
!tiidft                                        threshold=almo_scf_env%eps_filter, &
!tiidft                                        order=almo_scf_env%order_lanczos, &
!tiidft                                        !order=0, &
!tiidft                                        eps_lanczos=almo_scf_env%eps_lanczos, &
!tiidft                                        max_iter_lanczos=almo_scf_env%max_iter_lanczos)
!tiidft      ELSE IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_dm_sign) THEN
!tiidft         CALL invert_Hotelling(almo_scf_env%matrix_s_blk_inv(1), &
!tiidft                               almo_scf_env%matrix_s_blk(1), &
!tiidft                               threshold=almo_scf_env%eps_filter, &
!tiidft                               filter_eps=almo_scf_env%eps_filter)
!tiidft      ENDIF
!tiidft
!tiidft      CALL timestop(handle)
!tiidft
!tiidft   END SUBROUTINE almo_scf_init_ao_overlap
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief Selects the subroutine for the optimization of block-daigonal ALMOs.
!tiidft!>        Keep it short and clean.
!tiidft!> \param qs_env ...
!tiidft!> \param almo_scf_env ...
!tiidft!> \par History
!tiidft!>       2011.11 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE almo_scf_main(qs_env, almo_scf_env)
!tiidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_main', routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: handle, ispin, unit_nr
!tiidft      TYPE(cp_logger_type), POINTER                      :: logger
!tiidft
!tiidft      CALL timeset(routineN, handle)
!tiidft
!tiidft      ! get a useful output_unit
!tiidft      logger => cp_get_default_logger()
!tiidft      IF (logger%para_env%ionode) THEN
!tiidft         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
!tiidft      ELSE
!tiidft         unit_nr = -1
!tiidft      ENDIF
!tiidft
!tiidft      SELECT CASE (almo_scf_env%almo_update_algorithm)
!tiidft      CASE (almo_scf_pcg)
!tiidft
!tiidft         ! ALMO PCG optimizer as a special case of XALMO PCG
!tiidft         CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
!tiidft                                 almo_scf_env=almo_scf_env, &
!tiidft                                 optimizer=almo_scf_env%opt_block_diag_pcg, &
!tiidft                                 quench_t=almo_scf_env%quench_t_blk, &
!tiidft                                 matrix_t_in=almo_scf_env%matrix_t_blk, &
!tiidft                                 matrix_t_out=almo_scf_env%matrix_t_blk, &
!tiidft                                 assume_t0_q0x=.FALSE., &
!tiidft                                 perturbation_only=.FALSE., &
!tiidft                                 special_case=xalmo_case_block_diag)
!tiidft
!tiidft         DO ispin = 1, almo_scf_env%nspins
!tiidft            CALL orthogonalize_mos(ket=almo_scf_env%matrix_t_blk(ispin), &
!tiidft                                   overlap=almo_scf_env%matrix_sigma_blk(ispin), &
!tiidft                                   metric=almo_scf_env%matrix_s_blk(1), &
!tiidft                                   retain_locality=.TRUE., &
!tiidft                                   only_normalize=.FALSE., &
!tiidft                                   nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
!tiidft                                   eps_filter=almo_scf_env%eps_filter, &
!tiidft                                   order_lanczos=almo_scf_env%order_lanczos, &
!tiidft                                   eps_lanczos=almo_scf_env%eps_lanczos, &
!tiidft                                   max_iter_lanczos=almo_scf_env%max_iter_lanczos)
!tiidft         ENDDO
!tiidft         !CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env)
!tiidft
!tiidft      CASE (almo_scf_diag)
!tiidft
!tiidft         ! mixing/DIIS optimizer
!tiidft         CALL almo_scf_block_diagonal(qs_env, almo_scf_env, &
!tiidft                                      almo_scf_env%opt_block_diag_diis)
!tiidft
!tiidft      CASE (almo_scf_skip)
!tiidft
!tiidft         DO ispin = 1, almo_scf_env%nspins
!tiidft            CALL orthogonalize_mos(ket=almo_scf_env%matrix_t_blk(ispin), &
!tiidft                                   overlap=almo_scf_env%matrix_sigma_blk(ispin), &
!tiidft                                   metric=almo_scf_env%matrix_s_blk(1), &
!tiidft                                   retain_locality=.TRUE., &
!tiidft                                   only_normalize=.FALSE., &
!tiidft                                   nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
!tiidft                                   eps_filter=almo_scf_env%eps_filter, &
!tiidft                                   order_lanczos=almo_scf_env%order_lanczos, &
!tiidft                                   eps_lanczos=almo_scf_env%eps_lanczos, &
!tiidft                                   max_iter_lanczos=almo_scf_env%max_iter_lanczos)
!tiidft         ENDDO
!tiidft         !CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env)
!tiidft
!tiidft      END SELECT
!tiidft
!tiidft      ! we might need a copy of the converged KS and sigma_inv
!tiidft      DO ispin = 1, almo_scf_env%nspins
!tiidft         CALL dbcsr_copy(almo_scf_env%matrix_ks_0deloc(ispin), &
!tiidft                         almo_scf_env%matrix_ks(ispin))
!tiidft         CALL dbcsr_copy(almo_scf_env%matrix_sigma_inv_0deloc(ispin), &
!tiidft                         almo_scf_env%matrix_sigma_inv(ispin))
!tiidft      ENDDO
!tiidft
!tiidft      CALL timestop(handle)
!tiidft
!tiidft   END SUBROUTINE almo_scf_main
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief selects various post scf routines
!tiidft!> \param qs_env ...
!tiidft!> \param almo_scf_env ...
!tiidft!> \par History
!tiidft!>       2011.06 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE almo_scf_delocalization(qs_env, almo_scf_env)
!tiidft
!tiidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_delocalization', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: col, handle, hold, iblock_col, &
!tiidft                                                            iblock_row, ispin, mynode, &
!tiidft                                                            nblkcols_tot, nblkrows_tot, row, &
!tiidft                                                            unit_nr
!tiidft      LOGICAL                                            :: almo_experimental, tr
!tiidft      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_new_block
!tiidft      TYPE(cp_logger_type), POINTER                      :: logger
!tiidft      TYPE(dbcsr_distribution_type)                      :: dist
!tiidft      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: no_quench
!tiidft      TYPE(optimizer_options_type)                       :: arbitrary_optimizer
!tiidft
!tiidft      CALL timeset(routineN, handle)
!tiidft
!tiidft      ! get a useful output_unit
!tiidft      logger => cp_get_default_logger()
!tiidft      IF (logger%para_env%ionode) THEN
!tiidft         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
!tiidft      ELSE
!tiidft         unit_nr = -1
!tiidft      ENDIF
!tiidft
!tiidft      ! create a local optimizer that handles XALMO DIIS
!tiidft      ! the options of this optimizer are arbitrary because
!tiidft      ! XALMO DIIS SCF does not converge for yet unknown reasons and
!tiidft      ! currently used in the code to get perturbative estimates only
!tiidft      arbitrary_optimizer%optimizer_type = optimizer_diis
!tiidft      arbitrary_optimizer%max_iter = 3
!tiidft      arbitrary_optimizer%eps_error = 1.0E-6_dp
!tiidft      arbitrary_optimizer%ndiis = 2
!tiidft
!tiidft      SELECT CASE (almo_scf_env%deloc_method)
!tiidft      CASE (almo_deloc_x, almo_deloc_scf, almo_deloc_x_then_scf)
!tiidft
!tiidft         ! RZK-warning hack into the quenched routine:
!tiidft         ! create a quench matrix with all-all-all blocks 1.0
!tiidft         ! it is a waste of memory but since matrices are distributed
!tiidft         ! we can tolerate it for now
!tiidft         ALLOCATE (no_quench(almo_scf_env%nspins))
!tiidft         CALL dbcsr_create(no_quench(1), &
!tiidft                           template=almo_scf_env%matrix_t(1), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_get_info(no_quench(1), distribution=dist)
!tiidft         CALL dbcsr_distribution_get(dist, mynode=mynode)
!tiidft         CALL dbcsr_work_create(no_quench(1), &
!tiidft                                work_mutable=.TRUE.)
!tiidft         nblkrows_tot = dbcsr_nblkrows_total(no_quench(1))
!tiidft         nblkcols_tot = dbcsr_nblkcols_total(no_quench(1))
!tiidft         ! RZK-warning: is it a quadratic-scaling routine?
!tiidft         ! As a matter of fact it is! But this block treats
!tiidft         ! fully delocalized MOs. So it is unavoidable.
!tiidft         ! C'est la vie
!tiidft         DO row = 1, nblkrows_tot
!tiidft            DO col = 1, nblkcols_tot
!tiidft               tr = .FALSE.
!tiidft               iblock_row = row
!tiidft               iblock_col = col
!tiidft               CALL dbcsr_get_stored_coordinates(no_quench(1), &
!tiidft                                                 iblock_row, iblock_col, hold)
!tiidft               IF (hold .EQ. mynode) THEN
!tiidft                  NULLIFY (p_new_block)
!tiidft                  CALL dbcsr_reserve_block2d(no_quench(1), &
!tiidft                                             iblock_row, iblock_col, p_new_block)
!tiidft                  CPASSERT(ASSOCIATED(p_new_block))
!tiidft                  p_new_block(:, :) = 1.0_dp
!tiidft               ENDIF
!tiidft            ENDDO
!tiidft         ENDDO
!tiidft         CALL dbcsr_finalize(no_quench(1))
!tiidft         IF (almo_scf_env%nspins .GT. 1) THEN
!tiidft            DO ispin = 2, almo_scf_env%nspins
!tiidft               CALL dbcsr_create(no_quench(ispin), &
!tiidft                                 template=almo_scf_env%matrix_t(1), &
!tiidft                                 matrix_type=dbcsr_type_no_symmetry)
!tiidft               CALL dbcsr_copy(no_quench(ispin), no_quench(1))
!tiidft            ENDDO
!tiidft         ENDIF
!tiidft
!tiidft      END SELECT
!tiidft
!tiidft      SELECT CASE (almo_scf_env%deloc_method)
!tiidft      CASE (almo_deloc_none, almo_deloc_scf)
!tiidft
!tiidft         DO ispin = 1, almo_scf_env%nspins
!tiidft            CALL dbcsr_copy(almo_scf_env%matrix_t(ispin), &
!tiidft                            almo_scf_env%matrix_t_blk(ispin))
!tiidft         ENDDO
!tiidft
!tiidft      CASE (almo_deloc_x, almo_deloc_xk, almo_deloc_x_then_scf)
!tiidft
!tiidft         !!!! RZK-warning a whole class of delocalization methods
!tiidft         !!!! are commented out at the moment because some of their
!tiidft         !!!! routines have not been thoroughly tested.
!tiidft
!tiidft         !!!! if we have virtuals pre-optimize and truncate them
!tiidft         !!!IF (almo_scf_env%need_virtuals) THEN
!tiidft         !!!   SELECT CASE (almo_scf_env%deloc_truncate_virt)
!tiidft         !!!   CASE (virt_full)
!tiidft         !!!      ! simply copy virtual orbitals from matrix_v_full_blk to matrix_v_blk
!tiidft         !!!      DO ispin=1,almo_scf_env%nspins
!tiidft         !!!         CALL dbcsr_copy(almo_scf_env%matrix_v_blk(ispin),&
!tiidft         !!!                 almo_scf_env%matrix_v_full_blk(ispin))
!tiidft         !!!      ENDDO
!tiidft         !!!   CASE (virt_number,virt_occ_size)
!tiidft         !!!      CALL split_v_blk(almo_scf_env)
!tiidft         !!!      !CALL truncate_subspace_v_blk(qs_env,almo_scf_env)
!tiidft         !!!   CASE DEFAULT
!tiidft         !!!      CPErrorMessage(cp_failure_level,routineP,"illegal method for virtual space truncation")
!tiidft         !!!      CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
!tiidft         !!!   END SELECT
!tiidft         !!!ENDIF
!tiidft         !!!CALL harris_foulkes_correction(qs_env,almo_scf_env)
!tiidft
!tiidft         IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_pcg) THEN
!tiidft
!tiidft            CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
!tiidft                                    almo_scf_env=almo_scf_env, &
!tiidft                                    optimizer=almo_scf_env%opt_xalmo_pcg, &
!tiidft                                    quench_t=no_quench, &
!tiidft                                    matrix_t_in=almo_scf_env%matrix_t_blk, &
!tiidft                                    matrix_t_out=almo_scf_env%matrix_t, &
!tiidft                                    assume_t0_q0x=(almo_scf_env%xalmo_trial_wf .EQ. xalmo_trial_r0_out), &
!tiidft                                    perturbation_only=.TRUE., &
!tiidft                                    special_case=xalmo_case_fully_deloc)
!tiidft
!tiidft         ELSE IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_trustr) THEN
!tiidft
!tiidft            CALL almo_scf_xalmo_trustr(qs_env=qs_env, &
!tiidft                                       almo_scf_env=almo_scf_env, &
!tiidft                                       optimizer=almo_scf_env%opt_xalmo_trustr, &
!tiidft                                       quench_t=no_quench, &
!tiidft                                       matrix_t_in=almo_scf_env%matrix_t_blk, &
!tiidft                                       matrix_t_out=almo_scf_env%matrix_t, &
!tiidft                                       perturbation_only=.TRUE., &
!tiidft                                       special_case=xalmo_case_fully_deloc)
!tiidft
!tiidft         ELSE
!tiidft
!tiidft            CPABORT("Other algorithms do not exist")
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft      CASE (almo_deloc_xalmo_1diag)
!tiidft
!tiidft         IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_diag) THEN
!tiidft
!tiidft            almo_scf_env%perturbative_delocalization = .TRUE.
!tiidft            DO ispin = 1, almo_scf_env%nspins
!tiidft               CALL dbcsr_copy(almo_scf_env%matrix_t(ispin), &
!tiidft                               almo_scf_env%matrix_t_blk(ispin))
!tiidft            ENDDO
!tiidft            CALL almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, &
!tiidft                                            arbitrary_optimizer)
!tiidft
!tiidft         ELSE
!tiidft
!tiidft            CPABORT("Other algorithms do not exist")
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft      CASE (almo_deloc_xalmo_x)
!tiidft
!tiidft         IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_pcg) THEN
!tiidft
!tiidft            CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
!tiidft                                    almo_scf_env=almo_scf_env, &
!tiidft                                    optimizer=almo_scf_env%opt_xalmo_pcg, &
!tiidft                                    quench_t=almo_scf_env%quench_t, &
!tiidft                                    matrix_t_in=almo_scf_env%matrix_t_blk, &
!tiidft                                    matrix_t_out=almo_scf_env%matrix_t, &
!tiidft                                    assume_t0_q0x=(almo_scf_env%xalmo_trial_wf .EQ. xalmo_trial_r0_out), &
!tiidft                                    perturbation_only=.TRUE., &
!tiidft                                    special_case=xalmo_case_normal)
!tiidft
!tiidft         ELSE IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_trustr) THEN
!tiidft
!tiidft            CALL almo_scf_xalmo_trustr(qs_env=qs_env, &
!tiidft                                       almo_scf_env=almo_scf_env, &
!tiidft                                       optimizer=almo_scf_env%opt_xalmo_trustr, &
!tiidft                                       quench_t=almo_scf_env%quench_t, &
!tiidft                                       matrix_t_in=almo_scf_env%matrix_t_blk, &
!tiidft                                       matrix_t_out=almo_scf_env%matrix_t, &
!tiidft                                       perturbation_only=.TRUE., &
!tiidft                                       special_case=xalmo_case_normal)
!tiidft
!tiidft         ELSE
!tiidft
!tiidft            CPABORT("Other algorithms do not exist")
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft      CASE (almo_deloc_xalmo_scf)
!tiidft
!tiidft         IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_diag) THEN
!tiidft
!tiidft            CPABORT("Should not be here: convergence will fail!")
!tiidft
!tiidft            almo_scf_env%perturbative_delocalization = .FALSE.
!tiidft            DO ispin = 1, almo_scf_env%nspins
!tiidft               CALL dbcsr_copy(almo_scf_env%matrix_t(ispin), &
!tiidft                               almo_scf_env%matrix_t_blk(ispin))
!tiidft            ENDDO
!tiidft            CALL almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, &
!tiidft                                            arbitrary_optimizer)
!tiidft
!tiidft         ELSE IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_pcg) THEN
!tiidft
!tiidft            CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
!tiidft                                    almo_scf_env=almo_scf_env, &
!tiidft                                    optimizer=almo_scf_env%opt_xalmo_pcg, &
!tiidft                                    quench_t=almo_scf_env%quench_t, &
!tiidft                                    matrix_t_in=almo_scf_env%matrix_t_blk, &
!tiidft                                    matrix_t_out=almo_scf_env%matrix_t, &
!tiidft                                    assume_t0_q0x=(almo_scf_env%xalmo_trial_wf .EQ. xalmo_trial_r0_out), &
!tiidft                                    perturbation_only=.FALSE., &
!tiidft                                    special_case=xalmo_case_normal)
!tiidft
!tiidft            ! RZK-warning THIS IS A HACK TO GET ORBITAL ENERGIES
!tiidft            almo_experimental = .FALSE.
!tiidft            IF (almo_experimental) THEN
!tiidft               almo_scf_env%perturbative_delocalization = .TRUE.
!tiidft               !DO ispin=1,almo_scf_env%nspins
!tiidft               !   CALL dbcsr_copy(almo_scf_env%matrix_t(ispin),&
!tiidft               !           almo_scf_env%matrix_t_blk(ispin))
!tiidft               !ENDDO
!tiidft               CALL almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, &
!tiidft                                               arbitrary_optimizer)
!tiidft            ENDIF ! experimental
!tiidft
!tiidft         ELSE IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_trustr) THEN
!tiidft
!tiidft            CALL almo_scf_xalmo_trustr(qs_env=qs_env, &
!tiidft                                       almo_scf_env=almo_scf_env, &
!tiidft                                       optimizer=almo_scf_env%opt_xalmo_trustr, &
!tiidft                                       quench_t=almo_scf_env%quench_t, &
!tiidft                                       matrix_t_in=almo_scf_env%matrix_t_blk, &
!tiidft                                       matrix_t_out=almo_scf_env%matrix_t, &
!tiidft                                       perturbation_only=.FALSE., &
!tiidft                                       special_case=xalmo_case_normal)
!tiidft
!tiidft         ELSE
!tiidft
!tiidft            CPABORT("Other algorithms do not exist")
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft      CASE DEFAULT
!tiidft
!tiidft         CPABORT("Illegal delocalization method")
!tiidft
!tiidft      END SELECT
!tiidft
!tiidft      SELECT CASE (almo_scf_env%deloc_method)
!tiidft      CASE (almo_deloc_scf, almo_deloc_x_then_scf)
!tiidft
!tiidft         IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN
!tiidft            CPABORT("full scf is NYI for truncated virtual space")
!tiidft         ENDIF
!tiidft
!tiidft         IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_pcg) THEN
!tiidft
!tiidft            CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
!tiidft                                    almo_scf_env=almo_scf_env, &
!tiidft                                    optimizer=almo_scf_env%opt_xalmo_pcg, &
!tiidft                                    quench_t=no_quench, &
!tiidft                                    matrix_t_in=almo_scf_env%matrix_t, &
!tiidft                                    matrix_t_out=almo_scf_env%matrix_t, &
!tiidft                                    assume_t0_q0x=.FALSE., &
!tiidft                                    perturbation_only=.FALSE., &
!tiidft                                    special_case=xalmo_case_fully_deloc)
!tiidft
!tiidft         ELSE IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_trustr) THEN
!tiidft
!tiidft            CALL almo_scf_xalmo_trustr(qs_env=qs_env, &
!tiidft                                       almo_scf_env=almo_scf_env, &
!tiidft                                       optimizer=almo_scf_env%opt_xalmo_trustr, &
!tiidft                                       quench_t=no_quench, &
!tiidft                                       matrix_t_in=almo_scf_env%matrix_t, &
!tiidft                                       matrix_t_out=almo_scf_env%matrix_t, &
!tiidft                                       perturbation_only=.FALSE., &
!tiidft                                       special_case=xalmo_case_fully_deloc)
!tiidft
!tiidft         ELSE
!tiidft
!tiidft            CPABORT("Other algorithms do not exist")
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft      END SELECT
!tiidft
!tiidft      ! clean up
!tiidft      SELECT CASE (almo_scf_env%deloc_method)
!tiidft      CASE (almo_deloc_x, almo_deloc_scf, almo_deloc_x_then_scf)
!tiidft         DO ispin = 1, almo_scf_env%nspins
!tiidft            CALL dbcsr_release(no_quench(ispin))
!tiidft         ENDDO
!tiidft         DEALLOCATE (no_quench)
!tiidft      END SELECT
!tiidft
!tiidft      CALL timestop(handle)
!tiidft
!tiidft   END SUBROUTINE almo_scf_delocalization
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief orbital localization
!tiidft!> \param qs_env ...
!tiidft!> \param almo_scf_env ...
!tiidft!> \par History
!tiidft!>       2018.09 created [Ziling Luo]
!tiidft!> \author Ziling Luo
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE construct_nlmos(qs_env, almo_scf_env)
!tiidft
!tiidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'construct_nlmos', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: ispin
!tiidft
!tiidft      IF (almo_scf_env%construct_nlmos) THEN
!tiidft
!tiidft         DO ispin = 1, almo_scf_env%nspins
!tiidft
!tiidft            CALL orthogonalize_mos(ket=almo_scf_env%matrix_t(ispin), &
!tiidft                                   overlap=almo_scf_env%matrix_sigma(ispin), &
!tiidft                                   metric=almo_scf_env%matrix_s(1), &
!tiidft                                   retain_locality=.FALSE., &
!tiidft                                   only_normalize=.FALSE., &
!tiidft                                   nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
!tiidft                                   eps_filter=almo_scf_env%eps_filter, &
!tiidft                                   order_lanczos=almo_scf_env%order_lanczos, &
!tiidft                                   eps_lanczos=almo_scf_env%eps_lanczos, &
!tiidft                                   max_iter_lanczos=almo_scf_env%max_iter_lanczos)
!tiidft         ENDDO
!tiidft
!tiidft         CALL construct_nlmos_wrapper(qs_env, almo_scf_env, virtuals=.FALSE.)
!tiidft
!tiidft         IF (almo_scf_env%opt_nlmo_pcg%opt_penalty%virtual_nlmos) THEN
!tiidft            CALL construct_virtuals(almo_scf_env)
!tiidft            CALL construct_nlmos_wrapper(qs_env, almo_scf_env, virtuals=.TRUE.)
!tiidft         ENDIF
!tiidft
!tiidft         IF (almo_scf_env%opt_nlmo_pcg%opt_penalty%compactification_filter_start .GT. 0.0_dp) THEN
!tiidft            CALL nlmo_compactification(qs_env, almo_scf_env, almo_scf_env%matrix_t)
!tiidft         ENDIF
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft   END SUBROUTINE construct_nlmos
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief Calls NLMO optimization
!tiidft!> \param qs_env ...
!tiidft!> \param almo_scf_env ...
!tiidft!> \param virtuals ...
!tiidft!> \par History
!tiidft!>       2019.10 created [Ziling Luo]
!tiidft!> \author Ziling Luo
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE construct_nlmos_wrapper(qs_env, almo_scf_env, virtuals)
!tiidft
!tiidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft      LOGICAL, INTENT(IN)                                :: virtuals
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'construct_nlmos_wrapper', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      REAL(KIND=dp)                                      :: det_diff, prev_determinant
!tiidft
!tiidft      almo_scf_env%overlap_determinant = 1.0
!tiidft      ! KEEP: initial_vol_coeff = almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength
!tiidft      almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength = &
!tiidft         -1.0_dp*almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength !NEW1
!tiidft
!tiidft      ! loop over the strength of the orthogonalization penalty
!tiidft      prev_determinant = 10.0_dp
!tiidft      DO WHILE (almo_scf_env%overlap_determinant .GT. almo_scf_env%opt_nlmo_pcg%opt_penalty%final_determinant)
!tiidft
!tiidft         IF (.NOT. virtuals) THEN
!tiidft            CALL almo_scf_construct_nlmos(qs_env=qs_env, &
!tiidft                                          optimizer=almo_scf_env%opt_nlmo_pcg, &
!tiidft                                          matrix_s=almo_scf_env%matrix_s(1), &
!tiidft                                          matrix_mo_in=almo_scf_env%matrix_t, &
!tiidft                                          matrix_mo_out=almo_scf_env%matrix_t, &
!tiidft                                          template_matrix_sigma=almo_scf_env%matrix_sigma_inv, &
!tiidft                                          overlap_determinant=almo_scf_env%overlap_determinant, &
!tiidft                                          mat_distr_aos=almo_scf_env%mat_distr_aos, &
!tiidft                                          virtuals=virtuals, &
!tiidft                                          eps_filter=almo_scf_env%eps_filter)
!tiidft         ELSE
!tiidft            CALL almo_scf_construct_nlmos(qs_env=qs_env, &
!tiidft                                          optimizer=almo_scf_env%opt_nlmo_pcg, &
!tiidft                                          matrix_s=almo_scf_env%matrix_s(1), &
!tiidft                                          matrix_mo_in=almo_scf_env%matrix_v, &
!tiidft                                          matrix_mo_out=almo_scf_env%matrix_v, &
!tiidft                                          template_matrix_sigma=almo_scf_env%matrix_sigma_vv, &
!tiidft                                          overlap_determinant=almo_scf_env%overlap_determinant, &
!tiidft                                          mat_distr_aos=almo_scf_env%mat_distr_aos, &
!tiidft                                          virtuals=virtuals, &
!tiidft                                          eps_filter=almo_scf_env%eps_filter)
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft         det_diff = prev_determinant - almo_scf_env%overlap_determinant
!tiidft         almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength = &
!tiidft            almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength/ &
!tiidft            ABS(almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength_dec_factor)
!tiidft
!tiidft         IF (det_diff < almo_scf_env%opt_nlmo_pcg%opt_penalty%determinant_tolerance) THEN
!tiidft            EXIT
!tiidft         ENDIF
!tiidft         prev_determinant = almo_scf_env%overlap_determinant
!tiidft
!tiidft      ENDDO
!tiidft
!tiidft   END SUBROUTINE construct_nlmos_wrapper
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief Construct virtual orbitals
!tiidft!> \param almo_scf_env ...
!tiidft!> \par History
!tiidft!>       2019.10 created [Ziling Luo]
!tiidft!> \author Ziling Luo
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE construct_virtuals(almo_scf_env)
!tiidft
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'construct_virtuals', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: ispin, n
!tiidft      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
!tiidft      TYPE(dbcsr_type)                                   :: tempNV1, tempVOcc1, tempVOcc2, tempVV1, &
!tiidft                                                            tempVV2
!tiidft
!tiidft      DO ispin = 1, almo_scf_env%nspins
!tiidft
!tiidft         CALL dbcsr_create(tempNV1, &
!tiidft                           template=almo_scf_env%matrix_v(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(tempVOcc1, &
!tiidft                           template=almo_scf_env%matrix_vo(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(tempVOcc2, &
!tiidft                           template=almo_scf_env%matrix_vo(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(tempVV1, &
!tiidft                           template=almo_scf_env%matrix_sigma_vv(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(tempVV2, &
!tiidft                           template=almo_scf_env%matrix_sigma_vv(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft
!tiidft         ! Generate random virtual matrix
!tiidft         CALL dbcsr_init_random(almo_scf_env%matrix_v(ispin), &
!tiidft                                keep_sparsity=.FALSE.)
!tiidft
!tiidft         ! Project the orbital subspace out
!tiidft         CALL dbcsr_multiply("N", "N", 1.0_dp, &
!tiidft                             almo_scf_env%matrix_s(1), &
!tiidft                             almo_scf_env%matrix_v(ispin), &
!tiidft                             0.0_dp, tempNV1, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL dbcsr_multiply("T", "N", 1.0_dp, &
!tiidft                             tempNV1, &
!tiidft                             almo_scf_env%matrix_t(ispin), &
!tiidft                             0.0_dp, tempVOcc1, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL dbcsr_multiply("N", "N", 1.0_dp, &
!tiidft                             tempVOcc1, &
!tiidft                             almo_scf_env%matrix_sigma_inv(ispin), &
!tiidft                             0.0_dp, tempVOcc2, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL dbcsr_multiply("N", "T", 1.0_dp, &
!tiidft                             almo_scf_env%matrix_t(ispin), &
!tiidft                             tempVOcc2, &
!tiidft                             0.0_dp, tempNV1, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL dbcsr_add(almo_scf_env%matrix_v(ispin), tempNV1, 1.0_dp, -1.0_dp)
!tiidft
!tiidft         ! compute VxV overlap
!tiidft         CALL dbcsr_multiply("N", "N", 1.0_dp, &
!tiidft                             almo_scf_env%matrix_s(1), &
!tiidft                             almo_scf_env%matrix_v(ispin), &
!tiidft                             0.0_dp, tempNV1, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL dbcsr_multiply("T", "N", 1.0_dp, &
!tiidft                             almo_scf_env%matrix_v(ispin), &
!tiidft                             tempNV1, &
!tiidft                             0.0_dp, tempVV1, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL orthogonalize_mos(ket=almo_scf_env%matrix_v(ispin), &
!tiidft                                overlap=tempVV1, &
!tiidft                                metric=almo_scf_env%matrix_s(1), &
!tiidft                                retain_locality=.FALSE., &
!tiidft                                only_normalize=.FALSE., &
!tiidft                                nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
!tiidft                                eps_filter=almo_scf_env%eps_filter, &
!tiidft                                order_lanczos=almo_scf_env%order_lanczos, &
!tiidft                                eps_lanczos=almo_scf_env%eps_lanczos, &
!tiidft                                max_iter_lanczos=almo_scf_env%max_iter_lanczos)
!tiidft
!tiidft         ! compute VxV block of the KS matrix
!tiidft         CALL dbcsr_multiply("N", "N", 1.0_dp, &
!tiidft                             almo_scf_env%matrix_ks(ispin), &
!tiidft                             almo_scf_env%matrix_v(ispin), &
!tiidft                             0.0_dp, tempNV1, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL dbcsr_multiply("T", "N", 1.0_dp, &
!tiidft                             almo_scf_env%matrix_v(ispin), &
!tiidft                             tempNV1, &
!tiidft                             0.0_dp, tempVV1, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL dbcsr_get_info(tempVV1, nfullrows_total=n)
!tiidft         ALLOCATE (eigenvalues(n))
!tiidft         CALL cp_dbcsr_syevd(tempVV1, tempVV2, &
!tiidft                             eigenvalues, &
!tiidft                             para_env=almo_scf_env%para_env, &
!tiidft                             blacs_env=almo_scf_env%blacs_env)
!tiidft         DEALLOCATE (eigenvalues)
!tiidft
!tiidft         CALL dbcsr_multiply("N", "N", 1.0_dp, &
!tiidft                             almo_scf_env%matrix_v(ispin), &
!tiidft                             tempVV2, &
!tiidft                             0.0_dp, tempNV1, &
!tiidft                             filter_eps=almo_scf_env%eps_filter)
!tiidft
!tiidft         CALL dbcsr_copy(almo_scf_env%matrix_v(ispin), tempNV1)
!tiidft
!tiidft         CALL dbcsr_release(tempNV1)
!tiidft         CALL dbcsr_release(tempVOcc1)
!tiidft         CALL dbcsr_release(tempVOcc2)
!tiidft         CALL dbcsr_release(tempVV1)
!tiidft         CALL dbcsr_release(tempVV2)
!tiidft
!tiidft      ENDDO
!tiidft
!tiidft   END SUBROUTINE construct_virtuals
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief Compactify (set small blocks to zero) orbitals
!tiidft!> \param qs_env ...
!tiidft!> \param almo_scf_env ...
!tiidft!> \param matrix ...
!tiidft!> \par History
!tiidft!>       2019.10 created [Ziling Luo]
!tiidft!> \author Ziling Luo
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE nlmo_compactification(qs_env, almo_scf_env, matrix)
!tiidft
!tiidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:), &
!tiidft         INTENT(IN)                                      :: matrix
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'nlmo_compactification', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: iblock_col, iblock_col_size, iblock_row, &
!tiidft                                                            iblock_row_size, icol, irow, ispin, &
!tiidft                                                            Ncols, Nrows, nspins, para_group, &
!tiidft                                                            unit_nr
!tiidft      LOGICAL                                            :: element_by_element
!tiidft      REAL(KIND=dp)                                      :: energy, eps_local, eps_start, &
!tiidft                                                            max_element, spin_factor
!tiidft      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occ, retained
!tiidft      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p
!tiidft      TYPE(cp_logger_type), POINTER                      :: logger
!tiidft      TYPE(dbcsr_iterator_type)                          :: iter
!tiidft      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_p_tmp, matrix_t_tmp
!tiidft
!tiidft      ! define the output_unit
!tiidft      logger => cp_get_default_logger()
!tiidft      IF (logger%para_env%ionode) THEN
!tiidft         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
!tiidft      ELSE
!tiidft         unit_nr = -1
!tiidft      ENDIF
!tiidft
!tiidft      nspins = SIZE(matrix)
!tiidft      element_by_element = .FALSE.
!tiidft
!tiidft      IF (nspins .EQ. 1) THEN
!tiidft         spin_factor = 2.0_dp
!tiidft      ELSE
!tiidft         spin_factor = 1.0_dp
!tiidft      ENDIF
!tiidft
!tiidft      ALLOCATE (matrix_t_tmp(nspins))
!tiidft      ALLOCATE (matrix_p_tmp(nspins))
!tiidft      ALLOCATE (retained(nspins))
!tiidft      ALLOCATE (occ(2))
!tiidft
!tiidft      DO ispin = 1, nspins
!tiidft
!tiidft         ! init temporary storage
!tiidft         CALL dbcsr_create(matrix_t_tmp(ispin), &
!tiidft                           template=matrix(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_copy(matrix_t_tmp(ispin), matrix(ispin))
!tiidft
!tiidft         CALL dbcsr_create(matrix_p_tmp(ispin), &
!tiidft                           template=almo_scf_env%matrix_p(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_copy(matrix_p_tmp(ispin), almo_scf_env%matrix_p(ispin))
!tiidft
!tiidft      ENDDO
!tiidft
!tiidft      IF (unit_nr > 0) THEN
!tiidft         WRITE (unit_nr, *)
!tiidft         WRITE (unit_nr, '(T2,A)') &
!tiidft            "Energy dependence on the (block-by-block) filtering of the NLMO coefficients"
!tiidft         IF (unit_nr > 0) WRITE (unit_nr, '(T2,A13,A20,A20,A25)') &
!tiidft            "EPS filter", "Occupation Alpha", "Occupation Beta", "Energy"
!tiidft      ENDIF
!tiidft
!tiidft      eps_start = almo_scf_env%opt_nlmo_pcg%opt_penalty%compactification_filter_start
!tiidft      eps_local = MAX(eps_start, 10E-14_dp)
!tiidft
!tiidft      DO
!tiidft
!tiidft         IF (eps_local > 0.11_dp) EXIT
!tiidft
!tiidft         DO ispin = 1, nspins
!tiidft
!tiidft            retained(ispin) = 0
!tiidft            CALL dbcsr_work_create(matrix_t_tmp(ispin), work_mutable=.TRUE.)
!tiidft            CALL dbcsr_iterator_start(iter, matrix_t_tmp(ispin))
!tiidft            DO WHILE (dbcsr_iterator_blocks_left(iter))
!tiidft               CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p, &
!tiidft                                              row_size=iblock_row_size, col_size=iblock_col_size)
!tiidft               DO icol = 1, iblock_col_size
!tiidft
!tiidft                  IF (element_by_element) THEN
!tiidft
!tiidft                     DO irow = 1, iblock_row_size
!tiidft                        IF (ABS(data_p(irow, icol)) .LT. eps_local) THEN
!tiidft                           data_p(irow, icol) = 0.0_dp
!tiidft                        ELSE
!tiidft                           retained(ispin) = retained(ispin) + 1
!tiidft                        ENDIF
!tiidft                     ENDDO
!tiidft
!tiidft                  ELSE ! rows are blocked
!tiidft
!tiidft                     max_element = 0.0_dp
!tiidft                     DO irow = 1, iblock_row_size
!tiidft                        IF (ABS(data_p(irow, icol)) .GT. max_element) THEN
!tiidft                           max_element = ABS(data_p(irow, icol))
!tiidft                        ENDIF
!tiidft                     ENDDO
!tiidft                     IF (max_element .LT. eps_local) THEN
!tiidft                        DO irow = 1, iblock_row_size
!tiidft                           data_p(irow, icol) = 0.0_dp
!tiidft                        ENDDO
!tiidft                     ELSE
!tiidft                        retained(ispin) = retained(ispin) + iblock_row_size
!tiidft                     ENDIF
!tiidft
!tiidft                  ENDIF ! block rows?
!tiidft               ENDDO ! icol
!tiidft
!tiidft            ENDDO ! iterator
!tiidft            CALL dbcsr_iterator_stop(iter)
!tiidft            CALL dbcsr_finalize(matrix_t_tmp(ispin))
!tiidft            CALL dbcsr_filter(matrix_t_tmp(ispin), eps_local)
!tiidft
!tiidft            CALL dbcsr_get_info(matrix_t_tmp(ispin), group=para_group, &
!tiidft                                nfullrows_total=Nrows, &
!tiidft                                nfullcols_total=Ncols)
!tiidft            CALL mp_sum(retained(ispin), para_group)
!tiidft
!tiidft            !devide by the total no. elements
!tiidft            occ(ispin) = retained(ispin)/Nrows/Ncols
!tiidft
!tiidft            ! compute the global projectors (for the density matrix)
!tiidft            CALL almo_scf_t_to_proj( &
!tiidft               t=matrix_t_tmp(ispin), &
!tiidft               p=matrix_p_tmp(ispin), &
!tiidft               eps_filter=almo_scf_env%eps_filter, &
!tiidft               orthog_orbs=.FALSE., &
!tiidft               nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
!tiidft               s=almo_scf_env%matrix_s(1), &
!tiidft               sigma=almo_scf_env%matrix_sigma(ispin), &
!tiidft               sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
!tiidft               use_guess=.FALSE., &
!tiidft               algorithm=almo_scf_env%sigma_inv_algorithm, &
!tiidft               inv_eps_factor=almo_scf_env%matrix_iter_eps_error_factor, &
!tiidft               inverse_accelerator=almo_scf_env%order_lanczos, &
!tiidft               eps_lanczos=almo_scf_env%eps_lanczos, &
!tiidft               max_iter_lanczos=almo_scf_env%max_iter_lanczos, &
!tiidft               para_env=almo_scf_env%para_env, &
!tiidft               blacs_env=almo_scf_env%blacs_env)
!tiidft
!tiidft            ! compute dm from the projector(s)
!tiidft            CALL dbcsr_scale(matrix_p_tmp(ispin), spin_factor)
!tiidft
!tiidft         ENDDO
!tiidft
!tiidft         ! the KS matrix is updated outside the spin loop
!tiidft         CALL almo_dm_to_almo_ks(qs_env, &
!tiidft                                 matrix_p_tmp, &
!tiidft                                 almo_scf_env%matrix_ks, &
!tiidft                                 energy, &
!tiidft                                 almo_scf_env%eps_filter, &
!tiidft                                 almo_scf_env%mat_distr_aos)
!tiidft
!tiidft         IF (nspins .LT. 2) occ(2) = occ(1)
!tiidft         IF (unit_nr > 0) WRITE (unit_nr, '(T2,E13.3,F20.10,F20.10,F25.15)') &
!tiidft            eps_local, occ(1), occ(2), energy
!tiidft
!tiidft         eps_local = 2.0_dp*eps_local
!tiidft
!tiidft      ENDDO
!tiidft
!tiidft      DO ispin = 1, nspins
!tiidft
!tiidft         CALL dbcsr_release(matrix_t_tmp(ispin))
!tiidft         CALL dbcsr_release(matrix_p_tmp(ispin))
!tiidft
!tiidft      ENDDO
!tiidft
!tiidft      DEALLOCATE (matrix_t_tmp)
!tiidft      DEALLOCATE (matrix_p_tmp)
!tiidft      DEALLOCATE (occ)
!tiidft      DEALLOCATE (retained)
!tiidft
!tiidft   END SUBROUTINE nlmo_compactification
!tiidft
!tiidft! *****************************************************************************
!tiidft!> \brief after SCF we have the final density and KS matrices compute various
!tiidft!>        post-scf quantities
!tiidft!> \param qs_env ...
!tiidft!> \param almo_scf_env ...
!tiidft!> \par History
!tiidft!>       2015.03 created  [Rustam Z. Khaliullin]
!tiidft!> \author Rustam Z. Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE almo_scf_post(qs_env, almo_scf_env)
!tiidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tiidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_post', routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: handle, ispin
!tiidft      TYPE(cp_fm_type), POINTER                          :: mo_coeff
!tiidft      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_w
!tiidft      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_t_processed
!tiidft      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
!tiidft      TYPE(qs_scf_env_type), POINTER                     :: scf_env
!tiidft
!tiidft      CALL timeset(routineN, handle)
!tiidft
!tiidft      ! store matrices to speed up the next scf run
!tiidft      CALL almo_scf_store_extrapolation_data(almo_scf_env)
!tiidft
!tiidft      ! orthogonalize orbitals before returning them to QS
!tiidft      ALLOCATE (matrix_t_processed(almo_scf_env%nspins))
!tiidft      !ALLOCATE (matrix_v_processed(almo_scf_env%nspins))
!tiidft
!tiidft      DO ispin = 1, almo_scf_env%nspins
!tiidft
!tiidft         CALL dbcsr_create(matrix_t_processed(ispin), &
!tiidft                           template=almo_scf_env%matrix_t(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft
!tiidft         CALL dbcsr_copy(matrix_t_processed(ispin), &
!tiidft                         almo_scf_env%matrix_t(ispin))
!tiidft
!tiidft         !CALL dbcsr_create(matrix_v_processed(ispin), &
!tiidft         !                  template=almo_scf_env%matrix_v(ispin), &
!tiidft         !                  matrix_type=dbcsr_type_no_symmetry)
!tiidft
!tiidft         !CALL dbcsr_copy(matrix_v_processed(ispin), &
!tiidft         !                almo_scf_env%matrix_v(ispin))
!tiidft
!tiidft         IF (almo_scf_env%return_orthogonalized_mos) THEN
!tiidft
!tiidft            CALL orthogonalize_mos(ket=matrix_t_processed(ispin), &
!tiidft                                   overlap=almo_scf_env%matrix_sigma(ispin), &
!tiidft                                   metric=almo_scf_env%matrix_s(1), &
!tiidft                                   retain_locality=.FALSE., &
!tiidft                                   only_normalize=.FALSE., &
!tiidft                                   nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
!tiidft                                   eps_filter=almo_scf_env%eps_filter, &
!tiidft                                   order_lanczos=almo_scf_env%order_lanczos, &
!tiidft                                   eps_lanczos=almo_scf_env%eps_lanczos, &
!tiidft                                   max_iter_lanczos=almo_scf_env%max_iter_lanczos, &
!tiidft                                   smear=almo_scf_env%smear)
!tiidft         ENDIF
!tiidft
!tiidft      ENDDO
!tiidft
!tiidft      !! RS-WARNING: If smearing ALMO is requested, rescaled fully-occupied orbitals are returned to QS
!tiidft      !! RS-WARNING: Beware that QS will not be informed about electronic entropy.
!tiidft      !!             If you want a quick and dirty transfer to QS energy, uncomment the following hack:
!tiidft      !! IF (almo_scf_env%smear) THEN
!tiidft      !!    qs_env%energy%kTS = 0.0_dp
!tiidft      !!    DO ispin = 1, almo_scf_env%nspins
!tiidft      !!       qs_env%energy%kTS = qs_env%energy%kTS + almo_scf_env%kTS(ispin)
!tiidft      !!    END DO
!tiidft      !! END IF
!tiidft
!tiidft      ! return orbitals to QS
!tiidft      NULLIFY (mos, mo_coeff, scf_env)
!tiidft
!tiidft      CALL get_qs_env(qs_env, mos=mos, scf_env=scf_env)
!tiidft
!tiidft      DO ispin = 1, almo_scf_env%nspins
!tiidft         ! Currently only fm version of mo_set is usable.
!tiidft         ! First transform the matrix_t to fm version
!tiidft         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
!tiidft         CALL copy_dbcsr_to_fm(matrix_t_processed(ispin), mo_coeff)
!tiidft         CALL dbcsr_release(matrix_t_processed(ispin))
!tiidft      ENDDO
!tiidft      DO ispin = 1, almo_scf_env%nspins
!tiidft         CALL dbcsr_release(matrix_t_processed(ispin))
!tiidft      ENDDO
!tiidft      DEALLOCATE (matrix_t_processed)
!tiidft
!tiidft      ! calculate post scf properties
!tiidft      ! CALL almo_post_scf_compute_properties(qs_env, almo_scf_env)
!tiidft      CALL almo_post_scf_compute_properties(qs_env)
!tiidft
!tiidft      ! In normal scf procedure this is done by scf,
!tiidft      ! but it's not called by ALMO
!tiidft      ! best to release it somewhere else
!tiidft      IF (ASSOCIATED(scf_env)) THEN
!tiidft         CALL scf_env_release(scf_env)
!tiidft      ENDIF
!tiidft
!tiidft      ! compute the W matrix if associated
!tiidft      IF (almo_scf_env%calc_forces) THEN
!tiidft         CALL get_qs_env(qs_env, matrix_w=matrix_w)
!tiidft         IF (ASSOCIATED(matrix_w)) THEN
!tiidft            CALL calculate_w_matrix_almo(matrix_w, almo_scf_env)
!tiidft         ELSE
!tiidft            CPABORT("Matrix W is needed but not associated")
!tiidft         ENDIF
!tiidft      ENDIF
!tiidft
!tiidft      CALL timestop(handle)
!tiidft
!tiidft   END SUBROUTINE almo_scf_post
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief Do post scf calculations with ALMO
!tiidft!>        WARNING: ALMO post scf calculation may not work for certain quantities,
!tiidft!>        like forces, since ALMO wave function is only 'partially' optimized
!tiidft!> \param qs_env ...
!tiidft!> \par History
!tiidft!>       2016.12 created [Yifei Shi]
!tiidft!> \author Yifei Shi
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE almo_post_scf_compute_properties(qs_env)
!tiidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'almo_post_scf_compute_properties', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      CALL qs_scf_compute_properties(qs_env)
!tiidft
!tiidft   END SUBROUTINE almo_post_scf_compute_properties

END MODULE tidft

