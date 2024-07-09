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
   USE population_analyses,             ONLY: lowdin_population_analysis,&
                                              mulliken_population_analysis
   USE tidft_qs,                        ONLY: tidft_dm_to_qs_env
   USE tidft_post_scf,                  ONLY: tidft_postprocess
   USE almo_scf_methods,                ONLY: almo_scf_p_blk_to_t_blk,&
                                              almo_scf_t_rescaling,&
                                              almo_scf_t_to_proj,&
                                              distribute_domains,&
                                              orthogonalize_mos,&
                                              get_overlap
   USE physcon,                         ONLY: evolt
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
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE tidft_types,                     ONLY: tidft_env_type
   USE tidft_optimizer_methods,         ONLY: tidft_variational_optimization
   USE tidft_methods,                   ONLY: stitch_columns, &
                                              split_columns, &
                                              determinant_local, &
                                              matrix_invert_local, &
                                              gradVAR_scale_split_increment, &
                                              tidft_mainvar_orthonormalize, &
                                              dbcsr_set_element
   USE tidft_env_methods,               ONLY: tidft_env_deallocate,&
                                              tidft_init
!                                             tidft_env_deallocate_remaining
   USE tidft_qs,                        ONLY: tidft_dm_to_tidft_ks,&
                                              matrix_layout
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Khaliullin2013,&
                                              Kolafa2004,&
                                              Scheiber2018,&
                                              Staub2019,&
                                              cite_reference
   USE qs_operators_ao,                 ONLY: rRc_xyz_ao
   USE moments_utils,                   ONLY: get_reference_point
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_retain
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm, &
                                              dbcsr_allocate_matrix_set, &
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_retain
   USE cp_para_types,                   ONLY: cp_para_env_type
! CINDY: avoid using dbcsr_local
!   USE dbcsr_local,                     ONLY: dbcsr_init_random
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_binary_read, dbcsr_checksum, dbcsr_copy, dbcsr_create, &
        dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_info, dbcsr_get_stored_coordinates, dbcsr_init_p, dbcsr_dot, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_nblkcols_total, &
        dbcsr_nblkrows_total, dbcsr_p_type, dbcsr_release, dbcsr_reserve_block2d, dbcsr_scale, &
        dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_work_create, &
        dbcsr_print, dbcsr_get_diag, dbcsr_set_diag, dbcsr_transposed, dbcsr_binary_write, dbcsr_init_random
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
   USE qs_initial_guess,                ONLY: calculate_mopac_dm !, calculate_atomic_block_dm
   USE qs_atomic_block,                 ONLY: calculate_atomic_block_dm                                      
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              scf_env_release
   USE string_utilities,                ONLY: integer_to_string
   USE util,                            ONLY: sort
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

      ! create the initial guess for the MOs of all states
      CALL tidft_mainvar_initial_guess(qs_env, tidft_env)

      ! optimize MOs of the ground and excited states
      CALL tidft_variational_optimization(qs_env, tidft_env)

      !!! perform post-scf processing
      CALL tidft_postprocess(qs_env, tidft_env)

      CALL tidft_env_deallocate(tidft_env)

      ! CINDY: deallocate memory for selective state if optimize them state by state
!      IF (tidft_env%state_by_state) THEN
!         CALL tidft_env_deallocate_remaining(tidft_env)
!      ELSE
!         CALL tidft_env_deallocate(tidft_env)
!      ENDIF
      CALL timestop(handle)

   END SUBROUTINE tidft_entry

! **************************************************************************************************
!> \brief create the initial guess for main variables of all states
!> \param qs_env ...
!> \param tidft_env ...
!> \par History
!>       2021.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_mainvar_initial_guess(qs_env, tidft_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_mainvar_initial_guess', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: handle, iaspc, ispin, istore, naspc, &
                                                            nspins, unit_nr, nmos_stateI, &
                                                            stateI, ivar, nstates, idet, &
                                                            ndets, slaterI, stateJ, &
                                                            nmos_stateJ
!                                                           naos, alpha_pos, beta_pos
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
!                                                           overlap_element
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao 
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_type)                                   :: tempNO1, mat_proj_var1, temp_var
      TYPE(dbcsr_type)                                   :: temp2, mat_var23, mat_var23_sigma
      INTEGER                                            :: group, nsize
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diagonal
      TYPE(qs_energy_type), POINTER                      :: qs_energy

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
                      rho=rho)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      nspins = tidft_env%nspins
      nstates = tidft_env%nstates
      ndets=1
      IF (tidft_env%triplet_correction) ndets = 2

      ! create an initial guess
      SELECT CASE (tidft_env%initial_guess_type)
      CASE (atomic_guess)

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

         ! Create matrices of independent variables
         DO stateI = 1, nstates

            CALL dbcsr_create(mat_proj_var1, template=tidft_env%mat_s)

            DO ivar = 1, 3

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

               ! unrestricted case: no. of electron pairs is zero
               IF ( nmos_stateI .GT. 0 ) THEN

                  ! Generate random matrix
                  CALL dbcsr_init_random(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                          mini_seed=stateI+ivar*2+tidft_env%initial_guess_random_seed)

                  IF (stateI .EQ. 1) THEN ! ground state is treated differently

                     !CALL dbcsr_print(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))

                     ! Converts guess DM into matrices of independent variables
                     ! RZK-WARNING: It is important to put unpaired orbitals into ivar=2
                     ! This procedure does not guarantee this at all
                     ! As a result, it will be impossible to generate
                     ! properly unpaired guess for even-electron
                     ! open-shell systems
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

                     ! the percent of the random component in the guess
                     random_percentile = tidft_env%initial_guess_random_fraction
                     ! Obtain ES by slightly perturbing generated GS orbitals
                     CALL dbcsr_add(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), &
                        tidft_env%e_states(1)%slater_det(1)%mat_var(ivar), &
                        random_percentile, 1.0_dp - random_percentile)
                  ENDIF ! MOs

               ENDIF ! non-zero columns

            ENDDO ! end ivar
            CALL dbcsr_release(mat_proj_var1)

         ENDDO  ! end stateI

      CASE (restart_guess)
      ! Read from file

         project_name = logger%iter_info%project_name

         DO stateI = 1, nstates

            DO ivar = 1, 3

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

               IF ( nmos_stateI .GT. 0 ) THEN

                  !WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_", ispin, "_RESTART.dm"
                  WRITE (file_name, '(A,I0,A,I0,A)') TRIM(project_name)//"_TIDFT_STATE_", stateI, "_VAR_", ivar, "_RESTART.mo"
                  CALL dbcsr_create(temp_var, template=tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))
                  !CALL dbcsr_get_info(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), distribution=dist)
                  CALL dbcsr_get_info(temp_var, distribution=dist)
                  CALL dbcsr_binary_read(file_name, distribution=dist, &
                          matrix_new=tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))
                  cs_pos = dbcsr_checksum(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), pos=.TRUE.)
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A,E20.8)') "Read restart MOs "//TRIM(file_name)//" with checksum: ", cs_pos
                  ENDIF
                  CALL dbcsr_release(temp_var)

               ENDIF ! non-zero columns

            ENDDO ! end ivar

         ENDDO  ! end stateI

      CASE DEFAULT

         CPABORT("Illegal initial guess")

      END SELECT

      ! RZK-WARNING: uncomment MO write subroutine if appropriate conditions
      !CALL tidft_mainvar_write(tidft_env)
      CALL tidft_mainvar_orthonormalize(tidft_env)

      CALL timestop(handle)

   END SUBROUTINE tidft_mainvar_initial_guess

! **************************************************************************************************
!> \brief write main variables (MOs essentially) to disk
!> \param tidft_env ...
!> \par History
!>       2022.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_mainvar_write(tidft_env)
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_mainvar_write', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: handle, &
                                                            unit_nr, nmos_stateI, &
                                                            stateI, ivar, nstates
      REAL(KIND=dp)                                      :: cs_pos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist
      LOGICAL                                            :: write_mos

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nstates = tidft_env%nstates

      write_mos = .TRUE.
      IF (write_mos) THEN

         project_name = logger%iter_info%project_name

         DO stateI = 1, nstates

            DO ivar = 1, 3

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

               IF ( nmos_stateI .GT. 0 ) THEN

                  WRITE (file_name, '(A,I0,A,I0,A)') TRIM(project_name)//"_TIDFT_STATE_", stateI, "_VAR_", ivar, "_RESTART.mo"
                  cs_pos = dbcsr_checksum(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), pos=.TRUE.)
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A,E20.8)') "Writing restart MOs "//TRIM(file_name)//" with checksum: ", cs_pos
                  ENDIF
                  CALL dbcsr_binary_write(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar), file_name)

               ENDIF ! non-zero columns

            ENDDO ! end ivar

         ENDDO  ! end stateI

      ENDIF ! write mos

      CALL timestop(handle)

   END SUBROUTINE tidft_mainvar_write

END MODULE tidft

