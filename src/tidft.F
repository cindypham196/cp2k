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
      CALL tidft_post_scf(qs_env, tidft_env)

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

! **************************************************************************************************
!> \brief Print the final TIDFT excitation energies and oscillator strengths.
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_post_scf(qs_env, tidft_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tidft_post_scf'

      REAL(kind=dp), DIMENSION(3)                        :: reference_point
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipole_ao_xyz
      INTEGER                                            :: unit_nr
      CHARACTER(len=2)                                   :: lsd_str
      CHARACTER(len=20)                                  :: mult_str
      INTEGER                                            :: handle, ideriv, ispin, &
                                                            istate, nspins, &
                                                            nstates, idet, &
                                                            indexGS, indexES, &
                                                            nderivs, ndets
      REAL(kind=dp)                                      :: osc_strength, &
                                                            scaling_factor, &
                                                            excitationE
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trans_dipoles_from_GS
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: excitation_number
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: index_true, &
                                                            fake_e_per_block
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: energies, energies_temp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: tempNO1, tempNO2
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_proj

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nderivs = 3
      nspins = tidft_env%nspins
      nstates = tidft_env%nstates
      ndets=1
      IF (tidft_env%triplet_correction) ndets = 2

      ! compute components of the dipole operator in the length form
      NULLIFY (dipole_ao_xyz)
      CALL dbcsr_allocate_matrix_set(dipole_ao_xyz, nderivs)

      DO ideriv = 1, nderivs
         CALL dbcsr_init_p(dipole_ao_xyz(ideriv)%matrix)
         CALL dbcsr_copy(dipole_ao_xyz(ideriv)%matrix, tidft_env%mat_s)
      END DO

      !CALL get_reference_point( &
      !        rpoint = reference_point, &
      !        qs_env=qs_env, &
      !        reference=tidft_env%dipole_reference_type, &
      !        ref_point=tidft_env%dipole_reference_point)

      ! Compute the dipole matrix in the AO basis set
      CALL rRc_xyz_ao( &
              op=dipole_ao_xyz, &
              qs_env=qs_env, &
              !rc=reference_point, &
              rc=(/0.0_dp, 0.0_dp, 0.0_dp/), &
              order=1, &
              minimum_image=.FALSE., &
              soft=.FALSE.)

      ALLOCATE(index_true(nstates))
      ALLOCATE(energies(nstates))
      ALLOCATE(energies_temp(nstates))

      ! RZK-WARNING: which slater/state energy to use - theory is not developed
      DO istate = 1, nstates
         energies(istate) = tidft_env%e_states(istate)%state_energy
      ENDDO

      energies_temp(:) = energies(:)
      CALL sort(energies_temp, nstates, index_true)
      indexGS = index_true(1)

      ! transition dipole moment
      ALLOCATE (trans_dipoles_from_GS(nstates, nderivs, nspins))
      trans_dipoles_from_GS(:, :, :) = 0.0_dp

      ! excitation numbers
      ALLOCATE (excitation_number(nstates, 2)) ! second dim is spin
      excitation_number(: , :) = 0.0_dp

      idet = 1
      IF (tidft_env%triplet_correction) THEN
         CPWARN("Lowest-energy open-shell singlet determinant is used as the ground state.")
      ENDIF

      ! RZK-WARNING: If interstate_component is FALSE then ADJ will not be available
      DO ideriv = 1, nderivs
         DO ispin = 1, nspins

            ! Check if there are electrons of this spin
            IF (tidft_env%e_states(indexGS)%slater_det(idet)%nelectrons_ab(ispin) .GT. 0) THEN

               ! tempNO1 = D.T_GS
               CALL dbcsr_create(tempNO1, &
                       template=tidft_env%e_states(indexGS)%slater_det(idet)%mat_t(ispin))
               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                       dipole_ao_xyz(ideriv)%matrix, &
                       tidft_env%e_states(indexGS)%slater_det(idet)%mat_t(ispin), &
                       0.0_dp, tempNO1, &
                       filter_eps=tidft_env%eps_filter)

               DO istate = 2, nstates

                  indexES = index_true(istate)

                  CALL dbcsr_create(tempNO2, &
                       template=tidft_env%e_states(indexES)%slater_det(idet)%mat_t(ispin))

                  ! take the normalization into account
                  scaling_factor = 1.0_dp / SQRT( tidft_env%mat_phi_spin(indexGS,indexGS,ispin) * &
                     tidft_env%mat_phi_spin(indexES,indexES,ispin) )

                  ! tempNO2 = tempNO1.adj(sigma_GS,ES)
                  ! RZK-WARNING: Or perhaps tempNO2 = tempNO1.adj(sigma_JI)
                  IF (indexGS .LE. indexES) THEN
                     CALL dbcsr_multiply("N", "N", scaling_factor, &
                             tempNO1, &
                             tidft_env%mat_adj_sigma(indexGS, indexES, idet, ispin), &
                             0.0_dp, tempNO2, &
                             filter_eps=tidft_env%eps_filter)
                  ELSE
                     ! Note the transpose
                     CALL dbcsr_multiply("N", "T", scaling_factor, &
                             tempNO1, &
                             tidft_env%mat_adj_sigma(indexES, indexGS, idet, ispin), &
                             0.0_dp, tempNO2, &
                             filter_eps=tidft_env%eps_filter)
                  ENDIF

                  ! trans_dipoles_from_GS(ES,ideriv,ispin) = Trace [tempNO2 . T_ES^t]
                  CALL dbcsr_dot( &
                          tempNO2, &
                          tidft_env%e_states(indexES)%slater_det(idet)%mat_t(ispin), &
                          trans_dipoles_from_GS(indexES, ideriv, ispin))
                  CALL dbcsr_release(tempNO2)

               ENDDO ! excited states

               CALL dbcsr_release(tempNO1)

            ENDIF ! non-zero electrons

         ENDDO
      ENDDO

      ! RZK-WARNING: RZK disagrees that spins must be squared before their addition
      IF (nspins == 1) THEN
         !trans_dipoles_from_GS(:, :, 1) = SQRT(2.0_dp) * trans_dipoles_from_GS(:, :, 1)
         trans_dipoles_from_GS(:, :, 1) = 2.0_dp * trans_dipoles_from_GS(:, :, 1)
      ELSE
         !trans_dipoles_from_GS(:, :, 1) = SQRT(trans_dipoles_from_GS(:, :, 1)**2 + trans_dipoles_from_GS(:, :, 2)**2)
         trans_dipoles_from_GS(:, :, 1) = trans_dipoles_from_GS(:, :, 1) + trans_dipoles_from_GS(:, :, 2)
      END IF

      ! Compute the number of excited electrons
      DO ispin = 1, nspins

         ! Check if there are electrons of this spin
         IF (tidft_env%e_states(indexGS)%slater_det(idet)%nelectrons_ab(ispin) .GT. 0) THEN

            DO istate = 2, nstates

               indexES = index_true(istate)

               ! tempNO1 = D.T_GS
               CALL dbcsr_create(tempNO1, &
                       template=tidft_env%mat_sigma(indexGS, indexGS, idet, ispin))
               CALL dbcsr_create(tempNO2, &
                       template=tidft_env%mat_sigma(indexGS, indexGS, idet, ispin))

               ! tempNO1 = sigma_(0I) . siginv_(II)
               IF (indexGS .LE. indexES) THEN
                  CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          tidft_env%mat_sigma_inv(indexGS, idet, ispin), &
                          tidft_env%mat_sigma(indexGS, indexES, idet, ispin), &
                          0.0_dp, tempNO1, &
                          filter_eps=tidft_env%eps_filter)
                  CALL dbcsr_multiply("T", "N", 1.0_dp, &
                          tidft_env%mat_sigma(indexGS, indexES, idet, ispin), &
                          tempNO1, &
                          0.0_dp, tempNO2, &
                          filter_eps=tidft_env%eps_filter)
               ELSE
                  ! Note the transpose
                  CALL dbcsr_multiply("N", "T", 1.0_dp, &
                          tidft_env%mat_sigma_inv(indexGS, idet, ispin), &
                          tidft_env%mat_sigma(indexES, indexGS, idet, ispin), &
                          0.0_dp, tempNO1, &
                          filter_eps=tidft_env%eps_filter)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          tidft_env%mat_sigma(indexES, indexGS, idet, ispin), &
                          tempNO1, &
                          0.0_dp, tempNO2, &
                          filter_eps=tidft_env%eps_filter)
               ENDIF

               ! exc_no(ES,ispin) = Trace [ siginv_II . sigma_I0 . siginv_00 . sigma_0I] =
               !                  = Trace [ Rho_I . Rho_0 ]
               CALL dbcsr_dot( &
                       tidft_env%mat_sigma_inv(indexES, idet, ispin), &
                       tempNO2, &
                       excitation_number(indexES, ispin))

               excitation_number(indexES, ispin) = &
                  tidft_env%e_states(indexGS)%slater_det(idet)%nelectrons_ab(ispin) - &
                  excitation_number(indexES, ispin)

               CALL dbcsr_release(tempNO2)
               CALL dbcsr_release(tempNO1)

            ENDDO ! excited states

         ENDIF ! non-zero electrons

      ENDDO

      IF (nspins == 1) THEN
         excitation_number(:, 2) = excitation_number(:, 1)
      END IF

      ! Print header
      IF (tidft_env%restricted) THEN
         IF (tidft_env%unpairings .EQ. 0) THEN
            IF (tidft_env%desired_multip.EQ.1) THEN
               lsd_str = 'R'
            ELSE
               lsd_str = 'RO'
            ENDIF
         ELSE
            lsd_str = 'RO'
         ENDIF
      ELSE
         lsd_str = 'U'
      ENDIF

      IF (unit_nr > 0) THEN
         CALL integer_to_string(tidft_env%desired_multip, mult_str)

         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
         WRITE (unit_nr, '(1X,A,A,A)') lsd_str, "-TIDFT states of multiplicity ", TRIM(mult_str)
         WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)

         WRITE (unit_nr, '(1X,6X,A5,A5,A20,A7,A27,A9)') &
            "Exct", &
            "SCF", &
            "SCF energy", &
            "Exct", &
            "Transition dipole (a.u.)", &
            "Oscl str"
         WRITE (unit_nr, '(1X,6X,A5,A5,A20,A7,3A9,A9)') &
            "no.", &
            "no.", &
            "(a.u.)", &
            "(eV)", &
            "x", "y", "z", &
            "(a.u.)"
         WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
      ENDIF

      ! Include the GS into the print out
      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(1X,A6,I5,I5,F20.10,F7.2,3ES9.1E2,ES9.2E2)') &
            "TIDFT|", &
            0, &
            indexGS, &
            tidft_env%e_states(indexGS)%state_energy, &
            0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp

      END IF

      ! Print excited-state info state-by-state
      DO istate = 2, nstates
   
         indexES = index_true(istate) ! istate - index of increasing energy, indexES - index if TIDFT arrays

         IF (unit_nr > 0) THEN

            excitationE = energies(indexES) - energies(indexGS)

            osc_strength = 2.0_dp/3.0_dp * excitationE * &
               accurate_dot_product(trans_dipoles_from_GS(indexES, :, 1), trans_dipoles_from_GS(indexES, :, 1))

            WRITE (unit_nr, '(1X,A6,I5,I5,F20.10,F7.2,3ES9.1E2,ES9.2E2)') &
               "TIDFT|", &
               istate-1, &
               indexES, &
               tidft_env%e_states(indexES)%state_energy, &
               excitationE*evolt, &
               trans_dipoles_from_GS(indexES, 1:nderivs, 1), &
               osc_strength

         END IF
   
      ENDDO
   
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
      ENDIF

      ! Second table
      IF (unit_nr > 0) THEN
         CALL integer_to_string(tidft_env%desired_multip, mult_str)

         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(1X,6X,A5,A5,A42)') &
            "Exct", &
            "SCF", &
            "No. of excited electrons"
         WRITE (unit_nr, '(1X,6X,A5,A5,3A14)') &
            "no.", &
            "no.", &
            "Alpha", &
            "Beta", &
            "Total"
         WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
      ENDIF

      ! Include the GS into the print out
      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(1X,A6,I5,I5,3F14.10)') &
            "TIDFT|", &
            0, &
            indexGS, &
            0.0_dp, &
            0.0_dp, &
            0.0_dp

      END IF

      ! Print excited-state info state-by-state
      DO istate = 2, nstates
   
         indexES = index_true(istate) ! istate - index of increasing energy, indexES - index if TIDFT arrays

         IF (unit_nr > 0) THEN

            WRITE (unit_nr, '(1X,A6,I5,I5,3F14.10)') &
               "TIDFT|", &
               istate-1, &
               indexES, &
               excitation_number(indexES, 1), &
               excitation_number(indexES, 2), &
               SUM(excitation_number(indexES,:))

         END IF
   
      ENDDO
   
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
      ENDIF

      ! ==================================
      ! Print population analysis
      ALLOCATE(mat_proj(nspins))
      ALLOCATE(fake_e_per_block(1))
      fake_e_per_block(:) = -1

      !idet = 1

      DO istate = 1, nstates
   
         indexES = index_true(istate) ! istate - index of increasing energy, indexES - index if TIDFT arrays

         DO ispin = 1, nspins
   
            CALL dbcsr_create(mat_proj(ispin), template=tidft_env%mat_s)
   
            ! Create projectors only if there are electrons of this spin
            IF (tidft_env%e_states(indexES)%slater_det(idet)%nelectrons_ab(ispin).GT.0) THEN
   
               ! compute projectors from electronic states
               CALL almo_scf_t_to_proj(t=tidft_env%e_states(indexES)%slater_det(idet)%mat_t(ispin), &
                                       p=mat_proj(ispin), &
                                       eps_filter=tidft_env%eps_filter, &
                                       orthog_orbs=.FALSE., &
                                       nocc_of_domain=fake_e_per_block, &
                                       s=tidft_env%mat_s, &
                                       ! RZK-WARNING: sigma and sigma_inv will be overwritten
                                       sigma=tidft_env%mat_sigma(indexES,indexES,idet,ispin), &
                                       sigma_inv=tidft_env%mat_sigma_inv(indexES,idet,ispin), &
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
            CALL dbcsr_scale(mat_proj(1), 2.0_dp)
         ENDIF
   
         ! copy DM/projector into qs_env
         CALL tidft_dm_to_qs_env(qs_env, mat_proj)
         ! perform population analysis
         CALL mulliken_population_analysis(qs_env, &
                 output_unit = unit_nr, &
                 print_level = 2)

         ! re-compute projector from the dm
         IF (nspins == 1) THEN
            CALL dbcsr_scale(mat_proj(1), 0.5_dp)
         ENDIF

         DO ispin = 1, nspins
            CALL dbcsr_release(mat_proj(ispin))
         ENDDO
   
      ENDDO ! istate
   
      DEALLOCATE(mat_proj)
      DEALLOCATE(fake_e_per_block)
      ! ==================================

      CALL dbcsr_deallocate_matrix_set(dipole_ao_xyz)
      DEALLOCATE(excitation_number)
      DEALLOCATE(trans_dipoles_from_GS)
      DEALLOCATE(energies) 
      DEALLOCATE(energies_temp) 
      DEALLOCATE(index_true)

      CALL timestop(handle)

   END SUBROUTINE tidft_post_scf

   
END MODULE tidft

