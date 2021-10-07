!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Environment for TIDFT methods
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_env_methods

   USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                              almo_mat_dim_occ
   USE tidft_types,                     ONLY: tidft_env_type
   USE tidft_qs,                        ONLY: matrix_layout
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_retain
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_retain
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_complete_redistribute, dbcsr_copy, dbcsr_copy_into_existing, dbcsr_create, &
        dbcsr_desymmetrize, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_get_stored_coordinates, &
        dbcsr_multiply, dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_block2d, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_work_create, &
        dbcsr_type_symmetric
   USE input_constants,                 ONLY: &
        almo_constraint_distance, almo_deloc_none, almo_deloc_xalmo_1diag, &
        almo_domain_layout_atomic, almo_domain_layout_molecular, almo_frz_crystal, &
        almo_mat_distr_molecular, almo_scf_diag, almo_scf_skip, almo_scf_trustr, cg_hager_zhang, &
        do_bondparm_vdw, molecular_guess, tensor_orthogonal, virt_full, virt_minimal, virt_number, &
        xalmo_trial_r0_out, optimizer_bfgs
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_env_methods'

   PUBLIC :: tidft_env_create, tidft_env_deallocate

CONTAINS

! **************************************************************************************************
!> \brief Creation and basic initialization of the tidft environment
!> \param qs_env ...
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_env_create(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nallocate
      TYPE(tidft_env_type), POINTER                      :: tidft_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ALLOCATE (tidft_env)

      ! get basic quantities from the qs_env
      CALL get_qs_env(qs_env, input=input, dft_control=dft_control)

      ! parse the input file and set appropriate quantities
      CALL tidft_read_write_input(input, tidft_env)

!tiidft      ! set up the buffer for the history of matrices
!tiidft      almo_scf_env%almo_history%nstore = almo_scf_env%almo_extrapolation_order
!tiidft      almo_scf_env%almo_history%istore = 0
!tiidft      ! do not allocate zero
!tiidft      nallocate = MAX(1, almo_scf_env%almo_extrapolation_order)
!tiidft      ALLOCATE (almo_scf_env%almo_history%matrix_p_up_down(almo_scf_env%nspins, nallocate))
!tiidft      ALLOCATE (almo_scf_env%almo_history%matrix_t(almo_scf_env%nspins))
!tiidft      almo_scf_env%xalmo_history%nstore = almo_scf_env%xalmo_extrapolation_order
!tiidft      almo_scf_env%xalmo_history%istore = 0
!tiidft      nallocate = MAX(1, almo_scf_env%xalmo_extrapolation_order)
!tiidft      ALLOCATE (almo_scf_env%xalmo_history%matrix_p_up_down(almo_scf_env%nspins, nallocate))
!tiidft      !ALLOCATE (almo_scf_env%xalmo_history%matrix_x(almo_scf_env%nspins, nallocate))
!tiidft      ALLOCATE (almo_scf_env%xalmo_history%matrix_t(almo_scf_env%nspins))

      CALL tidft_init(qs_env, tidft_env)

      ! let qs_env know about the newly allocalted TIDFT environment
      CALL set_qs_env(qs_env, tidft_env=tidft_env)

      CALL timestop(handle)

   END SUBROUTINE tidft_env_create

! **************************************************************************************************
!> \brief Parses the ALMO input section
!> \param input ...
!> \param almo_scf_env ...
!> \par History
!>       2021.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_read_write_input(input, tidft_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_read_write_input', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, n_rep, &
                                                            jj, nmultips, GS_multips
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nstates 
      TYPE(section_vals_type), POINTER :: almo_analysis_section, almo_opt_diis_section, &
         almo_opt_pcg_section, tidft_section, matrix_iterate_section, nlmo_opt_pcg_section, &
         penalty_section, xalmo_opt_newton_pcg_section, xalmo_opt_pcg_section, &
         xalmo_opt_trustr_section

      CALL timeset(routineN, handle)

      tidft_section => section_vals_get_subs_vals(input, "DFT%TIDFT")
!tiidft      almo_opt_diis_section => section_vals_get_subs_vals(almo_scf_section, &
!tiidft                                                          "ALMO_OPTIMIZER_DIIS")

      CALL section_vals_val_get(tidft_section, "EPS_FILTER", &
              r_val=tidft_env%eps_filter)
      CALL section_vals_val_get(tidft_section, "N_EXCITED_STATES", &
              i_val=tidft_env%nstates)
      ! add the ground state
      tidft_env%nstates = tidft_env%nstates + 1
      ! Read in the number of excited states of different multiplicities
      !CALL section_vals_val_get(tidft_section,"N_EXCITED_STATES",i_rep_val=1, i_vals=tmplist)
      !nmultips = SIZE(tmplist)

      CALL section_vals_val_get(tidft_section, "SPIN_PURE", &
              l_val=tidft_env%spin_pure)
      CALL section_vals_val_get(tidft_section, "MAX_UNPAIR", &
              i_val=tidft_env%unpairings)

!tiidft      ! optimizers
!tiidft      CALL section_vals_val_get(almo_opt_diis_section, "EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_block_diag_diis%eps_error)
!tiidft      CALL section_vals_val_get(almo_opt_diis_section, "MAX_ITER", &
!tiidft                                i_val=almo_scf_env%opt_block_diag_diis%max_iter)
!tiidft      CALL section_vals_val_get(almo_opt_diis_section, "EPS_ERROR_EARLY", &
!tiidft                                r_val=almo_scf_env%opt_block_diag_diis%eps_error_early)
!tiidft      CALL section_vals_val_get(almo_opt_diis_section, "MAX_ITER_EARLY", &
!tiidft                                i_val=almo_scf_env%opt_block_diag_diis%max_iter_early)
!tiidft      CALL section_vals_val_get(almo_opt_diis_section, "N_DIIS", &
!tiidft                                i_val=almo_scf_env%opt_block_diag_diis%ndiis)
!tiidft
!tiidft      ! check for conflicts between options
!tiidft      IF (almo_scf_env%xalmo_trial_wf .EQ. xalmo_trial_r0_out .AND. &
!tiidft          almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_trustr) THEN
!tiidft         CPABORT("Trust region algorithm cannot optimize projected XALMOs")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%deloc_method .EQ. almo_deloc_xalmo_1diag .AND. &
!tiidft          almo_scf_env%xalmo_update_algorithm .NE. almo_scf_diag) THEN
!tiidft         CPABORT("1-step delocalization correction requires a different algorithm")
!tiidft      ENDIF

      CALL timestop(handle)

   END SUBROUTINE tidft_read_write_input

! **************************************************************************************************
!> \brief Initialization of tidft_env
!> \param qs_env ...
!> \param tidft_env ...
!> \param calc_forces ...
!> \par History
!>       2020.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_init(qs_env, tidft_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: ao, handle, i, iao, idomain, ispin, &
                                                            multip, naos, natoms, ndomains, nelec, &
                                                            nelec_a, nelec_b, nmols, nspins, &
                                                            unit_nr, nstates, ivar, nmultips, &
                                                            nslater, slaterI, stateI
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ! define the output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! set optimizers' types
      tidft_env%opt_bfgs%optimizer_type = optimizer_bfgs

      ! get info from the qs_env
      CALL get_qs_env(qs_env, &
                      nelectron_total=tidft_env%nelectrons, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control, &
                      input=input, &
                      has_unit_metric=tidft_env%orthogonal_basis, &
                      para_env=tidft_env%para_env, &
                      blacs_env=tidft_env%blacs_env)
      CALL cp_para_env_retain(tidft_env%para_env)
      CALL cp_blacs_env_retain(tidft_env%blacs_env)

      ! copy basic quantities
      IF (dft_control%nspins == 1 .OR. dft_control%nspins == 3) THEN
         tidft_env%restricted = .TRUE.
      ELSE
         tidft_env%restricted = .FALSE.
         IF (tidft_env%spin_pure) THEN
            CPABORT("Impossbile to obtain spin-pure determinants with unrestricted methods")
         ENDIF
      ENDIF

      tidft_env%desired_multip = 
      nstates =  SUM(tidft_env%nstates_multip)

      IF (tidft_env%restricted) THEN
         IF (tidft_env%unpairings .EQ. 0) THEN
            IF (tidft_env%desired_multip.EQ.1) THEN
               tidft_env%nspins = 1
               tidft_env%nvars = 1
            ELSE
               tidft_env%nspins = 2
               tidft_env%nvars = 2
            ENDIF
         ELSE
               tidft_env%nspins = 2
               tidft_env%nvars = 3
         ENDIF
      ELSE
         tidft_env%nspins = 2
         tidft_env%nvars = 2
      ENDIF

      IF (tidft_env%spin_pure) THEN
         IF ( .NOT. ( tidft_env%desired_multip.EQ.1 .AND. &
              tidft_env%unpairings.EQ.1 ) ) THEN
            CPABORT("Spin pure states cannot be enforced")
         ENDIF
      ENDIF

      ALLOCATE(tidft_env%e_states(nstates))
      nslater = 1
      IF (tidft_env%spin_pure) nslater = 2
      DO stateI = 1, nstates
         tidft_env%e_states(stateI)%nslater = nslater
         ALLOCATE(tidft_env%e_states(stateI)%slater_det(nslater))
         tidft_env%e_states(stateI)%slater_det(1)%desired_multip = tidft_env%desired_multip
         IF (tidft_env%spin_pure) tidft_env%e_states(stateI)%slater_det(2)%desired_multip = tidft_env%desired_multip + 2
         ! calculate the number of alpha and beta occupied orbitals from
         ! the number of electrons and multiplicity
         ! Na + Nb = Ne
         ! Na - Nb = Mult - 1 (assume Na > Nb as we do not have more info from get_molecule_set_info)
         ! Na = (Ne + Mult - 1)/2
         ! Nb = (Ne - Mult + 1)/2
         DO slaterI = 1, nslater
            tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(1) = &
               ( tidft_env%nelectrons + tidft_env%e_states(stateI)%slater_det(slaterI)%desired_multip - 1 ) / 2
            tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(2) = &
               ( tidft_env%nelectrons - tidft_env%e_states(stateI)%slater_det(slaterI)%desired_multip + 1 ) / 2
            tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(:)=0
            DO ivar = 1, tidft_env%nvars
               ! 1st entry: doubly occupied orbitals: (ne-(mult-1))/2-unpair
               IF (ivar .EQ. 1) tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(1) = &
                  (tidft_env%nelectrons - tidft_env%e_states(stateI)%slater_det(slaterI)%desired_multip + 1) / 2 - &
                  tidft_env%unpairings
               ! 2nd entry: singly occupied alpha orbitals: mult-1+unpair
               IF (ivar .EQ. 2) tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(2) = &
                  tidft_env%e_states(stateI)%slater_det(slaterI)%desired_multip - 1 + tidft_env%unpairings
               ! 3rd entry: singly occupied beta orbitals: unpair
               IF (ivar .EQ. 3) tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(3) = &
                  tidft_env%unpairings
            ENDDO
         ENDDO
      ENDDO

      CALL tidft_env_create_matrices(tidft_env, matrix_s(1)%matrix)

!tiidft      CALL tidft_init_ao_overlap(matrix_s(1)%matrix, tidft_env)
!tiidft      ! initialization of the KS matrix
!tiidft      CALL init_almo_ks_matrix_via_qs(qs_env, &
!tiidft                                      tidft_env%matrix_ks, &
!tiidft                                      tidft_env%mat_distr_aos, &
!tiidft                                      tidft_env%eps_filter)
!tiidft      CALL construct_qs_mos(qs_env, tidft_env)

!tiidft
!tiidft      ! FINISH setting job parameters here, print out job info
!tiidft      CALL tidft_print_job_info(tidft_env, unit_nr)
!tiidft

      CALL timestop(handle)

   END SUBROUTINE tidft_init

! **************************************************************************************************
!> \brief create various matrices
!> \param tidft_env ...
!> \param matrix_s0 ...
!> \par History
!>       2021.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_env_create_matrices(tidft_env, matrix_s0)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s0

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_create_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, &
                                                            nstates, nvars, &
                                                            nspins, ivar, &
                                                            stateI, &
                                                            nmos_stateI, &
                                                            nmos_stateI_flip, &
                                                            stateJ, &
                                                            nmos_stateJ

      CALL timeset(routineN, handle)

      CALL matrix_layout( &
              matrix_ao_template = matrix_s0, &
              matrix_new = tidft_env%mat_s , &
              name_new = 'm_s', &
              orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
              symmetry_new = dbcsr_type_symmetric)

      nstates = tidft_env%nstates
      nspins = tidft_env%nspins
      nvars = tidft_env%nvars

      ALLOCATE(tidft_env%mat_sigma(nstates,nstates,nspins))
      DO stateI = 1, nstates
         ALLOCATE(tidft_env%e_states(stateI)%slater_det(1)%mat_t(nspins))
         IF (tidft_env%spin_pure) ALLOCATE(tidft_env%e_states(stateI)%slater_det(2)%mat_t(nspins))
         ALLOCATE(tidft_env%e_states(stateI)%slater_det(1)%mat_var(nvars))
         ! It is not necessary to allocate mat_var for the other determinant(s)
      ENDDO

      DO ispin = 1, nspins

         DO stateI = 1, nstates

            nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%nelectrons_ab(ispin)

            CALL matrix_layout( &
                    matrix_ao_template = matrix_s0, &
                    matrix_new = tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin) , &
                    name_new = 'm_T_Ia', &
                    orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                    symmetry_new = dbcsr_type_no_symmetry, &
                    nmos = (/-1, nmos_stateI/) )

            IF (tidft_env%spin_pure) THEN

               nmos_stateI_flip = tidft_env%e_states(stateI)%slater_det(2)%nelectrons_ab(ispin)

               CALL matrix_layout( &
                       matrix_ao_template = matrix_s0, &
                       matrix_new = tidft_env%e_states(stateI)%slater_det(2)%mat_t(ispin) , &
                       name_new = 'm_Tflip_Ia', &
                       orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                       symmetry_new = dbcsr_type_no_symmetry, &
                       nmos = (/-1, nmos_stateI_flip/) )

            ENDIF

            DO stateJ = stateI, nstates

               nmos_stateJ = tidft_env%e_states(stateJ)%slater_det(1)%nelectrons_ab(ispin)

               CALL matrix_layout( &
                       matrix_ao_template = matrix_s0, &
                       matrix_new = tidft_env%mat_sigma(stateI, stateJ, ispin) , &
                       name_new = 'm_sig_IJab', &
                       orbital_type_new = (/almo_mat_dim_occ, almo_mat_dim_occ/), &
                       symmetry_new = dbcsr_type_no_symmetry, &
                       nmos = (/nmos_stateI, nmos_stateJ/) )

            ENDDO ! stateJ

         ENDDO ! stateI

      ENDDO ! ispin

      DO ivar = 1, nvars

         DO stateI = 1, nstates

            nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

            CALL matrix_layout( &
                    matrix_ao_template = matrix_s0, &
                    matrix_new = tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar) , &
                    name_new = 'm_X_Ia', &
                    orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                    symmetry_new = dbcsr_type_no_symmetry, &
                    nmos = (/-1, nmos_stateI/) )

         ENDDO ! stateI

      ENDDO ! ivar

      CALL timestop(handle)

   END SUBROUTINE tidft_env_create_matrices

! **************************************************************************************************
!> \brief deallocate tidft envirnoment data
!> \param tidft_env ...
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_env_deallocate(tidft_env)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_deallocate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: handle, ispin, &
                                                            nstates, nvars, &
                                                            nspins, ivar, &
                                                            stateI, &
                                                            nmos_stateI, &
                                                            nmos_stateI_flip, &
                                                            stateJ, &
                                                            nmos_stateJ

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nstates = tidft_env%nstates
      nspins = tidft_env%nspins
      nvars = tidft_env%nvars

      CALL dbcsr_release(tidft_env%mat_s)

      DO ispin = 1, nspins

         DO stateI = 1, nstates

            CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin))
            IF (tidft_env%spin_pure) THEN

               CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(2)%mat_t(ispin))

            ENDIF

            DO stateJ = stateI, nstates

               CALL dbcsr_release(tidft_env%mat_sigma(stateI, stateJ, ispin))

            ENDDO ! stateJ

         ENDDO ! stateI

      ENDDO ! ispin

      DO ivar = 1, nvars

         DO stateI = 1, nstates

            CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))

         ENDDO ! stateI

      ENDDO ! ivar

      DEALLOCATE(tidft_env%mat_sigma)
      DO stateI = 1, nstates
         DEALLOCATE(tidft_env%e_states(stateI)%slater_det(1)%mat_t)
         IF (tidft_env%spin_pure) DEALLOCATE(tidft_env%e_states(stateI)%slater_det(2)%mat_t)
         DEALLOCATE(tidft_env%e_states(stateI)%slater_det(1)%mat_var)
      ENDDO

      CALL cp_para_env_release(tidft_env%para_env)
      CALL cp_blacs_env_release(tidft_env%blacs_env)

      CALL timestop(handle)

   END SUBROUTINE tidft_env_deallocate

END MODULE tidft_env_methods

