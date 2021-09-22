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

   USE tidft_types,                     ONLY: tidft_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: &
        almo_constraint_distance, almo_deloc_none, almo_deloc_xalmo_1diag, &
        almo_domain_layout_atomic, almo_domain_layout_molecular, almo_frz_crystal, &
        almo_mat_distr_molecular, almo_scf_diag, almo_scf_skip, almo_scf_trustr, cg_hager_zhang, &
        do_bondparm_vdw, molecular_guess, tensor_orthogonal, virt_full, virt_minimal, virt_number, &
        xalmo_trial_r0_out
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

   PUBLIC :: tidft_env_create

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
!tiidft      almo_scf_env%nspins = dft_control%nspins
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
                                                            counter_i, k, jj, &
                                                            nmultips
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
!tiidft      almo_opt_pcg_section => section_vals_get_subs_vals(almo_scf_section, &
!tiidft                                                         "ALMO_OPTIMIZER_PCG")
!tiidft      xalmo_opt_pcg_section => section_vals_get_subs_vals(almo_scf_section, &
!tiidft                                                          "XALMO_OPTIMIZER_PCG")
!tiidft      xalmo_opt_trustr_section => section_vals_get_subs_vals(almo_scf_section, &
!tiidft                                                             "XALMO_OPTIMIZER_TRUSTR")
!tiidft      nlmo_opt_pcg_section => section_vals_get_subs_vals(almo_scf_section, &
!tiidft                                                         "NLMO_OPTIMIZER_PCG")
!tiidft      almo_analysis_section => section_vals_get_subs_vals(almo_scf_section, "ANALYSIS")
!tiidft      xalmo_opt_newton_pcg_section => section_vals_get_subs_vals(xalmo_opt_pcg_section, &
!tiidft                                                                 "XALMO_NEWTON_PCG_SOLVER")
!tiidft      matrix_iterate_section => section_vals_get_subs_vals(almo_scf_section, &
!tiidft                                                           "MATRIX_ITERATE")
!tiidft

      CALL section_vals_val_get(tidft_section, "EPS_FILTER", &
                                r_val=tidft_env%eps_filter)
      !CALL section_vals_val_get(tidft_section, "NSTATES", &
      !                          i_val=tidft_env%nstates)

      ! Read in the number of excited state
      CALL section_vals_val_get(tidft_section,"N_EXCITED_STATES",&
             i_rep_val=1, i_vals=tmplist)
      nmultips = SIZE(tmplist)
      ALLOCATE(nstates(nmultips))
      DO jj = 1, nmultips
         !WRITE(*,*) "Multiplicity: ", jj, "nstates: ", tmplist(jj)
         nstates(jj)=tmplist(jj)
      ENDDO

      WRITE(*,*) "In TIDFT we got EXCITED STATES: ", nstates(:)
      jj = 0
      ! CINDY: remove states with incompatible multiplicities
      !GS_multip = 2
      !tidft_env%nmultips=nmultips
      !tidft_env%nstates=nstates
      !CPABORT("Your error message")
      ! end CINDY

      DEALLOCATE(nstates)
!tiidft      CALL section_vals_val_get(almo_scf_section, "ALMO_ALGORITHM", &
!tiidft                                i_val=almo_scf_env%almo_update_algorithm)
!tiidft      CALL section_vals_val_get(almo_scf_section, "XALMO_ALGORITHM", &
!tiidft                                i_val=almo_scf_env%xalmo_update_algorithm)
!tiidft      CALL section_vals_val_get(almo_scf_section, "XALMO_TRIAL_WF", &
!tiidft                                i_val=almo_scf_env%xalmo_trial_wf)
!tiidft      CALL section_vals_val_get(almo_scf_section, "MO_OVERLAP_INV_ALG", &
!tiidft                                i_val=almo_scf_env%sigma_inv_algorithm)
!tiidft      CALL section_vals_val_get(almo_scf_section, "DELOCALIZE_METHOD", &
!tiidft                                i_val=almo_scf_env%deloc_method)
!tiidft      CALL section_vals_val_get(almo_scf_section, "XALMO_R_CUTOFF_FACTOR", &
!tiidft                                r_val=almo_scf_env%quencher_r0_factor)
!tiidft      CALL section_vals_val_get(almo_scf_section, "ALMO_EXTRAPOLATION_ORDER", &
!tiidft                                i_val=almo_scf_env%almo_extrapolation_order)
!tiidft      almo_scf_env%almo_extrapolation_order = MAX(0, almo_scf_env%almo_extrapolation_order)
!tiidft      CALL section_vals_val_get(almo_scf_section, "XALMO_EXTRAPOLATION_ORDER", &
!tiidft                                i_val=almo_scf_env%xalmo_extrapolation_order)
!tiidft      almo_scf_env%xalmo_extrapolation_order = MAX(0, almo_scf_env%xalmo_extrapolation_order)
!tiidft      CALL section_vals_val_get(almo_scf_section, "RETURN_ORTHOGONALIZED_MOS", &
!tiidft                                l_val=almo_scf_env%return_orthogonalized_mos)
!tiidft      CALL section_vals_val_get(almo_scf_section, "CONSTRUCT_NLMOS", &
!tiidft                                l_val=almo_scf_env%construct_nlmos)
!tiidft
!tiidft      CALL section_vals_val_get(matrix_iterate_section, "EPS_LANCZOS", &
!tiidft                                r_val=almo_scf_env%eps_lanczos)
!tiidft      CALL section_vals_val_get(matrix_iterate_section, "ORDER_LANCZOS", &
!tiidft                                i_val=almo_scf_env%order_lanczos)
!tiidft      CALL section_vals_val_get(matrix_iterate_section, "MAX_ITER_LANCZOS", &
!tiidft                                i_val=almo_scf_env%max_iter_lanczos)
!tiidft      CALL section_vals_val_get(matrix_iterate_section, "EPS_TARGET_FACTOR", &
!tiidft                                r_val=almo_scf_env%matrix_iter_eps_error_factor)
!tiidft
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
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_block_diag_pcg%eps_error)
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "MAX_ITER", &
!tiidft                                i_val=almo_scf_env%opt_block_diag_pcg%max_iter)
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "EPS_ERROR_EARLY", &
!tiidft                                r_val=almo_scf_env%opt_block_diag_pcg%eps_error_early)
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "MAX_ITER_EARLY", &
!tiidft                                i_val=almo_scf_env%opt_block_diag_pcg%max_iter_early)
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "MAX_ITER_OUTER_LOOP", &
!tiidft                                i_val=almo_scf_env%opt_block_diag_pcg%max_iter_outer_loop)
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "LIN_SEARCH_EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_block_diag_pcg%lin_search_eps_error)
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "LIN_SEARCH_STEP_SIZE_GUESS", &
!tiidft                                r_val=almo_scf_env%opt_block_diag_pcg%lin_search_step_size_guess)
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "CONJUGATOR", &
!tiidft                                i_val=almo_scf_env%opt_block_diag_pcg%conjugator)
!tiidft      CALL section_vals_val_get(almo_opt_pcg_section, "PRECONDITIONER", &
!tiidft                                i_val=almo_scf_env%opt_block_diag_pcg%preconditioner)
!tiidft
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_trustr%eps_error)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "MAX_ITER", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_trustr%max_iter)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "ALGORITHM", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_trustr%trustr_algorithm)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "EPS_ERROR_EARLY", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_trustr%eps_error_early)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "MAX_ITER_EARLY", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_trustr%max_iter_early)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "MAX_ITER_OUTER_LOOP", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_trustr%max_iter_outer_loop)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "ETA", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_trustr%rho_do_not_update)
!tiidft      almo_scf_env%opt_xalmo_trustr%rho_do_not_update = &
!tiidft         MIN(MAX(almo_scf_env%opt_xalmo_trustr%rho_do_not_update, 0.0_dp), 0.25_dp)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "MODEL_GRAD_NORM_RATIO", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_trustr%model_grad_norm_ratio)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "INITIAL_TRUST_RADIUS", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_trustr%initial_trust_radius)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "MAX_TRUST_RADIUS", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_trustr%max_trust_radius)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "CONJUGATOR", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_trustr%conjugator)
!tiidft      CALL section_vals_val_get(xalmo_opt_trustr_section, "PRECONDITIONER", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_trustr%preconditioner)
!tiidft
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_pcg%eps_error)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "MAX_ITER", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_pcg%max_iter)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "EPS_ERROR_EARLY", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_pcg%eps_error_early)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "MAX_ITER_EARLY", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_pcg%max_iter_early)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "MAX_ITER_OUTER_LOOP", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_pcg%max_iter_outer_loop)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "LIN_SEARCH_EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_pcg%lin_search_eps_error)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "LIN_SEARCH_STEP_SIZE_GUESS", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_pcg%lin_search_step_size_guess)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "PRECOND_FILTER_THRESHOLD", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_pcg%neglect_threshold)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "CONJUGATOR", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_pcg%conjugator)
!tiidft      CALL section_vals_val_get(xalmo_opt_pcg_section, "PRECONDITIONER", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_pcg%preconditioner)
!tiidft
!tiidft      penalty_section => section_vals_get_subs_vals(nlmo_opt_pcg_section, "PENALTY")
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%eps_error)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "MAX_ITER", &
!tiidft                                i_val=almo_scf_env%opt_nlmo_pcg%max_iter)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "EPS_ERROR_EARLY", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%eps_error_early)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "MAX_ITER_EARLY", &
!tiidft                                i_val=almo_scf_env%opt_nlmo_pcg%max_iter_early)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "MAX_ITER_OUTER_LOOP", &
!tiidft                                i_val=almo_scf_env%opt_nlmo_pcg%max_iter_outer_loop)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "LIN_SEARCH_EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%lin_search_eps_error)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "LIN_SEARCH_STEP_SIZE_GUESS", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%lin_search_step_size_guess)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "PRECOND_FILTER_THRESHOLD", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%neglect_threshold)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "CONJUGATOR", &
!tiidft                                i_val=almo_scf_env%opt_nlmo_pcg%conjugator)
!tiidft      CALL section_vals_val_get(nlmo_opt_pcg_section, "PRECONDITIONER", &
!tiidft                                i_val=almo_scf_env%opt_nlmo_pcg%preconditioner)
!tiidft      CALL section_vals_val_get(penalty_section, &
!tiidft                                "OPERATOR", &
!tiidft                                i_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%operator_type)
!tiidft      CALL section_vals_val_get(penalty_section, &
!tiidft                                "PENALTY_STRENGTH", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength)
!tiidft      CALL section_vals_val_get(penalty_section, &
!tiidft                                "PENALTY_STRENGTH_DECREASE_FACTOR", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength_dec_factor)
!tiidft      CALL section_vals_val_get(penalty_section, &
!tiidft                                "DETERMINANT_TOLERANCE", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%determinant_tolerance)
!tiidft      CALL section_vals_val_get(penalty_section, &
!tiidft                                "FINAL_DETERMINANT", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%final_determinant)
!tiidft      CALL section_vals_val_get(penalty_section, &
!tiidft                                "COMPACTIFICATION_FILTER_START", &
!tiidft                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%compactification_filter_start)
!tiidft      CALL section_vals_val_get(penalty_section, &
!tiidft                                "VIRTUAL_NLMOS", &
!tiidft                                l_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%virtual_nlmos)
!tiidft
!tiidft      CALL section_vals_val_get(xalmo_opt_newton_pcg_section, "EPS_ERROR", &
!tiidft                                r_val=almo_scf_env%opt_xalmo_newton_pcg_solver%eps_error)
!tiidft      CALL section_vals_val_get(xalmo_opt_newton_pcg_section, "MAX_ITER", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_newton_pcg_solver%max_iter)
!tiidft      CALL section_vals_val_get(xalmo_opt_newton_pcg_section, "MAX_ITER_OUTER_LOOP", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_newton_pcg_solver%max_iter_outer_loop)
!tiidft      CALL section_vals_val_get(xalmo_opt_newton_pcg_section, "PRECONDITIONER", &
!tiidft                                i_val=almo_scf_env%opt_xalmo_newton_pcg_solver%preconditioner)
!tiidft
!tiidft      CALL section_vals_val_get(almo_analysis_section, "_SECTION_PARAMETERS_", &
!tiidft                                l_val=almo_scf_env%almo_analysis%do_analysis)
!tiidft      CALL section_vals_val_get(almo_analysis_section, "FROZEN_MO_ENERGY_TERM", &
!tiidft                                i_val=almo_scf_env%almo_analysis%frozen_mo_energy_term)
!tiidft
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DOMAIN_LAYOUT_AOS",&
!tiidft      !                          i_val=almo_scf_env%domain_layout_aos)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DOMAIN_LAYOUT_MOS",&
!tiidft      !                          i_val=almo_scf_env%domain_layout_mos)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"MATRIX_CLUSTERING_AOS",&
!tiidft      !                          i_val=almo_scf_env%mat_distr_aos)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"MATRIX_CLUSTERING_MOS",&
!tiidft      !                          i_val=almo_scf_env%mat_distr_mos)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"CONSTRAINT_TYPE",&
!tiidft      !                          i_val=almo_scf_env%constraint_type)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"MU",&
!tiidft      !                          r_val=almo_scf_env%mu)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"FIXED_MU",&
!tiidft      !                          l_val=almo_scf_env%fixed_mu)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"EPS_USE_PREV_AS_GUESS",&
!tiidft      !                          r_val=almo_scf_env%eps_prev_guess)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"MIXING_FRACTION",&
!tiidft      !                          r_val=almo_scf_env%mixing_fraction)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_TENSOR_TYPE",&
!tiidft      !                          i_val=almo_scf_env%deloc_cayley_tensor_type)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_CONJUGATOR",&
!tiidft      !                          i_val=almo_scf_env%deloc_cayley_conjugator)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_MAX_ITER",&
!tiidft      !                          i_val=almo_scf_env%deloc_cayley_max_iter)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_USE_OCC_ORBS",&
!tiidft      !                          l_val=almo_scf_env%deloc_use_occ_orbs)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_USE_VIRT_ORBS",&
!tiidft      !                          l_val=almo_scf_env%deloc_cayley_use_virt_orbs)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_LINEAR",&
!tiidft      !                          l_val=almo_scf_env%deloc_cayley_linear)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_EPS_CONVERGENCE",&
!tiidft      !                          r_val=almo_scf_env%deloc_cayley_eps_convergence)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_OCC_PRECOND",&
!tiidft      !                          l_val=almo_scf_env%deloc_cayley_occ_precond)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_VIR_PRECOND",&
!tiidft      !                          l_val=almo_scf_env%deloc_cayley_vir_precond)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"ALMO_UPDATE_ALGORITHM_BD",&
!tiidft      !                          i_val=almo_scf_env%almo_update_algorithm)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_TRUNCATE_VIRTUALS",&
!tiidft      !                          i_val=almo_scf_env%deloc_truncate_virt)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"DELOC_VIRT_PER_DOMAIN",&
!tiidft      !                          i_val=almo_scf_env%deloc_virt_per_domain)
!tiidft      !
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_EPS_CONVERGENCE",&
!tiidft      !                          r_val=almo_scf_env%opt_k_eps_convergence)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_MAX_ITER",&
!tiidft      !                          i_val=almo_scf_env%opt_k_max_iter)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_OUTER_MAX_ITER",&
!tiidft      !                          i_val=almo_scf_env%opt_k_outer_max_iter)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_TRIAL_STEP_SIZE",&
!tiidft      !                          r_val=almo_scf_env%opt_k_trial_step_size)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJUGATOR",&
!tiidft      !                          i_val=almo_scf_env%opt_k_conjugator)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_TRIAL_STEP_SIZE_MULTIPLIER",&
!tiidft      !                          r_val=almo_scf_env%opt_k_trial_step_size_multiplier)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJ_ITER_START",&
!tiidft      !                          i_val=almo_scf_env%opt_k_conj_iter_start)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_PREC_ITER_START",&
!tiidft      !                          i_val=almo_scf_env%opt_k_prec_iter_start)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJ_ITER_FREQ_RESET",&
!tiidft      !                          i_val=almo_scf_env%opt_k_conj_iter_freq)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"OPT_K_PREC_ITER_FREQ_UPDATE",&
!tiidft      !                          i_val=almo_scf_env%opt_k_prec_iter_freq)
!tiidft      !
!tiidft      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_RADIUS_TYPE",&
!tiidft      !                          i_val=almo_scf_env%quencher_radius_type)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_R0_FACTOR",&
!tiidft      !                          r_val=almo_scf_env%quencher_r0_factor)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_R1_FACTOR",&
!tiidft      !                          r_val=almo_scf_env%quencher_r1_factor)
!tiidft      !!CALL section_vals_val_get(almo_scf_section,"QUENCHER_R0_SHIFT",&
!tiidft      !!                          r_val=almo_scf_env%quencher_r0_shift)
!tiidft      !!
!tiidft      !!CALL section_vals_val_get(almo_scf_section,"QUENCHER_R1_SHIFT",&
!tiidft      !!                          r_val=almo_scf_env%quencher_r1_shift)
!tiidft      !!
!tiidft      !!almo_scf_env%quencher_r0_shift = cp_unit_to_cp2k(&
!tiidft      !!   almo_scf_env%quencher_r0_shift,"angstrom")
!tiidft      !!almo_scf_env%quencher_r1_shift = cp_unit_to_cp2k(&
!tiidft      !!   almo_scf_env%quencher_r1_shift,"angstrom")
!tiidft      !
!tiidft      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_AO_OVERLAP_0",&
!tiidft      !                          r_val=almo_scf_env%quencher_s0)
!tiidft      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_AO_OVERLAP_1",&
!tiidft      !                          r_val=almo_scf_env%quencher_s1)
!tiidft
!tiidft      !CALL section_vals_val_get(almo_scf_section,"ENVELOPE_AMPLITUDE",&
!tiidft      !                          r_val=almo_scf_env%envelope_amplitude)
!tiidft
!tiidft      almo_scf_env%domain_layout_aos = almo_domain_layout_molecular
!tiidft      almo_scf_env%domain_layout_mos = almo_domain_layout_molecular
!tiidft      almo_scf_env%mat_distr_aos = almo_mat_distr_molecular
!tiidft      almo_scf_env%mat_distr_mos = almo_mat_distr_molecular
!tiidft
!tiidft      almo_scf_env%constraint_type = almo_constraint_distance
!tiidft      almo_scf_env%mu = -0.1_dp
!tiidft      almo_scf_env%fixed_mu = .FALSE.
!tiidft      almo_scf_env%mixing_fraction = 0.45_dp
!tiidft      almo_scf_env%eps_prev_guess = almo_scf_env%eps_filter/1000.0_dp
!tiidft
!tiidft      almo_scf_env%deloc_cayley_tensor_type = tensor_orthogonal
!tiidft      almo_scf_env%deloc_cayley_conjugator = cg_hager_zhang
!tiidft      almo_scf_env%deloc_cayley_max_iter = 100
!tiidft      almo_scf_env%deloc_use_occ_orbs = .TRUE.
!tiidft      almo_scf_env%deloc_cayley_use_virt_orbs = .FALSE.
!tiidft      almo_scf_env%deloc_cayley_linear = .FALSE.
!tiidft      almo_scf_env%deloc_cayley_eps_convergence = 1.0E-6_dp
!tiidft      almo_scf_env%deloc_cayley_occ_precond = .TRUE.
!tiidft      almo_scf_env%deloc_cayley_vir_precond = .TRUE.
!tiidft      almo_scf_env%deloc_truncate_virt = virt_full
!tiidft      almo_scf_env%deloc_virt_per_domain = -1
!tiidft
!tiidft      almo_scf_env%opt_k_eps_convergence = 1.0E-5_dp
!tiidft      almo_scf_env%opt_k_max_iter = 100
!tiidft      almo_scf_env%opt_k_outer_max_iter = 1
!tiidft      almo_scf_env%opt_k_trial_step_size = 0.05_dp
!tiidft      almo_scf_env%opt_k_conjugator = cg_hager_zhang
!tiidft      almo_scf_env%opt_k_trial_step_size_multiplier = 1.05_dp
!tiidft      almo_scf_env%opt_k_conj_iter_start = 0
!tiidft      almo_scf_env%opt_k_prec_iter_start = 0
!tiidft      almo_scf_env%opt_k_conj_iter_freq = 10000000
!tiidft      almo_scf_env%opt_k_prec_iter_freq = 10000000
!tiidft
!tiidft      almo_scf_env%quencher_radius_type = do_bondparm_vdw
!tiidft      almo_scf_env%quencher_r1_factor = almo_scf_env%quencher_r0_factor
!tiidft      !almo_scf_env%quencher_r0_shift=0.0_dp
!tiidft      !almo_scf_env%quencher_r1_shift=0.0_dp
!tiidft      !almo_scf_env%quencher_r0_shift = cp_unit_to_cp2k(&
!tiidft      !   almo_scf_env%quencher_r0_shift,"angstrom")
!tiidft      !almo_scf_env%quencher_r1_shift = cp_unit_to_cp2k(&
!tiidft      !   almo_scf_env%quencher_r1_shift,"angstrom")
!tiidft
!tiidft      almo_scf_env%quencher_s0 = 1.0E-4_dp
!tiidft      almo_scf_env%quencher_s1 = 1.0E-6_dp
!tiidft
!tiidft      almo_scf_env%envelope_amplitude = 1.0_dp
!tiidft
!tiidft      almo_scf_env%logical01 = .FALSE. ! md in eDOF space
!tiidft      almo_scf_env%logical02 = .TRUE. ! not used
!tiidft      almo_scf_env%logical03 = .TRUE. ! not used
!tiidft      almo_scf_env%logical04 = .TRUE. ! use preconditioner
!tiidft      almo_scf_env%logical05 = .FALSE. ! optimize theta
!tiidft
!tiidft      almo_scf_env%real01 = almo_scf_env%eps_filter/10.0_dp ! skip gradients
!tiidft      almo_scf_env%real02 = 0.0_dp ! not used
!tiidft      almo_scf_env%real03 = 0.0_dp ! not used
!tiidft      almo_scf_env%real04 = 0.5_dp ! mixing s-f precond
!tiidft
!tiidft      almo_scf_env%integer01 = 10 ! start eDOF-md
!tiidft      almo_scf_env%integer02 = 4 ! preconditioner type
!tiidft      almo_scf_env%integer03 = 0 ! not used
!tiidft      almo_scf_env%integer04 = 0 ! fixed number of line searches (no grad)
!tiidft      almo_scf_env%integer05 = 0 ! not used
!tiidft
!tiidft      ! check for conflicts between options
!tiidft      IF (almo_scf_env%xalmo_trial_wf .EQ. xalmo_trial_r0_out .AND. &
!tiidft          almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_trustr) THEN
!tiidft         CPABORT("Trust region algorithm cannot optimize projected XALMOs")
!tiidft      ENDIF
!tiidft
!tiidft      CALL section_vals_val_get(almo_scf_section, "XALMO_ALGORITHM", &
!tiidft                                i_val=almo_scf_env%xalmo_update_algorithm)
!tiidft      CALL section_vals_val_get(almo_scf_section, "XALMO_TRIAL_WF", &
!tiidft                                i_val=almo_scf_env%xalmo_trial_wf)
!tiidft      IF (almo_scf_env%deloc_method .EQ. almo_deloc_xalmo_1diag .AND. &
!tiidft          almo_scf_env%xalmo_update_algorithm .NE. almo_scf_diag) THEN
!tiidft         CPABORT("1-step delocalization correction requires a different algorithm")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%xalmo_trial_wf .EQ. xalmo_trial_r0_out .AND. &
!tiidft          almo_scf_env%almo_update_algorithm .EQ. almo_scf_skip .AND. &
!tiidft          almo_scf_env%almo_scf_guess .NE. molecular_guess) THEN
!tiidft         CPABORT("R0 projector requires optimized ALMOs")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%deloc_method .EQ. almo_deloc_none .AND. &
!tiidft          almo_scf_env%almo_update_algorithm .EQ. almo_scf_skip) THEN
!tiidft         CPABORT("No optimization requested")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%deloc_truncate_virt .EQ. virt_number .AND. &
!tiidft          almo_scf_env%deloc_virt_per_domain .LE. 0) THEN
!tiidft         CPABORT("specify a positive number of virtual orbitals")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%deloc_truncate_virt .EQ. virt_minimal) THEN
!tiidft         CPABORT("VIRT TRUNCATION TO MINIMAL BASIS IS NIY")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%domain_layout_mos .NE. almo_domain_layout_molecular) THEN
!tiidft         CPABORT("use MOLECULAR domains")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%domain_layout_aos .NE. almo_domain_layout_molecular) THEN
!tiidft         CPABORT("use MOLECULAR domains")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%mat_distr_mos .NE. almo_mat_distr_molecular) THEN
!tiidft         CPABORT("use MOLECULAR distr for MOs")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%mat_distr_aos == almo_mat_distr_molecular .AND. &
!tiidft          almo_scf_env%domain_layout_aos == almo_domain_layout_atomic) THEN
!tiidft         CPABORT("AO blocks cannot be larger than domains")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%mat_distr_mos == almo_mat_distr_molecular .AND. &
!tiidft          almo_scf_env%domain_layout_mos == almo_domain_layout_atomic) THEN
!tiidft         CPABORT("MO blocks cannot be larger than domains")
!tiidft      ENDIF
!tiidft
!tiidft      IF (almo_scf_env%quencher_r1_factor .GT. almo_max_cutoff_multiplier) THEN
!tiidft         CALL cp_abort(__LOCATION__, &
!tiidft                       "XALMO_R_CUTOFF_FACTOR is larger than almo_max_cutoff_multiplier. "// &
!tiidft                       "Increase the hard-coded almo_max_cutoff_multiplier")
!tiidft      ENDIF
!tiidft
!tiidft      ! check analysis settings
!tiidft      IF (almo_scf_env%almo_analysis%do_analysis) THEN
!tiidft
!tiidft         IF (almo_scf_env%almo_analysis%frozen_mo_energy_term == almo_frz_crystal &
!tiidft             .AND. almo_scf_env%almo_scf_guess .NE. molecular_guess) THEN
!tiidft            CPABORT("To compute frozen-MO energy term set ALMO_SCF_GUESS MOLECULAR")
!tiidft         ENDIF
!tiidft
!tiidft      ENDIF ! end analysis settings

      CALL timestop(handle)

   END SUBROUTINE tidft_read_write_input

END MODULE tidft_env_methods

