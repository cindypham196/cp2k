!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Variable-metric TIDFT optimization of molecular orbitals
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_optimizer_methods
   USE population_analyses,             ONLY: lowdin_population_analysis,&
                                              mulliken_population_analysis
   USE almo_scf_methods,                ONLY: almo_scf_t_to_proj,&
                                              get_overlap
   USE almo_scf_types,                  ONLY: optimizer_options_type
   USE cp_dbcsr_singular,               ONLY: cp_dbcsr_pseudoinverse, &
                                              cp_dbcsr_determinant_generalized, &
                                              cp_dbcsr_adjugate
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE tidft_guess_methods,             ONLY: construct_canonical_orbitals
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_filter, dbcsr_multiply, dbcsr_norm, &
        dbcsr_norm_maxabsnorm, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry, dbcsr_distribution_type, dbcsr_get_info, dbcsr_get_diag, &
        dbcsr_transposed, dbcsr_print, dbcsr_set_diag, dbcsr_dot, dbcsr_desymmetrize
   !USE dbcsr_api,                       ONLY: &
   !     dbcsr_add, dbcsr_add_on_diag, dbcsr_binary_read, dbcsr_checksum, dbcsr_copy, dbcsr_create, &
   !     dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
   !     dbcsr_get_info, dbcsr_get_stored_coordinates, &
   !     dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
   !     dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_nblkcols_total, &
   !     dbcsr_nblkrows_total, dbcsr_p_type, dbcsr_release, dbcsr_reserve_block2d, dbcsr_scale, &
   !     dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_work_create, &
   !     dbcsr_print, dbcsr_get_diag, dbcsr_set_diag, dbcsr_transposed
   USE input_constants,                 ONLY: &
        cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, &
        cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero, &
        tidft_prec_zero, tidft_prec_s, tidft_prec_h_s, tidft_prec_qhq_s, tidft_prec_qhq_qsq
   !USE input_constants,                 ONLY: &
   !     almo_deloc_none, almo_deloc_scf, almo_deloc_x, almo_deloc_x_then_scf, &
   !     almo_deloc_xalmo_1diag, almo_deloc_xalmo_scf, almo_deloc_xalmo_x, almo_deloc_xk, &
   !     almo_domain_layout_molecular, almo_mat_distr_atomic, almo_mat_distr_molecular, &
   !     almo_scf_diag, almo_scf_dm_sign, almo_scf_pcg, almo_scf_skip, almo_scf_trustr, &
   !     atomic_guess, molecular_guess, optimizer_diis, optimizer_lin_eq_pcg, optimizer_pcg, &
   !     optimizer_trustr, restart_guess, smear_fermi_dirac, virt_full, virt_number, virt_occ_size, &
   !     xalmo_case_block_diag, xalmo_case_fully_deloc, xalmo_case_normal, xalmo_trial_r0_out
   USE iterate_matrix,                  ONLY: determinant
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_sum
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE tidft_types,                     ONLY: tidft_env_type
   USE tidft_optimizer_types,           ONLY: tidft_optimizer_storage_type
   USE tidft_methods,                   ONLY: tidft_optimizer_mainvar_dot, &
                                              tidft_optimizer_mainvar_add, &
                                              tidft_optimizer_mainvar_multiply, &
                                              tidft_optimizer_mainvar_multiply_diff, &
                                              tidft_optimizer_mainvar_copy, &
                                              tidft_optimizer_gradient_norm, &
                                              tidft_optimizer_mainvar_print, &
                                              tidft_optimizer_mainvar_set, &
                                              tidft_mainvar_orthonormalize, &
                                              dbcsr_set_element, &
                                              stitch_columns, &
                                              split_columns, &
                                              determinant_local, &
                                              matrix_invert_local, &
                                              gradMO_to_gradVAR_increment, &
                                              precMO_to_precVAR_increment
   USE tidft_qs,                        ONLY: tidft_dm_to_tidft_ks
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_optimizer_methods'

   PUBLIC :: tidft_variational_optimization

   LOGICAL, PARAMETER :: debug_mode = .TRUE.
   LOGICAL, PARAMETER :: safe_mode = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Decides on whether we do state-by-state optimization or optimize all states simultaneously
!> \param qs_env ...
!> \param tidft_env ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_variational_optimization(qs_env, tidft_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_variational_optimization', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=59)                                  :: print_string
      INTEGER                                            :: unit_nr, &
                                                            nstates, &
                                                            current_no_of_states
      REAL(KIND=dp)                                      :: prev_determinant
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_energy_type), POINTER                      :: qs_energy
      REAL(KIND=dp)                                      :: save_initial_intra
      REAL(KIND=dp)                                      :: save_initial_inter
      LOGICAL, DIMENSION(:), ALLOCATABLE                 :: save_component_flags

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ALLOCATE(save_component_flags(SIZE(tidft_env%loss_function_components)))

      nstates = tidft_env%nstates

      IF (tidft_env%state_by_state) THEN

         save_initial_intra = tidft_env%intrastate_nonorthogonality_penalty_coeff
         save_initial_inter = tidft_env%interstate_nonorthogonality_penalty_coeff

         DO current_no_of_states = 1, nstates
   
            tidft_env%nstates = current_no_of_states

            IF (current_no_of_states .GE. 2 .AND. tidft_env%pre_orthogonalize_states) THEN
            ! For excited states, first inter-orthogonalize the current state to 
            ! those that are already converged

               IF (unit_nr > 0) THEN

                  WRITE (unit_nr, '()')
                  WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
                  WRITE (unit_nr, '(T2,A59,I6)') &
                     "INTER-ORTHOGONALIZATION OF STATE", &
                     !"OPTIMIZATION OF STATE", &
                     current_no_of_states
                  WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
                  WRITE (unit_nr, '()')

               ENDIF

               save_component_flags(:)=tidft_env%loss_function_components(:)
               tidft_env%loss_function_components(1) = .FALSE. ! energy is off
               tidft_env%loss_function_components(2:3) = .TRUE. ! nonorthogonality penalty is on

               CALL penalty_strength_variation(qs_env, tidft_env)

               tidft_env%loss_function_components(:)=save_component_flags(:)

               tidft_env%intrastate_nonorthogonality_penalty_coeff = &
                  save_initial_intra
               tidft_env%interstate_nonorthogonality_penalty_coeff = &
                  save_initial_inter

            ENDIF

            IF (unit_nr > 0) THEN

               WRITE (unit_nr, '()')
               WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
               WRITE (unit_nr, '(T2,A59,I6)') &
                  "OPTIMIZATION OF STATE", &
                  current_no_of_states
               WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
               WRITE (unit_nr, '()')

            ENDIF

            IF (current_no_of_states .EQ. 2) THEN
            ! Done with the ground state, 
            ! generate a better guess for excited states
               !CALL construct_canonical_orbitals(tidft_env, 1)
               ! Argument 1 indicates that the GS WF is used to generate the guess
            ENDIF

            IF (current_no_of_states .EQ. 1) THEN
            ! For the ground state, turn off the interstate penalty for faster optimization

               save_component_flags(:)=tidft_env%loss_function_components(:)
               tidft_env%loss_function_components(3) = .FALSE.

            ENDIF

            CALL penalty_strength_variation(qs_env, tidft_env)

            IF (current_no_of_states .EQ. 1) THEN
               tidft_env%loss_function_components(:)=save_component_flags(:)
            ENDIF

            tidft_env%intrastate_nonorthogonality_penalty_coeff = &
               save_initial_intra
            tidft_env%interstate_nonorthogonality_penalty_coeff = &
               save_initial_inter

         ENDDO
   
      ELSE

         IF (unit_nr > 0) THEN

            WRITE (unit_nr, '()')
            WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
            WRITE (unit_nr, '(T2,A59)') &
               "SIMULTANEOUS OPTIMIZATION OF ALL STATES"
            WRITE (unit_nr, '(T2,A)') REPEAT("=", 79)
            WRITE (unit_nr, '()')

         ENDIF

         CALL penalty_strength_variation(qs_env, tidft_env)

      ENDIF

      DEALLOCATE(save_component_flags)

   END SUBROUTINE tidft_variational_optimization

! **************************************************************************************************
!> \brief Optimizes the penalty strenght and calls fixed-strength TIDFT optimization
!> \param qs_env ...
!> \param tidft_env ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE penalty_strength_variation(qs_env, tidft_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'penalty_strength_variation', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=59)                                  :: print_string
      INTEGER                                            :: unit_nr
      REAL(KIND=dp)                                      :: prev_determinant
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_energy_type), POINTER                      :: qs_energy

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! loop over the orthogonalization penalty strength
      ! vary both intra and inter strengths together
      !prev_determinant = HUGE(1.0_dp)
      DO

         IF (unit_nr > 0) THEN

            WRITE (unit_nr, '()')
            WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
            print_string = "Penalty strength:"
            WRITE (unit_nr, '(T2,A59,F20.10)') print_string, &
               tidft_env%intrastate_nonorthogonality_penalty_coeff
            WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
            WRITE (unit_nr, '()')

            ! print the optimization header
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 26), &
               " TIDFT optimization of MOs ", REPEAT("-", 26)
            WRITE (unit_nr, '()')

            WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

         ENDIF

         !SELECT CASE (almo_scf_env%construct_nlmos)
         !CASE (almo_scf_pcg)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A13,A6,A23,A14,A14,A9)') "Method", "Iter", &
               "Objective Function", "Change", "Convergence", "Time"
            WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
         ENDIF
         CALL tidft_optimize_pcg( &
                 qs_env=qs_env, &
                 tidft_env=tidft_env, &
                 optimizer=tidft_env%tidft_optimizer_pcg)
         !CASE (almo_scf_trustr)
         !   CALL trust_r_report(unit_nr, &
         !                       iter_type=0, & ! print header, all values are ignored
         !                       iteration=0, &
         !                       radius=0.0_dp, &
         !                       loss=0.0_dp, &
         !                       delta_loss=0.0_dp, &
         !                       grad_norm=0.0_dp, &
         !                       predicted_reduction=0.0_dp, &
         !                       rho=0.0_dp, &
         !                       new=.TRUE., &
         !                       time=0.0_dp)
         !   CALL construct_nlmos_trustr(qs_env=qs_env, &
         !                               nlmo_env=nlmo_env, &
         !                               optimizer=almo_scf_env%opt_nlmo_trustr, &
         !                               matrix_s=almo_scf_env%matrix_s(1), &
         !                               matrix_mo_inout=almo_scf_env%matrix_t, &
         !                               matrix_templateOO=almo_scf_env%matrix_sigma_inv)
         !CASE DEFAULT
         !   CPABORT("Illegal optmizer")
         !END SELECT

         !det_diff = prev_determinant - tidft_env%det_mat_phi

         ! Check whether the determinant is sufficiently close to 1
         IF ( ABS(1.0_dp - tidft_env%det_mat_phi) .LE. tidft_env%allowed_determinant_deviaion &
             .OR. .NOT.tidft_env%loss_function_components(3)) THEN
            EXIT
         ENDIF

         !prev_determinant = tidft_env%det_mat_phi

         tidft_env%intrastate_nonorthogonality_penalty_coeff = &
            tidft_env%intrastate_nonorthogonality_penalty_coeff * &
            ABS(tidft_env%penalty_strength_inc_factor)

         tidft_env%interstate_nonorthogonality_penalty_coeff = &
            tidft_env%interstate_nonorthogonality_penalty_coeff * &
            ABS(tidft_env%penalty_strength_inc_factor)

      ENDDO

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

         print_string = "Final energy:"
         CALL get_qs_env(qs_env, energy=qs_energy)
         WRITE (unit_nr, '(T2,A59,F20.10)') &
            print_string, qs_energy%total

         print_string = "Final interstate determinant:"
         WRITE (unit_nr, '(T2,A59,F20.10)') &
            print_string, tidft_env%det_mat_phi

         print_string = "Final intrastate penalty strength:"
         WRITE (unit_nr, '(T2,A59,F20.10)') &
            print_string, tidft_env%intrastate_nonorthogonality_penalty_coeff

         print_string = "Final interstate penalty strength:"
         WRITE (unit_nr, '(T2,A59,F20.10)') &
            print_string, tidft_env%interstate_nonorthogonality_penalty_coeff

         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

      ENDIF

   END SUBROUTINE penalty_strength_variation

! **************************************************************************************************
!> \brief Optimization of MOs using PCG minimizers
!> \param qs_env ...
!> \param tidft_env ...
!> \param optimizer - controls the optimization algorithm
!> \par History
!>       2018.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimize_pcg(qs_env, tidft_env, optimizer)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env
      TYPE(optimizer_options_type), INTENT(INOUT)        :: optimizer

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_optimize_pcg', &
         routineP = moduleN//':'//routineN

      TYPE(tidft_optimizer_storage_type)                 :: optimizer_storage
      CHARACTER(LEN=30)                                  :: iter_type
      INTEGER :: cg_iteration, handle,  iteration, line_search_iteration, &
         linear_search_type, max_iter, outer_iteration, outer_max_iter, unit_nr
      LOGICAL                                            :: converged, just_started, line_search, &
                                                            outer_prepare_to_exit, &
                                                            prepare_to_exit, reset_conjugator
      REAL(KIND=dp) :: appr_sec_der, beta, denom, denom2, e0, e1, g0, g0sign, g1, g1sign, &
         grad_norm, line_search_error, energy_function, next_step_size_guess, &
         obj_function_ispin, objf_diff, objf_new, objf_old, penalty_func_new, step_size, t1, t2
      REAL(KIND=dp), DIMENSION(100)                      :: loss_components
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

! STEP 1: init
      CALL tidft_optimizer_allocate_storage( &
              data_storage = optimizer_storage, &
              nstates = tidft_env%nstates, &
              mat_var_template = tidft_env%e_states(1)%slater_det(1)%mat_var, &
              mat_nn_template = tidft_env%mat_s, &
              ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)

! STEP 2: Start iterations

      ! start the outer SCF loop
      outer_max_iter = optimizer%max_iter_outer_loop
      outer_prepare_to_exit = .FALSE.
      outer_iteration = 0
      grad_norm = 0.0_dp
      penalty_func_new = 0.0_dp
      linear_search_type = 1 ! quadratic search
      !linear_search_type = 2 ! safe restart, no quadratic assumption, takes more steps
      energy_function = 0.0_dp

      DO

         ! start the inner SCF loop
         max_iter = optimizer%max_iter
         prepare_to_exit = .FALSE.
         line_search = .FALSE.
         converged = .FALSE.
         iteration = 0
         cg_iteration = 0
         line_search_iteration = 0
         obj_function_ispin = 0.0_dp
         objf_new = 0.0_dp
         objf_old = 0.0_dp
         objf_diff = 0.0_dp
         line_search_error = 0.0_dp
         t1 = m_walltime()
         next_step_size_guess = 0.0_dp
         beta = 0.0_dp

         DO

            just_started = (iteration .EQ. 0) .AND. (outer_iteration .EQ. 0)

            ! save the previous gradient to compute beta
            ! do it only if the previous grad was computed
            ! for .NOT.line_search
            IF (line_search_iteration .EQ. 0 .AND. iteration .NE. 0) THEN
               CALL tidft_optimizer_mainvar_copy( &
                       optimizer_storage%prev_grad, &
                       optimizer_storage%grad, &
                       ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)
            ENDIF

! STEP 3: MAJOR Theta to MOs
! STEP 4: Obj function
! STEP 5: Gradient
            ! use the main variable to compute aux matrices
            ! compute objective (loss) function
            ! compute the gradient
            loss_components(:) = 0.0_dp
            IF (tidft_env%numerical_gradient) THEN
               CALL tidft_optimizer_numerical_gradient( &
                       qs_env = qs_env, &
                       tidft_env = tidft_env, &
                       dx = tidft_env%numerical_gradient_dx, &
                       loss_components_out = loss_components, &
                       gradient_out = optimizer_storage%grad )
            ELSE
               CALL mainvar_to_loss_and_gradient( &
                       qs_env = qs_env, &
                       tidft_env = tidft_env, &
                       loss_components_out = loss_components, &
                       gradient_out = optimizer_storage%grad )
            ENDIF

            energy_function = loss_components(1)
            penalty_func_new = loss_components(2) + loss_components(3)
            objf_new = SUM(loss_components)

! STEP 6: Convergence
            ! check convergence and other exit criteria
            CALL tidft_optimizer_gradient_norm( &
                    optimizer_storage%grad, &
                    optimizer%which_gnorm, &
                    grad_norm, &
                    tidft_env%e_states(1)%slater_det(1)%ncols_var, &
                    metric_inv = tidft_env%mat_s_inv)

            converged = (grad_norm .LE. optimizer%eps_error)
            IF (converged .OR. (iteration .GE. max_iter)) THEN
               prepare_to_exit = .TRUE.
            ENDIF

            ! it is not time to exit just yet
            IF (.NOT. prepare_to_exit) THEN

               ! check the gradient along the step direction
               ! and decide whether to switch to the line-search mode
               ! do not do this in the first iteration
               IF (iteration .NE. 0) THEN

                  ! enforce at least one line search
                  ! without even checking the error
                  IF (.NOT. line_search) THEN

                     line_search = .TRUE.
                     line_search_iteration = line_search_iteration + 1

                  ELSE

                     ! check the line-search error and decide whether to
                     ! change the direction
                     line_search_error = 0.0_dp
                     denom = 0.0_dp
                     denom2 = 0.0_dp

! STEP 7: Line-search error
                     CALL tidft_optimizer_mainvar_dot( &
                             mat_left = optimizer_storage%grad, &
                             mat_right = optimizer_storage%step, &
                             mat_left_covariant = .TRUE., &
                             mat_right_covariant = .FALSE., &
                             dot_product = line_search_error, &
                             ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)

                     CALL tidft_optimizer_mainvar_dot( &
                             mat_left = optimizer_storage%grad, &
                             mat_right = optimizer_storage%grad, &
                             mat_left_covariant = .TRUE., &
                             mat_right_covariant = .TRUE., &
                             metric_inv = tidft_env%mat_s_inv, &
                             dot_product = denom, &
                             ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var, &
                             eps_filter = tidft_env%eps_filter)

                     CALL tidft_optimizer_mainvar_dot( &
                             mat_left = optimizer_storage%step, &
                             mat_right = optimizer_storage%step, &
                             mat_left_covariant = .FALSE., &
                             mat_right_covariant = .FALSE., &
                             metric = tidft_env%mat_s, &
                             dot_product = denom2, &
                             ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var, &
                             eps_filter = tidft_env%eps_filter)

                     !DO ispin = 1, nspins

                     !   CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), tempreal)
                     !   line_search_error = line_search_error + tempreal
                     !   CALL dbcsr_dot(nlmo_env%grad(ispin), nlmo_env%grad(ispin), tempreal)
                     !   denom = denom + tempreal
                     !   CALL dbcsr_dot(step(ispin), step(ispin), tempreal)
                     !   denom2 = denom2 + tempreal

                     !ENDDO ! ispin

                     ! cosine of the angle between the step and grad
                     ! (must be close to zero at convergence)
                     line_search_error = line_search_error/SQRT(denom)/SQRT(denom2)

                     IF (ABS(line_search_error) .GT. optimizer%lin_search_eps_error) THEN
                        line_search = .TRUE.
                        line_search_iteration = line_search_iteration + 1
                     ELSE
                        line_search = .FALSE.
                        line_search_iteration = 0
                     ENDIF

                  ENDIF

               ENDIF ! iteration.ne.0

               IF (line_search) THEN
                  objf_diff = 0.0_dp
               ELSE
                  objf_diff = objf_new - objf_old
                  objf_old = objf_new
               ENDIF

               ! update the step direction
               IF (.NOT. line_search) THEN

                  cg_iteration = cg_iteration + 1

                  ! save the previous step
                  CALL tidft_optimizer_mainvar_copy( &
                          optimizer_storage%prev_step, &
                          optimizer_storage%step, &
                          ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)

! STEP 8: Grad to step
                  ! compute the new step
                  CALL tidft_optimizer_grad_to_step( &
                          tidft_env = tidft_env, &
                          precond = optimizer%preconditioner, &
                          grad = optimizer_storage%grad, &
                          step = optimizer_storage%step, &
                          prec_inv = optimizer_storage%prec_inv, &
                          prev_grad = optimizer_storage%prev_grad, &
                          prev_theta = optimizer_storage%prev_theta, &
                          iteration = iteration)

!tiidft                  CALL nlmo_env_loss_curvature_grad_to_step(nlmo_env=nlmo_env, &
!tiidft                                                            grad=nlmo_env%grad, &
!tiidft                                                            step=step, &
!tiidft                                                            prev_grad=prev_grad, &
!tiidft                                                            prev_m_theta=prev_m_theta, &
!tiidft                                                            iteration=iteration)

                  ! check whether we need to reset conjugate directions
                  IF (iteration .EQ. 0) THEN
                     reset_conjugator = .TRUE.
                  ENDIF
! STEP 9: Compute beta
                  ! compute the conjugation coefficient - beta
                  IF (.NOT. reset_conjugator) THEN
!RZK-WARNING: some fancy tensorially inconsistent conjugators should be eliminated
!             erase their aux matrices as well
!tiidft                     CALL compute_cg_beta( &
!tiidft                        beta=beta, &
!tiidft                        reset_conjugator=reset_conjugator, &
!tiidft                        conjugator=optimizer%conjugator, &
!tiidft                        grad=nlmo_env%grad, &
!tiidft                        prev_grad=prev_grad, &
!tiidft                        step=step, &
!tiidft                        prev_step=prev_step, &
!tiidft                        prev_minus_prec_grad=prev_minus_prec_grad)
                     CALL compute_cg_beta( &
                             beta = beta, &
                             reset_conjugator = reset_conjugator, &
                             conjugator = optimizer%conjugator, &
                             grad = optimizer_storage%grad, &
                             prev_grad = optimizer_storage%prev_grad, &
                             step = optimizer_storage%step, &
                             prev_step = optimizer_storage%prev_step, &
                             prev_minus_prec_grad = optimizer_storage%prev_minus_prec_grad, &
                             ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)

                  ENDIF

                  ! Do we need to reset the conjugator?
                  ! One test shows that resetting is absolutely necessary with Fletcher
                  ! and is also slightly better than the zero conjugator
                  IF (reset_conjugator) THEN

                     beta = 0.0_dp
                     IF (unit_nr > 0 .AND. (.NOT. just_started)) THEN
                        WRITE (unit_nr, '(T2,A35)') "Re-setting conjugator to zero"
                     ENDIF
                     reset_conjugator = .FALSE.

                  ENDIF

! STEP 10: save the preconditioned gradient 
                  ! save the preconditioned gradient (useful for beta)
!tiidft                  DO ispin = 1, nspins
!tiidft
!tiidft                     CALL dbcsr_copy(prev_minus_prec_grad(ispin), step(ispin))
!tiidft
!tiidft                     ! conjugate the step direction
!tiidft                     CALL dbcsr_add(step(ispin), prev_step(ispin), 1.0_dp, beta)
!tiidft
!tiidft                  ENDDO ! ispin
                  CALL tidft_optimizer_mainvar_copy( &
                          optimizer_storage%prev_minus_prec_grad, &
                          optimizer_storage%step, &
                          ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)
                  CALL tidft_optimizer_mainvar_add( &
                          optimizer_storage%step, &
                          optimizer_storage%prev_step, &
                          1.0_dp, beta, &
                          ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)

               ENDIF ! update the step direction
! STEP 11: Line search code
               ! estimate the step size
               IF (.NOT. line_search) THEN
                  ! we just changed the direction and
                  ! we have only E and grad from the current step
                  ! it is not enough to compute step_size - just guess it
                  e0 = objf_new
                  g0 = 0.0_dp
                  CALL tidft_optimizer_mainvar_dot( &
                          mat_left = optimizer_storage%grad, &
                          mat_right = optimizer_storage%step, &
                          mat_left_covariant = .TRUE., &
                          mat_right_covariant = .FALSE., &
                          dot_product = g0, &
                          ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)
                  !DO ispin = 1, nspins
                  !   CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), tempreal)
                  !   g0 = g0 + tempreal
                  !ENDDO ! ispin
                  g0sign = SIGN(1.0_dp, g0) ! sign of g0
                  IF (linear_search_type .EQ. 1) THEN ! this is quadratic LS
                     IF (iteration .EQ. 0) THEN
                        step_size = optimizer%lin_search_step_size_guess
                     ELSE
                        IF (next_step_size_guess .LE. 0.0_dp) THEN
                           step_size = optimizer%lin_search_step_size_guess
                        ELSE
                           ! take the last value
                           ! CINDY modification : trust previous step, confident to make more aggressive step
                           !step_size = next_step_size_guess*1.05_dp
                           ! CHANGE BACK TO PREVIOUS SETTING BELOW AFTER THE O3 TEST
                           step_size = optimizer%lin_search_step_size_guess  
                        ENDIF
                     ENDIF
                  ELSE IF (linear_search_type .EQ. 2) THEN ! this is cautious LS
                     ! this LS type is designed not to trust quadratic appr
                     ! so it always restarts from a safe step size
                     step_size = optimizer%lin_search_step_size_guess
                  ENDIF
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A12,A5,4A14)') " ", " ", "Line position", "Line function", &
                        "Line grad", "Next step"
                     WRITE (unit_nr, '(T2,A12,I5,4F14.4)') "Line search", cg_iteration, &
                        0.0_dp, e0, g0, step_size
                  ENDIF
                  next_step_size_guess = step_size
               ELSE ! this is not the first line search
                  e1 = objf_new
                  g1 = 0.0_dp
                  CALL tidft_optimizer_mainvar_dot( &
                          mat_left = optimizer_storage%grad, &
                          mat_right = optimizer_storage%step, &
                          mat_left_covariant = .TRUE., &
                          mat_right_covariant = .FALSE., &
                          dot_product = g1, &
                          ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)
                  !DO ispin = 1, nspins
                  !   CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), tempreal)
                  !   g1 = g1 + tempreal
                  !ENDDO ! ispin
                  g1sign = SIGN(1.0_dp, g1) ! sign of g1
                  IF (linear_search_type .EQ. 1) THEN
                     ! we have accumulated some points along this direction
                     ! use only the most recent g0 (quadratic approximation)
                     appr_sec_der = (g1 - g0)/step_size
                     !IF (unit_nr > 0) THEN
                     !   WRITE (unit_nr, '(A2,7F12.5)') &
                     !      "DT", e0, e1, g0, g1, appr_sec_der, step_size, -g1/appr_sec_der
                     !ENDIF
                     step_size = -g1/appr_sec_der
                  ELSE IF (linear_search_type .EQ. 2) THEN
                     ! alternative method for finding step size
                     ! do not use quadratic approximation, only gradient signs
                     IF (g1sign .NE. g0sign) THEN
                        step_size = -step_size/2.0; 
                     ELSE
                        step_size = step_size*1.5; 
                     ENDIF
                  ENDIF
                  ! end alternative LS types
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A12,A5,4A14)') " ", " ", "Line position", "Line function", &
                        "Line grad", "Next step"
                     WRITE (unit_nr, '(T2,A12,I5,4F14.4)') "Line search", cg_iteration, &
                        next_step_size_guess, e1, g1, step_size
                  ENDIF
                  e0 = e1
                  g0 = g1
                  g0sign = g1sign
                  next_step_size_guess = next_step_size_guess + step_size
               ENDIF

! STEP 12: Update theta
               ! update theta
               IF (.NOT. line_search) THEN ! we prepared to perform the first line search
               ! "previous" refers to the previous CG step, not the previous LS step
                  CALL tidft_optimizer_mainvar_copy_from_env( &
                          optimizer_storage%prev_theta, &
                          tidft_env)
               ENDIF

               CALL tidft_optimizer_mainvar_add_to_tidft_env( &
                       tidft_env, &
                       optimizer_storage%step, &
                       1.0_dp, step_size)

            ENDIF ! not.prepare_to_exit

            IF (line_search) THEN
               iter_type = "LS"
            ELSE
               iter_type = "CG"
            ENDIF

            t2 = m_walltime()
            IF (unit_nr > 0) THEN
               iter_type = TRIM("TIDFT OPT "//iter_type)
               WRITE (unit_nr, '(T2,A13,I6,F23.10,E14.5,F14.9,F9.2)') &
                  iter_type, iteration, &
                  objf_new, objf_diff, grad_norm, &
                  t2 - t1
               !WRITE (unit_nr, '(T2,A19,F23.10)') &
               !   "Energy:", energy_function
               !WRITE (unit_nr, '(T2,A19,F23.10)') &
               !   "Orthogonalization:", penalty_func_new
            ENDIF
            t1 = m_walltime()

            iteration = iteration + 1
            IF (prepare_to_exit) EXIT

         ENDDO ! inner loop

         IF (tidft_env%state_by_state) THEN
            CALL tidft_mainvar_orthonormalize(tidft_env, &
                    only_state = tidft_env%nstates)
         ELSE
            CALL tidft_mainvar_orthonormalize(tidft_env)
         ENDIF

         IF (converged .OR. (outer_iteration .GE. outer_max_iter)) THEN
            outer_prepare_to_exit = .TRUE.
         ENDIF

         outer_iteration = outer_iteration + 1
         IF (outer_prepare_to_exit) EXIT

      ENDDO ! outer loop

      ! clean up
      CALL tidft_optimizer_deallocate_storage( &
              data_storage = optimizer_storage, &
              nstates = tidft_env%nstates, &
              ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)

      IF (.NOT. converged) THEN
         CPABORT("Optimization not converged! ")
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE tidft_optimize_pcg

! **************************************************************************************************
!> \brief Use main variable to compute MOs, DMs, KSs, energy and energy gradients
!> \param qs_env ...
!> \param tidft_env ...
!> \par History
!>       2022.05 created as a separate subroutine [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE mainvar_to_loss_and_gradient(qs_env, tidft_env, &
      loss_components_out, gradient_out, silent)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env
      REAL(KIND=dp), DIMENSION(100), INTENT(OUT)         :: loss_components_out
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: gradient_out
      !TYPE(dbcsr_type), DIMENSION(:,:,:), OPTIONAL, &
      !   INTENT(INOUT)                                   :: gradient_components_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: silent

      CHARACTER(len=*), PARAMETER :: routineN = 'mainvar_to_loss_and_gradient', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: handle, iaspc, ispin, istore, naspc, &
                                                            nspins, unit_nr, nmos_stateI, &
                                                            stateI, ivar, nstates, idet, &
                                                            ndets, slaterI, stateJ, &
                                                            nmos_stateJ, first_state
      LOGICAL                                            :: aspc_guess, my_silent
      REAL(KIND=dp)                                      :: alpha, cs_pos, energy, kTS_sum, &
                                                            random_percentile, &
                                                            energy_singlet, energy_triplet, &
                                                            state_energy, obj_energy, &
                                                            interstate_nonorthogonality_penalty, &
                                                            intrastate_nonorthogonality_penalty, &
                                                            DetDiagonalSigma, penalty_stateI, &
                                                            determinant_spd
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_type)                                   :: tempNO1
      TYPE(dbcsr_type)                                   :: temp2, mat_var23, mat_var23_sigma
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: mat_proj
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: fake_e_per_block
      INTEGER                                            :: group, nsize
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diagonal
      TYPE(qs_energy_type), POINTER                      :: qs_energy
      LOGICAL                                            :: energy_component, &
                                                            interstate_penalty_component, &
                                                            intrastate_penalty_component, &
                                                            not_evaluated

      ! GRADIENT
      TYPE(dbcsr_type)                                   :: tempNO2, tempOO1, &
                                                            tempOO2, ST, STsiginv
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:,:)      :: mat_g, &
                                                            mat_g_intra, &
                                                            mat_g_inter
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:,:,:)    :: mat_g_energy
      INTEGER                                            :: opposite_spin

      CALL timeset(routineN, handle)

      my_silent=.FALSE.
      IF (PRESENT(silent)) my_silent=silent

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! supress output if necessary
      IF (my_silent) unit_nr = -1

      ! to re-use some almo_scf subroutines set 
      ! nocc_of_domain to the negative number
      ALLOCATE(fake_e_per_block(1))
      fake_e_per_block(:) = -1

      nspins = tidft_env%nspins
      nstates = tidft_env%nstates
      ndets=1
      IF (tidft_env%triplet_correction) ndets = 2

      first_state = 1
      IF (tidft_env%state_by_state) first_state = nstates

      energy_component=tidft_env%loss_function_components(1)
      intrastate_penalty_component=tidft_env%loss_function_components(2)
      interstate_penalty_component=tidft_env%loss_function_components(3)

!      ! CINDY's code starts here
!      ! Combine X1-X3 to generate MOs mat_t(alpha) and mat_t(beta) for ROKS and UKS
!      DO stateI = 1, nstates
!         DO slaterI = ndets, 1, -1 ! do T before S to avoid overwrite in the X1 U X2 step
!            DO ispin = 1, nspins
!               ! Do we have paired e-? X1
!               IF (tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(1) .NE. 0) THEN
!                  WRITE(*,*) "State, determinant, sipn: ", stateI, slaterI, ispin
!                  IF (tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(ispin+1) .NE. 0) THEN
!                     ! combine paired and unpaired electrons
!                     IF (slaterI .EQ. 2)  THEN ! Triplet
!                        IF (ispin .EQ. 1) THEN ! T_alpha
!                              ! Use the singlet matrix to store intermediate X1 U X2 matrix
!                              CALL stitch_columns(tidft_env%e_states(stateI)%slater_det(1)%mat_var(1), & ! X1 U X2
!                                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(2), &
!                                                   tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin))
!                              CALL stitch_columns(tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin), & ! X1 U X2 U X3
!                                                   tidft_env%e_states(stateI)%slater_det(1)%mat_var(3), &
!                                                   tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin))
!                        ELSE ! T_beta = X1, copy mat_var(1) into mat_t
!                              CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin), &
!                                              tidft_env%e_states(stateI)%slater_det(1)%mat_var(1))
!                        ENDIF
!                     ELSE ! Singlet state, here can do for both beta and alpha
!                           CALL stitch_columns(tidft_env%e_states(stateI)%slater_det(1)%mat_var(1), & ! X1 U X2 and X1 U X3
!                                                tidft_env%e_states(stateI)%slater_det(1)%mat_var(ispin+1), &
!                                                tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin))
!                     ENDIF
!                  ELSE !copy mat_var(1) into mat_t
!                        CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin), &
!                                        tidft_env%e_states(stateI)%slater_det(1)%mat_var(1))
!                  ENDIF
!
!               ELSE ! No paired e-. Do we have unpaired e-? n_alpha and n_beta are different
!                  IF (tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(ispin+1) .EQ. 0) THEN
!                     IF (ispin .EQ. 1) CPABORT("Number of alpha electrons cannot be zero!")
!                  ELSE ! No paired electrons, only unpaired. This also takes care of UKS
!
!
!                     !RZK-WARNING: Should deal with the triplet state here!
!                     IF (slaterI .EQ. 2)  THEN ! Triplet
!                        IF (ispin .EQ. 1) THEN ! T_alpha
!                        ! Use the singlet matrix to store intermediate X1 U X2 matrix
!                           CALL stitch_columns(tidft_env%e_states(stateI)%slater_det(1)%mat_var(2), & ! X1 U X2
!                                             tidft_env%e_states(stateI)%slater_det(1)%mat_var(3), &
!                                             tidft_env%e_states(stateI)%slater_det(slaterI)%mat_t(ispin))
!                        ELSE ! no T_beta for the triplet state, do nothing
!                        ENDIF
!                     ELSE ! Singlet state, here can do for both beta and alpha
!                        CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(1)%mat_t(ispin), &
!                            tidft_env%e_states(stateI)%slater_det(1)%mat_var(ispin+1))
!                     ENDIF
!
!                     
!                  ENDIF
!               ENDIF
!            ENDDO ! end ispin
!         ENDDO ! end slaterI
!      ENDDO  ! end stateI
      
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

      !IF (unit_nr > 0) THEN
      !   WRITE (unit_nr, *) "NStates", nstates
      !   WRITE (unit_nr, *) "NSpins", nspins
      !   WRITE (unit_nr, *) "NDets", ndets
      !ENDIF

      ! This code converts independent variables into MO matrices
      ! Should work for unrestricted case
      ! Reverse order of determinants allows to use empty singlet determinants
      ! to temporarily store the triplet determinant
      DO idet = ndets, 1, -1
         DO stateI = first_state, nstates
            DO ispin = 1, nspins

               !IF (unit_nr > 0) THEN
               !   WRITE (unit_nr, *) "idet, istate, ispin", idet, stateI, ispin
               !   WRITE (unit_nr, *) "  doubly-, singly-occupied, total orbitals: ", &
               !      tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(1), &
               !      tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(ispin+1), &
               !      tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)
               !ENDIF

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
                        !WRITE(*,*) "STITCH 4: ispin, istate, idet: ", ispin, stateI, idet
                        !WRITE(*,*) "L,R,C: ", &
                        ! tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(1), &
                        ! tidft_env%e_states(stateI)%slater_det(idet)%ncols_var(ispin+1), &
                        ! tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)
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

      obj_energy = 0.0_dp
      DO stateI = first_state, nstates
         tidft_env%e_states(stateI)%state_energy = 0.0_dp
         DO idet = 1, ndets
            tidft_env%e_states(stateI)%slater_det(idet)%slater_energy = 0.0_dp
         ENDDO
      ENDDO

      IF (energy_component) THEN
      ! Use MOs to compute the DMs, then KS matrices and then the energy

         ALLOCATE(mat_proj(nspins))

         DO stateI = first_state, nstates
   
            DO idet = 1, ndets
   
               DO ispin = 1, nspins
   
                  CALL dbcsr_create(mat_proj(ispin), template=tidft_env%mat_s)
   
                  ! Create projectors only if there are electrons of this spin
                  IF (tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin).GT.0) THEN
   
   
                     ! compute projectors from electronic states
                     CALL almo_scf_t_to_proj(t=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                             p=mat_proj(ispin), &
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
                  CALL dbcsr_scale(mat_proj(1), 2.0_dp)
               ENDIF
   
               ! RZK-WARNING: What happens if there are no beta electrons? And thus no projector
               CALL tidft_dm_to_tidft_ks(qs_env, &
                                         mat_proj, &
                                         tidft_env%e_states(stateI)%slater_det(idet)%mat_ks, &
                                         energy, &
                                         tidft_env%eps_filter, &
                                         smear=tidft_env%smear, &
                                         kTS_sum=kTS_sum)
   
! perform population analysis
!CALL mulliken_population_analysis(qs_env, &
!        output_unit = unit_nr, &
!        print_level = 2)

               ! re-compute projector from the dm
               IF (nspins == 1) THEN
                  CALL dbcsr_scale(mat_proj(1), 0.5_dp)
               ENDIF
   
               !IF (unit_nr > 0 .AND. ndets.EQ.2) THEN
               !   WRITE (unit_nr, '(T2,A28,I10,F40.10)') "Determinant", idet, energy
               !ENDIF
   
               IF (idet.EQ.1) THEN
                  energy_singlet = energy
                  energy_triplet = energy
               ENDIF
               IF (idet.EQ.2) energy_triplet = energy
               tidft_env%e_states(stateI)%slater_det(idet)%slater_energy = energy

               DO ispin = 1, nspins
                  CALL dbcsr_release(mat_proj(ispin))
               ENDDO
   
            ENDDO ! idet
   
            state_energy = 2*energy_singlet-energy_triplet
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A28,I10,F40.10)') "State ", stateI, state_energy
            ENDIF
            tidft_env%e_states(stateI)%state_energy = state_energy
   
            obj_energy = obj_energy + state_energy
   
         ENDDO ! stateI
   
         !IF (unit_nr > 0) THEN
         !   WRITE (unit_nr, '(T2,A38,F40.10)') "Sum of energies of all states", obj_energy
         !   WRITE (unit_nr, '()')
         !ENDIF
   
         DEALLOCATE(mat_proj)
   
      ENDIF ! energy_component

      ! Return the sum of all energies to QS
      CALL get_qs_env(qs_env, energy=qs_energy)
      qs_energy%total = obj_energy

      ! Pass the result up
      loss_components_out(1) = obj_energy 

      DEALLOCATE(fake_e_per_block)

      ! During state-by-state optimization
      ! nullifying mat_phi_spin erases info on fixed states
      !tidft_env%mat_phi_spin(:,:,:) = 0.0_dp
      tidft_env%mat_phi(:,:) = 0.0_dp
      tidft_env%mat_phi_inv(:,:) = 0.0_dp

      intrastate_nonorthogonality_penalty = 0.0_dp
      IF (intrastate_penalty_component) THEN 
      !! compute intrastate nonorthogonality penalty

         DO stateI = first_state, nstates

            IF (tidft_env%triplet_correction) THEN
               idet = 2
            ELSE
               idet = 1
            ENDIF

            !DO idet = 1, ndets

            !   IF (tidft_env%triplet_correction .AND. idet.EQ.1) cycle ! only triplet determinant is used

            DO ispin = 1, nspins

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)

               IF ( nmos_stateI .GT. 0) THEN

                  ! compute sigma unless it is already computed with the energy
                  IF (.NOT.energy_component) THEN
                     CALL get_overlap(bra=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                      ket=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                      overlap=tidft_env%mat_sigma(stateI,stateI,idet,ispin), &
                                      metric=tidft_env%mat_s, &
                                      retain_overlap_sparsity=.FALSE., &
                                      eps_filter=tidft_env%eps_filter)
                  ENDIF
                  ! compute the determinant and save it to the Phi matrix
                  CALL cp_dbcsr_determinant_generalized( &
                          matrix = tidft_env%mat_sigma(stateI,stateI,idet,ispin), &
                          det = tidft_env%mat_phi_spin(stateI,stateI,ispin), &
                          para_env = tidft_env%para_env, &
                          blacs_env = tidft_env%blacs_env)
                  !IF (unit_nr > 0) THEN
                  !   WRITE (unit_nr, *) "StateI, Slater, Spin: ", stateI, idet, ispin, &
                  !                      "Determinant: ", tidft_env%mat_phi_spin(stateI,stateI,ispin)
                  !ENDIF

                  !IF (debug_mode) THEN
                  !   CALL determinant(tidft_env%mat_sigma(stateI,stateI,idet,ispin), &
                  !           determinant_spd, &
                  !           tidft_env%eps_filter)
                  !   IF (unit_nr > 0) THEN
                  !      WRITE (unit_nr, *) "StateI, Slater, Spin: ", stateI, idet, ispin, &
                  !                         "SPD determi: ", determinant_spd
                  !   ENDIF
                  !ENDIF

                  ! compute the product of the diagonal elements
                  CALL dbcsr_get_info(tidft_env%mat_sigma(stateI,stateI,idet,ispin), &
                          nfullrows_total=nsize, group=group)
                  ALLOCATE (diagonal(nsize))
                  CALL dbcsr_get_diag(tidft_env%mat_sigma(stateI,stateI,idet,ispin), diagonal)
                  CALL mp_sum(diagonal, group)
                  DetDiagonalSigma = PRODUCT(diagonal)
                  DEALLOCATE(diagonal)

                  IF (unit_nr > 0) THEN
                     !WRITE (unit_nr, *) "StateI, Slater, Spin: ", stateI, idet, ispin, &
                     !                   "Det. diagon: ", DetDiagonalSigma
                     WRITE (unit_nr, '(T2,A23,I5,2I2,A25,F21.16)') &
                        "State, Slater, Spin: ", stateI, idet, ispin, &
                        "Normalized determinant: ", &
                        tidft_env%mat_phi_spin(stateI,stateI,ispin) / DetDiagonalSigma
                  ENDIF

                  penalty_stateI = - tidft_env%intrastate_nonorthogonality_penalty_coeff * &
                     ( LOG(tidft_env%mat_phi_spin(stateI,stateI,ispin)) - LOG(DetDiagonalSigma) )

                  !IF (unit_nr > 0) THEN
                  !   WRITE (unit_nr, '(T2,A28,I10,F40.10)') "State penalty ", stateI, penalty_stateI
                  !ENDIF

                  intrastate_nonorthogonality_penalty = intrastate_nonorthogonality_penalty + &
                     penalty_stateI

               ENDIF ! zero MOs for state I

            ENDDO ! ispin

            !ENDDO ! idet

         ENDDO ! stateI

      ENDIF ! intrastate_penalty_component

      !IF (unit_nr > 0) THEN
      !   WRITE (unit_nr, '(T2,A38,F40.10)') "Intrastate nonorthogonality penalty", &
      !      intrastate_nonorthogonality_penalty
      !   WRITE (unit_nr, '()')
      !ENDIF

      ! Pass the result up
      loss_components_out(2) = intrastate_nonorthogonality_penalty

      interstate_nonorthogonality_penalty = 0.0_dp
      IF (interstate_penalty_component) THEN 
      !! compute interstate nonorthogonality penalty

         DO stateI = 1, nstates
            DO idet = 1, ndets

               IF (tidft_env%triplet_correction .AND. idet.EQ.1) cycle ! only triplet determinant is used

               DO ispin = 1, nspins

                  nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)
                  IF ( nmos_stateI .GT. 0) THEN

                     DO stateJ = stateI, nstates
                        nmos_stateJ = tidft_env%e_states(stateJ)%slater_det(idet)%nelectrons_ab(ispin)

                        IF ( nmos_stateJ .GT. 0) THEN

                           ! skip the (I,I) overlap if it has been computed with energy or intrastate
                           not_evaluated = .NOT. (energy_component .OR. intrastate_penalty_component)
                           IF ( stateI.NE.stateJ .OR. not_evaluated ) THEN
                              CALL get_overlap(bra=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                               ket=tidft_env%e_states(stateJ)%slater_det(idet)%mat_t(ispin), &
                                               overlap=tidft_env%mat_sigma(stateI,stateJ,idet,ispin), &
                                               metric=tidft_env%mat_s, &
                                               retain_overlap_sparsity=.FALSE., &
                                               eps_filter=tidft_env%eps_filter)
                              !CALL dbcsr_print(tidft_env%mat_sigma(stateI,stateJ,idet,ispin))
                           ENDIF
                           ! compute the determinant and save it to the Phi matrix
                           not_evaluated = .NOT. intrastate_penalty_component
                           IF ( stateI.NE.stateJ .OR. not_evaluated ) THEN
                              CALL cp_dbcsr_determinant_generalized( &
                                      matrix = tidft_env%mat_sigma(stateI,stateJ,idet,ispin), &
                                      det = tidft_env%mat_phi_spin(stateI,stateJ,ispin), &
                                      para_env = tidft_env%para_env, &
                                      blacs_env = tidft_env%blacs_env)
                              tidft_env%mat_phi_spin(stateJ,stateI,ispin)=tidft_env%mat_phi_spin(stateI,stateJ,ispin)
                           ENDIF
                           !IF (unit_nr > 0) THEN
                           !   WRITE (unit_nr, *) "StateI, StateJ, Slater, Spin: ", stateI, stateJ, idet, ispin, &
                           !                      "Determinant: ", tidft_env%mat_phi_spin(stateI,stateJ,ispin)
                           !ENDIF
   
                        ENDIF ! zero MOs for state J

                     ENDDO ! stateJ

                  ENDIF ! zero MOs for state I

               ENDDO ! ispin

            ENDDO ! idet
         ENDDO ! stateI

         ! Print phi matrices
         !DO ispin = 1, nspins
         !   WRITE (*,*) "Phi matrix, spin ", ispin
         !   DO ivar = 1, nstates
         !      WRITE (*,*) tidft_env%mat_phi_spin(ivar, 1:nstates, ispin)
         !   ENDDO
         !   WRITE (*,*)
         !ENDDO
   
         ! Hadamard product of spin matrices to obtain full phi
         IF (nspins .EQ. 2) THEN
            tidft_env%mat_phi(:,:) = tidft_env%mat_phi_spin(:,:,1) * tidft_env%mat_phi_spin(:,:,2)
         ELSE
            tidft_env%mat_phi(:,:) = tidft_env%mat_phi_spin(:,:,1) * tidft_env%mat_phi_spin(:,:,1)
         ENDIF
   
         ! Print phi matrix
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, *) "Phi matrix"
            DO ivar = 1, nstates
               WRITE (unit_nr,'(T2,15ES9.1E2)') tidft_env%mat_phi(ivar, 1:nstates)
            ENDDO
         ENDIF

   
         ! compute the product of the diagonal elements
         DetDiagonalSigma = 1.0_dp
         DO ivar = 1, nstates
            DetDiagonalSigma = DetDiagonalSigma * tidft_env%mat_phi(ivar,ivar)
         ENDDO
         CALL determinant_local(tidft_env%mat_phi, determinant_spd, nstates)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A38,F40.10)') "Determinant of Diagonal(PHI): ", &
               DetDiagonalSigma
         ENDIF

         tidft_env%det_mat_phi = determinant_spd/DetDiagonalSigma

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A38,F40.10)') "Normalized determinant of PHI: ", &
               tidft_env%det_mat_phi
         ENDIF
   
         interstate_nonorthogonality_penalty = &
            - tidft_env%interstate_nonorthogonality_penalty_coeff * &
            ( LOG(tidft_env%det_mat_phi) )


      ENDIF ! interstate_penalty_component

      ! Pass the result up
      loss_components_out(3) = interstate_nonorthogonality_penalty

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A38,F40.10)') "Sum of energies of all states", obj_energy
         WRITE (unit_nr, '(T2,A38,F40.10)') "Intrastate nonorthogonality penalty", &
            intrastate_nonorthogonality_penalty
         WRITE (unit_nr, '(T2,A38,F40.10)') "Interstate nonorthogonality penalty", &
            interstate_nonorthogonality_penalty
         WRITE (unit_nr, '(T2,A38,F40.10)') "Objective function", &
            obj_energy + &
            intrastate_nonorthogonality_penalty + &
            interstate_nonorthogonality_penalty
         WRITE (unit_nr, '()')
      ENDIF

      ! GRADIENT CALCULATIONS

      ! Nullify the gradient before accumulation
      DO stateI = 1, nstates
         DO ivar = 1, 3
            IF (tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar).NE.0) THEN
               CALL dbcsr_set(gradient_out(stateI, ivar), 0.0_dp)
               !IF (gradient_components) THEN
               !   CALL dbcsr_set(gradient_components_out(stateI, ivar, 1), 0.0_dp)
               !   CALL dbcsr_set(gradient_components_out(stateI, ivar, 2), 0.0_dp)
               !   CALL dbcsr_set(gradient_components_out(stateI, ivar, 3), 0.0_dp)
               !ENDIF
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

      ! RZK-WARNING: not necessary to store T-gradients
      ! REmove after debugging
      ALLOCATE(mat_g_energy(nstates, nspins, ndets))
      ALLOCATE(mat_g_intra(nstates, nspins))
      ALLOCATE(mat_g_inter(nstates, nspins))

      DO ispin = 1, nspins
         DO stateI = first_state, nstates

            DO idet = 1, ndets
               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)
               IF ( nmos_stateI .GT. 0) THEN
                  CALL dbcsr_create(mat_g_energy(stateI, ispin, idet), &
                          template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
                  CALL dbcsr_set(mat_g_energy(stateI, ispin, idet), 0.0_dp)
               ENDIF ! zero MOs
            ENDDO ! idet

            IF (tidft_env%triplet_correction) THEN
               idet = 2
            ELSE
               idet = 1
            ENDIF
            CALL dbcsr_create(mat_g_intra(stateI, ispin), &
                    template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
            CALL dbcsr_create(mat_g_inter(stateI, ispin), &
                    template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
            CALL dbcsr_set(mat_g_intra(stateI, ispin), 0.0_dp)
            CALL dbcsr_set(mat_g_inter(stateI, ispin), 0.0_dp)

         ENDDO ! stateI
      ENDDO ! ispin

      IF (interstate_penalty_component) THEN
         tidft_env%mat_phi_inv(:,:)=tidft_env%mat_phi(:,:)
         CALL matrix_invert_local(tidft_env%mat_phi_inv, nstates)
         ! Print inverted phi
         !WRITE (*,*) "Phi inv matrix"
         !DO ivar = 1, nstates
         !   WRITE (*,*) tidft_env%mat_phi_inv(ivar, 1:nstates)
         !ENDDO
         !WRITE (*,*)
      ENDIF

      DO stateI = first_state, nstates
         DO idet = 1, ndets
            DO ispin = 1, nspins

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)

               IF ( nmos_stateI .GT. 0) THEN

                  !IF (unit_nr > 0) THEN
                  !   WRITE(unit_nr, '(T2,A28,3I4,A15,I4,A25,3I4)') &
                  !      "Convert (state,spin,det): ", stateI, ispin, idet, &
                  !      " containing ", &
                  !      nmos_stateI, "MOs into VARs containing", &
                  !      tidft_env%e_states(stateI)%slater_det(1)%ncols_var
                  !ENDIF

                  CALL dbcsr_create(tempNO1, &
                          template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
                  CALL dbcsr_create(tempOO1, &
                          template=tidft_env%mat_sigma_inv(stateI,idet,ispin))

                  CALL dbcsr_create(ST, &
                          template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
                  CALL dbcsr_create(STsiginv, &
                          template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))

                  ! evaluate STsiginv only for the penalty term, which means
                  ! Do NOT compute for triplet_correction and the first determinant
                  ! This means compute for triplet_correction==.FALSE. .OR. idet.EQ.2
                  IF (energy_component .OR. (.NOT.tidft_env%triplet_correction) .OR. idet.EQ.2 ) THEN

                     IF (.NOT.energy_component) THEN
                     !invert using cholesky for dense matrices
                        CALL dbcsr_copy(tidft_env%mat_sigma_inv(stateI,idet,ispin), &
                                tidft_env%mat_sigma(stateI,stateI,idet,ispin))
                        CALL cp_dbcsr_cholesky_decompose(tidft_env%mat_sigma_inv(stateI,idet,ispin), &
                                para_env=tidft_env%para_env, &
                                blacs_env=tidft_env%blacs_env)
                        CALL cp_dbcsr_cholesky_invert(tidft_env%mat_sigma_inv(stateI,idet,ispin), &
                                para_env=tidft_env%para_env, &
                                blacs_env=tidft_env%blacs_env, &
                                upper_to_full=.TRUE.)
                        CALL dbcsr_filter(tidft_env%mat_sigma_inv(stateI,idet,ispin), &
                                tidft_env%eps_filter)
                     ENDIF

                     ! T.siginv
                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                             tidft_env%mat_sigma_inv(stateI,idet,ispin), &
                             0.0_dp, tempNO1, &
                             filter_eps=tidft_env%eps_filter)

                     ! S.Tsiginv
                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             tidft_env%mat_s, &
                             tempNO1, &
                             0.0_dp, STsiginv, &
                             filter_eps=tidft_env%eps_filter)

                  ENDIF

                  ! Compute the energy component
                  IF (energy_component) THEN

                     ! F.Tsiginv
                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             tidft_env%e_states(stateI)%slater_det(idet)%mat_ks(ispin), &
                             tempNO1, &
                             0.0_dp, mat_g_energy(stateI, ispin, idet), &
                             filter_eps=tidft_env%eps_filter)

                     ! Tt.FTsiginv
                     CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                             mat_g_energy(stateI, ispin, idet), &
                             0.0_dp, tempOO1, &
                             filter_eps=tidft_env%eps_filter)

                     ! FTsiginv - STsiginv.TtFTsiginv
                     CALL dbcsr_multiply("N", "N", -1.0_dp, &
                             STsiginv, &
                             tempOO1, &
                             1.0_dp, mat_g_energy(stateI, ispin, idet), &
                             filter_eps=tidft_env%eps_filter)
                     ! Scale by 2.0
                     CALL dbcsr_scale(mat_g_energy(stateI, ispin, idet), 2.0_dp)
               


                     CALL gradMO_to_gradVAR_increment( &
                             mat_g_out=gradient_out(stateI, :), &
                             ncols_var=tidft_env%e_states(stateI)%slater_det(1)%ncols_var, &
                             mat_g_in=mat_g_energy(stateI, ispin, idet), &
                             mat_s=tidft_env%mat_s, &
                             idet=idet, &
                             ispin=ispin, &
                             nspins=nspins, &
                             triplet_correction=tidft_env%triplet_correction, &
                             energy_term=.TRUE.)


                     !! RZK-WARNING: IF components are given, do not accumulate the total gradient
                     !! It's a waste of time
                     !CALL gradMO_to_gradVAR_increment( &
                     !        mat_g_out=gradient_components_out(stateI, :, 1), &
                     !        ncols_var=tidft_env%e_states(stateI)%slater_det(1)%ncols_var, &
                     !        mat_g_in=mat_g_energy(stateI, ispin, idet), &
                     !        mat_s=tidft_env%mat_s, &
                     !        idet=idet, &
                     !        ispin=ispin, &
                     !        nspins=nspins, &
                     !        triplet_correction=tidft_env%triplet_correction, &
                     !        energy_term=.TRUE.)

                  ENDIF ! energy_component

                  IF (intrastate_penalty_component) THEN

                     ! Do NOT compute for triplet_correction and the first determinant
                     ! This means compute for triplet_correction==.FALSE. .OR. idet.EQ.2
                     IF ( .NOT.(tidft_env%triplet_correction .AND. idet.EQ.1) ) THEN

                        ! Intrastate term
                        ! G = STsiginv
                        CALL dbcsr_copy(mat_g_intra(stateI, ispin), STsiginv)
                        ! d^(-1)
                        CALL dbcsr_get_info(tidft_env%mat_sigma(stateI, stateI, idet, ispin), &
                                nfullrows_total=nsize, group=group)
                        ALLOCATE (diagonal(nsize))
                        CALL dbcsr_get_diag(tidft_env%mat_sigma(stateI, stateI, idet, ispin), diagonal)
                        diagonal(:) = 1.0_dp/diagonal(:) ! invert the diagonal
                        CALL dbcsr_set(tempOO1, 0.0_dp)
                        CALL dbcsr_set_diag(tempOO1, diagonal)
                        CALL dbcsr_filter(tempOO1, tidft_env%eps_filter)
                        DEALLOCATE(diagonal)
                        ! S.T
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%mat_s, &
                                tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                0.0_dp, ST, &
                                filter_eps=tidft_env%eps_filter)
                        ! G = G - ST.d^(-1)
                        CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                ST, &
                                tempOO1, &
                                1.0_dp, mat_g_intra(stateI, ispin), &
                                filter_eps=tidft_env%eps_filter)
                        CALL dbcsr_scale(mat_g_intra(stateI, ispin), & 
                                -2.0_dp*tidft_env%intrastate_nonorthogonality_penalty_coeff)

                        CALL gradMO_to_gradVAR_increment( &
                                mat_g_out=gradient_out(stateI, :), &
                                ncols_var=tidft_env%e_states(stateI)%slater_det(1)%ncols_var, &
                                mat_g_in=mat_g_intra(stateI, ispin), &
                                mat_s=tidft_env%mat_s, &
                                idet=idet, &
                                ispin=ispin, &
                                nspins=nspins, &
                                triplet_correction=tidft_env%triplet_correction, &
                                energy_term=.FALSE.)

                     ENDIF ! only triplet determinant is used in the penalty

                  ENDIF ! intrastate_penalty_component

                  ! Interstate term
                  IF (interstate_penalty_component) THEN

                        DO stateJ = 1, nstates

                           nmos_stateJ = tidft_env%e_states(stateJ)%slater_det(idet)%nelectrons_ab(ispin)

                           IF ( nmos_stateJ .GT. 0) THEN

                              ! Compute ADJUGATE
                              IF (stateJ .GE. stateI) THEN
                                 CALL dbcsr_create(tempOO2, template=tidft_env%mat_sigma(stateI,stateJ,idet,ispin))
                                 CALL dbcsr_copy(tempOO2, tidft_env%mat_sigma(stateI,stateJ,idet,ispin))
                              ELSE
                                 CALL dbcsr_create(tempOO2, template=tidft_env%mat_sigma(stateJ,stateI,idet,ispin))
                                 CALL dbcsr_transposed(tempOO2, tidft_env%mat_sigma(stateJ,stateI,idet,ispin))
                              ENDIF

                              ! RZK-WARNING time can be saved by computing ADJ only for half of the matrices
                              CALL cp_dbcsr_adjugate(matrix=tempOO2, &
                                      det = determinant_spd, &
                                      para_env = tidft_env%para_env, &
                                      blacs_env = tidft_env%blacs_env)
                              !CALL dbcsr_create(temp2, template=tidft_env%mat_sigma(stateI,stateJ,idet,ispin))
                              !CALL dbcsr_multiply("N", "N", 1.0_dp, tempOO2, &
                              !        tidft_env%mat_sigma(stateI,stateJ,idet,ispin), &
                              !        0.0_dp, temp2, &
                              !        filter_eps=tidft_env%eps_filter)
                              !CALL dbcsr_scale(temp2, 1.0_dp/determinant_spd)
                              !CALL dbcsr_print(temp2)
                              !CALL dbcsr_release(temp2)

                              ! save ADJ matrix
                              IF (stateJ .GE. stateI) THEN
                                 CALL dbcsr_copy(tidft_env%mat_adj_sigma(stateI,stateJ,idet,ispin), tempOO2)
                              ELSE
                                 CALL dbcsr_transposed(tidft_env%mat_adj_sigma(stateJ,stateI,idet,ispin), tempOO2)
                              ENDIF

                              ! This block saves ADJ and reduces the number of ADJ evaluations
                              ! RZK-WARNING: It does not work for some reason
                              !IF (stateJ .GE. stateI) THEN
                              !   CALL dbcsr_create(tempOO2, template=tidft_env%mat_sigma(stateI,stateJ,idet,ispin))
                              !   CALL dbcsr_copy(tidft_env%mat_adj_sigma(stateI,stateJ,idet,ispin), &
                              !           tidft_env%mat_sigma(stateI,stateJ,idet,ispin))
                              !   CALL cp_dbcsr_adjugate( &
                              !           matrix=tidft_env%mat_adj_sigma(stateI,stateJ,idet,ispin), &
                              !           det = determinant_spd, &
                              !           para_env = tidft_env%para_env, &
                              !           blacs_env = tidft_env%blacs_env)
                              !   CALL dbcsr_copy(tempOO2, tidft_env%mat_adj_sigma(stateI,stateJ,idet,ispin))
                              !   !CALL dbcsr_create(temp2, template=tidft_env%mat_sigma(stateI,stateJ,idet,ispin))
                              !   !CALL dbcsr_multiply("N", "N", 1.0_dp, &
                              !   !        tidft_env%mat_adj_sigma(stateI,stateJ,idet,ispin), &
                              !   !        tidft_env%mat_sigma(stateI,stateJ,idet,ispin), &
                              !   !        0.0_dp, temp2, &
                              !   !        filter_eps=tidft_env%eps_filter)
                              !   !CALL dbcsr_scale(temp2, 1.0_dp/determinant_spd)
                              !   !CALL dbcsr_print(temp2)
                              !   !CALL dbcsr_release(temp2)
                              !ELSE
                              !   CALL dbcsr_create(tempOO2, template=tidft_env%mat_adj_sigma(stateJ,stateI,idet,ispin))
                              !   CALL dbcsr_transposed(tempOO2, tidft_env%mat_adj_sigma(stateJ,stateI,idet,ispin))
                              !ENDIF

                              ! Do NOT compute for triplet_correction and the first determinant
                              ! This means compute for triplet_correction==.FALSE. .OR. idet.EQ.2
                              IF ( .NOT.(tidft_env%triplet_correction .AND. idet.EQ.1) ) THEN

                                 IF (ispin .EQ. 1) THEN
                                    opposite_spin = 2
                                 ELSE
                                    opposite_spin = 1
                                 ENDIF
                                 CALL dbcsr_scale(tempOO2, &
                                         tidft_env%mat_phi_spin(stateJ,stateI,opposite_spin) * &
                                         tidft_env%mat_phi_inv(stateJ,stateI) )

                                 ! S.T_Jindex
                                 ! RZK-WARNING: optimize N^2 re-evaluations to N
                                 CALL dbcsr_create(tempNO2, &
                                         template=tidft_env%e_states(stateJ)%slater_det(idet)%mat_t(ispin))
                                 not_evaluated = .NOT. intrastate_penalty_component
                                 IF (stateI.NE.stateJ .OR. not_evaluated) THEN
                                    CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                            tidft_env%mat_s, &
                                            tidft_env%e_states(stateJ)%slater_det(idet)%mat_t(ispin), &
                                            0.0_dp, tempNO2, &
                                            filter_eps=tidft_env%eps_filter)
                                 ELSE
                                    CALL dbcsr_copy(tempNO2,ST)
                                 ENDIF

                                 CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                         tempNO2, &
                                         tempOO2, &
                                         1.0_dp, mat_g_inter(stateI, ispin), &
                                         filter_eps=tidft_env%eps_filter)
                                 CALL dbcsr_release(tempNO2)

                              ENDIF ! only triplet determinant is used in the penalty

                              CALL dbcsr_release(tempOO2)

                           ENDIF ! zero MOs for state J

                        ENDDO ! stateJ loop

                     ! Do NOT compute for triplet_correction and the first determinant
                     ! This means compute for triplet_correction==.FALSE. .OR. idet.EQ.2
                     IF ( .NOT.(tidft_env%triplet_correction .AND. idet.EQ.1) ) THEN

                        CALL dbcsr_add(mat_g_inter(stateI, ispin), STsiginv, 1.0_dp, -1.0_dp)
                        CALL dbcsr_scale(mat_g_inter(stateI, ispin), & 
                                -2.0_dp*tidft_env%interstate_nonorthogonality_penalty_coeff)

                        CALL gradMO_to_gradVAR_increment( &
                                mat_g_out=gradient_out(stateI, :), &
                                ncols_var=tidft_env%e_states(stateI)%slater_det(1)%ncols_var, &
                                mat_g_in=mat_g_inter(stateI, ispin), &
                                mat_s=tidft_env%mat_s, &
                                idet=idet, &
                                ispin=ispin, &
                                nspins=nspins, &
                                triplet_correction=tidft_env%triplet_correction, &
                                energy_term=.FALSE.)

                     ENDIF ! only triplet determinant is used in the penalty

                  ENDIF ! interstate_penalty_component

                  CALL dbcsr_release(tempNO1)
                  CALL dbcsr_release(tempOO1)
                  CALL dbcsr_release(ST)
                  CALL dbcsr_release(STsiginv)

               ENDIF ! zero MOs

            ENDDO ! ispin
         ENDDO ! idet
      ENDDO ! stateI

      ! Release derivatives wrt T
      DO ispin = 1, nspins
         DO stateI = first_state, nstates
            DO idet = 1, ndets
               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)
               IF ( nmos_stateI .GT. 0) THEN
                  !CALL dbcsr_print(mat_g_energy(stateI, ispin, idet))
                  CALL dbcsr_release(mat_g_energy(stateI, ispin, idet))
               ENDIF ! zero MOs
            ENDDO ! idet
            !CALL dbcsr_print(mat_g_intra(stateI, ispin))
            CALL dbcsr_release(mat_g_intra(stateI, ispin))
            !CALL dbcsr_print(mat_g_inter(stateI, ispin))
            CALL dbcsr_release(mat_g_inter(stateI, ispin))
         ENDDO ! stateI
      ENDDO ! ispin
      DEALLOCATE(mat_g_energy)
      DEALLOCATE(mat_g_intra)
      DEALLOCATE(mat_g_inter)

      CALL timestop(handle)

   END SUBROUTINE mainvar_to_loss_and_gradient

! **************************************************************************************************
!> \brief Copy matrices describing main variables from tidft_env into storage: 
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_numerical_gradient(qs_env, tidft_env, dx, &
      loss_components_out, gradient_out)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env
      REAL(KIND=dp)                                      :: dx
      REAL(KIND=dp), DIMENSION(100), INTENT(OUT)         :: loss_components_out
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: gradient_out

      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: mat_x0, mat_dx, &
                                                            mat_current_x, &
                                                            grad
      REAL(KIND=dp), DIMENSION(100)                      :: loss
      REAL(KIND=dp)                                      :: y0, current_y_minus, &
                                                            current_y_plus, &
                                                            current_grad_plus, &
                                                            current_grad_minus, &
                                                            current_grad, &
                                                            sanity_check
      INTEGER                                            :: istate, ivar, nmos, &
                                                            nstates, nblkrows, &
                                                            nblkcols, iblkcol, &
                                                            iblkrow, irow, icol, &
                                                            plusminus, unit_nr, &
                                                            first_state
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes, &
                                                            col_blk_sizes
      TYPE(cp_logger_type), POINTER                      :: logger

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nstates = tidft_env%nstates

      first_state = 1
      IF (tidft_env%state_by_state) first_state = nstates

      ! set variables
      ALLOCATE(mat_x0(nstates,3))
      ALLOCATE(mat_dx(nstates,3))
      ALLOCATE(mat_current_x(nstates,3))
      ALLOCATE(grad(nstates,3))
      DO istate = 1, nstates
         DO ivar = 1, 3
            IF (tidft_env%e_states(istate)%slater_det(1)%ncols_var(ivar).NE.0) THEN
               CALL dbcsr_create(mat_x0(istate, ivar), &
                       template=tidft_env%e_states(istate)%slater_det(1)%mat_var(ivar))
               CALL dbcsr_create(mat_dx(istate, ivar), &
                       template=tidft_env%e_states(istate)%slater_det(1)%mat_var(ivar))
               CALL dbcsr_create(mat_current_x(istate, ivar), &
                       template=tidft_env%e_states(istate)%slater_det(1)%mat_var(ivar))
               CALL dbcsr_create(grad(istate, ivar), &
                       template=tidft_env%e_states(istate)%slater_det(1)%mat_var(ivar))
               CALL dbcsr_set(mat_x0(istate, ivar), 0.0_dp)
               CALL dbcsr_set(mat_dx(istate, ivar), 0.0_dp)
               CALL dbcsr_set(mat_current_x(istate, ivar), 0.0_dp)
               CALL dbcsr_set(grad(istate, ivar), 0.0_dp)
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

      ! update theta
      CALL tidft_optimizer_mainvar_copy_from_env( &
              mat_x0, &
              tidft_env)

      loss_components_out(:) = 0.0_dp
      CALL mainvar_to_loss_and_gradient( &
              qs_env = qs_env, &
              tidft_env = tidft_env, &
              loss_components_out = loss_components_out, &
              gradient_out = grad )
      y0 = SUM(loss_components_out)

WRITE(*,*)
WRITE(*,*) "ANALYTICAL GRADIENT"
WRITE(*,*)
      CALL tidft_optimizer_mainvar_print(grad, &
              tidft_env%e_states(1)%slater_det(1)%ncols_var)

      DO istate = first_state, nstates
         DO ivar = 1, 3
            nmos = tidft_env%e_states(istate)%slater_det(1)%ncols_var(ivar)
            IF (nmos.GT.0) THEN

               CALL dbcsr_get_info(mat_x0(istate,ivar), &
                       row_blk_size=row_blk_sizes, &
                       col_blk_size=col_blk_sizes, &
                       nblkrows_total=nblkrows, &
                       nblkcols_total=nblkcols)

               DO iblkrow = 1, nblkrows
                  DO irow = 1, row_blk_sizes(iblkrow)
                     DO iblkcol = 1, nblkcols
                        DO icol = 1, col_blk_sizes(iblkcol)
   
                           CALL tidft_optimizer_mainvar_set( &
                                   mat_dx, 0.0_dp, &
                                   tidft_env%e_states(1)%slater_det(1)%ncols_var)

                           CALL dbcsr_set_element(mat_dx(istate,ivar), dx, &
                                   iblkrow, irow, iblkcol, icol)
                           !WRITE(*,*) "CURRENT dx: ", istate, ivar, iblkrow, irow, iblkcol, icol
                           !CALL dbcsr_print(mat_dx(istate,ivar))

                           DO plusminus = -1, 1, 2
   
                              CALL tidft_optimizer_mainvar_copy( &
                                      mat_current_x, &
                                      mat_x0, &
                                      ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)
   
                              CALL tidft_optimizer_mainvar_add( &
                                      mat_current_x, &
                                      mat_dx, &
                                      1.0_dp, 1.0_dp*plusminus, &
                                      ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)
   
                              CALL tidft_optimizer_mainvar_copy_to_env( &
                                      tidft_env, &
                                      mat_current_x)
   
                              loss(:) = 0.0_dp
                              CALL mainvar_to_loss_and_gradient( &
                                      qs_env = qs_env, &
                                      tidft_env = tidft_env, &
                                      loss_components_out = loss, &
                                      gradient_out = grad, &
                                      silent = .TRUE.)

                              IF (plusminus.EQ.-1) THEN
                                 current_y_minus = SUM(loss)
                              ELSE
                                 current_y_plus = SUM(loss)
                              ENDIF

                           ENDDO !plusminus

                           current_grad_plus = (current_y_plus-y0)/dx
                           current_grad_minus = (current_y_minus-y0)/(-dx)
                           current_grad = 0.5_dp * ( current_grad_plus + current_grad_minus )

                           sanity_check = 100.0_dp*ABS(current_grad_plus-current_grad_minus)/current_grad
                           IF (sanity_check .GT. 5.0_dp) THEN
                              IF (unit_nr > 0) THEN
                                 WRITE(unit_nr, '(T2,A38,E12.1,F27.1,A1)') &
                                    "Numerical gradient is not accurate: ", current_grad, sanity_check, "%"
                              ENDIF
                              CPWARN("Plus-minus gradient values disagree! Reduce dx.")
                           ENDIF


                           CALL dbcsr_set_element(gradient_out(istate,ivar), current_grad, &
                                   iblkrow, irow, iblkcol, icol)

                        ENDDO ! icol
                     ENDDO ! iblkcol
                  ENDDO ! irow
               ENDDO ! iblkrow

            ENDIF

         ENDDO ! ivar
      ENDDO ! stateI

WRITE(*,*)
WRITE(*,*) "NUMERICAL GRADIENT"
WRITE(*,*)
      CALL tidft_optimizer_mainvar_print(gradient_out, &
              tidft_env%e_states(1)%slater_det(1)%ncols_var)

      ! set variables
      DO istate = 1, nstates
         DO ivar = 1, 3
            IF (tidft_env%e_states(istate)%slater_det(1)%ncols_var(ivar).NE.0) THEN
               CALL dbcsr_release(mat_x0(istate, ivar))
               CALL dbcsr_release(mat_dx(istate, ivar))
               CALL dbcsr_release(mat_current_x(istate, ivar))
               CALL dbcsr_release(grad(istate, ivar))
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI
      DEALLOCATE(mat_x0)
      DEALLOCATE(mat_dx)
      DEALLOCATE(mat_current_x)
      DEALLOCATE(grad)

   END SUBROUTINE tidft_optimizer_numerical_gradient

! **************************************************************************************************
!> \brief Copy matrices describing main variables from storage into tidft_env
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_copy_to_env(tidft_env, mat_source)

      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env 
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: mat_source

      INTEGER                                            :: istate, ivar, nmos, &
                                                            nstates, nvars, &
                                                            nstates2, nvars2

      nstates = SIZE(mat_source, 1)
      nvars = SIZE(mat_source, 2)
      CPASSERT(nvars.EQ.3)

      nstates2 = tidft_env%nstates
      CPASSERT(nstates.EQ.nstates2)

      DO istate = 1, nstates
         DO ivar = 1, nvars
            nmos = tidft_env%e_states(istate)%slater_det(1)%ncols_var(ivar)
            IF (nmos.GT.0) THEN
               CALL dbcsr_copy( &
                       tidft_env%e_states(istate)%slater_det(1)%mat_var(ivar), &
                       mat_source(istate,ivar))
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

   END SUBROUTINE tidft_optimizer_mainvar_copy_to_env

! **************************************************************************************************
!> \brief Copy matrices describing main variables from tidft_env into storage: 
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_copy_from_env(mat_dest, tidft_env)

      TYPE(tidft_env_type), INTENT(IN)                   :: tidft_env 
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: mat_dest

      INTEGER                                            :: istate, ivar, nmos, &
                                                            nstates, nvars, &
                                                            nstates2, nvars2

      nstates = SIZE(mat_dest, 1)
      nvars = SIZE(mat_dest, 2)
      CPASSERT(nvars.EQ.3)

      nstates2 = tidft_env%nstates
      CPASSERT(nstates.EQ.nstates2)

      DO istate = 1, nstates
         DO ivar = 1, nvars
            nmos = tidft_env%e_states(istate)%slater_det(1)%ncols_var(ivar)
            IF (nmos.GT.0) THEN
               CALL dbcsr_copy( &
                       mat_dest(istate,ivar), &
                       tidft_env%e_states(istate)%slater_det(1)%mat_var(ivar))
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

   END SUBROUTINE tidft_optimizer_mainvar_copy_from_env

! **************************************************************************************************
!> \brief Add matrices describing main variables: 
!>        mat_env = alpha * mat_env + beta * mat_source
!> \param ncols_var ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_mainvar_add_to_tidft_env(tidft_env, &
      mat_source, alpha, beta)

      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env 
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: mat_source
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta

      INTEGER                                            :: istate, ivar, nmos, &
                                                            nstates, nvars, &
                                                            nstates2, nvars2

      nstates = SIZE(mat_source, 1)
      nvars = SIZE(mat_source, 2)
      CPASSERT(nvars.EQ.3)

      nstates2 = tidft_env%nstates
      CPASSERT(nstates.EQ.nstates2)

      DO istate = 1, nstates
         DO ivar = 1, nvars
            nmos = tidft_env%e_states(istate)%slater_det(1)%ncols_var(ivar)
            IF (nmos.GT.0) THEN
               CALL dbcsr_add( &
                       tidft_env%e_states(istate)%slater_det(1)%mat_var(ivar), &
                       mat_source(istate,ivar), &
                       alpha, beta)
            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

   END SUBROUTINE tidft_optimizer_mainvar_add_to_tidft_env

! *****************************************************************************
!> \brief Compute beta for conjugate gradient algorithms
!> \param beta ...
!> \param reset_conjugator ...
!> \param conjugator ...
!> \param grad ...
!> \param prev_grad ...
!> \param step ...
!> \param prev_step ...
!> \param prev_minus_prec_grad ...
!> \par History
!>       2015.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE compute_cg_beta(beta, reset_conjugator, conjugator, &
      grad, prev_grad, step, prev_step, prev_minus_prec_grad, ncols_var)

      REAL(KIND=dp), INTENT(INOUT)                       :: beta
      LOGICAL, INTENT(INOUT)                             :: reset_conjugator
      INTEGER, INTENT(IN)                                :: conjugator
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: grad, prev_grad, step, prev_step
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT), &
         OPTIONAL                                        :: prev_minus_prec_grad
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var

      CHARACTER(len=*), PARAMETER                        :: routineN = 'compute_cg_beta'

      INTEGER                                            :: handle, i, nsize, unit_nr
      REAL(KIND=dp)                                      :: den, kappa, my_denom, my_numer, &
                                                            my_numer2, my_numer3, num, num2, num3, &
                                                            tau
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: m_tmp_no_1

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

!tiidft      IF (.NOT. PRESENT(prev_minus_prec_grad)) THEN
!tiidft         IF (conjugator .EQ. cg_fletcher_reeves .OR. &
!tiidft             conjugator .EQ. cg_polak_ribiere .OR. &
!tiidft             conjugator .EQ. cg_hager_zhang) THEN
!tiidft            CPABORT("conjugator needs more input")
!tiidft         ENDIF
!tiidft      ENDIF
!tiidft
!tiidft      nsize = SIZE(grad)
!tiidft
!tiidft      my_numer = 0.0_dp
!tiidft      my_numer2 = 0.0_dp
!tiidft      my_numer3 = 0.0_dp
!tiidft      my_denom = 0.0_dp
!tiidft
!tiidft      DO i = 1, nsize
!tiidft
!tiidft         CALL dbcsr_create(m_tmp_no_1, &
!tiidft                           template=grad(i), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)

         SELECT CASE (conjugator)
         CASE (cg_hestenes_stiefel)
            CPABORT("This conjugator NYI")
!tiidft            CALL dbcsr_copy(m_tmp_no_1, grad(i))
!tiidft            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), &
!tiidft                           1.0_dp, -1.0_dp)
!tiidft            CALL dbcsr_dot(m_tmp_no_1, step(i), num)
!tiidft            CALL dbcsr_dot(m_tmp_no_1, prev_step(i), den)
         CASE (cg_fletcher_reeves)
            CPABORT("This conjugator NYI")
!tiidft            CALL dbcsr_dot(grad(i), step(i), num)
!tiidft            CALL dbcsr_dot(prev_grad(i), prev_minus_prec_grad(i), den)
         CASE (cg_polak_ribiere)
            CPABORT("This conjugator NYI")
!tiidft            CALL dbcsr_dot(prev_grad(i), prev_minus_prec_grad(i), den)
!tiidft            CALL dbcsr_copy(m_tmp_no_1, grad(i))
!tiidft            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), 1.0_dp, -1.0_dp)
!tiidft            CALL dbcsr_dot(m_tmp_no_1, step(i), num)
         CASE (cg_fletcher)
            CALL tidft_optimizer_mainvar_dot( &
                    mat_left = grad, &
                    mat_right = step, &
                    mat_left_covariant = .TRUE., &
                    mat_right_covariant = .FALSE., &
                    dot_product = my_numer, &
                    ncols_var = ncols_var)
            CALL tidft_optimizer_mainvar_dot( &
                    mat_left = prev_grad, &
                    mat_right = prev_step, &
                    mat_left_covariant = .TRUE., &
                    mat_right_covariant = .FALSE., &
                    dot_product = my_denom, &
                    ncols_var = ncols_var)
!tiidft            CALL dbcsr_dot(grad(i), step(i), num)
!tiidft            CALL dbcsr_dot(prev_grad(i), prev_step(i), den)
         CASE (cg_liu_storey)
            CPABORT("This conjugator NYI")
!tiidft            CALL dbcsr_dot(prev_grad(i), prev_step(i), den)
!tiidft            CALL dbcsr_copy(m_tmp_no_1, grad(i))
!tiidft            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), 1.0_dp, -1.0_dp)
!tiidft            CALL dbcsr_dot(m_tmp_no_1, step(i), num)
         CASE (cg_dai_yuan)
            CPABORT("This conjugator NYI")
!tiidft            CALL dbcsr_dot(grad(i), step(i), num)
!tiidft            CALL dbcsr_copy(m_tmp_no_1, grad(i))
!tiidft            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), 1.0_dp, -1.0_dp)
!tiidft            CALL dbcsr_dot(m_tmp_no_1, prev_step(i), den)
         CASE (cg_hager_zhang)
            CPABORT("This conjugator NYI")
!tiidft            CALL dbcsr_copy(m_tmp_no_1, grad(i))
!tiidft            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), 1.0_dp, -1.0_dp)
!tiidft            CALL dbcsr_dot(m_tmp_no_1, prev_step(i), den)
!tiidft            CALL dbcsr_dot(m_tmp_no_1, prev_minus_prec_grad(i), num)
!tiidft            CALL dbcsr_dot(m_tmp_no_1, step(i), num2)
!tiidft            CALL dbcsr_dot(prev_step(i), grad(i), num3)
!tiidft            my_numer2 = my_numer2 + num2
!tiidft            my_numer3 = my_numer3 + num3
         CASE (cg_zero)
            num = 0.0_dp
            den = 1.0_dp
         CASE DEFAULT
            CPABORT("illegal conjugator")
         END SELECT
!tiidft         my_numer = my_numer + num
!tiidft         my_denom = my_denom + den
!tiidft
!tiidft         CALL dbcsr_release(m_tmp_no_1)
!tiidft
!tiidft      ENDDO ! i - nsize

!tiidft      DO i = 1, nsize

         SELECT CASE (conjugator)
         CASE (cg_hestenes_stiefel, cg_dai_yuan)
            beta = -1.0_dp*my_numer/my_denom
         CASE (cg_fletcher_reeves, cg_polak_ribiere, cg_fletcher, cg_liu_storey)
            beta = my_numer/my_denom
         CASE (cg_hager_zhang)
            kappa = -2.0_dp*my_numer/my_denom
            tau = -1.0_dp*my_numer2/my_denom
            beta = tau - kappa*my_numer3/my_denom
         CASE (cg_zero)
            beta = 0.0_dp
         CASE DEFAULT
            CPABORT("illegal conjugator")
         END SELECT

!tiidft      ENDDO ! i - nsize

      IF (beta .LT. 0.0_dp) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) " Resetting conjugator because beta is negative: ", beta
         ENDIF
         reset_conjugator = .TRUE.
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE compute_cg_beta

! **************************************************************************************************
!> \brief Compute step direction from the gradient
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_grad_to_step(tidft_env, precond, grad, step, prec_inv, &
      prev_grad, prev_theta, iteration)

      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env
      INTEGER, INTENT(IN)                                :: precond
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: grad
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: step
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: prec_inv
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(IN)       :: prev_grad, prev_theta
      INTEGER                                            :: iteration

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_optimizer_grad_to_step', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, nspins
      LOGICAL                                            :: d_bfgs, l_bfgs
      LOGICAL                                            :: energy_component
      REAL(KIND=dp)                                      :: bfgs_rho, bfgs_sum
      TYPE(dbcsr_type)                                   :: bfgs_s, bfgs_y, tempOccOcc1, &
                                                            tempOccOcc2, tempOccOcc3

!tiidft      nspins = SIZE(nlmo_env%m_theta)

!tiidft      d_bfgs = nlmo_env%d_bfgs
!tiidft      l_bfgs = nlmo_env%l_bfgs
!tiidft
!tiidft      ! if available use second derivative info - bfgs, hessian, preconditioner
!tiidft      IF (nlmo_env%hessian_type .EQ. xalmo_prec_identity) THEN ! no second derivatives
!tiidft

      SELECT CASE (precond)
      CASE (tidft_prec_h_s, tidft_prec_qhq_s, tidft_prec_qhq_qsq)

         IF (iteration .EQ. 0) THEN
            CALL tidft_compute_invert_preconditioner(tidft_env, precond, prec_inv)
         ENDIF

         CALL tidft_optimizer_mainvar_multiply_diff( &
                 "N", "N", &
                 -1.0_dp, &
                 prec_inv, &
                 grad, &
                 0.0_dp, &
                 step, &
                 tidft_env%eps_filter, &
                 tidft_env%e_states(1)%slater_det(1)%ncols_var)

      CASE (tidft_prec_s)

         CALL tidft_optimizer_mainvar_multiply( &
                 "N", "N", &
                 -1.0_dp, &
                 tidft_env%mat_s_inv, &
                 grad, &
                 0.0_dp, &
                 step, &
                 tidft_env%eps_filter, &
                 tidft_env%e_states(1)%slater_det(1)%ncols_var)

      CASE (tidft_prec_zero)
      ! Do not use any preconditioners: use negative gradient

         CALL tidft_optimizer_mainvar_add( &
                 step, &
                 grad, &
                 0.0_dp, -1.0_dp, &
                 ncols_var = tidft_env%e_states(1)%slater_det(1)%ncols_var)

      CASE DEFAULT
         CPABORT("Illegal preconditioner")
      END SELECT

!tiidft
!tiidft      ELSE ! use second derivatives
!tiidft
!tiidft         ! compute and invert hessian/precond?
!tiidft         IF (iteration .EQ. 0) THEN
!tiidft
!tiidft            IF (d_bfgs) THEN
!tiidft
!tiidft               ! create matrix filled with 1.0 here
!tiidft               CALL fill_matrix_with_ones(nlmo_env%m_model_hessian_inv(1))
!tiidft               IF (nspins .GT. 1) THEN
!tiidft                  DO ispin = 2, nspins
!tiidft                     CALL dbcsr_copy(nlmo_env%m_model_hessian_inv(ispin), &
!tiidft                                     nlmo_env%m_model_hessian_inv(1))
!tiidft                  ENDDO
!tiidft               ENDIF
!tiidft
!tiidft            ELSE IF (l_bfgs) THEN
!tiidft
!tiidft               CALL lbfgs_seed(nlmo_env%nlmo_lbfgs_history, nlmo_env%m_theta, grad)
!tiidft               DO ispin = 1, nspins
!tiidft                  CALL dbcsr_copy(step(ispin), grad(ispin))
!tiidft                  CALL dbcsr_scale(step(ispin), -1.0_dp)
!tiidft               ENDDO ! ispin
!tiidft
!tiidft            ELSE
!tiidft
!tiidft               CPWARN("THIS HESSIAN IS IDENTITY MATRIX")
!tiidft               ! computing preconditioner
!tiidft               DO ispin = 1, nspins
!tiidft                  ! TODO: write preconditioner code later
!tiidft                  ! For now, create matrix filled with 1.0 here
!tiidft                  CALL fill_matrix_with_ones(nlmo_env%m_model_hessian_inv(ispin))
!tiidft               ENDDO ! ispin
!tiidft
!tiidft            ENDIF
!tiidft
!tiidft         ELSE ! not iteration zero
!tiidft
!tiidft            ! update approx inverse hessian
!tiidft            IF (d_bfgs) THEN ! diagonal BFGS
!tiidft
!tiidft               DO ispin = 1, nspins
!tiidft
!tiidft                  CALL dbcsr_create(bfgs_s, &
!tiidft                                    template=nlmo_env%m_theta(ispin), &
!tiidft                                    matrix_type=dbcsr_type_no_symmetry)
!tiidft                  CALL dbcsr_create(bfgs_y, &
!tiidft                                    template=nlmo_env%m_theta(ispin), &
!tiidft                                    matrix_type=dbcsr_type_no_symmetry)
!tiidft                  CALL dbcsr_create(tempOccOcc1, &
!tiidft                                    template=nlmo_env%m_theta(ispin), &
!tiidft                                    matrix_type=dbcsr_type_no_symmetry)
!tiidft                  CALL dbcsr_create(tempOccOcc2, &
!tiidft                                    template=nlmo_env%m_theta(ispin), &
!tiidft                                    matrix_type=dbcsr_type_no_symmetry)
!tiidft                  CALL dbcsr_create(tempOccOcc3, &
!tiidft                                    template=nlmo_env%m_theta(ispin), &
!tiidft                                    matrix_type=dbcsr_type_no_symmetry)
!tiidft
!tiidft                  ! compute s and y
!tiidft                  CALL dbcsr_copy(bfgs_y, grad(ispin))
!tiidft                  CALL dbcsr_add(bfgs_y, prev_grad(ispin), 1.0_dp, -1.0_dp)
!tiidft                  CALL dbcsr_copy(bfgs_s, nlmo_env%m_theta(ispin))
!tiidft                  CALL dbcsr_add(bfgs_s, prev_m_theta(ispin), 1.0_dp, -1.0_dp)
!tiidft
!tiidft                  ! compute rho
!tiidft                  CALL dbcsr_dot(grad(ispin), step(ispin), bfgs_rho)
!tiidft                  bfgs_rho = 1.0_dp/bfgs_rho
!tiidft
!tiidft                  ! compute the sum of the squared elements of bfgs_y
!tiidft                  CALL dbcsr_dot(bfgs_y, bfgs_y, bfgs_sum)
!tiidft
!tiidft                  ! first term: start collecting new inv hessian in this temp matrix
!tiidft                  CALL dbcsr_copy(tempOccOcc2, nlmo_env%m_model_hessian_inv(ispin))
!tiidft
!tiidft                  ! second term: + rho * s * s
!tiidft                  CALL dbcsr_hadamard_product(bfgs_s, bfgs_s, tempOccOcc1)
!tiidft                  CALL dbcsr_add(tempOccOcc2, tempOccOcc1, 1.0_dp, bfgs_rho)
!tiidft
!tiidft                  ! third term: + rho^2 * s * s * H * sum_(y * y)
!tiidft                  CALL dbcsr_hadamard_product(tempOccOcc1, &
!tiidft                                              nlmo_env%m_model_hessian_inv(ispin), tempOccOcc3)
!tiidft                  CALL dbcsr_add(tempOccOcc2, tempOccOcc3, &
!tiidft                                 1.0_dp, bfgs_rho*bfgs_rho*bfgs_sum)
!tiidft
!tiidft                  ! fourth term: - 2 * rho * s * y * H
!tiidft                  CALL dbcsr_hadamard_product(bfgs_y, &
!tiidft                                              nlmo_env%m_model_hessian_inv(ispin), tempOccOcc1)
!tiidft                  CALL dbcsr_hadamard_product(bfgs_s, tempOccOcc1, tempOccOcc3)
!tiidft                  CALL dbcsr_add(tempOccOcc2, tempOccOcc3, &
!tiidft                                 1.0_dp, -2.0_dp*bfgs_rho)
!tiidft
!tiidft                  CALL dbcsr_copy(nlmo_env%m_model_hessian_inv(ispin), tempOccOcc2)
!tiidft
!tiidft                  CALL dbcsr_release(tempOccOcc1)
!tiidft                  CALL dbcsr_release(tempOccOcc2)
!tiidft                  CALL dbcsr_release(tempOccOcc3)
!tiidft                  CALL dbcsr_release(bfgs_y)
!tiidft                  CALL dbcsr_release(bfgs_s)
!tiidft
!tiidft               ENDDO ! ispin
!tiidft
!tiidft            ELSE IF (l_bfgs) THEN
!tiidft
!tiidft               CALL lbfgs_get_direction(nlmo_env%nlmo_lbfgs_history, nlmo_env%m_theta, grad, step)
!tiidft
!tiidft            ENDIF ! which method?
!tiidft
!tiidft         ENDIF ! compute approximate inverse hessian
!tiidft
!tiidft         ! apply precomputed hessian
!tiidft         IF (.NOT. l_bfgs) THEN
!tiidft
!tiidft            DO ispin = 1, nspins
!tiidft
!tiidft               CALL dbcsr_hadamard_product(nlmo_env%m_model_hessian_inv(ispin), &
!tiidft                                           grad(ispin), step(ispin))
!tiidft               CALL dbcsr_scale(step(ispin), -1.0_dp)
!tiidft
!tiidft            ENDDO ! ispin
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft      ENDIF ! second derivative type fork

   END SUBROUTINE tidft_optimizer_grad_to_step

! **************************************************************************************************
!> \brief Allocate temporary storage for an optimizer
!> \param 
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_allocate_storage(data_storage, nstates, &
      mat_var_template, mat_nn_template, ncols_var)

      TYPE(tidft_optimizer_storage_type)                 :: data_storage
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: mat_var_template
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_nn_template
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var
      INTEGER, INTENT(IN)                                :: nstates

      INTEGER                                            :: istate, ivar, &
                                                            nvars

      nvars = SIZE(mat_var_template)
      CPASSERT(nvars.EQ.3)

      ALLOCATE (data_storage%prev_theta(nstates,nvars))
      ALLOCATE (data_storage%prev_grad(nstates,nvars))
      ALLOCATE (data_storage%prev_step(nstates,nvars))
      ALLOCATE (data_storage%step(nstates,nvars))
      ALLOCATE (data_storage%grad(nstates,nvars))
      ALLOCATE (data_storage%prev_minus_prec_grad(nstates,nvars))

      ALLOCATE (data_storage%prec_inv(nstates,nvars))

      DO istate = 1, nstates
         DO ivar = 1, nvars
            IF (ncols_var(ivar).GT.0) THEN

               ! init storage
               CALL dbcsr_create(data_storage%prev_theta(istate,ivar), &
                                 template=mat_var_template(ivar), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(data_storage%prev_grad(istate,ivar), &
                                 template=mat_var_template(ivar), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(data_storage%prev_step(istate,ivar), &
                                 template=mat_var_template(ivar), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(data_storage%step(istate,ivar), &
                                 template=mat_var_template(ivar), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(data_storage%prev_minus_prec_grad(istate,ivar), &
                                 template=mat_var_template(ivar), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(data_storage%grad(istate,ivar), &
                                 template=mat_var_template(ivar), &
                                 matrix_type=dbcsr_type_no_symmetry)

               CALL dbcsr_create(data_storage%prec_inv(istate,ivar), &
                                 template=mat_nn_template, &
                                 matrix_type=dbcsr_type_no_symmetry)

               CALL dbcsr_set(data_storage%prev_theta(istate,ivar), 0.0_dp)
               CALL dbcsr_set(data_storage%prev_grad(istate,ivar), 0.0_dp)
               CALL dbcsr_set(data_storage%prev_step(istate,ivar), 0.0_dp)
               CALL dbcsr_set(data_storage%step(istate,ivar), 0.0_dp)
               CALL dbcsr_set(data_storage%prev_minus_prec_grad(istate,ivar), 0.0_dp)
               CALL dbcsr_set(data_storage%grad(istate,ivar), 0.0_dp)

               CALL dbcsr_set(data_storage%prec_inv(istate,ivar), 0.0_dp)

            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

   END SUBROUTINE tidft_optimizer_allocate_storage

! **************************************************************************************************
!> \brief Allocate temporary storage for an optimizer
!> \param 
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_optimizer_deallocate_storage(data_storage, nstates, &
      ncols_var)

      TYPE(tidft_optimizer_storage_type)                 :: data_storage
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ncols_var
      INTEGER, INTENT(IN)                                :: nstates

      INTEGER                                            :: istate, ivar, &
                                                            nvars

      nvars = 3

      DO istate = 1, nstates
         DO ivar = 1, nvars
            IF (ncols_var(ivar).GT.0) THEN

               ! init storage
               CALL dbcsr_release(data_storage%prev_theta(istate,ivar))
               CALL dbcsr_release(data_storage%prev_grad(istate,ivar))
               CALL dbcsr_release(data_storage%prev_step(istate,ivar))
               CALL dbcsr_release(data_storage%step(istate,ivar))
               CALL dbcsr_release(data_storage%prev_minus_prec_grad(istate,ivar))
               CALL dbcsr_release(data_storage%grad(istate,ivar))

               CALL dbcsr_release(data_storage%prec_inv(istate,ivar))

            ENDIF
         ENDDO ! ivar
      ENDDO ! stateI

      DEALLOCATE (data_storage%prev_theta)
      DEALLOCATE (data_storage%prev_grad)
      DEALLOCATE (data_storage%prev_step)
      DEALLOCATE (data_storage%step)
      DEALLOCATE (data_storage%grad)
      DEALLOCATE (data_storage%prev_minus_prec_grad)

      DEALLOCATE (data_storage%prec_inv)

   END SUBROUTINE tidft_optimizer_deallocate_storage

! **************************************************************************************************
!> \brief Computes and inverts preconditioner
!> \param 
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_compute_invert_preconditioner(tidft_env, precond, prec_inv)

      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env
      INTEGER, INTENT(IN)                                :: precond
      TYPE(dbcsr_type), DIMENSION(:,:), INTENT(INOUT)    :: prec_inv

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_compute_invert_preconditioner'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: temp_NN1
      TYPE(dbcsr_type)                                   :: temp_NO1
      TYPE(dbcsr_type)                                   :: temp_ST
      TYPE(dbcsr_type)                                   :: temp_FTsiginv
      TYPE(dbcsr_type)                                   :: temp_siginvTFTsiginv
      TYPE(dbcsr_type)                                   :: temp_prec
      INTEGER                                            :: stateI
      INTEGER                                            :: idet
      INTEGER                                            :: ispin
      INTEGER                                            :: ivar
      INTEGER                                            :: nmos_stateI
      INTEGER                                            :: first_state
      INTEGER                                            :: nstates
      INTEGER                                            :: ndets
      INTEGER                                            :: nspins
      INTEGER                                            :: min_zeros
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp)                                      :: scaling_factor_slater_determinant
      REAL(KIND=dp)                                      :: scaling_factor_spin
      REAL(KIND=dp)                                      :: zero_eigenvalue_eps_factor
      LOGICAL                                            :: energy_component
      LOGICAL                                            :: hamiltonian_project_R_out
      LOGICAL                                            :: overlap_project_R_out
      CHARACTER                                          :: mtype

      CALL timeset(routineN, handle)

      alpha = 1.0_dp
      zero_eigenvalue_eps_factor = 1.0_dp

      SELECT CASE (precond)
      CASE (tidft_prec_h_s)

         overlap_project_R_out = .FALSE.
         hamiltonian_project_R_out = .FALSE.

      CASE (tidft_prec_qhq_s)

         overlap_project_R_out = .FALSE.
         hamiltonian_project_R_out = .TRUE.

      CASE (tidft_prec_qhq_qsq)

         overlap_project_R_out = .TRUE.
         hamiltonian_project_R_out = .TRUE.

      CASE DEFAULT

         CPABORT("Cannot process this preconditioner")

      END SELECT

      nspins = tidft_env%nspins
      nstates = tidft_env%nstates
      ndets=1
      IF (tidft_env%triplet_correction) ndets = 2

      first_state = 1
      IF (tidft_env%state_by_state) first_state = nstates

      energy_component=tidft_env%loss_function_components(1)

      ! Use the KS matrix only if the energy component is on
      IF (energy_component) THEN

         ! Clear all states, not only the current state
         DO stateI = 1, nstates
            DO ivar = 1, 3
               IF (tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar).GT.0) THEN
                  CALL dbcsr_set(prec_inv(stateI,ivar), 0.0_dp)
               ENDIF
            ENDDO ! ivar
         ENDDO ! stateI

         DO stateI = first_state, nstates

            DO idet = 1, ndets

               scaling_factor_slater_determinant = 1.0_dp
               IF (tidft_env%triplet_correction) THEN
                  IF (idet.EQ.1) scaling_factor_slater_determinant = 2.0_dp
                  IF (idet.EQ.2) scaling_factor_slater_determinant = -1.0_dp
               ENDIF

               DO ispin = 1, nspins

                  scaling_factor_spin = 1.0_dp
                  IF (nspins .EQ. 1) scaling_factor_spin = 2.0_dp

                  nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)

                  IF ( nmos_stateI .GT. 0) THEN

                     !IF (unit_nr > 0) THEN
                     !   WRITE(unit_nr, '(T2,A28,3I4,A15,I4,A25,3I4)') &
                     !      "Convert (state,spin,det): ", stateI, ispin, idet, &
                     !      " containing ", &
                     !      nmos_stateI, "MOs into VARs containing", &
                     !      tidft_env%e_states(stateI)%slater_det(1)%ncols_var
                     !ENDIF

                     CALL dbcsr_create(temp_NO1, &
                             template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
                     CALL dbcsr_create(temp_ST, &
                             template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
                     CALL dbcsr_create(temp_FTsiginv, &
                             template=tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
                     CALL dbcsr_create(temp_siginvTFTsiginv, &
                             template=tidft_env%mat_sigma_inv(stateI,idet,ispin))
                     CALL dbcsr_create(temp_prec, &
                             template=tidft_env%mat_s, matrix_type = dbcsr_type_no_symmetry)

                     IF (hamiltonian_project_R_out) THEN
                        ! T.siginv
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                tidft_env%mat_sigma_inv(stateI,idet,ispin), &
                                0.0_dp, temp_NO1, &
                                filter_eps=tidft_env%eps_filter)
                        ! F.Tsiginv
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(stateI)%slater_det(idet)%mat_ks(ispin), &
                                temp_NO1, &
                                0.0_dp, temp_FTsiginv, &
                                filter_eps=tidft_env%eps_filter)
                        ! siginvT.FTsiginv
                        CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                temp_NO1, &
                                temp_FTsiginv, &
                                0.0_dp, temp_siginvTFTsiginv, &
                                filter_eps=tidft_env%eps_filter)
                     ENDIF

                     IF (hamiltonian_project_R_out .OR. overlap_project_R_out) THEN
                        ! S.T
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%mat_s, &
                                tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                                0.0_dp, temp_ST, &
                                filter_eps=tidft_env%eps_filter)
                     ENDIF

                     ! prec = F
                     CALL dbcsr_desymmetrize( &
                             tidft_env%e_states(stateI)%slater_det(idet)%mat_ks(ispin), &
                             temp_prec)

                     IF (hamiltonian_project_R_out) THEN
                        ! prec = prec - FTsiginv.(ST)^t
                        CALL dbcsr_multiply("N", "T", -1.0_dp, &
                                temp_FTsiginv, &
                                temp_ST, &
                                1.0_dp, temp_prec, &
                                filter_eps=tidft_env%eps_filter)
                        ! prec = prec - ST.(FTsiginv)^t
                        CALL dbcsr_multiply("N", "T", -1.0_dp, &
                                temp_ST, &
                                temp_FTsiginv, &
                                1.0_dp, temp_prec, &
                                filter_eps=tidft_env%eps_filter)
                        ! ST.(siginvTFTsiginv)
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                temp_ST, &
                                temp_siginvTFTsiginv, &
                                0.0_dp, temp_NO1, &
                                filter_eps=tidft_env%eps_filter)
                        ! prec = prec + STsiginvTFTsiginv.(ST)^t
                        CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                temp_NO1, &
                                temp_ST, &
                                1.0_dp, temp_prec, &
                                filter_eps=tidft_env%eps_filter)
                     ENDIF

                     ! prec = (1-SR)F(1-RS) + alpha * S
                     CALL dbcsr_create(temp_NN1, &
                             template=tidft_env%mat_s, matrix_type = dbcsr_type_no_symmetry)
                     CALL dbcsr_desymmetrize( &
                             tidft_env%mat_s, &
                             temp_NN1)
                     CALL dbcsr_add(temp_prec, &
                             temp_NN1, 1.0_dp, alpha)
                     CALL dbcsr_release(temp_NN1)

                     IF (overlap_project_R_out) THEN
                        ! ST.siginv
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                temp_ST, &
                                tidft_env%mat_sigma_inv(stateI,idet,ispin), &
                                0.0_dp, temp_NO1, &
                                filter_eps=tidft_env%eps_filter)
                        ! prec = prec - alpha * ST.(STsiginv)^t
                        CALL dbcsr_multiply("N", "T", -alpha, &
                                temp_ST, &
                                temp_NO1, &
                                1.0_dp, temp_prec, &
                                filter_eps=tidft_env%eps_filter)
                     ENDIF

                     ! Scale by the required factors
                     CALL dbcsr_scale(temp_prec, &
                             2.0_dp * scaling_factor_slater_determinant * scaling_factor_spin)

                     CALL precMO_to_precVAR_increment( &
                             mat_prec_out = prec_inv(stateI, :), &
                             ncols_var = tidft_env%e_states(stateI)%slater_det(1)%ncols_var, &
                             mat_prec_in = temp_prec, &
                             idet = idet, &
                             ispin = ispin)

                     CALL dbcsr_release(temp_NO1)
                     CALL dbcsr_release(temp_ST)
                     CALL dbcsr_release(temp_FTsiginv)
                     CALL dbcsr_release(temp_siginvTFTsiginv)
                     CALL dbcsr_release(temp_prec)

                  ENDIF ! zero MOs

               ENDDO ! ispin

            ENDDO ! idet

            ! Invert or pseudoinvert
            DO ivar = 1, 3
               IF (tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar).GT.0) THEN

                  SELECT CASE (precond)
                  CASE (tidft_prec_h_s, tidft_prec_qhq_s)
                  ! Invert using cholesky for dense matrices

                     CALL cp_dbcsr_cholesky_decompose( &
                             prec_inv(stateI,ivar), &
                             para_env=tidft_env%para_env, &
                             blacs_env=tidft_env%blacs_env)
                     CALL cp_dbcsr_cholesky_invert( &
                             prec_inv(stateI,ivar), &
                             para_env=tidft_env%para_env, &
                             blacs_env=tidft_env%blacs_env, &
                             upper_to_full=.TRUE.)
                     CALL dbcsr_filter( &
                             prec_inv(stateI,ivar), &
                             tidft_env%eps_filter)

                  CASE (tidft_prec_qhq_qsq)
                  ! Pseudo-invert using diagonalization

                     ! Determine how many eigenvectors are expected to be zero
                     ! If ndets == 1 then
                     ! Var1: nvar1+nvar2, nvar1+nvar3 (var1 includes spins alpha and beta)
                     ! Var2: nvar1+nvar2              (var2 includes only alpha spins)
                     ! Var3: nvar1+nvar3              (var3 includes only beta spins)
                     ! If ndets == 2 then the following values are added to the list
                     ! Var1: nvar1+nvar2+nvar3, nvar1 (var1 includes spins alpha and beta for the triplet state)
                     ! Var2: nvar1+nvar2+nvar3        (var2 includes only alpha electrons for the triplet state)
                     ! Var3: nvar1                    (var3 includes only beta electrons for the triplet state)
                     ! Thus minimum expected number of zeros in the preconditioner is:
                     ! Var1: If ndets == 1: nvar1, else: nvar1
                     ! (the number of zeros might happen to be higher if singly-occupied alpha and beta are the same)
                     ! Var2: If ndets == 1: nvar1+nvar2, else: nvar1+nvar2
                     ! Var3: If ndets == 1: nvar1+nvar3, else: nvar1 (fewer beta orbitals in the triplet state)
                     IF (ivar.EQ.1) THEN
                        min_zeros = &
                           tidft_env%e_states(stateI)%slater_det(1)%ncols_var(1)
                     ELSE IF (ivar.EQ.2) THEN
                        min_zeros = &
                           tidft_env%e_states(stateI)%slater_det(1)%ncols_var(1) + &
                           tidft_env%e_states(stateI)%slater_det(1)%ncols_var(2)
                     ELSE IF (ivar.EQ.3) THEN
                        IF (ndets.EQ.1) THEN
                           min_zeros = &
                              tidft_env%e_states(stateI)%slater_det(1)%ncols_var(1) + &
                              tidft_env%e_states(stateI)%slater_det(1)%ncols_var(3)
                        ELSE
                           min_zeros = &
                              tidft_env%e_states(stateI)%slater_det(1)%ncols_var(1)
                        ENDIF
                     ENDIF
                     CALL cp_dbcsr_pseudoinverse( &
                             matrix=prec_inv(stateI,ivar), &
                             zero_eps=tidft_env%eps_filter*zero_eigenvalue_eps_factor, &
                             min_eigenvalues_to_nullify=min_zeros, &
                             spd_expected=.TRUE., &
                             para_env=tidft_env%para_env, &
                             blacs_env=tidft_env%blacs_env)
                     CALL dbcsr_filter( &
                             prec_inv(stateI,ivar), &
                             tidft_env%eps_filter)

                  CASE DEFAULT

                     CPABORT("Cannot process this preconditioner")

                  END SELECT

               ENDIF ! zero columns check
            ENDDO ! ivar

         ENDDO ! stateI

      ELSE  ! energy_component

         CPABORT("This preconditioner needs the energy component")

      ENDIF ! energy_component

      CALL timestop(handle)

   END SUBROUTINE tidft_compute_invert_preconditioner

! **************************************************************************************************
!> \brief Optimization of ALMOs using trust region minimizers
!> \param qs_env ...
!> \param nlmo_env ...
!> \param optimizer   controls the optimization algorithm
!> \param matrix_s ...
!> \param matrix_mo_inout ...
!> \param matrix_templateOO ...
!> \par History
!>       2020.01 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
!tiidft   SUBROUTINE construct_nlmos_trustr(qs_env, nlmo_env, optimizer, &
!tiidft                                     matrix_s, matrix_mo_inout, &
!tiidft                                     matrix_templateOO)
!tiidft
!tiidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tiidft      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
!tiidft      TYPE(optimizer_options_type), INTENT(INOUT)        :: optimizer
!tiidft      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s
!tiidft      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:), &
!tiidft         INTENT(INOUT)                                   :: matrix_mo_inout
!tiidft      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:), &
!tiidft         INTENT(IN)                                      :: matrix_templateOO
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'construct_nlmos_trustr', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: handle, ispin, iteration, &
!tiidft                                                            iteration_type_to_report, nspins, &
!tiidft                                                            outer_iteration, unit_nr
!tiidft      LOGICAL                                            :: border_reached, inner_loop_success, &
!tiidft                                                            reset_conjugator, same_position, &
!tiidft                                                            scf_converged
!tiidft      REAL(kind=dp) :: beta, energy_start, energy_trial, eta, expected_reduction, &
!tiidft         fake_step_size_to_report, grad_norm_ratio, grad_norm_ref, loss_change_to_report, &
!tiidft         loss_start, loss_trial, model_grad_norm, penalty_start, penalty_trial, radius_current, &
!tiidft         radius_max, real_temp, rho, spin_factor, step_norm, step_size, t1, t1outer, t2, t2outer, &
!tiidft         y_scalar
!tiidft      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: grad_norm_spin
!tiidft      REAL(KIND=dp), DIMENSION(2)                        :: loss_components, overlap_determinant
!tiidft      TYPE(cp_logger_type), POINTER                      :: logger
!tiidft      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: m_model_Bd, m_model_d, m_model_r, &
!tiidft                                                            m_model_r_prev, m_model_rt, &
!tiidft                                                            m_model_rt_prev, m_theta_saved, &
!tiidft                                                            prev_step, step
!tiidft
!tiidft      CALL timeset(routineN, handle)
!tiidft
!tiidft      t1outer = m_walltime()
!tiidft
!tiidft      ! get a useful output_unit
!tiidft      logger => cp_get_default_logger()
!tiidft      IF (logger%para_env%mepos == logger%para_env%source) THEN
!tiidft         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
!tiidft      ELSE
!tiidft         unit_nr = -1
!tiidft      ENDIF
!tiidft
!tiidft      nspins = SIZE(matrix_mo_inout)
!tiidft      CALL nlmo_env_set_flags(nlmo_env=nlmo_env, &
!tiidft                              optimizer=optimizer)
!tiidft      IF (nlmo_env%hessian_type .EQ. xalmo_prec_lbfgs .OR. &
!tiidft          nlmo_env%hessian_type .EQ. xalmo_prec_dbfgs) THEN
!tiidft         CPABORT("BFGS cannot be combined with BFGS")
!tiidft      ENDIF
!tiidft      CALL nlmo_env_allocate_all(nlmo_env=nlmo_env, &
!tiidft                                 templateOO=matrix_templateOO, &
!tiidft                                 templateNO=matrix_mo_inout)
!tiidft      CALL nlmo_env_mainvar_initial_guess(nlmo_env=nlmo_env, &
!tiidft                                          m_mo_initial=matrix_mo_inout)
!tiidft      ! compute the overlap of the initial orbitals
!tiidft      CALL nlmo_env_compute_sigma0_ij(nlmo_env=nlmo_env, &
!tiidft                                      matrix_s=matrix_s, &
!tiidft                                      invert=.TRUE.)
!tiidft      ! use initial orbitals to compute the BK matrix
!tiidft      CALL nlmo_env_compute_BK_ij(nlmo_env=nlmo_env, &
!tiidft                                  qs_env=qs_env, &
!tiidft                                  m_metricNN=matrix_s, &
!tiidft                                  m_templateOO=nlmo_env%m_sigma0)
!tiidft
!tiidft      nspins = SIZE(matrix_mo_inout)
!tiidft      IF (nspins == 1) THEN
!tiidft         spin_factor = 2.0_dp
!tiidft      ELSE
!tiidft         spin_factor = 1.0_dp
!tiidft      ENDIF
!tiidft
!tiidft      ALLOCATE (grad_norm_spin(nspins))
!tiidft
!tiidft      ALLOCATE (prev_step(nspins))
!tiidft      ALLOCATE (step(nspins))
!tiidft      ALLOCATE (m_model_r(nspins))
!tiidft      ALLOCATE (m_model_rt(nspins))
!tiidft      ALLOCATE (m_model_d(nspins))
!tiidft      ALLOCATE (m_model_Bd(nspins))
!tiidft      ALLOCATE (m_model_r_prev(nspins))
!tiidft      ALLOCATE (m_model_rt_prev(nspins))
!tiidft      ALLOCATE (m_theta_saved(nspins))
!tiidft
!tiidft      DO ispin = 1, nspins
!tiidft
!tiidft         CALL dbcsr_create(prev_step(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(step(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(m_model_r(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(m_model_rt(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(m_model_d(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(m_model_Bd(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(m_model_r_prev(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(m_model_rt_prev(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft         CALL dbcsr_create(m_theta_saved(ispin), &
!tiidft                           template=nlmo_env%m_sigma0(ispin), &
!tiidft                           matrix_type=dbcsr_type_no_symmetry)
!tiidft
!tiidft         CALL dbcsr_set(step(ispin), 0.0_dp)
!tiidft         CALL dbcsr_set(prev_step(ispin), 0.0_dp)
!tiidft
!tiidft      ENDDO ! ispin
!tiidft
!tiidft      radius_max = optimizer%max_trust_radius
!tiidft      radius_current = MIN(optimizer%initial_trust_radius, radius_max)
!tiidft      ! eta must be between 0 and 0.25
!tiidft      eta = MIN(MAX(optimizer%rho_do_not_update, 0.0_dp), 0.25_dp)
!tiidft      energy_start = 0.0_dp
!tiidft      energy_trial = 0.0_dp
!tiidft      penalty_start = 0.0_dp
!tiidft      penalty_trial = 0.0_dp
!tiidft      loss_start = 0.0_dp ! sum of the energy and penalty
!tiidft      loss_trial = 0.0_dp
!tiidft
!tiidft      same_position = .FALSE.
!tiidft
!tiidft      ! use the main variable to compute aux matrices
!tiidft      CALL nlmo_env_mainvar_to_aux(nlmo_env)
!tiidft
!tiidft      ! compute objective (loss) function
!tiidft      loss_components(:) = 0.0_dp
!tiidft      CALL nlmo_env_loss_function(nlmo_env, loss_components, overlap_determinant)
!tiidft      energy_start = loss_components(1)
!tiidft      penalty_start = loss_components(2)
!tiidft      loss_start = SUM(loss_components)
!tiidft
!tiidft      ! start the outer step-size-adjustment loop
!tiidft      scf_converged = .FALSE.
!tiidft      adjust_r_loop: DO outer_iteration = 1, optimizer%max_iter_outer_loop
!tiidft
!tiidft         ! start the inner fixed-radius loop
!tiidft         border_reached = .FALSE.
!tiidft
!tiidft         DO ispin = 1, nspins
!tiidft            CALL dbcsr_set(step(ispin), 0.0_dp)
!tiidft            CALL dbcsr_filter(step(ispin), nlmo_env%eps_filter)
!tiidft         ENDDO
!tiidft
!tiidft         IF (.NOT. same_position) THEN
!tiidft
!tiidft            IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Compute model gradient"
!tiidft            CALL nlmo_env_loss_gradient(nlmo_env)
!tiidft
!tiidft         ENDIF ! skip_grad
!tiidft
!tiidft         ! check convergence and other exit criteria
!tiidft         DO ispin = 1, nspins
!tiidft            CALL dbcsr_norm(nlmo_env%grad(ispin), dbcsr_norm_maxabsnorm, &
!tiidft                            norm_scalar=grad_norm_spin(ispin))
!tiidft         ENDDO ! ispin
!tiidft         grad_norm_ref = MAXVAL(grad_norm_spin)
!tiidft
!tiidft         t2outer = m_walltime()
!tiidft         CALL trust_r_report(unit_nr, &
!tiidft                             iter_type=1, & ! only some data is important
!tiidft                             iteration=outer_iteration, &
!tiidft                             loss=loss_start, &
!tiidft                             delta_loss=0.0_dp, &
!tiidft                             grad_norm=grad_norm_ref, &
!tiidft                             predicted_reduction=0.0_dp, &
!tiidft                             rho=0.0_dp, &
!tiidft                             radius=radius_current, &
!tiidft                             new=.NOT. same_position, &
!tiidft                             time=t2outer - t1outer)
!tiidft         t1outer = m_walltime()
!tiidft
!tiidft         IF (grad_norm_ref .LE. optimizer%eps_error) THEN
!tiidft            scf_converged = .TRUE.
!tiidft            border_reached = .FALSE.
!tiidft            expected_reduction = 0.0_dp
!tiidft            IF (.NOT. (optimizer%early_stopping_on .AND. outer_iteration .EQ. 1)) &
!tiidft               EXIT adjust_r_loop
!tiidft         ELSE
!tiidft            scf_converged = .FALSE.
!tiidft         ENDIF
!tiidft
!tiidft         DO ispin = 1, nspins
!tiidft
!tiidft            CALL dbcsr_copy(m_model_r(ispin), nlmo_env%grad(ispin))
!tiidft            CALL dbcsr_scale(m_model_r(ispin), -1.0_dp)
!tiidft
!tiidft            CALL dbcsr_multiply("N", "N", 1.0_dp, &
!tiidft                                nlmo_env%m_siginv0(ispin), &
!tiidft                                m_model_r(ispin), &
!tiidft                                0.0_dp, m_model_rt(ispin), &
!tiidft                                filter_eps=nlmo_env%eps_filter)
!tiidft
!tiidft            CALL dbcsr_copy(m_model_d(ispin), m_model_rt(ispin))
!tiidft
!tiidft         ENDDO ! ispin
!tiidft
!tiidft         ! compute model Hessian
!tiidft         IF (.NOT. same_position) THEN
!tiidft
!tiidft            ! RZK-critical: Compute model Hessian
!tiidft            CALL nlmo_env_loss_hessian(nlmo_env)
!tiidft
!tiidft         ENDIF  ! not same position
!tiidft
!tiidft         ! print the header (argument values are ignored)
!tiidft         CALL fixed_r_report(unit_nr, &
!tiidft                             iter_type=0, &
!tiidft                             iteration=0, &
!tiidft                             step_size=0.0_dp, &
!tiidft                             border_reached=.FALSE., &
!tiidft                             curvature=0.0_dp, &
!tiidft                             grad_norm_ratio=0.0_dp, &
!tiidft                             time=0.0_dp)
!tiidft
!tiidft         IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Start inner loop"
!tiidft
!tiidft         t1 = m_walltime()
!tiidft         inner_loop_success = .FALSE.
!tiidft         ! trustr_steihaug, trustr_cauchy, trustr_dogleg
!tiidft         fixed_r_loop: DO iteration = 1, optimizer%max_iter
!tiidft
!tiidft            ! Step 2. Get curvature. If negative, step to the border
!tiidft            ! Compute B.d
!tiidft            CALL nlmo_env_hessian_apply(nlmo_env=nlmo_env, m_in=m_model_d, m_out=m_model_Bd)
!tiidft            !RZK-critical: y_scalar should be separate for each spin?
!tiidft            y_scalar = 0.0_dp
!tiidft            DO ispin = 1, nspins
!tiidft
!tiidft               ! Get y=d^T.B.d
!tiidft               CALL dbcsr_dot(m_model_d(ispin), m_model_Bd(ispin), real_temp)
!tiidft
!tiidft               y_scalar = y_scalar + real_temp
!tiidft
!tiidft            ENDDO ! ispin
!tiidft            IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Curvature: ", y_scalar
!tiidft
!tiidft            ! step to the border
!tiidft            IF (y_scalar .LT. 0.0_dp) THEN
!tiidft
!tiidft               CALL step_size_to_border( &
!tiidft                  step_size_out=step_size, &
!tiidft                  metric_in=nlmo_env%m_sigma0, &
!tiidft                  position_in=step, &
!tiidft                  direction_in=m_model_d, &
!tiidft                  trust_radius_in=radius_current, &
!tiidft                  eps_filter_in=nlmo_env%eps_filter)
!tiidft
!tiidft               DO ispin = 1, nspins
!tiidft                  CALL dbcsr_add(step(ispin), m_model_d(ispin), 1.0_dp, step_size)
!tiidft               ENDDO
!tiidft
!tiidft               border_reached = .TRUE.
!tiidft               inner_loop_success = .TRUE.
!tiidft
!tiidft               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
!tiidft                                                 reduction_out=expected_reduction, &
!tiidft                                                 step_in=step)
!tiidft
!tiidft               t2 = m_walltime()
!tiidft               CALL fixed_r_report(unit_nr, &
!tiidft                                   iter_type=2, &
!tiidft                                   iteration=iteration, &
!tiidft                                   step_size=step_size, &
!tiidft                                   border_reached=border_reached, &
!tiidft                                   curvature=y_scalar, &
!tiidft                                   grad_norm_ratio=expected_reduction, &
!tiidft                                   time=t2 - t1)
!tiidft
!tiidft               EXIT fixed_r_loop ! the inner loop
!tiidft
!tiidft            ENDIF ! y is negative
!tiidft
!tiidft            ! Step 3. Compute the step size along the direction
!tiidft            step_size = 0.0_dp
!tiidft            DO ispin = 1, nspins
!tiidft               CALL dbcsr_dot(m_model_r(ispin), m_model_rt(ispin), real_temp)
!tiidft               step_size = step_size + real_temp
!tiidft            ENDDO ! ispin
!tiidft            step_size = step_size/y_scalar
!tiidft            IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Proposed step size: ", step_size
!tiidft
!tiidft            ! Update the step matrix
!tiidft            DO ispin = 1, nspins
!tiidft               CALL dbcsr_copy(prev_step(ispin), step(ispin))
!tiidft               CALL dbcsr_add(step(ispin), m_model_d(ispin), 1.0_dp, step_size)
!tiidft            ENDDO
!tiidft
!tiidft            ! Compute step norm
!tiidft            CALL contravariant_matrix_norm( &
!tiidft               norm_out=step_norm, &
!tiidft               matrix_in=step, &
!tiidft               metric_in=nlmo_env%m_sigma0, &
!tiidft               eps_filter_in=nlmo_env%eps_filter)
!tiidft            IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Step norm: ", step_norm
!tiidft
!tiidft            ! Do not step beyond the trust radius
!tiidft            IF (step_norm .GT. radius_current) THEN
!tiidft
!tiidft               IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Norm is too large"
!tiidft               CALL step_size_to_border( &
!tiidft                  step_size_out=step_size, &
!tiidft                  metric_in=nlmo_env%m_sigma0, &
!tiidft                  position_in=prev_step, &
!tiidft                  direction_in=m_model_d, &
!tiidft                  trust_radius_in=radius_current, &
!tiidft                  eps_filter_in=nlmo_env%eps_filter)
!tiidft               IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Step size to border: ", step_size
!tiidft
!tiidft               DO ispin = 1, nspins
!tiidft                  CALL dbcsr_copy(step(ispin), prev_step(ispin))
!tiidft                  CALL dbcsr_add(step(ispin), m_model_d(ispin), 1.0_dp, step_size)
!tiidft               ENDDO
!tiidft
!tiidft               IF (debug_mode) THEN
!tiidft                  ! Compute step norm
!tiidft                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Extra norm evaluation"
!tiidft                  CALL contravariant_matrix_norm( &
!tiidft                     norm_out=step_norm, &
!tiidft                     matrix_in=step, &
!tiidft                     metric_in=nlmo_env%m_sigma0, &
!tiidft                     eps_filter_in=nlmo_env%eps_filter)
!tiidft                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Step norm: ", step_norm
!tiidft                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Current radius: ", radius_current
!tiidft               ENDIF
!tiidft
!tiidft               border_reached = .TRUE.
!tiidft               inner_loop_success = .TRUE.
!tiidft
!tiidft               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
!tiidft                                                 reduction_out=expected_reduction, &
!tiidft                                                 step_in=step)
!tiidft
!tiidft               t2 = m_walltime()
!tiidft               CALL fixed_r_report(unit_nr, &
!tiidft                                   iter_type=3, &
!tiidft                                   iteration=iteration, &
!tiidft                                   step_size=step_size, &
!tiidft                                   border_reached=border_reached, &
!tiidft                                   curvature=y_scalar, &
!tiidft                                   grad_norm_ratio=expected_reduction, &
!tiidft                                   time=t2 - t1)
!tiidft
!tiidft               EXIT fixed_r_loop ! the inner loop
!tiidft
!tiidft            ENDIF
!tiidft
!tiidft            IF (optimizer%trustr_algorithm .EQ. trustr_cauchy) THEN
!tiidft               ! trustr_steihaug, trustr_cauchy, trustr_dogleg
!tiidft
!tiidft               border_reached = .FALSE.
!tiidft               inner_loop_success = .TRUE.
!tiidft
!tiidft               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
!tiidft                                                 reduction_out=expected_reduction, &
!tiidft                                                 step_in=step)
!tiidft
!tiidft               t2 = m_walltime()
!tiidft               CALL fixed_r_report(unit_nr, &
!tiidft                                   iter_type=5, & ! Cauchy point
!tiidft                                   iteration=iteration, &
!tiidft                                   step_size=step_size, &
!tiidft                                   border_reached=border_reached, &
!tiidft                                   curvature=y_scalar, &
!tiidft                                   grad_norm_ratio=expected_reduction, &
!tiidft                                   time=t2 - t1)
!tiidft
!tiidft               EXIT fixed_r_loop ! the inner loop
!tiidft
!tiidft            ELSE IF (optimizer%trustr_algorithm .EQ. trustr_dogleg) THEN
!tiidft
!tiidft               ! invert B
!tiidft               IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...(Pseudo-)invert model Hessian"
!tiidft               CALL nlmo_env_invert_hessian(nlmo_env)
!tiidft
!tiidft               ! get pB = Binv.m_model_r = -Binv.grad
!tiidft               ! RZK-critical: implement an nlmo_env method
!tiidft               DO ispin = 1, nspins
!tiidft
!tiidft                  CALL dbcsr_multiply("N", "N", 1.0_dp, &
!tiidft                                      nlmo_env%m_model_hessian_inv(ispin), &
!tiidft                                      m_model_r(ispin), &
!tiidft                                      0.0_dp, m_model_Bd(ispin), &
!tiidft                                      filter_eps=nlmo_env%eps_filter)
!tiidft
!tiidft               ENDDO ! ispin
!tiidft
!tiidft               ! Compute norm of pB
!tiidft               CALL contravariant_matrix_norm( &
!tiidft                  norm_out=step_norm, &
!tiidft                  matrix_in=m_model_Bd, &
!tiidft                  metric_in=nlmo_env%m_sigma0, &
!tiidft                  eps_filter_in=nlmo_env%eps_filter)
!tiidft               IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...pB norm: ", step_norm
!tiidft
!tiidft               ! Do not step beyond the trust radius
!tiidft               IF (step_norm .LE. radius_current) THEN
!tiidft
!tiidft                  IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Full dogleg"
!tiidft
!tiidft                  border_reached = .FALSE.
!tiidft
!tiidft                  DO ispin = 1, nspins
!tiidft                     CALL dbcsr_copy(step(ispin), m_model_Bd(ispin))
!tiidft                  ENDDO
!tiidft
!tiidft                  fake_step_size_to_report = 2.0_dp
!tiidft                  iteration_type_to_report = 6
!tiidft
!tiidft               ELSE ! take a shorter dogleg step
!tiidft
!tiidft                  IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...pB norm is too large"
!tiidft
!tiidft                  border_reached = .TRUE.
!tiidft
!tiidft                  ! compute the dogleg vector = pB - pU
!tiidft                  ! this destroys -Binv.grad content
!tiidft                  DO ispin = 1, nspins
!tiidft                     CALL dbcsr_add(m_model_Bd(ispin), step(ispin), 1.0_dp, -1.0_dp)
!tiidft                  ENDDO
!tiidft
!tiidft                  CALL step_size_to_border( &
!tiidft                     step_size_out=step_size, &
!tiidft                     metric_in=nlmo_env%m_sigma0, &
!tiidft                     position_in=step, &
!tiidft                     direction_in=m_model_Bd, &
!tiidft                     trust_radius_in=radius_current, &
!tiidft                     eps_filter_in=nlmo_env%eps_filter)
!tiidft                  IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Step size to border: ", step_size
!tiidft                  IF (step_size .GE. 1.0_dp .OR. step_size .LE. 0.0_dp) THEN
!tiidft                     CPABORT("Wrong dog-leg step, should not be here")
!tiidft                  ENDIF
!tiidft
!tiidft                  DO ispin = 1, nspins
!tiidft                     CALL dbcsr_add(step(ispin), m_model_Bd(ispin), 1.0_dp, step_size)
!tiidft                  ENDDO
!tiidft
!tiidft                  fake_step_size_to_report = 1.0_dp + step_size
!tiidft                  iteration_type_to_report = 7
!tiidft
!tiidft               ENDIF ! full or partial dogleg?
!tiidft
!tiidft               IF (debug_mode) THEN
!tiidft                  ! Compute step norm
!tiidft                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Extra norm evaluation"
!tiidft                  CALL contravariant_matrix_norm( &
!tiidft                     norm_out=step_norm, &
!tiidft                     matrix_in=step, &
!tiidft                     metric_in=nlmo_env%m_sigma0, &
!tiidft                     eps_filter_in=nlmo_env%eps_filter)
!tiidft                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Step norm: ", step_norm
!tiidft                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Current radius: ", radius_current
!tiidft               ENDIF
!tiidft
!tiidft               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
!tiidft                                                 reduction_out=expected_reduction, &
!tiidft                                                 step_in=step)
!tiidft
!tiidft               inner_loop_success = .TRUE.
!tiidft
!tiidft               t2 = m_walltime()
!tiidft               CALL fixed_r_report(unit_nr, &
!tiidft                                   iter_type=iteration_type_to_report, &
!tiidft                                   iteration=iteration, &
!tiidft                                   step_size=fake_step_size_to_report, &
!tiidft                                   border_reached=border_reached, &
!tiidft                                   curvature=y_scalar, &
!tiidft                                   grad_norm_ratio=expected_reduction, &
!tiidft                                   time=t2 - t1)
!tiidft
!tiidft               EXIT fixed_r_loop ! the inner loop
!tiidft
!tiidft            ENDIF ! Non-iterative subproblem methods are done by now
!tiidft
!tiidft            ! Step 4: update model gradient
!tiidft            DO ispin = 1, nspins
!tiidft               ! save previous data
!tiidft               CALL dbcsr_copy(m_model_r_prev(ispin), m_model_r(ispin))
!tiidft               CALL dbcsr_add(m_model_r(ispin), m_model_Bd(ispin), &
!tiidft                              1.0_dp, -step_size)
!tiidft            ENDDO ! ispin
!tiidft
!tiidft            ! Model grad norm
!tiidft            DO ispin = 1, nspins
!tiidft               CALL dbcsr_norm(m_model_r(ispin), dbcsr_norm_maxabsnorm, &
!tiidft                               norm_scalar=grad_norm_spin(ispin))
!tiidft               !grad_norm_frob = dbcsr_frobenius_norm(grad(ispin)) / &
!tiidft               !                 dbcsr_frobenius_norm(quench_t(ispin))
!tiidft            ENDDO ! ispin
!tiidft            model_grad_norm = MAXVAL(grad_norm_spin)
!tiidft
!tiidft            ! Check norm reduction
!tiidft            grad_norm_ratio = model_grad_norm/grad_norm_ref
!tiidft            IF (grad_norm_ratio .LT. optimizer%model_grad_norm_ratio) THEN
!tiidft
!tiidft               border_reached = .FALSE.
!tiidft               inner_loop_success = .TRUE.
!tiidft
!tiidft               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
!tiidft                                                 reduction_out=expected_reduction, &
!tiidft                                                 step_in=step)
!tiidft
!tiidft               t2 = m_walltime()
!tiidft               CALL fixed_r_report(unit_nr, &
!tiidft                                   iter_type=4, &
!tiidft                                   iteration=iteration, &
!tiidft                                   step_size=step_size, &
!tiidft                                   border_reached=border_reached, &
!tiidft                                   curvature=y_scalar, &
!tiidft                                   grad_norm_ratio=expected_reduction, &
!tiidft                                   time=t2 - t1)
!tiidft
!tiidft               EXIT fixed_r_loop ! the inner loop
!tiidft
!tiidft            ENDIF
!tiidft
!tiidft            ! Step 5: update model direction
!tiidft            DO ispin = 1, nspins
!tiidft               ! save previous data
!tiidft               CALL dbcsr_copy(m_model_rt_prev(ispin), m_model_rt(ispin))
!tiidft            ENDDO ! ispin
!tiidft
!tiidft            DO ispin = 1, nspins
!tiidft
!tiidft               CALL dbcsr_multiply("N", "N", 1.0_dp, &
!tiidft                                   nlmo_env%m_siginv0(ispin), &
!tiidft                                   m_model_r(ispin), &
!tiidft                                   0.0_dp, m_model_rt(ispin), &
!tiidft                                   filter_eps=nlmo_env%eps_filter)
!tiidft
!tiidft            ENDDO ! ispin
!tiidft
!tiidft            CALL compute_cg_beta( &
!tiidft               beta=beta, &
!tiidft               reset_conjugator=reset_conjugator, &
!tiidft               conjugator=optimizer%conjugator, &
!tiidft               grad=m_model_r, &
!tiidft               prev_grad=m_model_r_prev, &
!tiidft               step=m_model_rt, &
!tiidft               prev_step=m_model_rt_prev)
!tiidft
!tiidft            DO ispin = 1, nspins
!tiidft               ! update direction
!tiidft               CALL dbcsr_add(m_model_d(ispin), m_model_rt(ispin), beta, 1.0_dp)
!tiidft            ENDDO ! ispin
!tiidft
!tiidft            t2 = m_walltime()
!tiidft            CALL fixed_r_report(unit_nr, &
!tiidft                                iter_type=1, &
!tiidft                                iteration=iteration, &
!tiidft                                step_size=step_size, &
!tiidft                                border_reached=border_reached, &
!tiidft                                curvature=y_scalar, &
!tiidft                                grad_norm_ratio=grad_norm_ratio, &
!tiidft                                time=t2 - t1)
!tiidft            t1 = m_walltime()
!tiidft
!tiidft         ENDDO fixed_r_loop
!tiidft         !!!! done with the inner loop
!tiidft         ! the inner loop must return: step, predicted reduction,
!tiidft         ! whether it reached the border and completed successfully
!tiidft
!tiidft         IF (.NOT. inner_loop_success) THEN
!tiidft            CPABORT("Inner loop did not produce solution")
!tiidft         ENDIF
!tiidft
!tiidft         ! Replace theta with the trial variables and re-evaluate the env
!tiidft         DO ispin = 1, nspins
!tiidft            CALL dbcsr_copy(m_theta_saved(ispin), nlmo_env%m_theta(ispin))
!tiidft            CALL dbcsr_add(nlmo_env%m_theta(ispin), step(ispin), 1.0_dp, 1.0_dp)
!tiidft         ENDDO ! ispin
!tiidft         CALL nlmo_env_mainvar_to_aux(nlmo_env)
!tiidft         loss_components = 0.0_dp
!tiidft         CALL nlmo_env_loss_function(nlmo_env, loss_components, overlap_determinant)
!tiidft         energy_trial = loss_components(1)
!tiidft         penalty_trial = loss_components(2)
!tiidft         loss_trial = SUM(loss_components)
!tiidft
!tiidft         ! analyze performance of the new set of varibles
!tiidft         rho = (loss_trial - loss_start)/expected_reduction
!tiidft         loss_change_to_report = loss_trial - loss_start
!tiidft
!tiidft         IF (rho < 0.25_dp) THEN
!tiidft            radius_current = 0.25_dp*radius_current
!tiidft         ELSE
!tiidft            IF (rho > 0.75_dp .AND. border_reached) THEN
!tiidft               radius_current = MIN(2.0_dp*radius_current, radius_max)
!tiidft            ENDIF
!tiidft         ENDIF ! radius adjustment
!tiidft
!tiidft         IF (rho > eta) THEN
!tiidft            ! The variable is inconsistent with gradient and hessian
!tiidft            ! It is consistent with aux matrices and loss
!tiidft            same_position = .FALSE.
!tiidft            loss_start = loss_trial
!tiidft            energy_start = energy_trial
!tiidft            penalty_start = penalty_trial
!tiidft         ELSE
!tiidft            same_position = .TRUE.
!tiidft            ! The variable is inconsistent with aux matrices and loss
!tiidft            ! It is consistent with gradient and hessian
!tiidft            DO ispin = 1, nspins
!tiidft               CALL dbcsr_copy(nlmo_env%m_theta(ispin), m_theta_saved(ispin))
!tiidft            ENDDO ! ispin
!tiidft         ENDIF ! finalize step
!tiidft
!tiidft         t2outer = m_walltime()
!tiidft         CALL trust_r_report(unit_nr, &
!tiidft                             iter_type=2, &
!tiidft                             iteration=outer_iteration, &
!tiidft                             loss=loss_trial, &
!tiidft                             delta_loss=loss_change_to_report, &
!tiidft                             grad_norm=0.0_dp, &
!tiidft                             predicted_reduction=expected_reduction, &
!tiidft                             rho=rho, &
!tiidft                             radius=radius_current, &
!tiidft                             new=.NOT. same_position, &
!tiidft                             time=t2outer - t1outer)
!tiidft         t1outer = m_walltime()
!tiidft
!tiidft      ENDDO adjust_r_loop
!tiidft
!tiidft      DO ispin = 1, nspins
!tiidft         CALL dbcsr_release(prev_step(ispin))
!tiidft         CALL dbcsr_release(step(ispin))
!tiidft         CALL dbcsr_release(m_model_r(ispin))
!tiidft         CALL dbcsr_release(m_model_rt(ispin))
!tiidft         CALL dbcsr_release(m_model_d(ispin))
!tiidft         CALL dbcsr_release(m_model_Bd(ispin))
!tiidft         CALL dbcsr_release(m_model_r_prev(ispin))
!tiidft         CALL dbcsr_release(m_model_rt_prev(ispin))
!tiidft         CALL dbcsr_release(m_theta_saved(ispin))
!tiidft      ENDDO ! ispin
!tiidft
!tiidft      DEALLOCATE (prev_step)
!tiidft      DEALLOCATE (step)
!tiidft      DEALLOCATE (m_model_r)
!tiidft      DEALLOCATE (m_model_rt)
!tiidft      DEALLOCATE (m_model_d)
!tiidft      DEALLOCATE (m_model_Bd)
!tiidft      DEALLOCATE (m_model_r_prev)
!tiidft      DEALLOCATE (m_model_rt_prev)
!tiidft      DEALLOCATE (m_theta_saved)
!tiidft
!tiidft      DEALLOCATE (grad_norm_spin)
!tiidft
!tiidft      DO ispin = 1, nspins
!tiidft         CALL dbcsr_copy(matrix_mo_inout(ispin), nlmo_env%m_mo(ispin))
!tiidft      ENDDO ! ispin
!tiidft
!tiidft      CALL nlmo_env_deallocate_all(nlmo_env)
!tiidft
!tiidft      IF (.NOT. scf_converged) THEN
!tiidft         CPABORT("Optimization not converged! ")
!tiidft      ENDIF
!tiidft
!tiidft      CALL timestop(handle)
!tiidft
!tiidft   END SUBROUTINE construct_nlmos_trustr

END MODULE tidft_optimizer_methods

