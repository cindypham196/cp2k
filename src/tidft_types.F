!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for time-independent DFT methods
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_types
   USE almo_scf_types,                  ONLY: optimizer_options_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_release,&
                                              dbcsr_type
   USE domain_submatrix_types,          ONLY: domain_map_type,&
                                              domain_submatrix_type
   USE input_constants,                 ONLY: &
        cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, cg_hager_zhang, cg_hestenes_stiefel, &
        cg_liu_storey, cg_polak_ribiere, cg_zero, optimizer_diis, optimizer_pcg, optimizer_trustr, &
        trustr_cauchy, trustr_dogleg, trustr_steihaug, xalmo_prec_domain, xalmo_prec_full, &
        xalmo_prec_zero
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_types'

   PUBLIC :: tidft_env_type, &
             tidft_env_release

   ! container for penalty terms that are added to the energy functional
   TYPE penalty_type

      REAL(KIND=dp)      :: final_determinant, penalty_strength, &
                            determinant_tolerance, penalty_strength_dec_factor, &
                            compactification_filter_start

   END TYPE penalty_type

   ! This structure describes a single Slater determinant
   TYPE slater_det_type

      ! The number of spins can be:
      ! 1 (closed-shell), 2 (open-shell), 3 (restricted open-shell)
      INTEGER                                            :: nspins, nelectrons_total, nbasis
      ! Electrons of different spin (spin)
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nelectrons_spin

      ! Occupied MO coeff NxO (spinJ)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: matrix_t

   END TYPE slater_det_type

   ! This structure decribes an electronic state
   TYPE e_state_type

      ! Number of Slater determinants describing the electronic state
      INTEGER                                            :: nslater, multiplicity
      TYPE(slater_det_type), DIMENSION(:), ALLOCATABLE   :: slater_det

   END TYPE e_state_type

   ! the structure contains general info about the system
   TYPE tidft_env_type

      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      ! Electronic states
      ! Number of multiplicities
      INTEGER                                            :: nmultips
      ! Number of states of each multiplictiy
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nstates
      ! Descriptors of electronic states (example: singlets, triplets, ...)
      TYPE(e_state_type), DIMENSION(:), ALLOCATABLE      :: e_states

      REAL(KIND=dp)                                      :: eps_filter
      INTEGER                                            :: order_lanczos
      REAL(KIND=dp)                                      :: eps_lanczos
      INTEGER                                            :: max_iter_lanczos

      ! forces
      LOGICAL                                            :: calc_forces, &
                                                            orthogonal_basis

      ! Options for various optimizers
      TYPE(optimizer_options_type)                       :: opt_pcg
      TYPE(optimizer_options_type)                       :: opt_trustr

      !!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! MATRICES !!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!

      ! Basis fucntions overlap NxN
      TYPE(dbcsr_type)                                   :: matrix_s
      ! Kohn-Sham Hamiltonian NxN (slater, spin)
      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: matrix_ks

      ! MO overlap OxO (stateI, stateJ)
      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: matrix_sigma

      ! Overlap of electronic states (stateI, stateJ, multiplicity)
      ! RZK-warning better structure is required
      REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE       :: matrix_phi

   END TYPE tidft_env_type

CONTAINS

!tiidft! **************************************************************************************************
!tiidft!> \brief Prints out the options of an optimizer
!tiidft!> \param optimizer   options to print
!tiidft!> \param unit_nr   output stream
!tiidft!> \par History
!tiidft!>       2014.10 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE print_optimizer_options(optimizer, unit_nr)
!tiidft
!tiidft      TYPE(optimizer_options_type), INTENT(IN)           :: optimizer
!tiidft      INTEGER, INTENT(IN)                                :: unit_nr
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'print_optimizer_options', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      CHARACTER(33)                                      :: conj_string, prec_string, type_string
!tiidft
!tiidft      IF (unit_nr .GT. 0) THEN
!tiidft
!tiidft         SELECT CASE (optimizer%optimizer_type)
!tiidft         CASE (optimizer_diis)
!tiidft            type_string = "DIIS"
!tiidft         CASE (optimizer_pcg)
!tiidft            type_string = "PCG"
!tiidft         CASE (optimizer_trustr)
!tiidft            type_string = "TRUST REGION"
!tiidft         END SELECT
!tiidft
!tiidft         WRITE (unit_nr, '(T4,A,T48,A33)') "optimizer type:", TRIM(type_string)
!tiidft         WRITE (unit_nr, '(T4,A,T48,I33)') "maximum iterations:", optimizer%max_iter
!tiidft         WRITE (unit_nr, '(T4,A,T48,E33.3)') "target error:", optimizer%eps_error
!tiidft
!tiidft         IF (optimizer%optimizer_type .EQ. optimizer_diis) THEN
!tiidft
!tiidft            WRITE (unit_nr, '(T4,A,T48,I33)') "maximum DIIS history:", optimizer%ndiis
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft         IF (optimizer%optimizer_type .EQ. optimizer_trustr .OR. &
!tiidft             optimizer%optimizer_type .EQ. optimizer_pcg) THEN
!tiidft
!tiidft            WRITE (unit_nr, '(T4,A,T48,I33)') "maximum outer loop iterations:", &
!tiidft               optimizer%max_iter_outer_loop
!tiidft
!tiidft            SELECT CASE (optimizer%preconditioner)
!tiidft            CASE (xalmo_prec_zero)
!tiidft               prec_string = "NONE"
!tiidft            CASE (xalmo_prec_domain)
!tiidft               prec_string = "0.5 KS + 0.5 S, DOMAINS"
!tiidft            CASE (xalmo_prec_full)
!tiidft               prec_string = "0.5 KS + 0.5 S, FULL"
!tiidft            END SELECT
!tiidft            WRITE (unit_nr, '(T4,A,T48,A33)') "preconditioner:", TRIM(prec_string)
!tiidft
!tiidft            SELECT CASE (optimizer%conjugator)
!tiidft            CASE (cg_zero)
!tiidft               conj_string = "Steepest descent"
!tiidft            CASE (cg_polak_ribiere)
!tiidft               conj_string = "Polak-Ribiere"
!tiidft            CASE (cg_fletcher_reeves)
!tiidft               conj_string = "Fletcher-Reeves"
!tiidft            CASE (cg_hestenes_stiefel)
!tiidft               conj_string = "Hestenes-Stiefel"
!tiidft            CASE (cg_fletcher)
!tiidft               conj_string = "Fletcher"
!tiidft            CASE (cg_liu_storey)
!tiidft               conj_string = "Liu-Storey"
!tiidft            CASE (cg_dai_yuan)
!tiidft               conj_string = "Dai-Yuan"
!tiidft            CASE (cg_hager_zhang)
!tiidft               conj_string = "Hager-Zhang"
!tiidft            END SELECT
!tiidft            WRITE (unit_nr, '(T4,A,T48,A33)') "conjugator:", TRIM(conj_string)
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft         IF (optimizer%optimizer_type .EQ. optimizer_pcg) THEN
!tiidft
!tiidft            WRITE (unit_nr, '(T4,A,T48,E33.3)') "line search step size guess:", &
!tiidft               optimizer%lin_search_step_size_guess
!tiidft            WRITE (unit_nr, '(T4,A,T48,E33.3)') "line search target error:", &
!tiidft               optimizer%lin_search_eps_error
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft         IF (optimizer%optimizer_type .EQ. optimizer_trustr) THEN
!tiidft
!tiidft            SELECT CASE (optimizer%trustr_algorithm)
!tiidft            CASE (trustr_steihaug)
!tiidft               conj_string = "Steihaug's CG"
!tiidft            CASE (trustr_cauchy)
!tiidft               conj_string = "Cauchy point"
!tiidft            CASE (trustr_dogleg)
!tiidft               conj_string = "Dogleg"
!tiidft            END SELECT
!tiidft            WRITE (unit_nr, '(T4,A,T48,A33)') "Subproblem algorithm:", TRIM(conj_string)
!tiidft
!tiidft            WRITE (unit_nr, '(T4,A,T48,E33.3)') "gradient decrease accepted:", &
!tiidft               optimizer%model_grad_norm_ratio
!tiidft            WRITE (unit_nr, '(T4,A,T48,E33.3)') "initial trust radius:", &
!tiidft               optimizer%initial_trust_radius
!tiidft            WRITE (unit_nr, '(T4,A,T48,E33.3)') "max trust radius:", &
!tiidft               optimizer%max_trust_radius
!tiidft            WRITE (unit_nr, '(T4,A,T48,E33.3)') "rho of no update lies between .0 and .25:", &
!tiidft               optimizer%rho_do_not_update
!tiidft
!tiidft         ENDIF
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft   END SUBROUTINE print_optimizer_options

! **************************************************************************************************
!> \brief release tidft envirnoment
!> \param tidft_env ...
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_env_release(tidft_env)
      TYPE(tidft_env_type), POINTER                   :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istore

      CALL timeset(routineN, handle)

!tiidft      ! delete history
!tiidft      DO ispin = 1, SIZE(almo_scf_env%almo_history%matrix_t)
!tiidft         DO istore = 1, MIN(almo_scf_env%almo_history%istore, almo_scf_env%almo_history%nstore)
!tiidft            CALL dbcsr_release(almo_scf_env%almo_history%matrix_p_up_down(ispin, istore))
!tiidft         ENDDO
!tiidft         IF (almo_scf_env%almo_history%istore .GT. 0) &
!tiidft            CALL dbcsr_release(almo_scf_env%almo_history%matrix_t(ispin))
!tiidft      ENDDO
!tiidft      DEALLOCATE (almo_scf_env%almo_history%matrix_p_up_down)
!tiidft      DEALLOCATE (almo_scf_env%almo_history%matrix_t)
!tiidft      ! delete xalmo history
!tiidft      DO ispin = 1, SIZE(almo_scf_env%xalmo_history%matrix_t)
!tiidft         DO istore = 1, MIN(almo_scf_env%xalmo_history%istore, almo_scf_env%xalmo_history%nstore)
!tiidft            CALL dbcsr_release(almo_scf_env%xalmo_history%matrix_p_up_down(ispin, istore))
!tiidft            !CALL dbcsr_release(almo_scf_env%xalmo_history%matrix_x(ispin, istore))
!tiidft         ENDDO
!tiidft         IF (almo_scf_env%xalmo_history%istore .GT. 0) &
!tiidft            CALL dbcsr_release(almo_scf_env%xalmo_history%matrix_t(ispin))
!tiidft      ENDDO
!tiidft      DEALLOCATE (almo_scf_env%xalmo_history%matrix_p_up_down)
!tiidft      !DEALLOCATE (almo_scf_env%xalmo_history%matrix_x)
!tiidft      DEALLOCATE (almo_scf_env%xalmo_history%matrix_t)

      DEALLOCATE (tidft_env)

      CALL timestop(handle)

   END SUBROUTINE tidft_env_release

END MODULE tidft_types

