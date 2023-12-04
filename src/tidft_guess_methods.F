!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Variable-metric TIDFT: methods to generate initial guess
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_guess_methods
   USE tidft_guess_types,               ONLY: tidft_excitation_list_type, &
                                              tidft_excitation_type, &
                                              orbital_set_type
   USE almo_scf_methods,                ONLY: almo_scf_t_to_proj,&
                                              get_overlap,&
                                              orthogonalize_mos
   USE almo_scf_types,                  ONLY: optimizer_options_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_filter, dbcsr_multiply, dbcsr_norm, &
        dbcsr_norm_maxabsnorm, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry, dbcsr_distribution_type, dbcsr_get_info, dbcsr_get_diag, &
        dbcsr_transposed, dbcsr_print, dbcsr_set_diag, dbcsr_dot, dbcsr_desymmetrize, &
        dbcsr_init_random
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
        cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero
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
                                              gradVAR_scale_split_increment, &
                                              precMO_to_precVAR_increment
   USE tidft_qs,                        ONLY: tidft_dm_to_tidft_ks
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_guess_methods'

   PUBLIC :: construct_canonical_orbitals, &
             guess_excited_states_orbitals

   LOGICAL, PARAMETER :: debug_mode = .TRUE.
   LOGICAL, PARAMETER :: safe_mode = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Constructs initial guess for molecular oribatals of excited states
!> \param tidft_env ...
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE guess_excited_states_orbitals(tidft_env)

      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'guess_excited_states_orbitals', &
         routineP = moduleN//':'//routineN

      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(tidft_excitation_list_type)                   :: excitation_list
      INTEGER                                            :: unit_nr, nspins, &
                                                            ispin, n, ivirtual, &
                                                            idet, norbs, iorb, &
                                                            count_unpairings, &
                                                            nalpha, nbeta
      INTEGER                                            :: nexcitations, &
                                                            iexcitation, &
                                                            generated_nstates, &
                                                            istate, state_index
      INTEGER                                            :: max_order
      INTEGER, DIMENSION(:,:), ALLOCATABLE               :: start_orbital_indices
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: index_true
      LOGICAL                                            :: check_unpairings
      REAL(KIND=dp)                                      :: eps_diff_orb_energy
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: delta, energies_temp
      TYPE(orbital_set_type)                             :: orbital_set
      TYPE(tidft_excitation_type)                        :: current_excitation

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! RZK TODO:
      !   * Create an input section within &TIDFT called &WF_GUESS
      !     It should contain controls for this procedure:
      !       MAX_EXC_ENERGY - maximum allowed excitation energy for the guess (pre-SCF)
      !       MIN_EXC_ENERGY - minimum allowed excitation energy for the guess (pre-SCF)
      !       MAX_EXC_ORDER  - maximum allowed number of excitations from the ground state
      !     The corresponding variables for these input keywords are:
      !       tidft_env%guess_max_excitation_energy (to be created)
      !       tidft_env%guess_min_excitation_energy (already implemented)
      !       tidft_env%guess_max_excitation_order (to be created)
      !     More controls will be added later
      !   * Test this procedure for UKS and RKS for He crystal (fast, login node) to see whether
      !     excitation pairs are correctly generated. That is, the test run should only print
      !     occupied-virtual pairs generated for each excited state guess. It is not necessary to
      !     create the guess orbitals yet.
      !   * NExt test is to do the same for ROKS singlet and ROKS triplet (with triplet correction FALSE).
      !   * Next text: same system ROKS singlet with triplet correction TRUE.
      !   * Repeat the test for another simple system (fast, login node). F3B..NH3 with small cell can be used.

      ! MAJOR STEP: Obtain canonical orbital energies and info on orbital location in matrices
      ! RZK-WARNING: make mainvar consistent with canonical orbitals - alternatively:
      ! Do not modify the orbitals in tidft_env, only get the energy (orbitals are needed later here though)
      nspins = tidft_env%nspins
      idet = 1 ! Use the "mixed" state to produce canonical orbitals

      ! The first index is spin, the second is occ/virt
      ALLOCATE(orbital_set%orbital_list(nspins, 2))
      DO ispin = 1, nspins

         DO ivirtual = 1, 2

            IF (ivirtual .EQ. 1) THEN
               CALL dbcsr_get_info(tidft_env%mat_sigma_inv(1,idet,ispin), nfullrows_total=norbs)
            ELSE
               CALL dbcsr_get_info(tidft_env%mat_sigma_vv(1,idet,ispin), nfullrows_total=norbs)
            ENDIF

            orbital_set%orbital_list(ispin, ivirtual)%norbs = norbs
            orbital_set%orbital_list(ispin, ivirtual)%ispin = ispin
            orbital_set%orbital_list(ispin, ivirtual)%ivirtual = ivirtual - 1

            IF (norbs .GT. 0) THEN
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%orb_energy(norbs))
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%index_within_matrix(norbs))
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%block_within_matrix(norbs))
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%index_within_block(norbs))
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%orb_has_spin_double(norbs))
               IF (nspins .EQ. 1) THEN ! RKS for all states
                  orbital_set%orbital_list(ispin, ivirtual)%orb_has_spin_double(:) = .TRUE.
               ELSE ! UKS and ROKS
                  orbital_set%orbital_list(ispin, ivirtual)%orb_has_spin_double(:) = .FALSE.
               ENDIF
            ENDIF

         ENDDO

      ENDDO

      CALL construct_canonical_orbitals(tidft_env, estate = 1, orbital_set = orbital_set)

      ! MAJOR STEP: Find out what orbitals are doubly and singly occupied
      ! For virtual orbitals, find out whether they CAN contain two or one electron(s)
      ! RZK-WARNING: This primitive algorithm only compares energies of alpha-beta orbitals
      ! In the future, norms of Delta T = T_alpha - T_beta should be analyzed

      ! This procedure is irrelevant for RKS and UKS because enforcing
      ! the correct number of unpairings is not required.
      ! Perform only for the case of ROKS
      count_unpairings = -1
      check_unpairings = tidft_env%restricted .AND. tidft_env%unpairings.GT.0
      IF (check_unpairings) THEN

         eps_diff_orb_energy = SQRT(tidft_env%eps_filter)

         count_unpairings = 0
         DO ivirtual = 1, 2

            nalpha = orbital_set%orbital_list(1, ivirtual)%norbs
            nbeta = orbital_set%orbital_list(2, ivirtual)%norbs
            ! RZK-WARNING: Implement cases with unequal number of electrons later
            CPASSERT(nalpha.EQ.nbeta)

            ALLOCATE(delta(nalpha))
            delta(:) = ABS(orbital_set%orbital_list(1, ivirtual)%orb_energy(:) - &
                       orbital_set%orbital_list(2, ivirtual)%orb_energy(:))
            DO iorb = 1, nalpha
               IF ( delta(iorb) .LT. eps_diff_orb_energy ) THEN
                  orbital_set%orbital_list(1, ivirtual)%orb_has_spin_double(iorb) = .TRUE.
                  orbital_set%orbital_list(2, ivirtual)%orb_has_spin_double(iorb) = .TRUE.
               ELSE
                  orbital_set%orbital_list(1, ivirtual)%orb_has_spin_double(iorb) = .FALSE.
                  orbital_set%orbital_list(2, ivirtual)%orb_has_spin_double(iorb) = .FALSE.
                  IF (ivirtual .EQ. 1) count_unpairings = count_unpairings + 1
               ENDIF
            ENDDO ! iorb
            DEALLOCATE(delta)

         ENDDO ! ivirtual

         IF (count_unpairings .GT. tidft_env%unpairings) THEN
            CPABORT("WF analysis found too many unpaired electrons")
         ENDIF

         ! RZK-WARNING: for triple corrction number of unpairings must be exctly equal to 1
         ! Here and in excitation generation

      ENDIF ! ROKS

      ! MAJOR STEP: Select occupied and virtual orbital from the set using various conditions
      ! RZK-WARNING: NYI

      ! MAJOR STEP: Loop over select occupied and virtual orbitals, keep desired in the list
      ALLOCATE(start_orbital_indices(nspins, 2))
      DO ispin = 1, nspins
         ! occupied are looped from top to bottom
         start_orbital_indices(ispin, 1) = &
            orbital_set%orbital_list(ispin, 1)%norbs
         ! virtuals are looped from bottom to top
         start_orbital_indices(ispin, 2) = 1
      ENDDO ! ispin

      ! Init storage structures
      nexcitations = tidft_env%nstates - 1
      max_order = tidft_env%guess_max_excitation_order
      CPASSERT(nexcitations.GT.0)
      CPASSERT(max_order.GT.0)
      excitation_list%max_nexcitations = nexcitations
      excitation_list%max_order = max_order
      excitation_list%cur_nexcitations = 0
      ALLOCATE(excitation_list%excitations(nexcitations))

      DO iexcitation = 1, nexcitations
         excitation_list%excitations(iexcitation)%order = -1
         excitation_list%excitations(iexcitation)%excitation_energy = 0.0_dp
         ALLOCATE(excitation_list%excitations(iexcitation)%delta_epsilon(max_order))
         ALLOCATE(excitation_list%excitations(iexcitation)%occupied(max_order))
         ALLOCATE(excitation_list%excitations(iexcitation)%virtual(max_order))
         ALLOCATE(excitation_list%excitations(iexcitation)%spin(max_order))
         excitation_list%excitations(iexcitation)%delta_epsilon(:) = 0.0_dp
         excitation_list%excitations(iexcitation)%occupied(:) = -1
         excitation_list%excitations(iexcitation)%virtual(:) = -1
         excitation_list%excitations(iexcitation)%spin(:) = -1
      ENDDO

      ALLOCATE(current_excitation%delta_epsilon(max_order))
      ALLOCATE(current_excitation%occupied(max_order))
      ALLOCATE(current_excitation%virtual(max_order))
      ALLOCATE(current_excitation%spin(max_order))
      current_excitation%order = 0
      current_excitation%excitation_energy = 0.0_dp
      current_excitation%unpairings_count = count_unpairings
      current_excitation%delta_epsilon(:) = 0.0_dp
      current_excitation%occupied(:) = -1
      current_excitation%virtual(:) = -1
      current_excitation%spin(:) = -1

      ! recursive subroutine for generaing multiple single excitations
      CALL generate_single_excitations( &
              current_order = 1, &
              excitation_list = excitation_list, &
              current_excitation = current_excitation, &
              orbital_set = orbital_set, &
              start_orbital_indices = start_orbital_indices, &
              min_excitation_energy = tidft_env%guess_min_excitation_energy, &
              check_unpairings = check_unpairings, &
              max_allowed_unpairings = tidft_env%unpairings)

      DEALLOCATE(current_excitation%delta_epsilon)
      DEALLOCATE(current_excitation%occupied)
      DEALLOCATE(current_excitation%virtual)
      DEALLOCATE(current_excitation%spin)

      ! MAJOR STEP: Order the list of excitations
      generated_nstates = excitation_list%cur_nexcitations
      ALLOCATE(index_true(generated_nstates))
      ALLOCATE(energies_temp(generated_nstates))

      DO istate = 1, generated_nstates
         energies_temp(istate) = excitation_list%excitations(istate)%excitation_energy
      ENDDO
      CALL sort(energies_temp, generated_nstates, index_true)
      !indexGS = index_true(1)

      ! MAJOR STEP: Swap mainvar matrix columns
      ! RZK-WARNING: Change tidft_env%nstates
      DO istate = 1, generated_nstates

         state_index = index_true(istate)
         !Process data from excitation_list%excitations(state_index) to swap columns
         !Convert orbital numbers to mainvar column numbers
         !Swap

      ENDDO

      ! MAJOR STEP: Clean up
      DEALLOCATE(index_true)
      DEALLOCATE(energies_temp)

      DO iexcitation = 1, nexcitations
         DEALLOCATE(excitation_list%excitations(iexcitation)%delta_epsilon)
         DEALLOCATE(excitation_list%excitations(iexcitation)%occupied)
         DEALLOCATE(excitation_list%excitations(iexcitation)%virtual)
         DEALLOCATE(excitation_list%excitations(iexcitation)%spin)
      ENDDO

      DEALLOCATE(excitation_list%excitations)

      DO ispin = 1, nspins

         DO ivirtual = 1, 2

            IF (unit_nr > 0) THEN
               WRITE(unit_nr, *) "Orbital energies for spin, virt: ", ispin, ivirtual
               WRITE(unit_nr, *) orbital_set%orbital_list(ispin, ivirtual)%orb_energy(:)
            ENDIF

            IF ( orbital_set%orbital_list(ispin, ivirtual)%norbs .GT. 0 ) THEN
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%orb_energy)
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%index_within_matrix)
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%block_within_matrix)
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%index_within_block)
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%orb_has_spin_double)
            ENDIF

         ENDDO

      ENDDO
      DEALLOCATE(start_orbital_indices)
      DEALLOCATE(orbital_set%orbital_list)

      CPABORT("EXCITED STATE GUESS WF HAS BEEN GENERATED")

   END SUBROUTINE guess_excited_states_orbitals

   RECURSIVE SUBROUTINE generate_single_excitations(current_order, excitation_list, &
      current_excitation, orbital_set, start_orbital_indices, min_excitation_energy, &
      check_unpairings, max_allowed_unpairings)

      INTEGER, INTENT(IN)                                :: current_order
      TYPE(tidft_excitation_list_type), INTENT(INOUT)    :: excitation_list
      TYPE(tidft_excitation_type), INTENT(IN)            :: current_excitation
      TYPE(orbital_set_type), INTENT(IN)                 :: orbital_set
      INTEGER, DIMENSION(:,:), INTENT(IN)                :: start_orbital_indices
      REAL(KIND=dp), INTENT(IN)                          :: min_excitation_energy
      LOGICAL, INTENT(IN)                                :: check_unpairings
      INTEGER, INTENT(IN)                                :: max_allowed_unpairings

      INTEGER                                            :: nspins, max_order, iocc, &
                                                            ispin, opposite_spin, &
                                                            nvirtuals, avirt, &
                                                            incr_unpairings, &
                                                            iexcitation, &
                                                            insertion_index
      INTEGER, DIMENSION(:,:), ALLOCATABLE               :: next_level_start_indices
      TYPE(tidft_excitation_type)                        :: next_level_excitation
      REAL(KIND=dp)                                      :: spin_factor, max_energy
      LOGICAL                                            :: add_to_excitation_list, found, &
                                                            originates_from_doubly_occupied, &
                                                            goes_to_doubly_occupied_virtual, &
                                                            max_order_condition, &
                                                            deeper_level

      nspins = SIZE(orbital_set%orbital_list, 1)
      max_order = excitation_list%max_order

      ALLOCATE(next_level_start_indices(nspins, 2))

      ALLOCATE(next_level_excitation%delta_epsilon(max_order))
      ALLOCATE(next_level_excitation%occupied(max_order))
      ALLOCATE(next_level_excitation%virtual(max_order))
      ALLOCATE(next_level_excitation%spin(max_order))

      ! copy indices
      next_level_start_indices(:,:) = start_orbital_indices(:,:)

      IF (nspins .EQ. 2) THEN
         spin_factor = 1.0_dp
      ELSE
         spin_factor = 2.0_dp
      ENDIF

      DO ispin = 1, nspins

         opposite_spin = 1
         IF (nspins .EQ. 2) THEN
            IF (ispin .EQ. 1) THEN
               opposite_spin = 2
            ELSE
               opposite_spin = 1
            ENDIF
         ENDIF

         next_level_start_indices(opposite_spin, :) = &
            start_orbital_indices(opposite_spin, :)

         DO iocc = start_orbital_indices(ispin, 1), 1, -1

            nvirtuals = orbital_set%orbital_list(ispin, 2)%norbs

            DO avirt = start_orbital_indices(ispin, 2), nvirtuals, +1

               ! copy current excitation into the next level excitation
               CALL copy_excitation(dest = next_level_excitation, &
                       source = current_excitation)

               ! update the next level excitation
               next_level_excitation%order = current_order
               next_level_excitation%spin(current_order) = ispin
               next_level_excitation%occupied(current_order) = iocc
               next_level_excitation%virtual(current_order) = avirt
               next_level_excitation%delta_epsilon(current_order) = &
                  spin_factor * ( &
                  orbital_set%orbital_list(ispin, 2)%orb_energy(avirt) - &
                  orbital_set%orbital_list(ispin, 1)%orb_energy(iocc) )
               next_level_excitation%excitation_energy = &
                  next_level_excitation%excitation_energy + &
                  next_level_excitation%delta_epsilon(current_order)

               ! copy the newly created excitation to the list if it satisfies:
               !   energy filter requirememnts (delE should be more than the min excitation energy, no max)
               !   unpairings no. requirement (ROKS)
               !   RZK_WARNING: NYI - spin symmetry requirement (in some cases alpha excitations are equivalent to beta)
               add_to_excitation_list = .TRUE.

               ! energy filter requirement
               IF (SUM(next_level_excitation%delta_epsilon) .LT. min_excitation_energy) THEN
                  add_to_excitation_list = .FALSE.
               ENDIF

               ! No need to check if the previous ADD criteria failed
               IF (check_unpairings .AND. add_to_excitation_list) THEN

                  incr_unpairings = 0

                  ! Check if the excitation is from a doubly-occupied orbital
                  originates_from_doubly_occupied = .FALSE.
                  IF (orbital_set%orbital_list(ispin, 1)%orb_has_spin_double(iocc)) THEN
                     ! This is an occupied orbital that can hold two electrons
                     ! Is there another electron - of the opposite spin - in the orbital?
                     ! Check if the other electron is already excited
                     found = .FALSE.
                     DO iexcitation = 1, current_order - 1
                        IF (current_excitation%occupied(iexcitation) .EQ. iocc .AND. &
                            current_excitation%spin(iexcitation) .EQ. opposite_spin ) &
                            found = .TRUE.
                     ENDDO
                     IF (.NOT. found) originates_from_doubly_occupied = .TRUE.
                  ENDIF
                  IF (originates_from_doubly_occupied) &
                     incr_unpairings = incr_unpairings + 1

                  ! Check if the excitation is to a doubly-occupied orbital
                  goes_to_doubly_occupied_virtual = .FALSE.
                  IF (orbital_set%orbital_list(ispin, 2)%orb_has_spin_double(avirt)) THEN
                     ! This is a virtual orbital that can hold two electrons
                     ! Is there another electron - of the opposite spin - already in the orbital?
                     found = .FALSE.
                     DO iexcitation = 1, current_order - 1
                        IF (current_excitation%virtual(iexcitation) .EQ. avirt .AND. &
                            current_excitation%spin(iexcitation) .EQ. opposite_spin ) &
                            found = .TRUE.
                     ENDDO
                     IF (found) goes_to_doubly_occupied_virtual = .TRUE.
                  ENDIF
                  IF (goes_to_doubly_occupied_virtual) &
                     incr_unpairings = incr_unpairings - 1

                  next_level_excitation%unpairings_count = &
                     next_level_excitation%unpairings_count + incr_unpairings

                  IF (next_level_excitation%unpairings_count .GT. max_allowed_unpairings) &
                     add_to_excitation_list = .FALSE.

               ENDIF ! check unpairings

               IF (add_to_excitation_list) THEN

                  insertion_index = -1

                  ! If there are empty slots add to the next empty slot
                  IF (excitation_list%cur_nexcitations .LT. excitation_list%max_nexcitations) THEN
                     insertion_index = excitation_list%cur_nexcitations + 1
                     excitation_list%cur_nexcitations = excitation_list%cur_nexcitations + 1
                  ELSE
                  ! Find the location with the maximum excitation energy
                     max_energy = 0.0_dp
                     DO iexcitation = 1, excitation_list%cur_nexcitations
                        IF (excitation_list%excitations(iexcitation)%excitation_energy .GT. &
                           max_energy) THEN
                           insertion_index = iexcitation
                           max_energy = excitation_list%excitations(iexcitation)%excitation_energy
                        ENDIF
                     ENDDO ! excitation loop
                     IF (next_level_excitation%excitation_energy .GT. max_energy) &
                        insertion_index = -1
                  ENDIF

                  ! Copy excitation to the list
                  IF (insertion_index .GT. 0) THEN
                     CALL copy_excitation(dest = excitation_list%excitations(insertion_index), &
                             source = next_level_excitation)
                  ENDIF

               ENDIF

               ! prepare for the next level excitation
               next_level_start_indices(ispin, 1) = iocc - 1
               next_level_start_indices(ispin, 2) = avirt + 1

               max_order_condition = current_order .LT. max_order
               deeper_level = max_order_condition

               IF (deeper_level) THEN

                  CALL generate_single_excitations( &
                          current_order = current_order + 1, &
                          excitation_list = excitation_list, &
                          current_excitation = next_level_excitation, &
                          orbital_set = orbital_set, &
                          start_orbital_indices = next_level_start_indices, &
                          min_excitation_energy = min_excitation_energy, &
                          check_unpairings = check_unpairings, &
                          max_allowed_unpairings = max_allowed_unpairings)


               ENDIF

            ENDDO ! a virt

         ENDDO ! iocc

      ENDDO ! ispin

      DEALLOCATE(next_level_excitation%delta_epsilon)
      DEALLOCATE(next_level_excitation%occupied)
      DEALLOCATE(next_level_excitation%virtual)
      DEALLOCATE(next_level_excitation%spin)

      DEALLOCATE(next_level_start_indices)

   END SUBROUTINE generate_single_excitations

! **************************************************************************************************
!> \brief Constructs canonical orbitals. They are stored in mat_t of the environment
!> \param tidft_env ...
!> \param estate - electronic state for which the orbitals are constructed
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_canonical_orbitals(tidft_env, estate, orbital_set)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env
      INTEGER, INTENT(IN)                             :: estate
      TYPE(orbital_set_type), INTENT(INOUT)           :: orbital_set

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_canonical_orbitals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr, ispin, &
                                                            idet, ivirtual, &
                                                            norbs, ndets, nspins, &
                                                            nmos_stateI
      TYPE(dbcsr_type)                                   :: tempNV1, tempVOcc1, &
                                                            tempVOcc2, tempVV1, &
                                                            tempVV2, tempNO1, &
                                                            tempOO1, tempOO2
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: fake_e_per_block
      TYPE(cp_logger_type), POINTER                      :: logger

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ALLOCATE(fake_e_per_block(1))
      fake_e_per_block(:) = -1

      nspins = tidft_env%nspins
      ndets=1
      IF (tidft_env%triplet_correction) ndets = 2
      idet = 1 ! Use the singlet state to produce canonical orbitals

      DO ispin = 1, nspins

         nmos_stateI = tidft_env%e_states(estate)%slater_det(idet)%nelectrons_ab(ispin)

         IF ( nmos_stateI .GT. 0) THEN

            CALL dbcsr_create(tempNO1, &
                              template=tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempOO1, &
                              template=tidft_env%mat_sigma(estate,estate,idet,ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempOO2, &
                              template=tidft_env%mat_sigma(estate,estate,idet,ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            CALL orthogonalize_mos(ket=tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                   overlap=tidft_env%mat_sigma(estate,estate,idet,ispin), &
                                   metric=tidft_env%mat_s, &
                                   retain_locality=.FALSE., &
                                   only_normalize=.FALSE., &
                                   nocc_of_domain=fake_e_per_block, &
                                   eps_filter=tidft_env%eps_filter, &
                                   order_lanczos=tidft_env%order_lanczos, &
                                   eps_lanczos=tidft_env%eps_lanczos, &
                                   max_iter_lanczos=tidft_env%max_iter_lanczos)

            ! compute the OxO block of the KS matrix
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_ks(ispin), &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                0.0_dp, tempNO1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                tempNO1, &
                                0.0_dp, tempOO1, &
                                filter_eps=tidft_env%eps_filter)

            CALL cp_dbcsr_syevd(tempOO1, tempOO2, &
                                orbital_set%orbital_list(ispin, 1)%orb_energy, &
                                !o_eigenvalues(:,ispin), &
                                para_env=tidft_env%para_env, &
                                blacs_env=tidft_env%blacs_env)

            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                tempOO2, &
                                0.0_dp, tempNO1, &
                                filter_eps=tidft_env%eps_filter)

            ! RZK-WARNING: this makes T inconsistent with mainvar
            CALL dbcsr_copy(tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), tempNO1)

            ! Obtain info on the position of orbitals within the matrix
            CALL get_orbital_indices(tempNO1, &
                    orbital_set%orbital_list(ispin, 1)%block_within_matrix, &
                    orbital_set%orbital_list(ispin, 1)%index_within_block)

            CALL dbcsr_release(tempNO1)
            CALL dbcsr_release(tempOO1)
            CALL dbcsr_release(tempOO2)

            ! Virtual orbitals
            CALL dbcsr_create(tempNV1, &
                              template=tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempVOcc1, &
                              template=tidft_env%e_states(estate)%slater_det(idet)%mat_vo(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempVOcc2, &
                              template=tidft_env%e_states(estate)%slater_det(idet)%mat_vo(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempVV1, &
                              template=tidft_env%mat_sigma_vv(estate,idet,ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempVV2, &
                              template=tidft_env%mat_sigma_vv(estate,idet,ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            ! Generate random virtual matrix
            CALL dbcsr_init_random(tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                   keep_sparsity=.FALSE.)

            ! Project the orbital subspace out
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%mat_s, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                0.0_dp, tempNV1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                tempNV1, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                0.0_dp, tempVOcc1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tempVOcc1, &
                                tidft_env%mat_sigma_inv(estate,idet,ispin), &
                                0.0_dp, tempVOcc2, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                tempVOcc2, &
                                0.0_dp, tempNV1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_add(tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), tempNV1, 1.0_dp, -1.0_dp)

            ! compute VxV overlap
            !CALL dbcsr_multiply("N", "N", 1.0_dp, &
            !                    tidft_env%mat_s, &
            !                    tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
            !                    0.0_dp, tempNV1, &
            !                    filter_eps=tidft_env%eps_filter)

            !CALL dbcsr_multiply("T", "N", 1.0_dp, &
            !                    tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
            !                    tempNV1, &
            !                    0.0_dp, tempVV1, &
            !                    filter_eps=tidft_env%eps_filter)

            CALL orthogonalize_mos(ket=tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                   overlap=tempVV1, &
                                   metric=tidft_env%mat_s, &
                                   retain_locality=.FALSE., &
                                   only_normalize=.FALSE., &
                                   nocc_of_domain=fake_e_per_block, &
                                   eps_filter=tidft_env%eps_filter, &
                                   order_lanczos=tidft_env%order_lanczos, &
                                   eps_lanczos=tidft_env%eps_lanczos, &
                                   max_iter_lanczos=tidft_env%max_iter_lanczos)

            ! compute VxV block of the KS matrix
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_ks(ispin), &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                0.0_dp, tempNV1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                tempNV1, &
                                0.0_dp, tempVV1, &
                                filter_eps=tidft_env%eps_filter)

            !CALL dbcsr_get_info(tempVV1, nfullrows_total=n)
            CALL cp_dbcsr_syevd(tempVV1, tempVV2, &
                                orbital_set%orbital_list(ispin, 2)%orb_energy, &
                                !v_eigenvalues(:,ispin), &
                                para_env=tidft_env%para_env, &
                                blacs_env=tidft_env%blacs_env)

            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                tempVV2, &
                                0.0_dp, tempNV1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_copy(tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), tempNV1)

            ! Obtain info on the position of orbitals within the matrix
            CALL get_orbital_indices(tempNV1, &
                    orbital_set%orbital_list(ispin, 2)%block_within_matrix, &
                    orbital_set%orbital_list(ispin, 2)%index_within_block)

            CALL dbcsr_release(tempNV1)
            CALL dbcsr_release(tempVOcc1)
            CALL dbcsr_release(tempVOcc2)
            CALL dbcsr_release(tempVV1)
            CALL dbcsr_release(tempVV2)

         ENDIF

      ENDDO

      DEALLOCATE(fake_e_per_block)

   END SUBROUTINE construct_canonical_orbitals

   SUBROUTINE get_orbital_indices(matrix, block_within_matrix, index_within_block)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix 
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: block_within_matrix, &
                                                            index_within_block

      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes
      INTEGER                                            :: mo_start, iblock, imo

      CALL dbcsr_get_info(matrix, row_blk_size=blk_sizes)
      mo_start = 1
      DO iblock = 1, SIZE(blk_sizes)
         block_within_matrix(mo_start:mo_start+blk_sizes(iblock)-1) = iblock
         index_within_block(mo_start:mo_start+blk_sizes(iblock)-1) = (/ (imo, imo = 1, blk_sizes(iblock)) /)
         mo_start = mo_start + blk_sizes(iblock)
      ENDDO

   END SUBROUTINE get_orbital_indices

   SUBROUTINE copy_excitation(dest, source)

      TYPE(tidft_excitation_type), INTENT(IN)             :: source
      TYPE(tidft_excitation_type), INTENT(OUT)            :: dest

      CPASSERT(SIZE(source%delta_epsilon).EQ.SIZE(dest%delta_epsilon))

      ! copy data 
      dest%order = source%order
      dest%excitation_energy = source%excitation_energy
      dest%unpairings_count = source%unpairings_count
      dest%delta_epsilon(:) = source%delta_epsilon(:)
      dest%occupied(:) = source%occupied(:)
      dest%virtual(:) = source%virtual(:)
      dest%spin(:) = source%spin(:)

   END SUBROUTINE copy_excitation

END MODULE tidft_guess_methods

