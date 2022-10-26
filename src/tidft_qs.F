!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface between TIDFT and QS
!> \par History
!>       2021.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_qs
   USE tidft_types,                     ONLY: tidft_env_type
   USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                              almo_mat_dim_occ,&
                                              almo_mat_dim_virt,&
                                              almo_mat_dim_virt_disc,&
                                              almo_mat_dim_virt_full,&
                                              almo_scf_env_type
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE dbcsr_api,                       ONLY: &
        dbcsr_complete_redistribute, dbcsr_copy, dbcsr_copy_into_existing, dbcsr_create, &
        dbcsr_desymmetrize, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_get_stored_coordinates, &
        dbcsr_multiply, dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_block2d, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_work_create
   USE input_constants,                 ONLY: almo_constraint_ao_overlap,&
                                              almo_constraint_block_diagonal,&
                                              almo_constraint_distance,&
                                              almo_domain_layout_molecular,&
                                              almo_mat_distr_atomic,&
                                              almo_mat_distr_molecular,&
                                              do_bondparm_covalent,&
                                              do_bondparm_vdw
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_allgather
   USE molecule_types,                  ONLY: get_molecule_set_info,&
                                              molecule_type
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              scf_env_create
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_qs'

   PUBLIC :: matrix_layout, tidft_init_ks_matrices, tidft_dm_to_tidft_ks, &
             tidft_dm_to_qs_env

CONTAINS

! **************************************************************************************************
!> \brief create the ALMO matrix templates
!> \param matrix_new ...
!> \param matrix_qs ...
!> \param almo_scf_env ...
!> \param name_new ...
!> \param orbital_type ...
!> \param symmetry_new ...
!> \param spin_key ...
!> \param init_domains ...
!> \par History
!>       2021.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE matrix_layout(matrix_ao_template, matrix_new, name_new, &
      orbital_type_new, suggested_mo_blk_size, symmetry_new, nmos, verbose)

      TYPE(tidft_env_type)                               :: tidft_env
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_ao_template
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_new
      CHARACTER(len=*), INTENT(IN)                       :: name_new
      INTEGER, DIMENSION(2), INTENT(IN)                  :: orbital_type_new
      INTEGER, INTENT(IN)                                :: suggested_mo_blk_size
      CHARACTER, INTENT(IN)                              :: symmetry_new
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: nmos
      LOGICAL, OPTIONAL                                  :: verbose

      INTEGER                                            :: dimen, handle, hold, iatom, iblock_col, &
                                                            iblock_row, imol, mynode, natoms, &
                                                            nblkrows_tot, nlength, nmols, row, &
                                                            ncores1D, unit_nr
      INTEGER, DIMENSION(2)                              :: local_nmos
      INTEGER, DIMENSION(:), POINTER :: blk_distr, blk_sizes, block_sizes_new, col_distr_new, &
         col_sizes_new, distr_new_array, row_distr_new, row_sizes_new
      TYPE(dbcsr_distribution_type)                      :: dist_new, dist_qs
      LOGICAL                                            :: my_verbose
      TYPE(cp_logger_type), POINTER                      :: logger

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (PRESENT(nmos)) THEN
         local_nmos = nmos
      ELSE
         local_nmos = (/-1, -1/)
      ENDIF

      IF (PRESENT(verbose)) THEN
         my_verbose = verbose
      ELSE
         my_verbose = .FALSE.
      ENDIF

      ! obtain distribution from the qs matrix - it might be useful
      ! to get the structure of the AO dimensions
      CALL dbcsr_get_info(matrix_ao_template, distribution=dist_qs)

      DO dimen = 1, 2 ! 1 - row, 2 - column dimension

         IF (orbital_type_new(dimen) == almo_mat_dim_aobasis) THEN ! this dimension is AO

            ! structure of an AO dimension can be copied from matrix_ao_template
            IF (dimen == 1) THEN !rows
               CALL dbcsr_distribution_get(dist_qs, row_dist=blk_distr)
            ELSE !columns
               CALL dbcsr_distribution_get(dist_qs, col_dist=blk_distr)
            ENDIF

            CALL dbcsr_get_info(matrix_ao_template, row_blk_size=blk_sizes)

            ALLOCATE (block_sizes_new(SIZE(blk_sizes)))
            ALLOCATE (distr_new_array(SIZE(blk_distr)))
            block_sizes_new(:) = blk_sizes(:)
            distr_new_array(:) = blk_distr(:)

         ELSE ! this dimension is not AO

            IF (orbital_type_new(dimen) == almo_mat_dim_occ .OR. &
                orbital_type_new(dimen) == almo_mat_dim_virt) THEN

               IF (.NOT.PRESENT(nmos)) CPABORT("Implementation error: supply nmos")

               ! for MOs, sizes of blocks and their distibution is created from scratch
               IF (dimen == 1) THEN !rows
                  CALL dbcsr_distribution_get(dist_qs, nprows=ncores1D)
               ELSE !columns
                  CALL dbcsr_distribution_get(dist_qs, npcols=ncores1D)
               ENDIF

               CALL create_distribute_mo_blocks( &
                       norbs = local_nmos(dimen), &
                       ncores1D = ncores1D, &
                       suggested_blk_size = suggested_mo_blk_size, &
                       blk_sizes = block_sizes_new, &
                       blk_distr = distr_new_array)

               IF (my_verbose .AND. unit_nr > 0) THEN
                  WRITE(*,*) "NCORES(rows, then cols): ", ncores1D
               ENDIF

            ELSE

               CPABORT("Illegal dimension")

            ENDIF ! end choosing dim size (atomic, occ, virt)

         ENDIF ! end choosing dimension size (AOs vs .NOT.AOs)

         ! create final arrays
         IF (dimen == 1) THEN !rows
            ALLOCATE (row_sizes_new(SIZE(block_sizes_new)))
            ALLOCATE (row_distr_new(SIZE(distr_new_array)))
            row_sizes_new = block_sizes_new
            row_distr_new = distr_new_array
         ELSE !columns
            ALLOCATE (col_sizes_new(SIZE(block_sizes_new)))
            ALLOCATE (col_distr_new(SIZE(distr_new_array)))
            col_sizes_new = block_sizes_new
            col_distr_new = distr_new_array
         ENDIF

         IF (ASSOCIATED(block_sizes_new)) DEALLOCATE(block_sizes_new)
         IF (ASSOCIATED(distr_new_array)) DEALLOCATE(distr_new_array)

      ENDDO ! both rows and columns are done

      IF (my_verbose .AND. unit_nr > 0) THEN
         WRITE(*,*) "ROW: sizes"
         WRITE(*,*) row_sizes_new
         WRITE(*,*) "ROW: distr"
         WRITE(*,*) row_distr_new
         WRITE(*,*) "COL: sizes"
         WRITE(*,*) col_sizes_new
         WRITE(*,*) "COL: distr"
         WRITE(*,*) col_distr_new
      ENDIF

      ! Create the distribution
      CALL dbcsr_distribution_new(dist_new, template=dist_qs, &
                                  row_dist=row_distr_new, col_dist=col_distr_new, &
                                  reuse_arrays=.TRUE.)

      ! Create the matrix
      CALL dbcsr_create(matrix_new, name_new, &
                        dist_new, symmetry_new, &
                        row_sizes_new, col_sizes_new, reuse_arrays=.TRUE.)
      CALL dbcsr_finalize(matrix_new)

      CALL dbcsr_distribution_release(dist_new)

      IF (ASSOCIATED(row_sizes_new)) DEALLOCATE(row_sizes_new)
      IF (ASSOCIATED(row_distr_new)) DEALLOCATE(row_distr_new)
      IF (ASSOCIATED(col_sizes_new)) DEALLOCATE(col_sizes_new)
      IF (ASSOCIATED(col_distr_new)) DEALLOCATE(col_distr_new)

   END SUBROUTINE matrix_layout

! **************************************************************************************************
!> \brief Distribute matrix blocks
!> \param norbs ...
!> \param ncores1D ...
!> \param suggested_block_size: -1 means no suggestion
!> \param blk_sizes ...
!> \param blk_distr ...
!> \par History
!>       2021.09 created [Cindy Pham]
!> \author Cindy Pham
! **************************************************************************************************
   SUBROUTINE create_distribute_mo_blocks(norbs, ncores1D, suggested_blk_size, blk_sizes, blk_distr)


      INTEGER, INTENT(IN)                                :: norbs, ncores1D, suggested_blk_size

      INTEGER, DIMENSION(:), POINTER, INTENT(INOUT)      :: blk_sizes, blk_distr

      INTEGER, DIMENSION(:), POINTER                     :: trial_blk_sizes

      INTEGER                                            :: norbs_per_block, ndelta, &
                                                            igroup, best_diff, &
                                                            remainder, I, J, nblocks, &
                                                            ngroups, diff, best_nblocks, &
                                                            nblocks_suggested, &
                                                            nblocks_delta
      REAL                                               :: istd2, best_std2


!     ! Case 1: the number of orbitals is less than the number of cores
!      IF (norbs .LT. ncores1D) THEN
!
!      ALLOCATE(blk_sizes(norbs), blk_distr(norbs))
!      blk_sizes(:) = 1
!      blk_distr(:) = (/ ( I, I = 0, norbs-1, 1) /)
!     !! Case 2: given a suggested_blk_size, the number of distributed blks = ncores1D*ngroups
!     ELSE
!        ngroups = NINT(1.0_dp*norbs/suggested_blk_size/ncores1D)
!        IF (ngroups .EQ. 0) ngroups = 1
!        ! Vary ngroups within 10%
!        ndelta = CEILING(1.0_dp*ngroups/10)
!        ! assign a random large value
!        best_std2 = 1.0_dp*(norbs**2+1)
!        best_diff = norbs
!        ! Choose the best ngroups
!        DO igroup = ngroups-ndelta, ngroups+ndelta ! group_size
!           IF (igroup .EQ. 0) CYCLE
!           nblocks = ncores1D*igroup
!           norbs_per_block = FLOOR((1.0_dp*norbs)/nblocks)
!           ALLOCATE(trial_blk_sizes(nblocks))
!           trial_blk_sizes(:) = norbs_per_block
!           remainder = MOD(norbs, nblocks)
!           IF (remainder.GT.0) trial_blk_sizes(1:remainder) = trial_blk_sizes(1:remainder) + 1
!           diff = ABS(norbs_per_block - suggested_blk_size)
           ! Compute standard deviation square across all trial_blk_sizes
           !istd2 = 0.0_dp
           !DO I = 1, nblocks
            !  istd2 = istd2 + trial_blk_sizes(I)**2
           !ENDDO
           !istd2 = istd2 - (1.0_dp*norbs/nblocks)**2
           ! Compute standard deviation across all trial blk_sizes
!           istd2 = ABS(SUM(trial_blk_sizes(:)**2)/nblocks - (1.0_dp*norbs/nblocks)**2)
           ! Choose ngroups with better distribution
           !IF (istd2 .LT. best_std2 .OR. ABS(istd2 - best_std2) .LT. 1E-8) THEN
           !IF (istd2 .LT. best_std2) THEN
              !best_std2 = istd2
              !best_nblocks = nblocks
              !IF (ASSOCIATED(blk_sizes)) DEALLOCATE(blk_sizes)
              !ALLOCATE(blk_sizes(nblocks))
              !blk_sizes(:) = trial_blk_sizes(:)
              !IF (istd2 .LE. best_std2 .OR. diff .LE. best_diff) THEN
              !   best_diff = diff
              !   best_nblocks = nblocks
              !   ALLOCATE(blk_sizes(best_nblock))
              !   blk_sizes = trial_blk_sizes ! best option whose value is the closest to suggested_blk_size
              !ENDIF
           !ENDIF
            ! Choose ngroups with better distribution
!            IF (istd2 .LT. best_std2 .OR. ABS(istd2 - best_std2) .LT. 1E-8) THEN
!               best_std2 = istd2
!               IF (istd2 .LE. best_std2 .OR. diff .LE. best_diff) THEN
!                  best_diff = diff
!                  best_nblocks = nblocks
!                  ngroups = igroup !for checking purpose
!                  ALLOCATE(blk_sizes(best_nblocks))
!                  blk_sizes = trial_blk_sizes!copy ! best copy whose value is the closest to suggested_blk_size
!               ENDIF
!            ENDIF
!         END DO
!         DEALLOCATE(trial_blk_sizes)
        !END DO
        ! Find blk_distr()
!        ALLOCATE(blk_distr(best_nblocks))
!        I = 0
!        DO J=1, best_nblocks
!           IF (I .LT. ncores1D) THEN
!              blk_distr(J) = I
!              I = I + 1
!           ELSE
!              blk_distr(J) = 0
!              I = 1
!           ENDIF
!        END DO
!
!ENDIF

      nblocks_suggested = CEILING(1.0_dp*norbs/suggested_blk_size)
      IF (nblocks_suggested .LE. 0) nblocks_suggested = 1

      ! Vary the number of blocks within 10%
      nblocks_delta = CEILING(1.0_dp*nblocks_suggested/10)

      ! assign a random large value
      best_std2 = 1.0_dp*(norbs**2+1)
      best_diff = norbs

      ! Choose the best nblocks
      ! RZK-WARNING: Unlike the previous scheme,
      ! This loop is over nblocks, not over ngroups
      ! This might lead to worse load-balancing but does not introduce errors
      ! In the future an even better load-balancing scheme can be designed
      ! One idea is to loop over blocks but to make sure the interval includes nearby groups
      ! An important check should be added: no blocks containing zero orbitals
      DO nblocks = nblocks_suggested-nblocks_delta, nblocks_suggested+nblocks_delta

         IF (nblocks .LE. 0) CYCLE
         norbs_per_block = FLOOR((1.0_dp*norbs)/nblocks)

         ALLOCATE(trial_blk_sizes(nblocks))
         trial_blk_sizes(:) = norbs_per_block
         remainder = MOD(norbs, nblocks)
         IF (remainder.GT.0) trial_blk_sizes(1:remainder) = trial_blk_sizes(1:remainder) + 1

         ! Compute standard deviation square across all trial_blk_sizes
         istd2 = 0.0_dp
         DO I = 1, nblocks
            istd2 = istd2 + trial_blk_sizes(I)**2
         ENDDO
         istd2 = istd2/nblocks - (1.0_dp*norbs/nblocks)**2

         ! Choose nblocks with better distribution
         IF (istd2 .LT. best_std2) THEN
            best_std2 = istd2
            best_nblocks = nblocks
            IF (ASSOCIATED(blk_sizes)) DEALLOCATE(blk_sizes) 
            ALLOCATE(blk_sizes(nblocks))
            blk_sizes(:) = trial_blk_sizes(:)
         ENDIF

         DEALLOCATE(trial_blk_sizes)

      END DO

      ! Find blk_distr()
      ALLOCATE(blk_distr(best_nblocks))
      I = 0
      DO J=1, best_nblocks
         IF (I .LT. ncores1D) THEN
            blk_distr(J) = I
            I = I + 1
         ELSE
            blk_distr(J) = 0
            I = 1
         ENDIF
      END DO

   END SUBROUTINE create_distribute_mo_blocks

! **************************************************************************************************
!> \brief Initialization of the QS and TIDFT KS matrices
!> \param qs_env ...
!> \param tidft_env ...
!> \param mat_distr_aos ...
!> \param eps_filter ...
!> \par History
!>       2021.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_init_ks_matrices(qs_env, tidft_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type)                               :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_init_ks_matrices'

      INTEGER                                            :: handle, ispin, nspins, &
                                                            stateI
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_qs_ks, matrix_qs_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      NULLIFY (sab_orb)

      ! get basic quantities from the qs_env
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_qs_s, &
                      matrix_ks=matrix_qs_ks, &
                      ks_env=ks_env, &
                      sab_orb=sab_orb)

      nspins = tidft_env%nspins
      CPASSERT(nspins.EQ.dft_control%nspins)

      ! create matrix_ks in the QS env if necessary
      IF (.NOT. ASSOCIATED(matrix_qs_ks)) THEN
         CALL dbcsr_allocate_matrix_set(matrix_qs_ks, nspins)
         DO ispin = 1, nspins
            ALLOCATE (matrix_qs_ks(ispin)%matrix)
            CALL dbcsr_create(matrix_qs_ks(ispin)%matrix, &
                              template=matrix_qs_s(1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_qs_ks(ispin)%matrix, sab_orb)
            CALL dbcsr_set(matrix_qs_ks(ispin)%matrix, 0.0_dp)
         ENDDO
         CALL set_ks_env(ks_env, matrix_ks=matrix_qs_ks)
      ENDIF

      ! copy to tidft_env
      DO ispin = 1, tidft_env%nspins

         DO stateI = 1, tidft_env%nstates

            CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(1)%mat_ks(ispin), matrix_qs_ks(ispin)%matrix)
            CALL dbcsr_filter(tidft_env%e_states(stateI)%slater_det(1)%mat_ks(ispin), tidft_env%eps_filter)

            IF (tidft_env%triplet_correction) THEN

               CALL dbcsr_copy(tidft_env%e_states(stateI)%slater_det(2)%mat_ks(ispin), matrix_qs_ks(ispin)%matrix)
               CALL dbcsr_filter(tidft_env%e_states(stateI)%slater_det(2)%mat_ks(ispin), tidft_env%eps_filter)

            ENDIF

         ENDDO ! stateI

      ENDDO ! ispin

      CALL timestop(handle)

   END SUBROUTINE tidft_init_ks_matrices

!tidft! **************************************************************************************************
!tidft!> \brief Create MOs in the QS env to be able to return ALMOs to QS
!tidft!> \param qs_env ...
!tidft!> \param almo_scf_env ...
!tidft!> \par History
!tidft!>       2016.12 created [Yifei Shi]
!tidft!> \author Yifei Shi
!tidft! **************************************************************************************************
!tidft   SUBROUTINE construct_qs_mos(qs_env, almo_scf_env)
!tidft
!tidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tidft
!tidft      CHARACTER(len=*), PARAMETER                        :: routineN = 'construct_qs_mos'
!tidft
!tidft      INTEGER                                            :: handle, ispin, ncol_fm, nrow_fm
!tidft      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
!tidft      TYPE(cp_fm_type), POINTER                          :: mo_fm_copy
!tidft      TYPE(dft_control_type), POINTER                    :: dft_control
!tidft      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
!tidft      TYPE(qs_scf_env_type), POINTER                     :: scf_env
!tidft
!tidft      CALL timeset(routineN, handle)
!tidft
!tidft      ! create and init scf_env (this is necessary to return MOs to qs)
!tidft      NULLIFY (mos, mo_fm_copy, fm_struct_tmp, scf_env)
!tidft      CALL scf_env_create(scf_env)
!tidft
!tidft      !CALL qs_scf_env_initialize(qs_env, scf_env)
!tidft      CALL set_qs_env(qs_env, scf_env=scf_env)
!tidft      CALL get_qs_env(qs_env, dft_control=dft_control, mos=mos)
!tidft
!tidft      CALL dbcsr_get_info(almo_scf_env%matrix_t(1), nfullrows_total=nrow_fm, nfullcols_total=ncol_fm)
!tidft
!tidft      ! allocate and init mo_set
!tidft      DO ispin = 1, almo_scf_env%nspins
!tidft
!tidft         ! Currently only fm version of mo_set is usable.
!tidft         ! First transform the matrix_t to fm version
!tidft         CALL allocate_mo_set(mo_set=mos(ispin)%mo_set, &
!tidft                              nao=nrow_fm, &
!tidft                              nmo=ncol_fm, &
!tidft                              nelectron=almo_scf_env%nelectrons_total, &
!tidft                              n_el_f=REAL(almo_scf_env%nelectrons_total, dp), &
!tidft                              maxocc=2.0_dp, &
!tidft                              flexible_electron_count=dft_control%relax_multiplicity)
!tidft
!tidft         CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_fm, ncol_global=ncol_fm, &
!tidft                                  context=almo_scf_env%blacs_env, &
!tidft                                  para_env=almo_scf_env%para_env)
!tidft
!tidft         CALL cp_fm_create(mo_fm_copy, fm_struct_tmp, name="t_orthogonal_converted_to_fm")
!tidft         CALL cp_fm_struct_release(fm_struct_tmp)
!tidft         !CALL copy_dbcsr_to_fm(almo_scf_env%matrix_t(ispin), mo_fm_copy)
!tidft
!tidft         CALL init_mo_set(mos(ispin)%mo_set, fm_ref=mo_fm_copy, name='fm_mo')
!tidft
!tidft         CALL cp_fm_release(mo_fm_copy)
!tidft
!tidft      END DO
!tidft
!tidft      CALL timestop(handle)
!tidft
!tidft   END SUBROUTINE construct_qs_mos

! **************************************************************************************************
!> \brief supplies density matrix to the qs env
!> \param qs_env ...
!> \param matrix_p ...
!> \par History
!>       2021.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_dm_to_qs_env(qs_env, matrix_p)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), DIMENSION(:)                     :: matrix_p

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tidft_dm_to_qs_env'

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      NULLIFY (rho, rho_ao)
      nspins = SIZE(matrix_p)
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ! set the new density matrix
      DO ispin = 1, nspins
         CALL dbcsr_copy(rho_ao(ispin)%matrix, matrix_p(ispin))
      END DO
      CALL qs_rho_update_rho(rho, qs_env=qs_env)
      CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CALL timestop(handle)

   END SUBROUTINE tidft_dm_to_qs_env

! **************************************************************************************************
!> \brief uses the density matrix to compute QS KS matrix and the new energy
!> \param qs_env ...
!> \param matrix_p ...
!> \param energy_total ...
!> \param smear ...
!> \param kTS_sum ...
!> \par History
!>       2021.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_dm_to_qs_ks(qs_env, matrix_p, energy_total, smear, kTS_sum)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), DIMENSION(:)                     :: matrix_p
      REAL(KIND=dp)                                      :: energy_total
      LOGICAL, INTENT(IN), OPTIONAL                      :: smear
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: kTS_sum

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tidft_dm_to_qs_ks'

      INTEGER                                            :: handle
      LOGICAL                                            :: smearing
      REAL(KIND=dp)                                      :: entropic_term
      TYPE(qs_energy_type), POINTER                      :: energy

      CALL timeset(routineN, handle)

      IF (PRESENT(smear)) THEN
         smearing = smear
      ELSE
         smearing = .FALSE.
      ENDIF

      IF (PRESENT(kTS_sum)) THEN
         entropic_term = kTS_sum
      ELSE
         entropic_term = 0.0_dp
      ENDIF

      NULLIFY (energy)
      CALL get_qs_env(qs_env, energy=energy)
      CALL tidft_dm_to_qs_env(qs_env, matrix_p)
      CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., &
                               print_active=.TRUE.)

      !! Add electronic entropy contribution if smearing is requested
      !! Previous QS entropy is replaced by the sum of the entropy for each spin
      IF (smearing) THEN
         energy%total = energy%total - energy%kTS + entropic_term
      END IF

      energy_total = energy%total

      CALL timestop(handle)

   END SUBROUTINE tidft_dm_to_qs_ks

! **************************************************************************************************
!> \brief uses the density matrix to compute TIDFT KS matrix and the new energy
!> \param qs_env ...
!> \param matrix_p ...
!> \param matrix_ks ...
!> \param energy_total ...
!> \param eps_filter ...
!> \param smear ...
!> \param kTS_sum ...
!> \par History
!>       2021.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_dm_to_tidft_ks(qs_env, matrix_p, matrix_ks, energy_total, eps_filter, &
                                 smear, kTS_sum)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), DIMENSION(:)                     :: matrix_p, matrix_ks
      REAL(KIND=dp)                                      :: energy_total, eps_filter
      !INTEGER, INTENT(IN)                                :: mat_distr_aos
      LOGICAL, INTENT(IN), OPTIONAL                      :: smear
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: kTS_sum

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_dm_to_tidft_ks'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: smearing
      REAL(KIND=dp)                                      :: entropic_term
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_qs_ks

      CALL timeset(routineN, handle)

      IF (PRESENT(smear)) THEN
         smearing = smear
      ELSE
         smearing = .FALSE.
      ENDIF

      IF (PRESENT(kTS_sum)) THEN
         entropic_term = kTS_sum
      ELSE
         entropic_term = 0.0_dp
      ENDIF

      ! update KS matrix in the QS env
      CALL tidft_dm_to_qs_ks(qs_env, matrix_p, energy_total, & !mat_distr_aos, &
                            smear=smearing, &
                            kTS_sum=entropic_term)

      nspins = SIZE(matrix_ks)

      ! get KS matrix from the QS env and convert to the ALMO format
      CALL get_qs_env(qs_env, matrix_ks=matrix_qs_ks)
      DO ispin = 1, nspins
         CALL dbcsr_copy(matrix_ks(ispin), matrix_qs_ks(ispin)%matrix)
         !CALL matrix_qs_to_almo(matrix_qs_ks(ispin)%matrix, &
         !                       matrix_ks(ispin), &
         !                       mat_distr_aos, .FALSE.)
         CALL dbcsr_filter(matrix_ks(ispin), eps_filter)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE tidft_dm_to_tidft_ks

!tidft! **************************************************************************************************
!tidft!> \brief update qs_env total energy
!tidft!> \param qs_env ...
!tidft!> \param energy ...
!tidft!> \param energy_singles_corr ...
!tidft!> \par History
!tidft!>       2021.12 created [Rustam Z Khaliullin]
!tidft!> \author Rustam Z Khaliullin
!tidft! **************************************************************************************************
!tidft   SUBROUTINE almo_scf_update_ks_energy(qs_env, energy, energy_singles_corr)
!tidft
!tidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tidft      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: energy, energy_singles_corr
!tidft
!tidft      TYPE(qs_energy_type), POINTER                      :: qs_energy
!tidft
!tidft      CALL get_qs_env(qs_env, energy=qs_energy)
!tidft
!tidft      IF (PRESENT(energy_singles_corr)) THEN
!tidft         qs_energy%singles_corr = energy_singles_corr
!tidft      ELSE
!tidft         qs_energy%singles_corr = 0.0_dp
!tidft      ENDIF
!tidft
!tidft      IF (PRESENT(energy)) THEN
!tidft         qs_energy%total = energy
!tidft      ENDIF
!tidft
!tidft      qs_energy%total = qs_energy%total + qs_energy%singles_corr
!tidft
!tidft   END SUBROUTINE almo_scf_update_ks_energy

!tidft! **************************************************************************************************
!tidft!> \brief Creates the matrix that imposes absolute locality on MOs
!tidft!> \param qs_env ...
!tidft!> \param almo_scf_env ...
!tidft!> \par History
!tidft!>       2011.11 created [Rustam Z. Khaliullin]
!tidft!> \author Rustam Z. Khaliullin
!tidft! **************************************************************************************************
!tidft   SUBROUTINE almo_scf_construct_quencher(qs_env, almo_scf_env)
!tidft
!tidft      TYPE(qs_environment_type), POINTER                 :: qs_env
!tidft      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
!tidft
!tidft      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_construct_quencher'
!tidft
!tidft      CHARACTER                                          :: sym
!tidft      INTEGER :: col, contact_atom_1, contact_atom_2, domain_col, domain_map_local_entries, &
!tidft         domain_row, global_entries, global_list_length, grid1, GroupID, handle, hold, iatom, &
!tidft         iatom2, iblock_col, iblock_row, idomain, idomain2, ientry, igrid, ineig, ineighbor, &
!tidft         iNode, inode2, ipair, ispin, jatom, jatom2, jdomain2, local_list_length, &
!tidft         max_domain_neighbors, max_neig, mynode, nblkcols_tot, nblkrows_tot, nblks, ndomains, &
!tidft         neig_temp, nnode2, nNodes, row, unit_nr
!tidft      INTEGER, ALLOCATABLE, DIMENSION(:) :: current_number_neighbors, domain_entries_cpu, &
!tidft         domain_map_global, domain_map_local, first_atom_of_molecule, global_list, &
!tidft         last_atom_of_molecule, list_length_cpu, list_offset_cpu, local_list, offset_for_cpu
!tidft      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: domain_grid, domain_neighbor_list, &
!tidft                                                            domain_neighbor_list_excessive
!tidft      LOGICAL                                            :: already_listed, block_active, &
!tidft                                                            delayed_increment, found, &
!tidft                                                            max_neig_fails, tr
!tidft      REAL(KIND=dp)                                      :: contact1_radius, contact2_radius, &
!tidft                                                            distance, distance_squared, overlap, &
!tidft                                                            r0, r1, s0, s1, trial_distance_squared
!tidft      REAL(KIND=dp), DIMENSION(3)                        :: rab
!tidft      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_new_block
!tidft      TYPE(cell_type), POINTER                           :: cell
!tidft      TYPE(cp_logger_type), POINTER                      :: logger
!tidft      TYPE(dbcsr_distribution_type)                      :: dist
!tidft      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
!tidft      TYPE(dbcsr_type)                                   :: matrix_s_sym
!tidft      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
!tidft      TYPE(neighbor_list_iterator_p_type), &
!tidft         DIMENSION(:), POINTER                           :: nl_iterator, nl_iterator2
!tidft      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
!tidft         POINTER                                         :: sab_almo
!tidft      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
!tidft
!tidft      CALL timeset(routineN, handle)
!tidft
!tidft      ! get a useful output_unit
!tidft      logger => cp_get_default_logger()
!tidft      IF (logger%para_env%ionode) THEN
!tidft         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
!tidft      ELSE
!tidft         unit_nr = -1
!tidft      ENDIF
!tidft
!tidft      ndomains = almo_scf_env%ndomains
!tidft
!tidft      CALL get_qs_env(qs_env=qs_env, &
!tidft                      particle_set=particle_set, &
!tidft                      molecule_set=molecule_set, &
!tidft                      cell=cell, &
!tidft                      matrix_s=matrix_s, &
!tidft                      sab_almo=sab_almo)
!tidft
!tidft      ! if we are dealing with molecules get info about them
!tidft      IF (almo_scf_env%domain_layout_mos == almo_domain_layout_molecular .OR. &
!tidft          almo_scf_env%domain_layout_aos == almo_domain_layout_molecular) THEN
!tidft         ALLOCATE (first_atom_of_molecule(almo_scf_env%nmolecules))
!tidft         ALLOCATE (last_atom_of_molecule(almo_scf_env%nmolecules))
!tidft         CALL get_molecule_set_info(molecule_set, &
!tidft                                    mol_to_first_atom=first_atom_of_molecule, &
!tidft                                    mol_to_last_atom=last_atom_of_molecule)
!tidft      ENDIF
!tidft
!tidft      ! create a symmetrized copy of the ao overlap
!tidft      CALL dbcsr_create(matrix_s_sym, &
!tidft                        template=almo_scf_env%matrix_s(1), &
!tidft                        matrix_type=dbcsr_type_no_symmetry)
!tidft      CALL dbcsr_get_info(almo_scf_env%matrix_s(1), &
!tidft                          matrix_type=sym)
!tidft      IF (sym .EQ. dbcsr_type_no_symmetry) THEN
!tidft         CALL dbcsr_copy(matrix_s_sym, almo_scf_env%matrix_s(1))
!tidft      ELSE
!tidft         CALL dbcsr_desymmetrize(almo_scf_env%matrix_s(1), &
!tidft                                 matrix_s_sym)
!tidft      ENDIF
!tidft
!tidft      ALLOCATE (almo_scf_env%quench_t(almo_scf_env%nspins))
!tidft      ALLOCATE (almo_scf_env%domain_map(almo_scf_env%nspins))
!tidft
!tidft      !DO ispin=1,almo_scf_env%nspins
!tidft      ispin = 1
!tidft
!tidft      ! create the sparsity template for the occupied orbitals
!tidft      CALL matrix_almo_create(matrix_new=almo_scf_env%quench_t(ispin), &
!tidft                              matrix_qs=matrix_s(1)%matrix, &
!tidft                              almo_scf_env=almo_scf_env, &
!tidft                              name_new="T_QUENCHER", &
!tidft                              orbital_type=(/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
!tidft                              symmetry_new=dbcsr_type_no_symmetry, &
!tidft                              spin_key=ispin, &
!tidft                              init_domains=.FALSE.)
!tidft
!tidft      ! initialize distance quencher
!tidft      CALL dbcsr_work_create(almo_scf_env%quench_t(ispin), &
!tidft                             work_mutable=.TRUE.)
!tidft
!tidft      nblkrows_tot = dbcsr_nblkrows_total(almo_scf_env%quench_t(ispin))
!tidft      nblkcols_tot = dbcsr_nblkcols_total(almo_scf_env%quench_t(ispin))
!tidft
!tidft      CALL dbcsr_get_info(almo_scf_env%quench_t(ispin), distribution=dist)
!tidft      CALL dbcsr_distribution_get(dist, numnodes=nNodes, group=GroupID, mynode=mynode)
!tidft
!tidft      ! create global atom neighbor list from the local lists
!tidft      ! first, calculate number of local pairs
!tidft      local_list_length = 0
!tidft      CALL neighbor_list_iterator_create(nl_iterator, sab_almo)
!tidft      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
!tidft         ! nnode - total number of neighbors for iatom
!tidft         ! inode - current neighbor count
!tidft         CALL get_iterator_info(nl_iterator, &
!tidft                                iatom=iatom2, jatom=jatom2, inode=inode2, nnode=nnode2)
!tidft         !WRITE(*,*) "GET INFO: ",iatom2, jatom2, inode2, nnode2
!tidft         IF (inode2 == 1) THEN
!tidft            local_list_length = local_list_length + nnode2
!tidft         END IF
!tidft      END DO
!tidft      CALL neighbor_list_iterator_release(nl_iterator)
!tidft
!tidft      ! second, extract the local list to an array
!tidft      ALLOCATE (local_list(2*local_list_length))
!tidft      local_list(:) = 0
!tidft      local_list_length = 0
!tidft      CALL neighbor_list_iterator_create(nl_iterator2, sab_almo)
!tidft      DO WHILE (neighbor_list_iterate(nl_iterator2) == 0)
!tidft         CALL get_iterator_info(nl_iterator2, &
!tidft                                iatom=iatom2, jatom=jatom2)
!tidft         local_list(2*local_list_length + 1) = iatom2
!tidft         local_list(2*local_list_length + 2) = jatom2
!tidft         local_list_length = local_list_length + 1
!tidft      ENDDO ! end loop over pairs of atoms
!tidft      CALL neighbor_list_iterator_release(nl_iterator2)
!tidft
!tidft      ! third, communicate local length to the other nodes
!tidft      ALLOCATE (list_length_cpu(nNodes), list_offset_cpu(nNodes))
!tidft      CALL mp_allgather(2*local_list_length, list_length_cpu, GroupID)
!tidft
!tidft      ! fourth, create a global list
!tidft      list_offset_cpu(1) = 0
!tidft      DO iNode = 2, nNodes
!tidft         list_offset_cpu(iNode) = list_offset_cpu(iNode - 1) + &
!tidft                                  list_length_cpu(iNode - 1)
!tidft      ENDDO
!tidft      global_list_length = list_offset_cpu(nNodes) + list_length_cpu(nNodes)
!tidft
!tidft      ! fifth, communicate all list data
!tidft      ALLOCATE (global_list(global_list_length))
!tidft      CALL mp_allgather(local_list, global_list, &
!tidft                        list_length_cpu, list_offset_cpu, GroupID)
!tidft      DEALLOCATE (list_length_cpu, list_offset_cpu)
!tidft      DEALLOCATE (local_list)
!tidft
!tidft      ! calculate maximum number of atoms surrounding the domain
!tidft      ALLOCATE (current_number_neighbors(almo_scf_env%ndomains))
!tidft      current_number_neighbors(:) = 0
!tidft      global_list_length = global_list_length/2
!tidft      DO ipair = 1, global_list_length
!tidft         iatom2 = global_list(2*(ipair - 1) + 1)
!tidft         jatom2 = global_list(2*(ipair - 1) + 2)
!tidft         idomain2 = almo_scf_env%domain_index_of_atom(iatom2)
!tidft         jdomain2 = almo_scf_env%domain_index_of_atom(jatom2)
!tidft         ! add to the list
!tidft         current_number_neighbors(idomain2) = current_number_neighbors(idomain2) + 1
!tidft         ! add j,i with i,j
!tidft         IF (idomain2 .NE. jdomain2) THEN
!tidft            current_number_neighbors(jdomain2) = current_number_neighbors(jdomain2) + 1
!tidft         ENDIF
!tidft      ENDDO
!tidft      max_domain_neighbors = MAXVAL(current_number_neighbors)
!tidft
!tidft      ! use the global atom neighbor list to create a global domain neighbor list
!tidft      ALLOCATE (domain_neighbor_list_excessive(ndomains, max_domain_neighbors))
!tidft      current_number_neighbors(:) = 1
!tidft      DO ipair = 1, ndomains
!tidft         domain_neighbor_list_excessive(ipair, 1) = ipair
!tidft      ENDDO
!tidft      DO ipair = 1, global_list_length
!tidft         iatom2 = global_list(2*(ipair - 1) + 1)
!tidft         jatom2 = global_list(2*(ipair - 1) + 2)
!tidft         idomain2 = almo_scf_env%domain_index_of_atom(iatom2)
!tidft         jdomain2 = almo_scf_env%domain_index_of_atom(jatom2)
!tidft         already_listed = .FALSE.
!tidft         DO ineighbor = 1, current_number_neighbors(idomain2)
!tidft            IF (domain_neighbor_list_excessive(idomain2, ineighbor) .EQ. jdomain2) THEN
!tidft               already_listed = .TRUE.
!tidft               EXIT
!tidft            ENDIF
!tidft         ENDDO
!tidft         IF (.NOT. already_listed) THEN
!tidft            ! add to the list
!tidft            current_number_neighbors(idomain2) = current_number_neighbors(idomain2) + 1
!tidft            domain_neighbor_list_excessive(idomain2, current_number_neighbors(idomain2)) = jdomain2
!tidft            ! add j,i with i,j
!tidft            IF (idomain2 .NE. jdomain2) THEN
!tidft               current_number_neighbors(jdomain2) = current_number_neighbors(jdomain2) + 1
!tidft               domain_neighbor_list_excessive(jdomain2, current_number_neighbors(jdomain2)) = idomain2
!tidft            ENDIF
!tidft         ENDIF
!tidft      ENDDO ! end loop over pairs of atoms
!tidft      EALLOCATE (global_list)
!tidft
!tidft      max_domain_neighbors = MAXVAL(current_number_neighbors)
!tidft      ALLOCATE (domain_neighbor_list(ndomains, max_domain_neighbors))
!tidft      domain_neighbor_list(:, :) = 0
!tidft      domain_neighbor_list(:, :) = domain_neighbor_list_excessive(:, 1:max_domain_neighbors)
!tidft      DEALLOCATE (domain_neighbor_list_excessive)
!tidft
!tidft      ALLOCATE (almo_scf_env%domain_map(ispin)%index1(ndomains))
!tidft      ALLOCATE (almo_scf_env%domain_map(ispin)%pairs(max_domain_neighbors*ndomains, 2))
!tidft      almo_scf_env%domain_map(ispin)%pairs(:, :) = 0
!tidft      almo_scf_env%domain_map(ispin)%index1(:) = 0
!tidft      domain_map_local_entries = 0
!tidft
!tidft      ! RZK-warning intermediate [0,1] quencher values are ill-defined
!tidft      ! for molecules (not continuous and conceptually inadequate)
!tidft
!tidft      ! O(N) loop over domain pairs
!tidft      DO row = 1, nblkrows_tot
!tidft         DO col = 1, current_number_neighbors(row)
!tidft            tr = .FALSE.
!tidft            iblock_row = row
!tidft            iblock_col = domain_neighbor_list(row, col)
!tidft            CALL dbcsr_get_stored_coordinates(almo_scf_env%quench_t(ispin), &
!tidft                                              iblock_row, iblock_col, hold)
!tidft
!tidft            IF (hold .EQ. mynode) THEN
!tidft
!tidft               ! Translate indices of distribution blocks to indices of domain blocks
!tidft               ! Rows are AOs
!tidft               domain_row = almo_scf_env%domain_index_of_ao_block(iblock_row)
!tidft               ! Columns are electrons (i.e. MOs)
!tidft               domain_col = almo_scf_env%domain_index_of_mo_block(iblock_col)
!tidft
!tidft               SELECT CASE (almo_scf_env%constraint_type)
!tidft               CASE (almo_constraint_block_diagonal)
!tidft
!tidft                  block_active = .FALSE.
!tidft                  ! type of electron groups
!tidft                  IF (almo_scf_env%domain_layout_mos == almo_domain_layout_molecular) THEN
!tidft
!tidft                     ! type of ao domains
!tidft                     IF (almo_scf_env%domain_layout_aos == almo_domain_layout_molecular) THEN
!tidft
!tidft                        ! ao domains are molecular / electron groups are molecular
!tidft                        IF (domain_row == domain_col) THEN
!tidft                           block_active = .TRUE.
!tidft                        ENDIF
!tidft
!tidft                     ELSE ! ao domains are atomic
!tidft
!tidft                        ! ao domains are atomic / electron groups are molecular
!tidft                        CPABORT("Illegal: atomic domains and molecular groups")
!tidft
!tidft                     ENDIF
!tidft
!tidft                  ELSE ! electron groups are atomic
!tidft
!tidft                     ! type of ao domains
!tidft                     IF (almo_scf_env%domain_layout_aos == almo_domain_layout_molecular) THEN
!tidft
!tidft                        ! ao domains are molecular / electron groups are atomic
!tidft                        CPABORT("Illegal: molecular domains and atomic groups")
!tidft
!tidft                     ELSE
!tidft
!tidft                        ! ao domains are atomic / electron groups are atomic
!tidft                        IF (domain_row == domain_col) THEN
!tidft                           block_active = .TRUE.
!tidft                        ENDIF
!tidft
!tidft                     ENDIF
!tidft
!tidft                  ENDIF ! end type of electron groups
!tidft
!tidft                  IF (block_active) THEN
!tidft
!tidft                     NULLIFY (p_new_block)
!tidft                     CALL dbcsr_reserve_block2d(almo_scf_env%quench_t(ispin), &
!tidft                                                iblock_row, iblock_col, p_new_block)
!tidft                     CPASSERT(ASSOCIATED(p_new_block))
!tidft                     p_new_block(:, :) = 1.0_dp
!tidft
!tidft                     IF (domain_map_local_entries .GE. max_domain_neighbors*almo_scf_env%ndomains) THEN
!tidft                        CPABORT("weird... max_domain_neighbors is exceeded")
!tidft                     ENDIF
!tidft                     almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries + 1, 1) = iblock_row
!tidft                     almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries + 1, 2) = iblock_col
!tidft                     domain_map_local_entries = domain_map_local_entries + 1
!tidft
!tidft                  ENDIF
!tidft
!tidft               CASE (almo_constraint_ao_overlap)
!tidft
!tidft                  ! type of electron groups
!tidft                  IF (almo_scf_env%domain_layout_mos == almo_domain_layout_molecular) THEN
!tidft
!tidft                     ! type of ao domains
!tidft                     IF (almo_scf_env%domain_layout_aos == almo_domain_layout_molecular) THEN
!tidft
!tidft                        ! ao domains are molecular / electron groups are molecular
!tidft
!tidft                        ! compute the maximum overlap between the atoms of the two molecules
!tidft                        CALL dbcsr_get_block_p(matrix_s_sym, &
!tidft                                               iblock_row, iblock_col, p_new_block, found)
!tidft                        IF (found) THEN
!tidft                           !   CPErrorMessage(cp_failure_level,routineP,"S block not found")
!tidft                           !   CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
!tidft                           overlap = MAXVAL(ABS(p_new_block))
!tidft                        ELSE
!tidft                           overlap = 0.0_dp
!tidft                        ENDIF
!tidft
!tidft                     ELSE ! ao domains are atomic
!tidft
!tidft                        ! ao domains are atomic / electron groups are molecular
!tidft                        ! overlap_between_atom_and_molecule(atom=domain_row,molecule=domain_col)
!tidft                        CPABORT("atomic domains and molecular groups - NYI")
!tidft
!tidft                     ENDIF
!tidft
!tidft                  ELSE ! electron groups are atomic
!tidft
!tidft                     ! type of ao domains
!tidft                     IF (almo_scf_env%domain_layout_aos == almo_domain_layout_molecular) THEN
!tidft
!tidft                        ! ao domains are molecular / electron groups are atomic
!tidft                        ! overlap_between_atom_and_molecule(atom=domain_col,molecule=domain_row)
!tidft                        CPABORT("molecular domains and atomic groups - NYI")
!tidft
!tidft                     ELSE
!tidft
!tidft                        ! ao domains are atomic / electron groups are atomic
!tidft                        ! compute max overlap between atoms: domain_row and domain_col
!tidft                        CALL dbcsr_get_block_p(matrix_s_sym, &
!tidft                                               iblock_row, iblock_col, p_new_block, found)
!tidft                        IF (found) THEN
!tidft                           !   CPErrorMessage(cp_failure_level,routineP,"S block not found")
!tidft                           !   CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
!tidft                           overlap = MAXVAL(ABS(p_new_block))
!tidft                        ELSE
!tidft                           overlap = 0.0_dp
!tidft                        ENDIF
!tidft
!tidft                     ENDIF
!tidft
!tidft                  ENDIF ! end type of electron groups
!tidft
!tidft                  s0 = -LOG10(ABS(almo_scf_env%quencher_s0))
!tidft                  s1 = -LOG10(ABS(almo_scf_env%quencher_s1))
!tidft                  IF (overlap .EQ. 0.0_dp) THEN
!tidft                     overlap = -LOG10(ABS(almo_scf_env%eps_filter)) + 100.0_dp
!tidft                  ELSE
!tidft                     overlap = -LOG10(overlap)
!tidft                  ENDIF
!tidft                  IF (s0 .LT. 0.0_dp) THEN
!tidft                     CPABORT("S0 is less than zero")
!tidft                  ENDIF
!tidft                  IF (s1 .LE. 0.0_dp) THEN
!tidft                     CPABORT("S1 is less than or equal to zero")
!tidft                  ENDIF
!tidft                  IF (s0 .GE. s1) THEN
!tidft                     CPABORT("S0 is greater than or equal to S1")
!tidft                  ENDIF
!tidft
!tidft                  ! Fill in non-zero blocks if AOs are close to the electron center
!tidft                  IF (overlap .LT. s1) THEN
!tidft                     NULLIFY (p_new_block)
!tidft                     CALL dbcsr_reserve_block2d(almo_scf_env%quench_t(ispin), &
!tidft                                                iblock_row, iblock_col, p_new_block)
!tidft                     CPASSERT(ASSOCIATED(p_new_block))
!tidft                     IF (overlap .LE. s0) THEN
!tidft                        p_new_block(:, :) = 1.0_dp
!tidft                        !WRITE(*,'(A15,2I7,3F8.3,E11.3)') "INTRA-BLOCKS: ",&
!tidft                        !   iblock_col, iblock_row, s0, s1, overlap, p_new_block(1,1)
!tidft                     ELSE
!tidft                        p_new_block(:, :) = 1.0_dp/(1.0_dp + EXP(-(s0 - s1)/(s0 - overlap) - (s0 - s1)/(overlap - s1)))
!tidft                        !WRITE(*,'(A15,2I7,3F8.3,E11.3)') "INTER-BLOCKS: ",&
!tidft                        !   iblock_col, iblock_row, s0, s1, overlap, p_new_block(1,1)
!tidft                     ENDIF
!tidft
!tidft                     IF (ABS(p_new_block(1, 1)) .GT. ABS(almo_scf_env%eps_filter)) THEN
!tidft                        IF (domain_map_local_entries .GE. max_domain_neighbors*almo_scf_env%ndomains) THEN
!tidft                           CPABORT("weird... max_domain_neighbors is exceeded")
!tidft                        ENDIF
!tidft                        almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries + 1, 1) = iblock_row
!tidft                        almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries + 1, 2) = iblock_col
!tidft                        domain_map_local_entries = domain_map_local_entries + 1
!tidft                     ENDIF
!tidft
!tidft                  ENDIF
!tidft
!tidft               CASE (almo_constraint_distance)
!tidft
!tidft                  ! type of electron groups
!tidft                  IF (almo_scf_env%domain_layout_mos == almo_domain_layout_molecular) THEN
!tidft
!tidft                     ! type of ao domains
!tidft                     IF (almo_scf_env%domain_layout_aos == almo_domain_layout_molecular) THEN
!tidft
!tidft                        ! ao domains are molecular / electron groups are molecular
!tidft
!tidft                        ! compute distance between molecules: domain_row and domain_col
!tidft                        ! distance between molecules is defined as the smallest
!tidft                        ! distance among all atom pairs
!tidft                        IF (domain_row == domain_col) THEN
!tidft                           distance = 0.0_dp
!tidft                           contact_atom_1 = first_atom_of_molecule(domain_row)
!tidft                           contact_atom_2 = first_atom_of_molecule(domain_col)
!tidft                        ELSE
!tidft                           distance_squared = 1.0E+100_dp
!tidft                           contact_atom_1 = -1
!tidft                           contact_atom_2 = -1
!tidft                           DO iatom = first_atom_of_molecule(domain_row), last_atom_of_molecule(domain_row)
!tidft                              DO jatom = first_atom_of_molecule(domain_col), last_atom_of_molecule(domain_col)
!tidft                                 rab(:) = pbc(particle_set(iatom)%r(:), particle_set(jatom)%r(:), cell)
!tidft                                 trial_distance_squared = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
!tidft                                 IF (trial_distance_squared .LT. distance_squared) THEN
!tidft                                    distance_squared = trial_distance_squared
!tidft                                    contact_atom_1 = iatom
!tidft                                    contact_atom_2 = jatom
!tidft                                 ENDIF
!tidft                              ENDDO ! jatom
!tidft                           ENDDO ! iatom
!tidft                           CPASSERT(contact_atom_1 .GT. 0)
!tidft                           distance = SQRT(distance_squared)
!tidft                        ENDIF
!tidft
!tidft                     ELSE ! ao domains are atomic
!tidft
!tidft                        ! ao domains are atomic / electron groups are molecular
!tidft                        !distance_between_atom_and_molecule(atom=domain_row,molecule=domain_col)
!tidft                        CPABORT("atomic domains and molecular groups - NYI")
!tidft
!tidft                     ENDIF
!tidft
!tidft                  ELSE ! electron groups are atomic
!tidft
!tidft                     ! type of ao domains
!tidft                     IF (almo_scf_env%domain_layout_aos == almo_domain_layout_molecular) THEN
!tidft
!tidft                        ! ao domains are molecular / electron groups are atomic
!tidft                        !distance_between_atom_and_molecule(atom=domain_col,molecule=domain_row)
!tidft                        CPABORT("molecular domains and atomic groups - NYI")
!tidft
!tidft                     ELSE
!tidft
!tidft                        ! ao domains are atomic / electron groups are atomic
!tidft                        ! compute distance between atoms: domain_row and domain_col
!tidft                        rab(:) = pbc(particle_set(domain_row)%r(:), particle_set(domain_col)%r(:), cell)
!tidft                        distance = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
!tidft                        contact_atom_1 = domain_row
!tidft                        contact_atom_2 = domain_col
!tidft
!tidft                     ENDIF
!tidft
!tidft                  ENDIF ! end type of electron groups
!tidft
!tidft                  ! get atomic radii to compute distance cutoff threshold
!tidft                  IF (almo_scf_env%quencher_radius_type == do_bondparm_covalent) THEN
!tidft                     CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_1)%atomic_kind, &
!tidft                                          rcov=contact1_radius)
!tidft                     CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_2)%atomic_kind, &
!tidft                                          rcov=contact2_radius)
!tidft                  ELSE IF (almo_scf_env%quencher_radius_type == do_bondparm_vdw) THEN
!tidft                     CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_1)%atomic_kind, &
!tidft                                          rvdw=contact1_radius)
!tidft                     CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_2)%atomic_kind, &
!tidft                                          rvdw=contact2_radius)
!tidft                  ELSE
!tidft                     CPABORT("Illegal quencher_radius_type")
!tidft                  END IF
!tidft                  contact1_radius = cp_unit_to_cp2k(contact1_radius, "angstrom")
!tidft                  contact2_radius = cp_unit_to_cp2k(contact2_radius, "angstrom")
!tidft
!tidft                  !RZK-warning the procedure is faulty for molecules:
!tidft                  ! the closest contacts should be found using
!tidft                  ! the element specific radii
!tidft
!tidft                  ! compute inner and outer cutoff radii
!tidft                  r0 = almo_scf_env%quencher_r0_factor*(contact1_radius + contact2_radius)
!tidft                  !+almo_scf_env%quencher_r0_shift
!tidft                  r1 = almo_scf_env%quencher_r1_factor*(contact1_radius + contact2_radius)
!tidft                  !+almo_scf_env%quencher_r1_shift
!tidft
!tidft                  IF (r0 .LT. 0.0_dp) THEN
!tidft                     CPABORT("R0 is less than zero")
!tidft                  ENDIF
!tidft                  IF (r1 .LE. 0.0_dp) THEN
!tidft                     CPABORT("R1 is less than or equal to zero")
!tidft                  ENDIF
!tidft                  IF (r0 .GT. r1) THEN
!tidft                     CPABORT("R0 is greater than or equal to R1")
!tidft                  ENDIF
!tidft
!tidft                  ! Fill in non-zero blocks if AOs are close to the electron center
!tidft                  IF (distance .LT. r1) THEN
!tidft                     NULLIFY (p_new_block)
!tidft                     CALL dbcsr_reserve_block2d(almo_scf_env%quench_t(ispin), &
!tidft                                                iblock_row, iblock_col, p_new_block)
!tidft                     CPASSERT(ASSOCIATED(p_new_block))
!tidft                     IF (distance .LE. r0) THEN
!tidft                        p_new_block(:, :) = 1.0_dp
!tidft                        !WRITE(*,'(A15,2I7,5F8.3,E11.3)') "INTRA-BLOCKS: ",&
!tidft                        !   iblock_col, iblock_row, contact1_radius,&
!tidft                        !   contact2_radius, r0, r1, distance, p_new_block(1,1)
!tidft                     ELSE
!tidft                        ! remove the intermediate values from the quencher temporarily
!tidft                        CPABORT("")
!tidft                        p_new_block(:, :) = 1.0_dp/(1.0_dp + EXP((r1 - r0)/(r0 - distance) + (r1 - r0)/(r1 - distance)))
!tidft                        !WRITE(*,'(A15,2I7,5F8.3,E11.3)') "INTER-BLOCKS: ",&
!tidft                        !   iblock_col, iblock_row, contact1_radius,&
!tidft                        !   contact2_radius, r0, r1, distance, p_new_block(1,1)
!tidft                     ENDIF
!tidft
!tidft                     IF (ABS(p_new_block(1, 1)) .GT. ABS(almo_scf_env%eps_filter)) THEN
!tidft                        IF (domain_map_local_entries .GE. max_domain_neighbors*almo_scf_env%ndomains) THEN
!tidft                           CPABORT("weird... max_domain_neighbors is exceeded")
!tidft                        ENDIF
!tidft                        almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries + 1, 1) = iblock_row
!tidft                        almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries + 1, 2) = iblock_col
!tidft                        domain_map_local_entries = domain_map_local_entries + 1
!tidft                     ENDIF
!tidft
!tidft                  ENDIF
!tidft
!tidft               CASE DEFAULT
!tidft                  CPABORT("Illegal constraint type")
!tidft               END SELECT
!tidft
!tidft            ENDIF ! mynode
!tidft
!tidft         ENDDO
!tidft      ENDDO ! end O(N) loop over pairs
!tidft
!tidft      DEALLOCATE (domain_neighbor_list)
!tidft      DEALLOCATE (current_number_neighbors)
!tidft
!tidft      CALL dbcsr_finalize(almo_scf_env%quench_t(ispin))
!tidft      !CALL dbcsr_scale(almo_scf_env%quench_t(ispin),&
!tidft      !        almo_scf_env%envelope_amplitude)
!tidft      CALL dbcsr_filter(almo_scf_env%quench_t(ispin), &
!tidft                        almo_scf_env%eps_filter)
!tidft
!tidft      ! check that both domain_map and quench_t have the same number of entries
!tidft      nblks = dbcsr_get_num_blocks(almo_scf_env%quench_t(ispin))
!tidft      IF (nblks .NE. domain_map_local_entries) THEN
!tidft         CPABORT("number of blocks is wrong")
!tidft      ENDIF
!tidft
!tidft      ! first, communicate map sizes on the other nodes
!tidft      ALLOCATE (domain_entries_cpu(nNodes), offset_for_cpu(nNodes))
!tidft      CALL mp_allgather(2*domain_map_local_entries, domain_entries_cpu, GroupID)
!tidft
!tidft      ! second, create
!tidft      offset_for_cpu(1) = 0
!tidft      DO iNode = 2, nNodes
!tidft         offset_for_cpu(iNode) = offset_for_cpu(iNode - 1) + &
!tidft                                 domain_entries_cpu(iNode - 1)
!tidft      ENDDO
!tidft      global_entries = offset_for_cpu(nNodes) + domain_entries_cpu(nNodes)
!tidft
!tidft      ! communicate all entries
!tidft      ALLOCATE (domain_map_global(global_entries))
!tidft      ALLOCATE (domain_map_local(2*domain_map_local_entries))
!tidft      DO ientry = 1, domain_map_local_entries
!tidft         domain_map_local(2*(ientry - 1) + 1) = almo_scf_env%domain_map(ispin)%pairs(ientry, 1)
!tidft         domain_map_local(2*ientry) = almo_scf_env%domain_map(ispin)%pairs(ientry, 2)
!tidft      ENDDO
!tidft      CALL mp_allgather(domain_map_local, domain_map_global, &
!tidft                        domain_entries_cpu, offset_for_cpu, GroupID)
!tidft      DEALLOCATE (domain_entries_cpu, offset_for_cpu)
!tidft      DEALLOCATE (domain_map_local)
!tidft
!tidft      DEALLOCATE (almo_scf_env%domain_map(ispin)%index1)
!tidft      DEALLOCATE (almo_scf_env%domain_map(ispin)%pairs)
!tidft      ALLOCATE (almo_scf_env%domain_map(ispin)%index1(ndomains))
!tidft      ALLOCATE (almo_scf_env%domain_map(ispin)%pairs(global_entries/2, 2))
!tidft      almo_scf_env%domain_map(ispin)%pairs(:, :) = 0
!tidft      almo_scf_env%domain_map(ispin)%index1(:) = 0
!tidft
!tidft      ! unpack the received data into a local variable
!tidft      ! since we do not know the maximum global number of neighbors
!tidft      ! try one. if fails increase the maximum number and try again
!tidft      ! until it succeeds
!tidft      max_neig = max_domain_neighbors
!tidft      max_neig_fails = .TRUE.
!tidft      max_neig_loop: DO WHILE (max_neig_fails)
!tidft         ALLOCATE (domain_grid(almo_scf_env%ndomains, 0:max_neig))
!tidft         domain_grid(:, :) = 0
!tidft         ! init the number of collected neighbors
!tidft         domain_grid(:, 0) = 1
!tidft         ! loop over the records
!tidft         global_entries = global_entries/2
!tidft         DO ientry = 1, global_entries
!tidft            ! get the center
!tidft            grid1 = domain_map_global(2*ientry)
!tidft            ! get the neighbor
!tidft            ineig = domain_map_global(2*(ientry - 1) + 1)
!tidft            ! check boundaries
!tidft            IF (domain_grid(grid1, 0) .GT. max_neig) THEN
!tidft               ! this neighbor will overstep the boundaries
!tidft               ! stop the trial and increase the max number of neighbors
!tidft               DEALLOCATE (domain_grid)
!tidft               max_neig = max_neig*2
!tidft               CYCLE max_neig_loop
!tidft            ENDIF
!tidft            ! for the current center loop over the collected neighbors
!tidft            ! to insert the current record in a numerical order
!tidft            delayed_increment = .FALSE.
!tidft            DO igrid = 1, domain_grid(grid1, 0)
!tidft               ! compare the current neighbor with that already in the 'book'
!tidft               IF (ineig .LT. domain_grid(grid1, igrid)) THEN
!tidft                  ! if this one is smaller then insert it here and pick up the one
!tidft                  ! from the book to continue inserting
!tidft                  neig_temp = ineig
!tidft                  ineig = domain_grid(grid1, igrid)
!tidft                  domain_grid(grid1, igrid) = neig_temp
!tidft               ELSE
!tidft                  IF (domain_grid(grid1, igrid) .EQ. 0) THEN
!tidft                     ! got the empty slot now - insert the record
!tidft                     domain_grid(grid1, igrid) = ineig
!tidft                     ! increase the record counter but do it outside the loop
!tidft                     delayed_increment = .TRUE.
!tidft                  ENDIF
!tidft               ENDIF
!tidft            ENDDO
!tidft            IF (delayed_increment) THEN
!tidft               domain_grid(grid1, 0) = domain_grid(grid1, 0) + 1
!tidft            ELSE
!tidft               ! should not be here - all records must be inserted
!tidft               CPABORT("all records must be inserted")
!tidft            ENDIF
!tidft         ENDDO
!tidft         max_neig_fails = .FALSE.
!tidft      ENDDO max_neig_loop
!tidft      DEALLOCATE (domain_map_global)
!tidft
!tidft      ientry = 1
!tidft      DO idomain = 1, almo_scf_env%ndomains
!tidft         DO ineig = 1, domain_grid(idomain, 0) - 1
!tidft            almo_scf_env%domain_map(ispin)%pairs(ientry, 1) = domain_grid(idomain, ineig)
!tidft            almo_scf_env%domain_map(ispin)%pairs(ientry, 2) = idomain
!tidft            ientry = ientry + 1
!tidft         ENDDO
!tidft         almo_scf_env%domain_map(ispin)%index1(idomain) = ientry
!tidft      ENDDO
!tidft      DEALLOCATE (domain_grid)
!tidft
!tidft      !ENDDO ! ispin
!tidft      IF (almo_scf_env%nspins .EQ. 2) THEN
!tidft         CALL dbcsr_copy(almo_scf_env%quench_t(2), &
!tidft                         almo_scf_env%quench_t(1))
!tidft         almo_scf_env%domain_map(2)%pairs(:, :) = &
!tidft            almo_scf_env%domain_map(1)%pairs(:, :)
!tidft         almo_scf_env%domain_map(2)%index1(:) = &
!tidft            almo_scf_env%domain_map(1)%index1(:)
!tidft      ENDIF
!tidft
!tidft      CALL dbcsr_release(matrix_s_sym)
!tidft
!tidft      IF (almo_scf_env%domain_layout_mos == almo_domain_layout_molecular .OR. &
!tidft          almo_scf_env%domain_layout_aos == almo_domain_layout_molecular) THEN
!tidft         DEALLOCATE (first_atom_of_molecule)
!tidft         DEALLOCATE (last_atom_of_molecule)
!tidft      ENDIF
!tidft
!tidft      !mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
!tidft      !   dbcsr_distribution(almo_scf_env%quench_t(ispin))))
!tidft      !nblkrows_tot = dbcsr_nblkrows_total(almo_scf_env%quench_t(ispin))
!tidft      !nblkcols_tot = dbcsr_nblkcols_total(almo_scf_env%quench_t(ispin))
!tidft      !DO row = 1, nblkrows_tot
!tidft      !   DO col = 1, nblkcols_tot
!tidft      !      tr = .FALSE.
!tidft      !      iblock_row = row
!tidft      !      iblock_col = col
!tidft      !      CALL dbcsr_get_stored_coordinates(almo_scf_env%quench_t(ispin),&
!tidft      !              iblock_row, iblock_col, tr, hold)
!tidft      !      CALL dbcsr_get_block_p(almo_scf_env%quench_t(ispin),&
!tidft      !              row, col, p_new_block, found)
!tidft      !      write(*,*) "RST_NOTE:", mynode, row, col, hold, found
!tidft      !   ENDDO
!tidft      !ENDDO
!tidft
!tidft      CALL timestop(handle)
!tidft
!tidft   END SUBROUTINE almo_scf_construct_quencher
!tidft
!tidft! *****************************************************************************
!tidft!> \brief Compute matrix W (energy-weighted density matrix) that is needed
!tidft!>        for the evaluation of forces
!tidft!> \param matrix_w ...
!tidft!> \param almo_scf_env ...
!tidft!> \par History
!tidft!>       2015.03 created [Rustam Z. Khaliullin]
!tidft!> \author Rustam Z. Khaliullin
!tidft! **************************************************************************************************
!tidft   SUBROUTINE calculate_w_matrix_almo(matrix_w, almo_scf_env)
!tidft      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_w
!tidft      TYPE(almo_scf_env_type)                            :: almo_scf_env
!tidft
!tidft      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_almo'
!tidft
!tidft      INTEGER                                            :: handle, ispin
!tidft      REAL(KIND=dp)                                      :: scaling
!tidft      TYPE(dbcsr_type)                                   :: tmp_nn1, tmp_no1, tmp_oo1, tmp_oo2
!tidft
!tidft      CALL timeset(routineN, handle)
!tidft
!tidft      IF (almo_scf_env%nspins == 1) THEN
!tidft         scaling = 2.0_dp
!tidft      ELSE
!tidft         scaling = 1.0_dp
!tidft      ENDIF
!tidft
!tidft      DO ispin = 1, almo_scf_env%nspins
!tidft
!tidft         CALL dbcsr_create(tmp_nn1, template=almo_scf_env%matrix_s(1), &
!tidft                           matrix_type=dbcsr_type_no_symmetry)
!tidft         CALL dbcsr_create(tmp_no1, template=almo_scf_env%matrix_t(ispin), &
!tidft                           matrix_type=dbcsr_type_no_symmetry)
!tidft         CALL dbcsr_create(tmp_oo1, template=almo_scf_env%matrix_sigma_inv(ispin), &
!tidft                           matrix_type=dbcsr_type_no_symmetry)
!tidft         CALL dbcsr_create(tmp_oo2, template=almo_scf_env%matrix_sigma_inv(ispin), &
!tidft                           matrix_type=dbcsr_type_no_symmetry)
!tidft
!tidft         CALL dbcsr_copy(tmp_nn1, almo_scf_env%matrix_ks(ispin))
!tidft         ! 1. TMP_NO1=F.T
!tidft         CALL dbcsr_multiply("N", "N", scaling, tmp_nn1, almo_scf_env%matrix_t(ispin), &
!tidft                             0.0_dp, tmp_no1, filter_eps=almo_scf_env%eps_filter)
!tidft         ! 2. TMP_OO1=T^(tr).TMP_NO1=T^(tr).(FT)
!tidft         CALL dbcsr_multiply("T", "N", 1.0_dp, almo_scf_env%matrix_t(ispin), tmp_no1, &
!tidft                             0.0_dp, tmp_oo1, filter_eps=almo_scf_env%eps_filter)
!tidft         ! 3. TMP_OO2=TMP_OO1.siginv=(T^(tr)FT).siginv
!tidft         CALL dbcsr_multiply("N", "N", 1.0_dp, tmp_oo1, almo_scf_env%matrix_sigma_inv(ispin), &
!tidft                             0.0_dp, tmp_oo2, filter_eps=almo_scf_env%eps_filter)
!tidft         ! 4. TMP_OO1=siginv.TMP_OO2=siginv.(T^(tr)FTsiginv)
!tidft         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_sigma_inv(ispin), tmp_oo2, &
!tidft                             0.0_dp, tmp_oo1, filter_eps=almo_scf_env%eps_filter)
!tidft         ! 5. TMP_NO1=T.TMP_OO1.=T.(siginvT^(tr)FTsiginv)
!tidft         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_t(ispin), tmp_oo1, &
!tidft                             0.0_dp, tmp_no1, filter_eps=almo_scf_env%eps_filter)
!tidft         ! 6. TMP_NN1=TMP_NO1.T^(tr)=(TsiginvT^(tr)FTsiginv).T^(tr)=RFR
!tidft         CALL dbcsr_multiply("N", "T", 1.0_dp, tmp_no1, almo_scf_env%matrix_t(ispin), &
!tidft                             0.0_dp, tmp_nn1, filter_eps=almo_scf_env%eps_filter)
!tidft         CALL matrix_almo_to_qs(tmp_nn1, matrix_w(ispin)%matrix, almo_scf_env%mat_distr_aos)
!tidft
!tidft         CALL dbcsr_release(tmp_nn1)
!tidft         CALL dbcsr_release(tmp_no1)
!tidft         CALL dbcsr_release(tmp_oo1)
!tidft         CALL dbcsr_release(tmp_oo2)
!tidft
!tidft      ENDDO
!tidft
!tidft      CALL timestop(handle)
!tidft
!tidft   END SUBROUTINE calculate_w_matrix_almo

END MODULE tidft_qs

