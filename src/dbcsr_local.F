MODULE dbcsr_local

   !! Modified DBCSR routines not yet available in the official release
   USE dbcsr_blas_operations, ONLY: set_larnv_seed
   USE dbcsr_api,                       ONLY: &
        dbcsr_complete_redistribute, dbcsr_copy, dbcsr_copy_into_existing, dbcsr_create, &
        dbcsr_desymmetrize, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_get_stored_coordinates, &
        dbcsr_multiply, dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_block2d, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_work_create, &
        dbcsr_put_block
   USE dbcsr_kinds, ONLY: dp, &
                          int_4, &
                          int_8, &
                          real_4, &
                          real_8, &
                          sp

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_local'

   PUBLIC :: dbcsr_init_random

CONTAINS

   SUBROUTINE dbcsr_init_random(matrix, keep_sparsity, mini_seed)
      !! ... TODO : unify with other version which is generic in the data_type
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, OPTIONAL                                  :: keep_sparsity
      INTEGER, INTENT(IN), OPTIONAL                      :: mini_seed
      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_random'

      INTEGER                                            :: col, col_size, handle, hold, iseed(4), &
                                                            mynode, ncol, nrow, row, row_size, &
                                                            stored_col, stored_row, my_mini_seed
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      LOGICAL                                            :: found, my_keep_sparsity, tr
      REAL(real_8), ALLOCATABLE, DIMENSION(:)            :: rnd
      REAL(real_8), DIMENSION(:, :), POINTER             :: buff, data_d
      TYPE(dbcsr_distribution_type)                      :: dist

!   ---------------------------------------------------------------------------

      my_keep_sparsity = .FALSE.
      IF (PRESENT(keep_sparsity)) my_keep_sparsity = keep_sparsity

      my_mini_seed = 1
      IF (PRESENT(mini_seed)) my_mini_seed = mini_seed

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      !row_blk_size => array_data(matrix%row_blk_size)
      !col_blk_size => array_data(matrix%col_blk_size)
      CALL dbcsr_get_info(matrix, distribution=dist)
      CALL dbcsr_distribution_get(dist, mynode=mynode)
      CALL dbcsr_work_create(matrix, work_mutable=.TRUE.)

      ALLOCATE (rnd(MAXVAL(row_blk_size)*MAXVAL(col_blk_size)))
      nrow = dbcsr_nblkrows_total(matrix)
      ncol = dbcsr_nblkcols_total(matrix)
      DO row = 1, nrow
      DO col = 1, ncol
         row_size = row_blk_size(row)
         col_size = col_blk_size(col)
         tr = .FALSE.
         stored_row = row
         stored_col = col
         CALL dbcsr_get_stored_coordinates(matrix, stored_row, stored_col, hold)
         IF (hold .EQ. mynode) THEN
            CALL dbcsr_get_block_p(matrix, stored_row, stored_col, data_d, tr, found)
            IF (found .OR. (.NOT. my_keep_sparsity)) THEN
               ! set the seed for dlarnv, is here to guarantee same value of the random numbers
               ! for all layouts (and block distributions)
               CALL set_larnv_seed(row, nrow, col, ncol, my_mini_seed, iseed)
               CALL dlarnv(1, iseed, row_size*col_size, rnd(1))
            ENDIF
            IF (found) THEN
               CALL dcopy(row_size*col_size, rnd, 1, data_d, 1)
            ELSE
               IF (.NOT. my_keep_sparsity) THEN
                  ALLOCATE (buff(row_size, col_size))
                  CALL dcopy(row_size*col_size, rnd, 1, buff, 1)
                  CALL dbcsr_put_block(matrix, stored_row, stored_col, buff)
                  DEALLOCATE (buff)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      DEALLOCATE (rnd)

      CALL dbcsr_finalize(matrix)
      CALL timestop(handle)

   END SUBROUTINE dbcsr_init_random

END MODULE dbcsr_local

