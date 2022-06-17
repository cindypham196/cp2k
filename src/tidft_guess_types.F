!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for TIDFT initial guess
!> \par History
!>       2020.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_guess_types
   USE dbcsr_api,                       ONLY: dbcsr_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_guess_types'

   PUBLIC :: tidft_excitation_list_type, &
             tidft_excitation_type, &
             orbital_set_type

   TYPE tidft_excitation_list_type

      TYPE(tidft_excitation_type), DIMENSION(:), &
         ALLOCATABLE                                     :: excitations
      INTEGER                                            :: max_order
      INTEGER                                            :: max_nexcitations
      INTEGER                                            :: cur_nexcitations

   END TYPE tidft_excitation_list_type

   TYPE tidft_excitation_type

      INTEGER                                            :: order
      REAL(KIND=dp)                                      :: excitation_energy
      INTEGER                                            :: unpairings_count
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: delta_epsilon
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: occupied
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: virtual
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: spin

   END TYPE tidft_excitation_type

   TYPE orbital_list_type

      INTEGER                                            :: norbs
      INTEGER                                            :: ispin
      INTEGER                                            :: ivirtual ! 0 - occ, 1 - virtual
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: orb_energy
      LOGICAL, DIMENSION(:), ALLOCATABLE                 :: orb_has_spin_double
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: index_within_matrix ! position within matrix
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: block_within_matrix
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: index_within_block ! position witihn block

   END TYPE orbital_list_type

   TYPE orbital_set_type

      TYPE(orbital_list_type), DIMENSION(:,:), &
         ALLOCATABLE                                     :: orbital_list

   END TYPE orbital_set_type

END MODULE tidft_guess_types

