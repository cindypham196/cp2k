!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for time-independent DFT optimizers
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_optimizer_types
   USE dbcsr_api,                       ONLY: dbcsr_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_optimizer_types'

   PUBLIC :: tidft_optimizer_storage_type

   ! the structure contains all data necessery to perform optimization
   TYPE tidft_optimizer_storage_type

      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: grad
      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: prev_theta
      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: prev_grad
      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: prev_step
      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: step
      TYPE(dbcsr_type), DIMENSION(:,:), ALLOCATABLE      :: prev_minus_prec_grad

   END TYPE tidft_optimizer_storage_type

END MODULE tidft_optimizer_types

