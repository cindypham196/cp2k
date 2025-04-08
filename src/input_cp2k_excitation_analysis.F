!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Analysis of excitated states
!> \par History
!>      07.2024 factored out from input_cp2k_properties_dft [Rustam Khaliullin]
!> \author RZK
! **************************************************************************************************
MODULE input_cp2k_excitation_analysis

   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level,&
                                              low_print_level,&
                                              medium_print_level,&
                                              silent_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        current_gauge_atom, current_gauge_r, current_gauge_r_and_step_func, &
        current_orb_center_atom, current_orb_center_box, current_orb_center_common, &
        current_orb_center_wannier, do_et_ddapc, do_full_density, do_no_et, do_spin_density, &
        gto_cartesian, gto_spherical, oe_gllb, oe_lb, oe_none, oe_saop, oe_shift, &
        ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, scan_x, scan_xy, &
        scan_xyz, scan_xz, scan_y, scan_yz, scan_z, tddfpt_dipole_berry, tddfpt_dipole_length, &
        tddfpt_dipole_velocity, tddfpt_kernel_full, tddfpt_kernel_none, tddfpt_kernel_stda, &
        use_mom_ref_coac, use_mom_ref_com, use_mom_ref_user, use_mom_ref_zero
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              integer_t,&
                                              lchar_t,&
                                              logical_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_excitation_analysis'

   PUBLIC :: create_nto_analysis_section, &
             create_mos_molden_section

CONTAINS

! **************************************************************************************************
!> \brief Section controls NTO printing
!> \param section ...
!> \par History
!>       2024.07 factored out [Rustam Z Khaliullin]
! **************************************************************************************************
   SUBROUTINE create_nto_analysis_section(section)

      TYPE(section_type), POINTER                        :: section
      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL cp_print_key_section_create(section, __LOCATION__, name="NTO_ANALYSIS", &
                                       description="Perform a natural transition orbital analysis.", &
                                       print_level=medium_print_level)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="THRESHOLD", &
                          description="Threshold for sum of NTO eigenvalues considered", &
                          usage="Threshold 0.95", &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.975_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBE_FILES", &
                          description="Print NTOs on Cube Files", &
                          usage="CUBE_FILES {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTENSITY_THRESHOLD", &
                          description="Threshold for oscillator strength to screen states.", &
                          usage="Intensity_threshold 0.01", &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATE_LIST", &
                          description="Specifies a list of states for the NTO calculations.", &
                          usage="STATE_LIST {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). Provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_nto_analysis_section

! **************************************************************************************************
!> \brief Section controls MOLDEN printing
!> \param section ...
!> \par History
!>       2024.07 factored out [Rustam Z Khaliullin]
! **************************************************************************************************
   SUBROUTINE create_mos_molden_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL cp_print_key_section_create(section, __LOCATION__, "MOS_MOLDEN", &
                                       description="Write the NTO in Molden file format, for visualisation.", &
                                       print_level=debug_print_level + 1, add_last=add_last_numeric, filename="MOS")

      NULLIFY(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NDIGITS", &
                          description="Specifies the number of significant digits retained. 3 is OK for visualization.", &
                          usage="NDIGITS {int}", &
                          default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GTO_KIND", &
                          description="Representation of Gaussian-type orbitals", &
                          default_i_val=gto_spherical, &
                          enum_c_vals=s2a("CARTESIAN", "SPHERICAL"), &
                          enum_desc=s2a( &
                          "Cartesian Gaussian orbitals. Use with caution", &
                          "Spherical Gaussian orbitals. Incompatible with VMD"), &
                          enum_i_vals=(/gto_cartesian, gto_spherical/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mos_molden_section

END MODULE input_cp2k_excitation_analysis
