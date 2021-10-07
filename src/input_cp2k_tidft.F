!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief input for the TIDFT section
!> \author Rustam Khaliullin
! **************************************************************************************************
MODULE input_cp2k_tidft
   USE bibliography,                    ONLY: Khaliullin2007,&
                                              Khaliullin2008,&
                                              Khaliullin2013,&
                                              Scheiber2018
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              low_print_level
   USE input_constants,                 ONLY: &
        almo_deloc_none, almo_deloc_scf, almo_deloc_x, almo_deloc_x_then_scf, &
        almo_deloc_xalmo_1diag, almo_deloc_xalmo_scf, almo_deloc_xalmo_x, almo_frz_crystal, &
        almo_frz_none, almo_scf_diag, almo_scf_pcg, almo_scf_skip, almo_scf_trustr, atomic_guess, &
        cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, cg_hager_zhang, cg_hestenes_stiefel, &
        cg_liu_storey, cg_polak_ribiere, cg_zero, molecular_guess, op_loc_berry, op_loc_pipek, &
        optimizer_diis, optimizer_lin_eq_pcg, optimizer_pcg, optimizer_trustr, &
        spd_inversion_dense_cholesky, spd_inversion_ls_hotelling, spd_inversion_ls_taylor, &
        trustr_cauchy, trustr_dogleg, trustr_steihaug, xalmo_prec_domain, xalmo_prec_full, &
        xalmo_prec_zero, xalmo_trial_r0_out, xalmo_trial_simplex
   USE input_val_types,                 ONLY: char_t,&
                                              integer_t,&
                                              lchar_t,&
                                              logical_t,&
                                              real_t
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_tidft'

   INTEGER, PARAMETER, PRIVATE :: optimizer_block_diagonal_diis = 1
   INTEGER, PARAMETER, PRIVATE :: optimizer_block_diagonal_pcg = 2
   INTEGER, PARAMETER, PRIVATE :: optimizer_xalmo_pcg = 3
   INTEGER, PARAMETER, PRIVATE :: optimizer_xalmo_trustr = 4
   INTEGER, PARAMETER, PRIVATE :: optimizer_newton_pcg_solver = 5
   INTEGER, PARAMETER, PRIVATE :: optimizer_nlmo_pcg = 6

   PUBLIC :: create_tidft_section

CONTAINS

! **************************************************************************************************
!> \brief create the TIDFT section
!> \param section ...
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Cindy Pham, Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_tidft_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_tidft_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="TIDFT", &
                          description="Time-independent density functional theory.", &
                          n_keywords=30, n_subsections=10, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Matrix sparsity filter.", &
                          usage="EPS_FILTER 1.e-6", default_r_val=1.e-9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_EXCITED_STATES", &
                          description="Specifies the number of excited states (same multiplicity as the ground state).", &
                          usage="N_EXCITED_STATES 3", default_i_val=2)
      !CALL keyword_create(keyword, __LOCATION__, name="N_EXCITED_STATES", &
      !                    description="Specifies the number of excited states: singlet, doublet, triplet, ...", &
      !                    usage="N_EXCITED_STATES 3 1 4 0", &
      !                    n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_UNPAIR", &
                          description="Specifies the number of allowed electron unpairings in ROKS", &
                          usage="MAX_UNPAIR 0", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_PURE", &
                          description="Obtain spin pure electronic states.", &
                          usage="SPIN_PURE .TRUE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      !tiidftCALL keyword_create(keyword, __LOCATION__, name="ALMO_SCF_GUESS", &
      !tiidft                    description="The method to generate initial ALMOs.", &
      !tiidft                    usage="ALMO_SCF_GUESS MOLECULAR", &
      !tiidft                    default_i_val=atomic_guess, &
      !tiidft                    enum_c_vals=s2a("MOLECULAR", "ATOMIC"), &
      !tiidft                    enum_desc=s2a("SCF calculations on single molecules controlled by the regular SCF "// &
      !tiidft                                  "keywords outside ALMO options. This kind of calculation is expensive "// &
      !tiidft                                  "and only recommended if ALMO SCF does not converge from the ATOMIC guess.", &
      !tiidft                                  "Superposition of atomic densities."), &
      !tiidft                    enum_i_vals=(/molecular_guess, atomic_guess/))
      !tiidftCALL section_add_keyword(section, keyword)
      !tiidftCALL keyword_release(keyword)

      !tiidftCALL keyword_create(keyword, __LOCATION__, name="MO_OVERLAP_INV_ALG", &
      !tiidft                    description="Algorithm to invert MO overlap matrix.", &
      !tiidft                    usage="MO_OVERLAP_INV_ALG LS_HOTELLING", &
      !tiidft                    default_i_val=spd_inversion_ls_hotelling, &
      !tiidft                    enum_c_vals=s2a("LS_HOTELLING", "LS_TAYLOR", "DENSE_CHOLESKY"), &
      !tiidft                    enum_desc=s2a("Linear scaling iterative Hotelling algorithm. Fast for large sparse matrices.", &
      !tiidft                                  "Linear scaling algorithm based on Taylor expansion of (A+B)^(-1).", &
      !tiidft                                  "Stable but dense Cholesky algorithm. Cubically scaling."), &
      !tiidft                    enum_i_vals=(/spd_inversion_ls_hotelling, spd_inversion_ls_taylor, spd_inversion_dense_cholesky/))
      !tiidftCALL section_add_keyword(section, keyword)
      !tiidftCALL keyword_release(keyword)

      !tiidftCALL keyword_create(keyword, __LOCATION__, name="ALMO_EXTRAPOLATION_ORDER", &
      !tiidft                    description="Number of previous states used for the ASPC extrapolation of the ALMO "// &
      !tiidft                    "initial guess. 0 implies that the guess is given by ALMO_SCF_GUESS at each step.", &
      !tiidft                    usage="ALMO_EXTRAPOLATION_ORDER 3", default_i_val=3)
      !tiidftCALL section_add_keyword(section, keyword)
      !tiidftCALL keyword_release(keyword)

      !tiidftCALL keyword_create(keyword, __LOCATION__, name="XALMO_ALGORITHM", &
      !tiidft                    description="Specifies the algorithm to update ALMOs on eXtended domains (XALMOs).", &
      !tiidft                    usage="XALMO_ALGORITHM TRUST_REGION", &
      !tiidft                    default_i_val=almo_scf_pcg, &
      !tiidft                    enum_c_vals=s2a("DIAG", "PCG", "TRUST_REGION"), &
      !tiidft                    enum_desc=s2a("DIIS-accelerated diagonalization.", &
      !tiidft                                  "PCG algorithm controlled by XALMO_OPTIMIZER_PCG.", &
      !tiidft                                  "Trust-region optimizer controlled by XALMO_OPTIMIZER_TRUSTR"), &
      !tiidft                    enum_i_vals=(/almo_scf_diag, almo_scf_pcg, almo_scf_trustr/))
      !tiidftCALL section_add_keyword(section, keyword)
      !tiidftCALL keyword_release(keyword)

      !tiidftCALL keyword_create(keyword, __LOCATION__, name="RETURN_ORTHOGONALIZED_MOS", &
      !tiidft                    description="Orthogonalize Kohn-Sham orbitls before they are returned"// &
      !tiidft                    "to Quickstep (i.e. for calculation of properties)", &
      !tiidft                    usage="RETURN_ORTHOGONALIZED_MOS .TRUE.", default_l_val=.TRUE., &
      !tiidft                    lone_keyword_l_val=.TRUE.)
      !tiidftCALL section_add_keyword(section, keyword)
      !tiidftCALL keyword_release(keyword)

      !tiidftNULLIFY (subsection)
      !tiidftCALL create_optimizer_section(subsection, optimizer_block_diagonal_diis)
      !tiidftCALL section_add_subsection(section, subsection)
      !tiidftCALL section_release(subsection)

      !tiidftNULLIFY (subsection)
      !tiidftCALL create_optimizer_section(subsection, optimizer_block_diagonal_pcg)
      !tiidftCALL section_add_subsection(section, subsection)
      !tiidftCALL section_release(subsection)

      !tiidftNULLIFY (subsection)
      !tiidftCALL create_optimizer_section(subsection, optimizer_xalmo_pcg)
      !tiidftCALL section_add_subsection(section, subsection)
      !tiidftCALL section_release(subsection)

      !tiidftNULLIFY (subsection)
      !tiidftCALL create_optimizer_section(subsection, optimizer_xalmo_trustr)
      !tiidftCALL section_add_subsection(section, subsection)
      !tiidftCALL section_release(subsection)

      !tiidftNULLIFY (subsection)
      !tiidftCALL create_optimizer_section(subsection, optimizer_nlmo_pcg)
      !tiidftCALL section_add_subsection(section, subsection)
      !tiidftCALL section_release(subsection)

      !tiidftNULLIFY (subsection)
      !tiidftCALL create_matrix_iterate_section(subsection)
      !tiidftCALL section_add_subsection(section, subsection)
      !tiidftCALL section_release(subsection)

   END SUBROUTINE create_tidft_section

!tiidft! **************************************************************************************************
!tiidft!> \brief The optimizer section is a collection of keywords that are similar
!tiidft!>        to all optimization methods (e.g. target error, number of iterations)
!tiidft!> \param section ...
!tiidft!> \param optimizer_id   allows to adapt the standard section for specific needs
!tiidft!> \par History
!tiidft!>       2012.03 created [Rustam Z Khaliullin]
!tiidft!>       2014.10 fully integrated [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   RECURSIVE SUBROUTINE create_optimizer_section(section, optimizer_id)
!tiidft
!tiidft      TYPE(section_type), POINTER                        :: section
!tiidft      INTEGER, INTENT(IN)                                :: optimizer_id
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'create_optimizer_section', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      INTEGER                                            :: optimizer_type
!tiidft      TYPE(keyword_type), POINTER                        :: keyword
!tiidft      TYPE(section_type), POINTER                        :: subsection
!tiidft
!tiidft      CPASSERT(.NOT. ASSOCIATED(section))
!tiidft      NULLIFY (section)
!tiidft
!tiidft      ! choose the name of the section
!tiidft      SELECT CASE (optimizer_id)
!tiidft      CASE (optimizer_block_diagonal_diis)
!tiidft         CALL section_create(section, __LOCATION__, name="ALMO_OPTIMIZER_DIIS", &
!tiidft                             description="Controls the iterative DIIS-accelerated optimization of block-diagonal ALMOs.", &
!tiidft                             n_keywords=5, n_subsections=0, repeats=.FALSE.)
!tiidft         optimizer_type = optimizer_diis
!tiidft      CASE (optimizer_block_diagonal_pcg)
!tiidft         CALL section_create(section, __LOCATION__, name="ALMO_OPTIMIZER_PCG", &
!tiidft                             description="Controls the PCG optimization of block-diagonal ALMOs.", &
!tiidft                             n_keywords=9, n_subsections=1, repeats=.FALSE.)
!tiidft         optimizer_type = optimizer_pcg
!tiidft      CASE (optimizer_nlmo_pcg)
!tiidft         CALL section_create(section, __LOCATION__, name="NLMO_OPTIMIZER_PCG", &
!tiidft                             description="Controls the PCG optimization of nonorthogonal localized MOs.", &
!tiidft                             n_keywords=9, n_subsections=1, repeats=.FALSE.)
!tiidft         optimizer_type = optimizer_pcg
!tiidft         NULLIFY (subsection)
!tiidft         CALL create_penalty_section(subsection)
!tiidft         CALL section_add_subsection(section, subsection)
!tiidft         CALL section_release(subsection)
!tiidft      CASE (optimizer_xalmo_pcg)
!tiidft         CALL section_create(section, __LOCATION__, name="XALMO_OPTIMIZER_PCG", &
!tiidft                             description="Controls the PCG optimization of extended ALMOs.", &
!tiidft                             n_keywords=10, n_subsections=2, repeats=.FALSE.)
!tiidft         NULLIFY (subsection)
!tiidft         CALL create_optimizer_section(subsection, optimizer_newton_pcg_solver)
!tiidft         CALL section_add_subsection(section, subsection)
!tiidft         CALL section_release(subsection)
!tiidft         optimizer_type = optimizer_pcg
!tiidft      CASE (optimizer_xalmo_trustr)
!tiidft         CALL section_create(section, __LOCATION__, name="XALMO_OPTIMIZER_TRUSTR", &
!tiidft                             description="Controls the trust-region optimization of extended ALMOs. "// &
!tiidft                             "Trust radius is varied in the outer loop. Once the trust radius is "// &
!tiidft                             "chosen (and fixed) the model function can be minized using various "// &
!tiidft                             "approaches. Currently, an iterative conjugate-gradient approach is "// &
!tiidft                             "used and controlled by the inner loop", &
!tiidft                             n_keywords=10, n_subsections=0, repeats=.FALSE.)
!tiidft         optimizer_type = optimizer_trustr
!tiidft      CASE (optimizer_newton_pcg_solver)
!tiidft         CALL section_create(section, __LOCATION__, name="XALMO_NEWTON_PCG_SOLVER", &
!tiidft                             description="Controls an iterative solver of the Newton-Raphson linear equation.", &
!tiidft                             n_keywords=4, n_subsections=0, repeats=.FALSE.)
!tiidft         optimizer_type = optimizer_lin_eq_pcg
!tiidft      CASE DEFAULT
!tiidft         CPABORT("No default values allowed")
!tiidft      END SELECT
!tiidft
!tiidft      NULLIFY (keyword)
!tiidft
!tiidft      ! add common keywords
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
!tiidft                          description="Maximum number of iterations", &
!tiidft                          usage="MAX_ITER 100", default_i_val=20)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="EPS_ERROR", &
!tiidft                          description="Target value of the MAX norm of the error", &
!tiidft                          usage="EPS_ERROR 1.E-6", default_r_val=1.0E-5_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      ! add common keywords
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_EARLY", &
!tiidft                          description="Maximum number of iterations for truncated SCF "// &
!tiidft                          "(e.g. Langevin-corrected MD). Negative values mean that this keyword is not used.", &
!tiidft                          usage="MAX_ITER_EARLY 5", default_i_val=-1)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="EPS_ERROR_EARLY", &
!tiidft                          description="Target value of the MAX norm of the error for truncated SCF "// &
!tiidft                          "(e.g. Langevin-corrected MD). Negative values mean that this keyword is not used.", &
!tiidft                          usage="EPS_ERROR_EARLY 1.E-2", default_r_val=-1.0_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      ! add keywords specific to each type
!tiidft      IF (optimizer_type .EQ. optimizer_diis) THEN
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="N_DIIS", &
!tiidft                             description="Number of error vectors to be used in the DIIS "// &
!tiidft                             "optimization procedure", &
!tiidft                             usage="N_DIIS 5", default_i_val=6)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft      IF (optimizer_type .EQ. optimizer_pcg .OR. &
!tiidft          optimizer_type .EQ. optimizer_lin_eq_pcg .OR. &
!tiidft          optimizer_type .EQ. optimizer_trustr) THEN
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_OUTER_LOOP", &
!tiidft                             description="Maximum number of iterations in the outer loop. "// &
!tiidft                             "Use the outer loop to update the preconditioner and reset the conjugator. "// &
!tiidft                             "This can speed up convergence significantly.", &
!tiidft                             usage="MAX_ITER 10", default_i_val=0)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="PRECONDITIONER", &
!tiidft                             description="Select a preconditioner for the conjugate gradient optimization", &
!tiidft                             usage="PRECONDITIONER DOMAIN", &
!tiidft                             default_i_val=xalmo_prec_domain, &
!tiidft                             enum_c_vals=s2a("NONE", "DEFAULT", "DOMAIN", "FULL"), &
!tiidft                             enum_desc=s2a("Do not use preconditioner", &
!tiidft                                           "Same as DOMAIN preconditioner", &
!tiidft                                           "Invert preconditioner domain-by-domain."// &
!tiidft                                           " The main component of the linear scaling algorithm", &
!tiidft                                           "Solve linear equations step=-H.grad on the entire space"), &
!tiidft                             enum_i_vals=(/xalmo_prec_zero, xalmo_prec_domain, &
!tiidft                                           xalmo_prec_domain, xalmo_prec_full/))
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft      IF (optimizer_type .EQ. optimizer_pcg) THEN
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="LIN_SEARCH_EPS_ERROR", &
!tiidft                             description="Target value of the gradient norm during the linear search", &
!tiidft                             usage="LIN_SEARCH_EPS_ERROR 1.E-2", default_r_val=1.0E-3_dp)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="LIN_SEARCH_STEP_SIZE_GUESS", &
!tiidft                             description="The size of the first step in the linear search", &
!tiidft                             usage="LIN_SEARCH_STEP_SIZE_GUESS 0.1", default_r_val=1.0_dp)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="PRECOND_FILTER_THRESHOLD", &
!tiidft                             description="Select eigenvalues of the preconditioner "// &
!tiidft                             "that are smaller than the threshold and project out the "// &
!tiidft                             "corresponding eigenvectors from the gradient. No matter "// &
!tiidft                             "how large the threshold is the maximum number of projected "// &
!tiidft                             "eienvectors for a fragment equals to the number of occupied "// &
!tiidft                             "orbitals of fragment's neighbors.", &
!tiidft                             usage="PRECOND_FILTER_THRESHOLD 0.1", default_r_val=0.5_dp)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft      IF (optimizer_type .EQ. optimizer_pcg .OR. &
!tiidft          optimizer_type .EQ. optimizer_trustr) THEN
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="CONJUGATOR", &
!tiidft                             description="Various methods to compute step directions in the PCG optimization", &
!tiidft                             usage="CONJUGATOR POLAK_RIBIERE", &
!tiidft                             default_i_val=cg_hager_zhang, &
!tiidft                             enum_c_vals=s2a("ZERO", "POLAK_RIBIERE", "FLETCHER_REEVES", &
!tiidft                                             "HESTENES_STIEFEL", "FLETCHER", "LIU_STOREY", "DAI_YUAN", "HAGER_ZHANG"), &
!tiidft                             enum_desc=s2a("Steepest descent", "Polak and Ribiere", &
!tiidft                                           "Fletcher and Reeves", "Hestenes and Stiefel", &
!tiidft                                           "Fletcher (Conjugate descent)", "Liu and Storey", &
!tiidft                                           "Dai and Yuan", "Hager and Zhang"), &
!tiidft                             enum_i_vals=(/cg_zero, cg_polak_ribiere, cg_fletcher_reeves, &
!tiidft                                           cg_hestenes_stiefel, cg_fletcher, cg_liu_storey, &
!tiidft                                           cg_dai_yuan, cg_hager_zhang/))
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft      IF (optimizer_type .EQ. optimizer_trustr) THEN
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="ALGORITHM", &
!tiidft                             description="Selects an algorithm to solve the fixed-radius subproblem", &
!tiidft                             usage="ALGORITHM CG", &
!tiidft                             default_i_val=trustr_cauchy, &
!tiidft                             enum_c_vals=s2a("CG", "CAUCHY", "DOGLEG"), &
!tiidft                             enum_desc=s2a("Steihaug's iterative CG algorithm that does not invert model Hessian", &
!tiidft                                           "Compute simple Cauchy point", &
!tiidft                                           "Dogleg optimizer"), &
!tiidft                             enum_i_vals=(/trustr_steihaug, trustr_cauchy, trustr_dogleg/))
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="ETA", &
!tiidft                             description="Must be between 0.0 and 0.25. Rho value below which the "// &
!tiidft                             "optimization of the model function is not accepted and the "// &
!tiidft                             "optimization is restarted from the same point but decreased "// &
!tiidft                             "trust radius. Rho is the ratio of the actual over predicted "// &
!tiidft                             "change in the objective function", &
!tiidft                             usage="ETA 0.1", default_r_val=0.25_dp)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="MODEL_GRAD_NORM_RATIO", &
!tiidft                             description="Stop the fixed-trust-radius (inner) loop optimization "// &
!tiidft                             "once the ratio of the current norm of the model gradient over the "// &
!tiidft                             "initial norm drops below this threshold", &
!tiidft                             usage="MODEL_GRAD_NORM_RATIO 1.E-2", default_r_val=0.01_dp)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="INITIAL_TRUST_RADIUS", &
!tiidft                             description="Initial trust radius", &
!tiidft                             usage="INITIAL_TRUST_RADIUS 0.1", default_r_val=0.1_dp)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft         CALL keyword_create(keyword, __LOCATION__, name="MAX_TRUST_RADIUS", &
!tiidft                             description="Maximum allowed trust radius", &
!tiidft                             usage="MAX_TRUST_RADIUS 1.0", default_r_val=2.0_dp)
!tiidft         CALL section_add_keyword(section, keyword)
!tiidft         CALL keyword_release(keyword)
!tiidft
!tiidft      ENDIF
!tiidft
!tiidft   END SUBROUTINE create_optimizer_section
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief The section controls iterative matrix operations like SQRT or inverse
!tiidft!> \param section ...
!tiidft!> \par History
!tiidft!>       2017.05 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE create_matrix_iterate_section(section)
!tiidft
!tiidft      TYPE(section_type), POINTER                        :: section
!tiidft
!tiidft      TYPE(keyword_type), POINTER                        :: keyword
!tiidft
!tiidft      CPASSERT(.NOT. ASSOCIATED(section))
!tiidft      NULLIFY (section)
!tiidft
!tiidft      CALL section_create(section, __LOCATION__, name="MATRIX_ITERATE", &
!tiidft                          description="Controls linear scaling iterative procedure on matrices: inversion, sqrti, etc. "// &
!tiidft                          "High-order Lanczos accelerates convergence provided it can estimate the eigenspectrum correctly.", &
!tiidft                          n_keywords=4, n_subsections=0, repeats=.FALSE.)
!tiidft
!tiidft      NULLIFY (keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="EPS_TARGET_FACTOR", &
!tiidft                          description="Multiplication factor that determines acceptable error in the iterative procedure. "// &
!tiidft                          "Acceptable error = EPS_TARGET_FACTOR * EPS_FILTER", &
!tiidft                          usage="EPS_TARGET_FACTOR 100.0", default_r_val=10.0_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="EPS_LANCZOS", &
!tiidft                          description="Threshold for Lanczos eigenvalue estimation.", &
!tiidft                          usage="EPS_LANCZOS 1.0E-4", default_r_val=1.0E-3_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="ORDER_LANCZOS", &
!tiidft                          description="Order of the Lanczos estimator. Use 0 to turn off. Do not use 1.", &
!tiidft                          usage="ORDER_LANCZOS 5", default_i_val=3)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_LANCZOS", &
!tiidft                          description="Maximum number of Lanczos iterations.", &
!tiidft                          usage="MAX_ITER_LANCZOS 64", default_i_val=128)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft   END SUBROUTINE create_matrix_iterate_section
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief The section controls penalty methods
!tiidft!> \param section ...
!tiidft!> \par History
!tiidft!>       2018.01 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE create_penalty_section(section)
!tiidft
!tiidft      TYPE(section_type), POINTER                        :: section
!tiidft
!tiidft      TYPE(keyword_type), POINTER                        :: keyword
!tiidft
!tiidft      CPASSERT(.NOT. ASSOCIATED(section))
!tiidft      NULLIFY (section)
!tiidft
!tiidft      CALL section_create(section, __LOCATION__, name="PENALTY", &
!tiidft                          description="Add penalty terms to the energy functional.", &
!tiidft                          n_keywords=3, n_subsections=0, repeats=.FALSE.)
!tiidft
!tiidft      NULLIFY (keyword)
!tiidft
!tiidft      CALL keyword_create( &
!tiidft         keyword, __LOCATION__, name="OPERATOR", &
!tiidft         description="Type of opertator which defines the spread functional", &
!tiidft         usage="OPERATOR PIPEK", &
!tiidft         enum_c_vals=s2a("BERRY", "PIPEK"), &
!tiidft         enum_i_vals=(/op_loc_berry, op_loc_pipek/), &
!tiidft         default_i_val=op_loc_berry)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="PENALTY_STRENGTH", &
!tiidft                          description="Strength of the orthogonalization penalty", &
!tiidft                          usage="PENALTY_STRENGTH 1.1", default_r_val=1.1_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="PENALTY_STRENGTH_DECREASE_FACTOR", &
!tiidft                          description="Factor that decreases the strength of the orthogonalization penalty.", &
!tiidft                          usage="PENALTY_STRENGTH_DECREASE_FACTOR 1.1", default_r_val=1.1_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="DETERMINANT_TOLERANCE", &
!tiidft                          description="Stop the optimization of the penalty strength if the determinant of the overlap "// &
!tiidft                          "changes less than this tolerance threshold.", &
!tiidft                          usage="DETERMINANT_TOLERANCE 1.0E-4", default_r_val=1.0E-3_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="FINAL_DETERMINANT", &
!tiidft                          description="The final determinant that obtained after optimization.", &
!tiidft                          usage="FINAL_DETERMINANT 0.1", default_r_val=0.1_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="COMPACTIFICATION_FILTER_START", &
!tiidft                          description="Set orbital coefficients with absolute value smaller than this value to zero.", &
!tiidft                          usage="COMPACTIFICATION_FILTER_START 1.e-6", default_r_val=-1.0_dp)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="VIRTUAL_NLMOS", &
!tiidft                          description="Localize virtual oribtals", &
!tiidft                          usage="VIRTUAL_NLMOS .TRUE.", default_l_val=.FALSE.)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft   END SUBROUTINE create_penalty_section
!tiidft
!tiidft! **************************************************************************************************
!tiidft!> \brief The section controls electronic structure analysis based on ALMOs
!tiidft!> \param section ...
!tiidft!> \par History
!tiidft!>       2014.10 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE create_almo_analysis_section(section)
!tiidft
!tiidft      TYPE(section_type), POINTER                        :: section
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_analysis_section', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      TYPE(keyword_type), POINTER                        :: keyword
!tiidft      TYPE(section_type), POINTER                        :: subsection
!tiidft
!tiidft      CPASSERT(.NOT. ASSOCIATED(section))
!tiidft      NULLIFY (section)
!tiidft
!tiidft      CALL section_create(section, __LOCATION__, name="ANALYSIS", &
!tiidft                          description="Controls electronic structure analysis based on ALMOs and XALMOs.", &
!tiidft                          n_keywords=2, n_subsections=1, repeats=.FALSE., &
!tiidft                          citations=(/Khaliullin2007, Khaliullin2008/))
!tiidft
!tiidft      NULLIFY (keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
!tiidft                          description="Activation of ALMO-based electronic structure analysis.", &
!tiidft                          usage="&ANALYSIS T", &
!tiidft                          default_l_val=.FALSE., &
!tiidft                          lone_keyword_l_val=.TRUE.)
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      CALL keyword_create(keyword, __LOCATION__, name="FROZEN_MO_ENERGY_TERM", &
!tiidft                          description="Perform calculations on single molecules to compute the frozen density term", &
!tiidft                          usage="FROZEN_MO_ENERGY_TERM SUBLATTICE", default_i_val=almo_frz_none, &
!tiidft                          !enum_c_vals=s2a("SKIP", "ISOLATED", "SUBLATTICE"),&
!tiidft                          enum_c_vals=s2a("SKIP", "SUBLATTICE"), &
!tiidft                          enum_desc=s2a("Do not compute the frozen MO energy term.", &
!tiidft                                        !"Use isolated gas-phase molecules as the reference.",&
!tiidft                                        "Use energies of single molecules in their positions in the crystal "// &
!tiidft                                        "cell as the reference. "// &
!tiidft                                        "This term does not have an interpretation if fragmetns are charged."), &
!tiidft                          enum_i_vals=(/almo_frz_none, &
!tiidft                                        !almo_frz_isolated,&
!tiidft                                        almo_frz_crystal/))
!tiidft      CALL section_add_keyword(section, keyword)
!tiidft      CALL keyword_release(keyword)
!tiidft
!tiidft      NULLIFY (subsection)
!tiidft      CALL create_almo_print_section(subsection)
!tiidft      CALL section_add_subsection(section, subsection)
!tiidft      CALL section_release(subsection)
!tiidft
!tiidft   END SUBROUTINE create_almo_analysis_section
!tiidft
!tiidft! *****************************************************************************
!tiidft!> \brief Create a section with all print keys for ALMO methods
!tiidft!> \param section ...
!tiidft!> \par History
!tiidft!>       2016.05 created [Rustam Z Khaliullin]
!tiidft!> \author Rustam Z Khaliullin
!tiidft! **************************************************************************************************
!tiidft   SUBROUTINE create_almo_print_section(section)
!tiidft
!tiidft      TYPE(section_type), POINTER                        :: section
!tiidft
!tiidft      CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_print_section', &
!tiidft         routineP = moduleN//':'//routineN
!tiidft
!tiidft      TYPE(section_type), POINTER                        :: print_key
!tiidft
!tiidft      CPASSERT(.NOT. ASSOCIATED(section))
!tiidft      NULLIFY (section, print_key)
!tiidft
!tiidft      CALL section_create(section, __LOCATION__, name="PRINT", &
!tiidft                          description="Controls printing of detailed ALMO decomosition analysis results", &
!tiidft                          n_keywords=0, n_subsections=2, repeats=.TRUE.)
!tiidft
!tiidft      ! create print keys within this subsection
!tiidft      CALL cp_print_key_section_create(print_key, __LOCATION__, "ALMO_EDA_CT", &
!tiidft                                       description="Controls printing of the electron transfer terms in "// &
!tiidft                                       "ALMO energy decomposition analysis. "// &
!tiidft                                       "File format: first column is the index of the electron acceptor "// &
!tiidft                                       "fragment, second -- is the index of the electron donor fragment, "// &
!tiidft                                       "third - the energy change (a.u.) associated with the electron transfer.", &
!tiidft                                       print_level=low_print_level, common_iter_levels=1, &
!tiidft                                       filename="EDA")
!tiidft      CALL section_add_subsection(section, print_key)
!tiidft      CALL section_release(print_key)
!tiidft
!tiidft      CALL cp_print_key_section_create(print_key, __LOCATION__, "ALMO_CTA", &
!tiidft                                       description="Controls printing of the electron transfer terms in "// &
!tiidft                                       "ALMO charge transfer analysis. "// &
!tiidft                                       "File format: first column is the index of the electron acceptor "// &
!tiidft                                       "fragment, second -- is the index of the electron donor fragment, "// &
!tiidft                                       "third - the change (a.u.) transferred between the two fragments.", &
!tiidft                                       print_level=low_print_level, common_iter_levels=1, &
!tiidft                                       filename="CTA")
!tiidft      CALL section_add_subsection(section, print_key)
!tiidft      CALL section_release(print_key)
!tiidft
!tiidft   END SUBROUTINE create_almo_print_section

END MODULE input_cp2k_tidft
