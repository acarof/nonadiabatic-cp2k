!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for the MIXED environment
!> \par History
!>      10.2008 created [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
MODULE input_cp2k_mixed
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             low_print_level,&
                                             medium_print_level
  USE input_constants,                 ONLY: mix_coupled,&
                                             mix_generic,&
                                             mix_linear_combination,&
                                             mix_minimum,&
!                                            mix_restrained
!> CHANGE_AC
!> 15-12-11
!> To create a new way to combine the force
!> To get the adiabatic force 
                                             mix_restrained, &
                                             mix_adiab
!> End CHANGE_AC

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
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mixed'

  PUBLIC :: create_mix_section

CONTAINS

! *****************************************************************************
!> \brief Create the input section for MIXED.
!> \param section the section to create
!> \author fschiff
! *****************************************************************************
  SUBROUTINE create_mix_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mix_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: sub2section, sub3section, &
                                                subsection
    TYPE(section_type), POINTER              :: print_key


    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,name="MIXED",&
         description="This section contains all information to run with a hamiltonian "//&
                     "defined by a mixing of force_evals",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    NULLIFY(keyword, subsection)

    CALL keyword_create(keyword, name="MIXING_TYPE",&
         description="The type of mixing to be employed",&
         usage="MIXING_TYPE LINEAR_COMBINATION",&
         default_i_val=mix_linear_combination,&
         enum_c_vals=s2a("LINEAR_COMBINATION", &
                         "MINIMUM",&
                         "COUPLED",&
                         "RESTRAINT",&
!                        "GENMIX"),&
!> CHANGE_AC
!> 15-12-11
!> To create a new way to combine the force
!> To get the adiabatic force 
                        "GENMIX",&
                        "ADIABATIC"),&
!> End CHANGE_AC
         enum_desc=s2a("Linear combination of force envs (support only 2 force_evals)", &
                       "Use the force env with the minimum energy (support only 2 force_evals)",&
                       "Consider the force envs as a two state system with a given"//&
                       " coupling matrix element (support only 2 force_evals)",&
                       "Use the difference between the energy of the force envs as a"//&
                       " restraint on the first (support only 2 force_evals)",&
!                       "Defines a user-driven generica coupling (support for an unlimited number of force_eval)"),&
!         enum_i_vals=(/mix_linear_combination,mix_minimum,mix_coupled,mix_restrained,mix_generic/))
!> CHANGE_AC
!> 15-12-11
!> To create a new way to combine the force
!> To get the adiabatic force 
                       "Defines a user-driven generica coupling (support for an unlimited number of force_eval)",&
                       "Combines the diabatic forces to get the adiabatic forces"),&
         enum_i_vals=(/mix_linear_combination,mix_minimum,mix_coupled,mix_restrained,mix_generic, &
                       mix_adiab /))
!> End CHANGE_AC
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="GROUP_PARTITION",&
         description="gives the exact number of processors for each group."//&
         " If not specified processors allocated will be equally distributed for"//&
         " the specified subforce_eval, trying to build a number of groups equal to the"//&
         " number of subforce_eval specified.",&
         usage="group_partition  2 2 4 2 4 ", type_of_var=integer_t, n_var=-1)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NGROUPS",variants=(/"NGROUP"/),&
         description="Gives the wanted number of groups. If not specified the number"//&
         " of groups is set to the number of subforce_eval defined.",&
         usage="ngroups 4", type_of_var=integer_t)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! Double force_eval
    CALL section_create(subsection,name="LINEAR",&
         description="Linear combination between two force_eval:  F= lambda F1 + (1-lambda) F2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="LAMBDA",&
         description="Specify the mixing parameter lambda in the formula.",&
         usage="lambda <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL section_create(subsection,name="COUPLING",&
         description="Coupling between two force_eval: E=(E1+E2 - sqrt((E1-E2)**2+4*H12**2))/2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="COUPLING_PARAMETER",&
         description="Coupling parameter H12 used in the coupling",&
         usage="COUPLING_PARAMETER <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL section_create(subsection,name="RESTRAINT",&
         description="Restraint between two force_eval: E = E1 + k*(E1-E2-t)**2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="RESTRAINT_TARGET",&
         description="Target value of the restraint (t) ",&
         usage="RESTRAINT_TARGET <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RESTRAINT_STRENGTH",&
         description="Strength of the restraint (k) in "//&
         "k*(E1-E2-t)**2" ,&
         usage="RESTRAINT_STRENGTH <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    ! Adiabatic force eval
!> CHANGE_AC
!> 15-12-11
!> Create the input part for the Non Adiabatic Molecular Dynamics

    CALL section_create(subsection,name="ADIABATIC",&
         description="Defines the methods to obtain the adiabatic forces.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)

    CALL keyword_create(keyword, name="METHOD_PROPAGATION",&
         description="Defines the methods to obtain the adiabatic forces.",&
         usage="METHOD_PROPAGATION : CLASSICAL_PATH, FSSH",  type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)


    CALL keyword_create(keyword, name="ELECTR_PROPAGATION",&
         description="Defines basis and methods to propagate the wf.",&
         usage="DIABATIC_PROP ",  type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

!> 15-12-11
    CALL keyword_create(keyword, name="METHOD_COUPLING",&
         description="Defines the methods to obtain the adiabatic forces.",&
         usage="METHOD_COUPLING : AOM",  type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SURF_HOP_CHOICE",&
         description="Decide whether to use or not the trivial hopping correction",&
         usage="SURF_HOP_CHOICE: BASIC EXACT", type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="THRESHOLD_TRIVIAL_CROSS",&
     !    description="parameter used as a threshold to activate trivial correction", &
      !    usage="1e-03", type_of_var=real_t)
    !CALL section_add_keyword(subsection,keyword)
    !CALL keyword_release(keyword)
  
    CALL keyword_create(keyword, name="DECOHERENCE_CORRECTIONS",&
         description="Defines the methods to decohere the wavefunction if required.",&
         usage="DECOHERENCE_CORRECTIONS: BASIC EXACT", type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    
    CALL keyword_create(keyword, name="DECO_TIME",&
         description="Defines the decoherence time used by the damping",&
         usage="DECO_TIME: BASIC EXACT", type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    
    CALL keyword_create(keyword, name="SPURIOUS_TRANSFER_CORR",&
         description="Choose if you want to activate the spurious tranfer corre.",&
         usage="SPURIOUS_TRANSFER_CORR: True", default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.) 
        ! type_of_var=char_t,&
        ! n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EDC_C",&
         description="C parameter in the Energy-based deco. correction formula",&
         usage="1.0D0", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EDC_E0",&
         description="E0 parameter in the Energy-based deco. correction formula",&
         usage="0.1D0", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="THRESHOLD_TAU_FORCES",&
         description="threshold used in the forces based calculation",&
         usage="0.1D0", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TEMPERATURE_FG_WIDTH",&
         description="Temperature to define FG width, currectly 298 K",&
         usage="0.0009444571", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INITIALIZATION",&
         description="Decide how to initiate the system for wavefunction and surface.",&
         usage="DIABATIC or ADIABATIC", & 
         type_of_var=char_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)


    CALL keyword_create(keyword, name="RK_PROPAGATION",&
         description="Decide which RK scheme is used: Diab  or Adiab or HS",&
         usage="DIABATIC or ADIABATIC or HS", & 
         type_of_var=char_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="REPRESENTATION",&
         description="Decide which representation is used: Diab  or Adiab",&
         usage="DIABATIC_BASIS or ADIABATIC_BASIS", & 
         type_of_var=char_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)


    CALL keyword_create(keyword, name="FIRST_DIABAT",&
         description="starting molecule number",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="FIRST_ADIABAT",&
         description="starting molecule number",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword,name="SEED",&
         description="Initial seed for the global (pseudo)random number "//&
         "generator to create a stream of normally Gaussian "//&
         "distributed random numbers.",&
         usage="SEED <INTEGER>",default_i_val=2000)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    NULLIFY(sub2section)
    CALL section_create(sub2section,name="AOM",&
         description="Defines the methods to obtain the adiabatic forces.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)

    CALL keyword_create(keyword, name="NUMBER_DIABATIC_STATES",&
         description="Number of diabatic states to consider",&
         usage="2", type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NUMBER_ATOMS_PER_SITE",&
         description="Number of atoms per site",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NUMBER_ORBITALS",&
         description="Number of orbitals per sites",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CUT_OFF_BETWEEN_SITES",&
         description="Distance cut-off for H_ab calculation",&
         usage="6.0D0", type_of_var=real_t, &
         unit_str="angstrom")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CUT_OFF_CONNECTIVITY",&
         description="Distance cut-off for H_ab calculation",&
         usage="3.50D0", type_of_var=real_t, &
         unit_str="bohr")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)


    CALL keyword_create(keyword, name="SCALING_FACTOR",&
         description="Hab = scaling*Sab",&
         usage="0.50D0*0.065190D0, HARTREE", type_of_var=real_t,&
         unit_str="hartree")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="DECOHERENCE_CRITERION",&
         description="Criterion to decohere",&
         usage="1.0E-05", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CBAR",&
         description="Sab constant factor",&
         usage="0.5082 WITHOUT UNITS", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ATOMIC_OVERLAP_CUTOFF",&
         description="An atomic overlap cutoff",&
         usage="1.0D-17", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NUCLEAR_TIMESTEP",&
         description="Timestep for atomic motion, must be coherent",&
         usage="timestep 0.5", type_of_var=real_t, &
         unit_str = "fs")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ELECTRONIC_PARTIAL_STEP",&
         description="Number of electronic steps by atomic step",&
         usage="5", type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="COLLAPSE",&
         description="Decide if we use the decoherence procedure",&
         usage="True", & 
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RANDOM_DERIVATIVE",&
         description="Decide if we use forward and backward derivative "//&
                      "in a random way (True) or only forward (False).",&
         usage="True", &
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ANALYTICS",&
         description="Decide if we use the analytical equation rather than "//&
                      "the numerical expression. Only for 2x2 systems.",&
         usage="True", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CENTER_OF_MASS",&
         description="Decide if do all electronic propagation and velocities "//&
                     "rescaling into the center-of-mass frame.",&
         usage="True", &
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="METHOD_RESCALING",&
         description="Defines the methods to rescale the velocities if required.",&
         usage="METHOD_RESCALING: BASIC EXACT", type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="METHOD_ADIABATIC_NACV",&
         description="Defines the methods to compute the adiabatic NACV for the rescaling.",&
         usage="METHOD_RESCALING: BASIC EXACT", type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="METHOD_REVERSAL",&
         description="Defines the methods to reverse the velocities after a frustrated hop if required.",&
         usage="METHOD_RESCALING: BASIC EXACT", type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NACV_INCREMENT",&
         description="Spatial Increment to calculate the NACV with Finite Difference.",&
         usage="1E-3*1.889725989D0 (bohr)", & 
         default_r_val=0.00188972598_dp,    &
         unit_str="bohr",                   &
         type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PPI_FILE_NAME",&
         description="Specifies the filename that contains the lookup table",&
         usage="PPI_FILE_NAME <FILENAME>",type_of_var=lchar_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PSIGMA_FILE_NAME",&
         description="Specifies the filename that contains the lookup table",&
         usage="PPI_FILE_NAME <FILENAME>",type_of_var=lchar_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)



    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="The HOMO coefficient in the format:"//&
         "<p><tt>ATOMIC_KIND  LIST NO HOMOS HOMOP</tt></p>",&
         repeats=.TRUE., usage="{{String} {Int} {Int} {Real} {Real}}",&
         type_of_var=lchar_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)

    CALL section_create(sub2section,name="PRINT",&
         description="Section of possible print options in ADIABATIC env.",&
         n_keywords=0, n_subsections=1, repeats=.FALSE.)

    NULLIFY(print_key)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="Controls the printing of information during the evaluation of "//&
         "the mixed environment. ",&
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key_section=print_key,&
          name="COEFFICIENTS",&
          description="Control the output of the wavefunction coefficients (diabatic basis).",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)
    !----------PC
     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="PVECS",&
          description="Control the output of the p-vectors for PCs visualisation script.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)
     !-------------------

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="POPULATIONS",&
          description="Control the output of the diabatic populations.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="ADIABAT_ENERGIES",&
          description="Control the output of the adiabatic energies.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="HAMILTONIAN",&
          description="Control the output of the Hamiltonian.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="ACTIVE_STATE_INFO",&
          description="Print the info relative to the active state (state, energy, wavefunction.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="SH_INFO",&
          description="Print the info relative to the surface hopping event:"//&
          "Probability, Frustrated Hop, Decoherence, etc.", &
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)


     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="NACE",&
          description="Control the output of the NACE.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "fs^-1")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="NACV",&
          description="Control the output of the NACV.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "bohr^-1")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="HOP_NACV",&
          description="Control the output of the NACV.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "bohr^-1")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="DENSITY_MATRIX",&
          description="Control the output of the density matrix (diabatic basis).",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="STATE",&
          description="Control the output of the present state.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="EXACT_FORCES",&
          description="Control the output of the exact forces in case of ANALYTICS.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)
!> End CHANGE_AC


    ! Multiple force_eval
    CALL section_create(subsection,name="GENERIC",&
         description="User driven coupling between two or more force_eval.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="MIXING_FUNCTION",&
         description="Specifies the mixing functional form in mathematical notation.",&
         usage="MIXING_FUNCTION (E1+E2-LOG(E1/E2))",  type_of_var=lchar_t,&
         n_var=1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="VARIABLES",&
         description="Defines the variables of the functional form. To allow an efficient"//&
         " mapping the order of the energy variables will be considered identical to the"//&
         " order of the force_eval in the force_eval_order list.",&
         usage="VARIABLES x",  type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PARAMETERS",&
         description="Defines the parameters of the functional form",&
         usage="PARAMETERS a b D",  type_of_var=char_t,&
         n_var=-1, repeats=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="VALUES",&
         description="Defines the values of parameter of the functional form",&
         usage="VALUES ",  type_of_var=real_t,&
         n_var=-1, repeats=.TRUE., unit_str="internal_cp2k")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="UNITS",&
         description="Optionally, allows to define valid CP2K unit strings for each parameter value. "//&
                     "It is assumed that the corresponding parameter value is specified in this unit.",&
         usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K",  type_of_var=char_t,&
         n_var=-1, repeats=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="DX",&
         description="Parameter used for computing the derivative with the Ridders' method.",&
         usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ERROR_LIMIT",&
         description="Checks that the error in computing the derivative is not larger than "//&
         "the value set; in case error is larger a warning message is printed.",&
         usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    ! Mapping of atoms
    NULLIFY(sub2section, sub3section)
    CALL section_create(subsection,name="MAPPING",&
         description="Defines the mapping of atoms for the different force_eval with the mixed force_eval."//&
         " The default is to have a mapping 1-1 between atom index (i.e. all force_eval share the same"//&
         " geometrical structure). The mapping is based on defining fragments and the mapping the "//&
         " fragments between the several force_eval and the mixed force_eval",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    ! Mixed force_eval
    CALL section_create(sub2section,name="FORCE_EVAL_MIXED",&
         description="Defines the fragments for the mixed force_eval (reference)",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL section_create(sub3section,name="FRAGMENT",&
         description="Fragment definition",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the fragment defined",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Starting and ending atomic index defining one fragment must be provided",&
         usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.TRUE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(sub2section,sub3section)
    CALL section_release(sub3section)
    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)

    ! All other force_eval
    CALL section_create(sub2section,name="FORCE_EVAL",&
         description="Defines the fragments and the mapping for each force_eval (an integer index (ID) "//&
         "needs to be provided as parameter)",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="DEFINE_FRAGMENTS",&
         description="Specify the fragments definition of the force_eval through the fragments of the"//&
         " force_eval_mixed. This avoids the pedantic definition of the fragments for the force_eval,"//&
         " assuming the order of the fragments for the specified force_eval is the same as the sequence "//&
         " of integers provided. Easier to USE should be preferred to the specification of the single fragments.",&
         usage="DEFINE_FRAGMENTS <INTEGER> .. <INTEGER>", type_of_var=integer_t, n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the force_eval for which fragments and mappings are provided",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL section_create(sub3section,name="FRAGMENT",&
         description="Fragment definition",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the fragment defined",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Starting and ending atomic index defining one fragment must be provided",&
         usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.FALSE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAP",&
         description="Provides the index of the fragment of the MIXED force_eval mapped on the"//&
         " locally defined fragment.",&
         usage="MAP <INTEGER>", type_of_var=integer_t, n_var=1, repeats=.FALSE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(sub2section,sub3section)
    CALL section_release(sub3section)
    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)

    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL create_print_mix_section(subsection)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)
  END SUBROUTINE create_mix_section

! *****************************************************************************
!> \brief Create the print section for mixed
!> \param section the section to create
!> \author teo
! *****************************************************************************
  SUBROUTINE create_print_mix_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_mix_section', &
      routineP = moduleN//':'//routineN

    TYPE(section_type), POINTER              :: print_key

    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,name="print",&
         description="Section of possible print options in MIXED env.",&
         n_keywords=0, n_subsections=1, repeats=.FALSE.)

    NULLIFY(print_key)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="Controls the printing of information during the evaluation of "//&
         "the mixed environment. ",&
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key,"DIPOLE",&
         description="Controls the printing of dipole information. "//&
         "Requires the DIPOLE calculation be active for all subforce_eval.", &
         print_level=medium_print_level,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)


  END SUBROUTINE create_print_mix_section

END MODULE input_cp2k_mixed
