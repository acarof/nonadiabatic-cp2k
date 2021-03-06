!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input section for NEGF based quantum transport calculations 
!>        (integration with the quantum transport code OMEN)
!>
!> \par History
!>       07.2013 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! *****************************************************************************
MODULE input_cp2k_transport
  USE bibliography,                    ONLY: Bruck2014
  USE cp_dbcsr_interface,              ONLY: csr_dbcsr_blkrow_dist,&
                                             csr_eqrow_ceil_dist,&
                                             csr_eqrow_floor_dist
  USE input_constants,                 ONLY: &
       injmethod_beyn, injmethod_evp, injsign_negative, injsign_positive, &
       linsolver_mumps, linsolver_splitsolve, linsolver_superlu, &
       transport_diagscalapack, transport_localscf, transport_negf, &
       transport_transmission, transport_writeout
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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_transport'

  PUBLIC :: create_transport_section

CONTAINS

! *****************************************************************************
!> \brief creates the TRABSPORT section
!> \param[inout] section the section to be created
! *****************************************************************************
  SUBROUTINE create_transport_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_transport_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,"TRANSPORT",&
         description="Specifies the parameters for transport, sets parameters for the OMEN code, "//&
                     "see also http://www.nano-tcad.ethz.ch/ ",&
         citations=(/Bruck2014/),&
         n_keywords=19, n_subsections=1, repeats=.FALSE.)

    NULLIFY (keyword)

    CALL keyword_create(keyword, name="TRANSPORT_METHOD",&
         description="Preferred method for transport calculations.",&
         usage="TRANSPORT_METHOD <method>",&
         default_i_val=transport_negf,&
         enum_c_vals=s2a( "DIAG", "TRANSPORT", "LOCAL_SCF", "TRANSMISSION", "WRITEOUT"),&
         enum_desc=s2a("standard diagonalization using ScaLapack routines",&
         "self-consistent CP2K and OMEN transport calculations",&
         "CP2K valence Hamiltonian + OMEN self-consistent calculations on conduction electrons",&
         "self-consistent transmission calculations without applied bias voltage",&
         "write out the Hamiltonian and overlap matrices"),&
         enum_i_vals=(/transport_diagscalapack, transport_negf, transport_localscf, & 
                       transport_transmission, transport_writeout/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="N_ABSCISSAE",&
         description="The number of abscissae per integration interval on the real axis.",&
         usage="N_ABSCISSAE <integer>", default_i_val=64)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="N_KPOINTS",&
         description="The number of k points for determination of the singularities.",&
         usage="N_KPOINTS <integer>", default_i_val=64)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NUM_INTERVAL",&
         description="Max number of energy points per small interval.",&
         usage="NUM_INTERVAL <integer>", default_i_val=10)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TASKS_PER_POINT",&
         description="Minimum number of tasks per energy point.",&
         usage="TASKS_PER_POINT <integer>", default_i_val=1)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="COLZERO_THRESHOLD",&
         description="The smallest number that is not zero in the full diagonalization part.",&
         usage="COLZERO_THRESHOLD <real>", default_r_val=1.0E-12_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EPS_LIMIT",&
         description="The smallest eigenvalue that is kept.",&
         usage="EPS_LIMIT <real>", default_r_val=1.0E-4_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EPS_DECAY",&
         description="The smallest imaginary part that a decaying eigenvalue may have not to be "//&
         "considered as propagating.",&
         usage="EPS_DECAY <real>", default_r_val=1.0E-4_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EPS_SINGULARITY_CURVATURES",&
         description="Filter for degenerate bands in the bandstructure.",&
         usage="EPS_SINGULARITY_CURVATURES <real>", default_r_val=1.0E-12_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EPS_MU",&
         description="Accuracy to which the Fermi level should be determined.",&
         usage="EPS_MU <real>", default_r_val=1.0E-6_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EPS_EIGVAL_DEGEN",&
         description="Filter for degenerate bands in the injection vector.",&
         usage="EPS_EIGVAL_DEGEN <real>", default_r_val=1.0E-6_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ENERGY_INTERVAL",&
         description="Average distance for big intervals in energy vector.",&
         usage="ENERGY_INTERVAL <real>", default_r_val=1.0E-3_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MIN_INTERVAL",&
         description="Smallest enery distance in energy vector.",&
         usage="MIN_INTERVAL <real>", default_r_val=1.0E-4_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TEMPERATURE",&
         description="Temperature.",&
         usage="TEMPERATURE [K] 300.0", default_r_val=300.0_dp, unit_str="K")
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ROW_DISTRIBUTION",&
         description="How to distribute matrix rows over tasks.",&
         enum_c_vals=s2a("DBCSR","CEILING","FLOOR"),&
         enum_desc=s2a("Each task holds floor(N_BLOCK/TASKS_PER_POINT) DBCSR block rows for a total of N_BLOCK block rows",&
                       "Each task holds ceiling(N/TASKS_PER_POINT) rows for a total of N matrix rows",&
                       "Each task holds floor(N/TASKS_PER_POINT) rows for a total of N matrix rows"),&
         enum_i_vals=(/csr_dbcsr_blkrow_dist,csr_eqrow_ceil_dist,csr_eqrow_floor_dist/),&
         default_i_val=csr_eqrow_floor_dist)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CSR_SCREENING",&
         description="Whether distance screening should be applied to improve sparsity of CSR matrices.",&
         default_l_val=.TRUE., lone_keyword_l_val = .TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="LINEAR_SOLVER",&
         description="Solver to solve the linear system of equations.",&
         usage="LINEAR_SOLVER <solver>", &
         default_i_val=linsolver_superlu,&
         enum_c_vals=s2a( "SplitSolve", "SuperLU", "MUMPS"),&
         enum_i_vals=(/linsolver_splitsolve, linsolver_superlu, linsolver_mumps/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INJECTION_METHOD",&
         description="Method to solve the eigenvalue problem for the open boundary conditions.",&
         usage="INJECTION_METHOD <method>", &
         default_i_val=injmethod_beyn,&
         enum_c_vals=s2a( "EVP", "BEYN"),&
         enum_desc=s2a("Full eigenvalue solver.",&
                       "Beyn eigenvalue solver."),&
         enum_i_vals=(/injmethod_evp, injmethod_beyn/))
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CUTOUT", &
         description="Degrees of freedom at the beginning and the end of the structure where the "//&
         "density should not be changed.", &
         usage="CUTOUT <integer> <integer>", &
         n_var=2, default_i_vals=(/0, 0/))
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    NULLIFY(subsection)

    CALL create_contact_section(subsection)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

  END SUBROUTINE create_transport_section

! *****************************************************************************
!> \brief Creates the section for creating contacts.
!> \param[inout] section the section to be created
! *****************************************************************************
  SUBROUTINE create_contact_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_contact_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,name="CONTACT", &
         description="Parameters for defining device contacts.",&
         n_keywords=5, n_subsections=0, repeats=.TRUE.)

    NULLIFY(keyword)

    CALL keyword_create(keyword, name="BANDWIDTH",&
         description="The number of neighboring unit cells that one unit cell interacts with.",&
         usage="BANDWIDTH <integer>", default_i_val=2)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="START",&
         description="Index of the atom where the contact unit cell starts.",&
         usage="START <integer>", default_i_val=-1)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="N_ATOMS",&
         description="Number of atoms in the contact unit cell.",&
         usage="N_ATOMS <integer>", default_i_val=0)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INJECTION_SIGN",&
         description="Contact unit cell interacts with unit cells to the right (positive) or "//&
         "to the left (negative).",&
         usage="INJECTION_SIGN <integer>", &
         default_i_val=injsign_positive,&
         enum_c_vals=s2a( "POSITIVE", "NEGATIVE"),&
         enum_desc=s2a("When the contact unit cell is at the upper left corner of the Hamiltonian.",&
                       "When the contact unit cell is at the lower right corner of the Hamiltonian."),&
         enum_i_vals=(/injsign_positive, injsign_negative/))
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INJECTING_CONTACT",&
         description="whether or not the contact can inject electrons.",&
         default_l_val=.TRUE., lone_keyword_l_val = .TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

  END SUBROUTINE create_contact_section

END MODULE input_cp2k_transport
