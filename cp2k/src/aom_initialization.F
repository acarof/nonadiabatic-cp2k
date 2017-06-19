!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   !
!   Copyright (C) 2000 - 2015  CP2K developers group
!   !
!-----------------------------------------------------------------------------!

!!!!!!!!!! CHANGE_AC !!!!!!!!!!!!!

! *****************************************************************************
!> \brief Interface for the adiabatic force calculations
!> \par History
!>      carof, 2016.01.15: creation, put only ffsh_aom subroutines
! *****************************************************************************

! *****************************************************************************
!>      jacob: comment
! construct_ham_X_SH_full_degen
!
! This program is designed to propagate a single excess electronic charge
! along a molecular dynamics trajectory, either during the trajectory
! (to do surface hopping simulations) or afterwards (classical path approximation).
!
! This program takes molecular coordinates as inputs.
! From these coordinate inputs, it calculates the electronic Hamiltonian
! in a molecular orbital basis, using the overlapFINAL module to calculate
! the off-diagonals and either accepting the diagonals as inputs or calculating
! them in some other way.
!
! With the electronic Hamiltonian, the program then propagates the electron using the RK4 method.
! At the final electronic state, the probabilities for hopping between adiabatic surfaces are calculated.
! The program then decides what adiabatic surface the system is on and calculates
! the appropriate atomic forces.
!
! In the classical path approximation, the above steps are repeated until the
! end of the input molecular trajectory.  For surface-hopping, the data here is fed back into
! the molecular dynamics and the trajectory moves forward another time step.
! *****************************************************************************

MODULE aom_initialization 

     USE kinds,                        ONLY: default_string_length,&
                                             dp
     USE aom_overlapfinal,             ONLY: connect_list2
     USE aom_utils,                    ONLY: read_listcoeff
     USE aom_types,                    ONLY: aom_env_type
     USE input_section_types,          ONLY: section_vals_type,&
                                             section_vals_list_get,&
                                             section_vals_val_get
     USE cp_linked_list_val,           ONLY: cp_sll_val_type
     USE parallel_rng_types,           ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream
#include "./base/base_uses.f90"

     IMPLICIT NONE

     PRIVATE
 
     CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='aom_initialization'   

     PUBLIC :: create_aom_env

CONTAINS

  SUBROUTINE create_aom_env(nadiab, aom_atoms, norbitals, natom, &
                                  positions_tp,     &
                                  coupling_section, &
                                  adiab_section,    &
                                  aom_env)


     INTEGER, INTENT(IN)                            :: nadiab,  aom_atoms, norbitals, natom
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     TYPE(section_vals_type), POINTER               :: coupling_section
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(aom_env_type), INTENT(INOUT)              :: aom_env

     CHARACTER(len=*), PARAMETER :: routineN ='create_aom_env', &
     routineP = moduleN//':'//routineN

     integer :: i, j, mat, c, coord
     REAL(KIND=dp) :: cutoff_connect
     REAL(KIND=dp), dimension(:,:), allocatable :: crds
     REAL(KIND=dp), dimension(:,:), allocatable :: Atp
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pcoefftot, scoefftot
     LOGICAL, DIMENSION(:), ALLOCATABLE :: atomlist_mask
     integer, dimension(:), ALLOCATABLE           :: ext_elements
     TYPE(cp_sll_val_type), POINTER               :: listcoeff

     LOGICAL    :: do_random
     REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed
     INTEGER    :: input_seed

     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_CONNECTIVITY",&
             r_val=cutoff_connect)
     CALL section_vals_val_get(coupling_section,"RANDOM_DERIVATIVE", l_val=do_random)

     mat = norbitals*nadiab !> Size of the vectorial space where the wavefunctions are projected
     c = nadiab*aom_atoms !> Numbers of active atoms

     ALLOCATE(aom_env%population(mat))
     ALLOCATE(aom_env%ctin(mat))
     ALLOCATE(aom_env%hamiltonian(mat,mat))
     ALLOCATE(aom_env%eigenvector_matrix(mat,mat))
     ALLOCATE(aom_env%connlist1(4,aom_atoms))
     ALLOCATE(aom_env%neighbourlist(nadiab, nadiab,0:aom_atoms,0:aom_atoms))
     ALLOCATE(aom_env%Stt(mat,mat))
     ALLOCATE(aom_env%Stptp(mat,mat))
     ALLOCATE(aom_env%Stmt(mat,mat))
     ALLOCATE(aom_env%Sttm(mat,mat))
     ALLOCATE(aom_env%Stmtm(mat,mat))
     ALLOCATE(aom_env%NN(mat, mat))
     ALLOCATE(aom_env%rm_t(3, natom))
     ALLOCATE(aom_env%rm_tm(3, natom))
     ALLOCATE(aom_env%energies_t(nadiab))
     ALLOCATE(aom_env%adiab_energies_t(mat))
     ALLOCATE(aom_env%bneighbour(3,0:aom_atoms))
     ALLOCATE(aom_env%nacv(3,mat,mat,mat*aom_atoms))
     ALLOCATE(aom_env%diabatic_forces(3,natom, mat))
     ALLOCATE(aom_env%forces_adiab(3,natom))
     ALLOCATE(aom_env%forces_adiab_new(3,natom))
     allocate(aom_env%adiab_nace_tp(mat,mat))
     ALLOCATE(aom_env%pcoeffs(aom_atoms,norbitals))
     ALLOCATE(aom_env%scoeffs(aom_atoms))
     aom_env%futur_decohere = .FALSE.
     allocate(aom_env%elements(c))

     IF (do_random) THEN
        CALL section_vals_val_get(adiab_section,"SEED",i_val=input_seed)
        initial_seed(:,:) = REAL(input_seed,KIND=dp)
        CALL create_rng_stream(rng_stream=aom_env%rng_stream,&
                             name="aom fssh",&
                             distribution_type=UNIFORM, &
                             seed=initial_seed)
     ENDIF

     !> READ THE AOM COEFFICIENT
     CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff) !Read them in the input files
     CALL read_listcoeff(listcoeff, ext_elements, aom_env%ext_pres, &
                         scoefftot, pcoefftot, natom) !Transform them into four arrays
     aom_env%elements = PACK(ext_elements, aom_env%ext_pres) !Only select active atoms
     allocate(atomlist_mask(c))
     atomlist_mask =.NOT.((scoefftot.eq.0).AND.(pcoefftot.eq.0)) !Only select atoms with orbitals
     allocate(aom_env%atomlist(COUNT(atomlist_mask(1:aom_atoms))))
     aom_env%atomlist = PACK( (/ (i, i=1,aom_atoms) /), atomlist_mask(1:aom_atoms) ) !List atoms with orbitals
     allocate(aom_env%connlist2(4,size(aom_env%atomlist)))
     do i=1,aom_atoms
          do j=1,norbitals
                aom_env%pcoeffs(i,j) = pcoefftot( i + (j-1)*norbitals)
          end do
          aom_env%scoeffs(i) = scoefftot(i)
     end do

     !> CONSTRUCT THE CONNLIST
     !we fill out the connlist arrays before calculating the nacv
     allocate(crds(3, c))
     do coord=1,3
        crds(coord,:) = PACK(positions_tp(coord,:), aom_env%ext_pres)
     enddo
     ALLOCATE(Atp(4,aom_atoms))
     Atp(2:4,1:aom_atoms) = crds(1:3,1:aom_atoms)
     call connect_list2(Atp(2:4,:), aom_env%atomlist, aom_env%connlist1, aom_env%connlist2, cutoff_connect)

  END SUBROUTINE create_aom_env

END MODULE aom_initialization
