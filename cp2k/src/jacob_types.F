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
! (to do surface hopping simulations) or afterwards (classical path
! approximation).
!
! This program takes molecular coordinates as inputs.
! From these coordinate inputs, it calculates the electronic Hamiltonian
! in a molecular orbital basis, using the overlapFINAL module to calculate
! the off-diagonals and either accepting the diagonals as inputs or calculating
! them in some other way.
!
! With the electronic Hamiltonian, the program then propagates the electron
! using the RK4 method.
! At the final electronic state, the probabilities for hopping between adiabatic
! surfaces are calculated.
! The program then decides what adiabatic surface the system is on and
! calculates
! the appropriate atomic forces.
!
! In the classical path approximation, the above steps are repeated until the
! end of the input molecular trajectory.  For surface-hopping, the data here is
! fed back into
! the molecular dynamics and the trajectory moves forward another time step.
! *****************************************************************************

MODULE jacob_types
  USE kinds,                           ONLY: dp

IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'jacob_types'

   TYPE aom_env_type
!> To store the coordinates in step before and the one before
     REAL(KIND = dp), DIMENSION(:,:), POINTER         :: rm1,&             ! position step min 1
                                                         rm2,&             ! position step min 2
                                                         velocities_plus,& ! velocities step
                                                         velocities,&      ! velocities step min 1
                                                         forces_adiab      ! forces at step 0
!> To store energies in the Jacob's Method
     REAL(KIND = dp), DIMENSION(:), POINTER :: energies   
!> Old and new state number
     INTEGER                                :: new_state
     INTEGER                                :: old_state
!     LOGICAL                                :: inrun
!> To store the previous overlap matrix
     REAL(KIND=dp), dimension(:,:), POINTER :: Stt, Stmt, Sttm, Stmtm
!> To store some connectivity information
     integer, dimension(:,:), POINTER       :: connlist1, connlist2, bneighbour
!> To store the previous coefficients
     COMPLEX(KIND=dp), DIMENSION(:),POINTER              :: ctin
!> To print the populations
     REAL(KIND=dp), DIMENSION(:), POINTER                :: population
!> To store some connectivity information
     INTEGER, DIMENSION(:,:,:,:), POINTER                  :: neighbourlist
!> To print the Hamiltonian
     REAL(KIND=dp), DIMENSION(:,:), POINTER              :: hamiltonian
!> To print the NACEs
     REAL(KIND=dp), DIMENSION(:,:), POINTER              :: nace
!> To print the NACVs
     REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER          :: nacv
!> Store information of a previous hop- required for the decoherence
!procedure
     LOGICAL                                             :: futur_decohere
     REAL(KIND = dp)                                     :: factor
!> If hop
     LOGICAL                                     :: is_hop
   END TYPE aom_env_type


   TYPE sh_print_env_type
!> The probability array for the hop
     REAL(KIND = dp), DIMENSION(:), POINTER      :: proba
!> Random number for test hop
     REAL(KIND = dp)                             :: hop
!> If renormalization occurs
     LOGICAL                                     :: renorm
!> If renormalization occurs: the sum of probabilities
     REAL(KIND = dp)                             :: sumprob
     REAL(KIND = dp)                             :: factor
!> Old and new state number
     INTEGER                                     :: old_state
     INTEGER                                     :: new_state
!> Delta energy between old_state and putative new_state
     REAL(KIND = dp)                             :: delta_energy
!> If hop occurs
     LOGICAL                                     :: is_hop
!> If decoherence occurs
     LOGICAL                                     :: decohere
!> If it's the initial step
     LOGICAL                                     :: is_first
!> Random number for initial choice
     REAL(KIND = dp)                             :: initrandom
!> Initial population in the adiabatic states
     REAL(KIND = dp), DIMENSION(:), POINTER      :: pop_init
!> First adiabat after the random choice
     INTEGER                                     :: first_adiabat
   END TYPE sh_print_env_type

   PUBLIC :: aom_env_type, &
             sh_print_env_type


END MODULE jacob_types
