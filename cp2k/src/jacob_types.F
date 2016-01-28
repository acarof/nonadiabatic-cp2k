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

   TYPE aom_basket_type
!> To store the coordinates in step before and the one before
     REAL(KIND = dp), DIMENSION(:,:), POINTER         :: rm1,&       ! position step min 1
                                                         rm2         ! position step min 2
! To store energies in the Jacob's Method
     REAL(KIND = dp), DIMENSION(:), POINTER :: energies   
     INTEGER                                :: state
     LOGICAL                                :: inrun
     REAL(KIND=dp), dimension(:,:), POINTER :: Stt, Stmt, Sttm, Stmtm
     integer, dimension(:,:), POINTER       :: connlist1, connlist2, bneighbour
     COMPLEX(KIND=dp), DIMENSION(:),POINTER              :: ctin
     REAL(KIND=dp), DIMENSION(:), POINTER                :: population
     INTEGER, DIMENSION(:,:,:), POINTER         :: neighbourlist
     REAL(KIND=dp), DIMENSION(:,:), POINTER              :: hamiltonian
   END TYPE aom_basket_type

   PUBLIC :: aom_basket_type

END MODULE jacob_types
