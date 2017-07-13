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

MODULE aom_types
  USE kinds,                           ONLY: default_string_length, &
                                             dp
  USE input_section_types,             ONLY: section_vals_type

  USE parallel_rng_types,              ONLY: rng_stream_type


IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'aom_types'



   TYPE sh_print_env_type
!> General variable
     INTEGER                                     :: nadiab, norbitals
     REAL(KIND = dp)                             :: times, deltatMD
     TYPE(section_vals_type), POINTER            :: adiab_section
     REAL(KIND = dp), DIMENSION(:), POINTER      :: adiab_energies_t   
!> Variable for the first timestep
     REAL(KIND = dp), DIMENSION(:), POINTER      :: pop_init
     INTEGER                                     :: first_adiabat
     REAL(KIND = dp)                             :: initrandom
!> Variable for Tully's test
     LOGICAL                                     :: attempt_hop
     REAL(KIND = dp)                             :: sumprob
     REAL(KIND = dp)                             :: enermin, exprob, sumflux 
     LOGICAL                                     :: successful_hop
     LOGICAL                                     :: additional_prints
     INTEGER                                     :: old_state
     INTEGER                                     :: new_state
     LOGICAL                                     :: renorm
     REAL(KIND = dp), DIMENSION(:), POINTER      :: proba
     REAL(KIND = dp), DIMENSION(:), POINTER      :: tau_deco
     REAL(KIND = dp)                             :: hop
!Variable for the rescaling
     REAL(KIND = dp)                             :: delta_energy
     REAL(KIND = dp)                             :: kinetic_energy
     REAL(KIND = dp)                             :: A, B
     REAL(KIND = dp)                             :: rescale_factor
     REAL(KIND = dp), DIMENSION(3)               :: momentum_shift
!> Variable for the reversing
     LOGICAL                                     :: reverse
     CHARACTER(LEN=default_string_length)        :: reversal_keyword
     REAL(KIND = dp)                             :: product_pd, product_dfn, product_dfo
!> Variable for the decoherence
     LOGICAL                                     :: decohere
     LOGICAL                                     :: collapse
!> To store energies in the Jacob's Method
!> The probability array for the hop
!> Random number for test hop
!> If renormalization occurs
!> If renormalization occurs: the sum of probabilities
!> Old and new state number
!> Delta energy between old_state and putative new_state
!> Random number for initial choice
!> Initial population in the adiabatic states
!> First adiabat after the random choice
   END TYPE sh_print_env_type


   TYPE aom_env_type
!> To store the coordinates in step before and the one before
     REAL(KIND = dp), DIMENSION(:,:), POINTER         :: rm_t,&             ! position step min 1
                                                         rm_tm,&             ! position step min 2
                                                         velocities_tp,& ! velocities step
                                                         velocities_t,&      ! velocities step min 1
                                                         forces_adiab, &      ! forces at step 0
                                                         forces_adiab_new
!> To store energies in the Jacob's Method
     REAL(KIND = dp), DIMENSION(:), POINTER :: energies_t   
     REAL(KIND = dp), DIMENSION(:), POINTER :: adiab_energies_t 
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: F
     REAL(KIND=dp), dimension(:,:), POINTER :: At, Atm, Atp, Bt, Btm, & 
                                                         Btp, siteenergies
     REAL(KIND = dp)                        :: kinetic_energy
     REAL(KIND = dp)                        :: delta_energy, scaling, cbar
!> Old and new state number
     INTEGER                                :: new_state
     INTEGER                                :: active_state
     INTEGER                                :: old_state
!     LOGICAL                                :: inrun
     REAL(KIND=dp), dimension(:,:), POINTER :: coms
     REAL(KIND=dp), dimension(:,:,:), POINTER :: crds
     REAL(KIND=dp), dimension(:,:), POINTER :: dist
!> To store the previous overlap matrix
     REAL(KIND=dp), dimension(:,:), POINTER :: Stt, Stmt, Sttm, Stmtm, &
                                             Stptp, Sttp, Stpt
!     REAL(KIND=dp), dimension(:,:), POINTER :: Stt, Sttp, Stptp, Stpt, Stmt,
!Sttm, Stmtm
!> To store some connectivity information
     integer, dimension(:,:), POINTER       :: connlist1, connlist2, bneighbour
!> To store the previous coefficients
     COMPLEX(KIND=dp), DIMENSION(:),POINTER              :: ctin
!> To print the populations
     REAL(KIND=dp), DIMENSION(:), POINTER                :: population
!> To store some connectivity information
     INTEGER, DIMENSION(:,:,:,:), POINTER                  :: neighbourlist
     REAL(KIND=dp), dimension(:,:,:), allocatable          :: pcoeffsrun
!> To print the Hamiltonian
     COMPLEX(KIND=dp), DIMENSION(:,:), POINTER              :: H_t, H_tp
     REAL(KIND=dp), DIMENSION(:,:), POINTER              :: eigenvector_matrix
     INTEGER, DIMENSION(:,:), POINTER                    :: NN
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER            :: diabatic_forces
     LOGICAL, DIMENSION(:), POINTER                      :: ext_pres
     INTEGER, DIMENSION(:), POINTER                      :: elements, atomlist
!> To print the NACEs
     REAL(KIND=dp), DIMENSION(:,:), POINTER              :: nace
!> To print the NACVs
     REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER          :: nacv
!> Store information of a previous hop- required for the decoherence
!procedure
     LOGICAL                                             :: futur_decohere
     REAL(KIND = dp)                                     :: factor
!> If hop
     LOGICAL                                     :: attempt_hop
     LOGICAL                                     :: successful_hop
!> Adiabatic nace at time t and t-dt
     !REAL(KIND=dp), dimension(:,:), POINTER      :: adiab_nace_tp, adiab_nace_t
     REAL(KIND=dp), dimension(:,:), POINTER      :: adiab_nace_tp, adiab_nace_t
     TYPE(sh_print_env_type)                     :: sh_print_env
     TYPE(rng_stream_type), POINTER              :: rng_stream
     REAL(KIND=dp), DIMENSION(:,:), POINTER      :: pcoeffs
     REAL(KIND=dp), DIMENSION(:), POINTER         :: scoeffs
     REAL(KIND = dp)                              :: energy_new
     INTEGER, dimension(:), POINTER          :: class_forbidden_hop_t

   END TYPE aom_env_type

   PUBLIC :: aom_env_type, &
             sh_print_env_type


END MODULE aom_types
