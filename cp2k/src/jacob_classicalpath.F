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

MODULE jacob_classicalpath 
 
     USE  kinds,                           ONLY: default_string_length,&
                                               dp


     USE jacob_overlapfinal,                 ONLY:  connect_list2, &
                                             calc_sab, &
                                             sab_screen_new, &
                                             calc_dRSab

     USE jacob_utils,                        ONLY: makeX

     USE jacob_types,                     ONLY: aom_basket_type
     USE input_section_types,             ONLY: section_vals_get,&
                                                section_vals_get_subs_vals,&
                                                section_vals_type,&
                                                section_vals_list_get,&
                                                section_vals_val_get
     USE jacob_print,                     ONLY: adiab_print_coeff, &
                                                adiab_print_pop,   &
                                                adiab_print_matrix_density


     !USE jacob_types,                     ONLY: aom_basket_type

     !   USE jacob_utils
     !   USE jacob_types,                     ONLY: aom_basket_type

#include "./base/base_uses.f90"


     implicit none

        PRIVATE
 
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='jacob_construct_ham_X_SH_full_degen_classicalpath'   

        PUBLIC :: construct_ham_X_SH_full_degen_classicalpath

CONTAINS


        SUBROUTINE construct_ham_X_SH_full_degen_classicalpath(nadiab, &
                   ext_atoms, norbitals, ext_cutoff, ext_scaling, &
                   ext_init, ext_cbar, ext_minao, ext_timestep, ext_el_steps, &
                   positions_plus, energies_plus, aom_basket, ext_elements, &
                   ext_pres, ext_scoeff, ext_pcoeff, itimes, adiab_section)

        integer          :: nadiab,  norbitals, ext_init
        integer          :: ext_atoms, ext_el_steps
        REAL(KIND=dp) :: ext_cutoff, ext_scaling, ext_cbar, ext_minao, ext_timestep
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_plus
        REAL(kind=dp), DIMENSION(:),  POINTER    :: energies_plus
        TYPE(aom_basket_type)                    :: aom_basket
        integer, dimension(:), ALLOCATABLE, INTENT(IN)        :: ext_elements
        LOGICAL, dimension(:), ALLOCATABLE,  INTENT(IN)       :: ext_pres
        REAL(KIND=dp), dimension(:), ALLOCATABLE, INTENT(IN)  :: ext_scoeff,ext_pcoeff
        INTEGER, POINTER                                       :: itimes
        TYPE(section_vals_type), POINTER         :: adiab_section


      !  REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions, positions_plus
      !  REAL(kind=dp), DIMENSION(:),  POINTER    :: energies,  energies_plus

        CHARACTER(len=*), PARAMETER :: routineN ='construct_ham_X_SH_full_degen_classicalpath', &
        routineP = moduleN//':'//routineN

     !logical :: inrun, debug, decohere, collapse
!     logical :: inrun
        integer :: coord
     !integer :: mol, i, j, k, l, m, d, stat, atoms, times, c, t, t2, init, inc, filelength, state, mat, orb, bondsno
     integer :: mol, i, j, k, l, m, d, atoms, c, t, t2, init, inc, old_state, mat, orb
!     integer, dimension(:), allocatable :: atomlist, pfiles, elements
     integer, dimension(:), allocatable :: atomlist, elements
     integer, dimension(:,:), allocatable :: NN, connlist1, connlist2, aneighbour, bneighbour
     integer, dimension(:,:,:), allocatable :: neighbourlist
     !REAL(KIND=dp) :: cutoff, scaling, deltatMD, deltatE, lambdai, distance
     REAL(KIND=dp) :: cutoff, scaling, deltatMD, deltatE
     REAL(KIND=dp) :: S
     !COMPLEX(KIND=dp), dimension(2,2) :: CT, DCT
     !character(len=100) :: mdcrd, mdcrdm, mdcrdp, bondlistfile, filler, scoefficients, mdcrdtype, atomlistfile, elementfile
     !character(len=100) ::  scoefficients
!     character(len=100), dimension(:), allocatable :: pcoefficients
     !character(len=1) :: zgemm_type, input, progtype
     character(len=1) :: zgemm_type
     !REAL(KIND=dp) :: ctmag, bohrang, mev, kcal, offset, cbar, nsave, minao
     REAL(KIND=dp) :: ctmag, bohrang, mev, kcal, cbar, minao
     !REAL(KIND=dp), dimension(:), allocatable :: scoeffs, E1, E2
     REAL(KIND=dp), dimension(:), allocatable :: scoeffs
     REAL(KIND=dp), dimension(:,:,:,:), allocatable :: nacv
     !REAL(KIND=dp), dimension(:,:,:), allocatable :: crds, pcoeffsrun, F
     REAL(KIND=dp), dimension(:,:,:), allocatable :: crds, pcoeffsrun
     !REAL(KIND=dp), dimension(:,:), allocatable :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm, Imat, pcoeffs, tclforces
     REAL(KIND=dp), dimension(:,:), allocatable :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm, Imat, pcoeffs
     !REAL(KIND=dp), dimension(:,:), allocatable :: Stt2, Sttp2, Stptp2, Stpt2, Stmt2, Sttm2, Stmtm2
     !REAL(KIND=dp), dimension(:,:), allocatable :: Sttp2, Stmtm2
     REAL(KIND=dp), dimension(:,:), allocatable :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     !REAL(KIND=dp), dimension(:,:), allocatable :: C1, C2, dist
     REAL(KIND=dp), dimension(:,:), allocatable :: dist
     !REAL(KIND=dp) :: progtime
     !COMPLEX(KIND=dp) :: alpha, beta, ctcomp, DH11, DH22, DH12
     COMPLEX(KIND=dp) :: alpha, beta, ctcomp
     COMPLEX(KIND=dp), dimension(:), allocatable :: ctin, ctout, ctotal
     !COMPLEX(KIND=dp), dimension(:,:), allocatable :: H1, H2, H12, Hdiff, Hsave, X1, X2, X12, deltaA
     COMPLEX(KIND=dp), dimension(:,:), allocatable :: H1, H2, H12, Hdiff, Hsave, X1, X2, X12, deltaA, H2save
     COMPLEX(KIND=dp), dimension(:,:), allocatable :: Ktot, K2, K3, K4
     !REAL(KIND=dp), dimension(:,:), allocatable :: At, Atm, Atp, Atpp, Bt, Btm, Btp, coms, siteenergies, Atal
     REAL(KIND=dp), dimension(:,:), allocatable :: At, Atm, Atp, Bt, Btm, Btp, coms, siteenergies
!     integer, dimension(:,:), allocatable :: bonds

     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pcoefftot, scoefftot
     INTEGER, DIMENSION(:), ALLOCATABLE :: elementtot
     LOGICAL, DIMENSION(:), ALLOCATABLE :: atomlist_mask
!     INTEGER                            :: my_itimes

     REAL(KIND=dp)                      :: dR,times

!print*, energies_plus(:), "kkkk"

     !We begin by reading in the input parameters
     !There are currently two modes this program can run in
     !In Classical Path Approximation mode, the program reads in a pre-computed trajectory
     !In Surface Hopping mode, the program runs alongside non-adiabatic MD
     !This variable is 'S' for Surface Hopping and 'C' for Classical Path
!     progtype = 'S'
     !name of coordinate file(s)
     !technically only mdcrd matters for progtype == C
!     mdcrdm = './coordinates_minus.txt'
!     mdcrd = './coordinates.txt'
!     mdcrdp = './coordinates_plus.txt'
     !'type' of file so that program can interpret it
     !currently supported:
     ! - "AMBER" for AMBER .mdcrd file
     ! - "NAMD" for NAMD .dcd file
     ! - "XYZ" for .xyz file
     ! - "SH" is format printed out by NAMD surfhop code and should be default for progtype == S
!     if (progtype == 'S') then
!          mdcrdtype = 'SH'
!     else
!          mdcrdtype = 'AMBER'
!     end if
     !number of molecules
!     mol = 10
     mol = nadiab ! CHANGE_AC
     !atoms per molecule
!     atoms = 6
        atoms = ext_atoms ! CHANGE_AC
     !for Classical Path Approximation, this is the number of timesteps in pre-computed trajectory
     !otherwise, 2
!     if (progtype == 'S') then
!          times = 2
!     elseif (progtype == 'C') then
          !times = 1000
!          times = 2 ! CHANGE_AC
!     end if
     !bonds per molecule
!     bondsno = 5
     !name of file containing bondlist
!     bondlistfile = './bondlistethylene'
     !name of file containing list of atomic number of each atom in molecule
!     elementfile = './atomnolist'
     !orbitals per molecule
     !orb = 1
     orb = norbitals ! CHANGE_AC
     !name of s-coefficient file
!     scoefficients = './coeffjhomos'
     !intermolecular cutoff
     !cutoff = 6.0D0
     cutoff = ext_cutoff ! CHANGE_AC
     !scaling factor, Hab = scaling*Sab
     !scaling = 0.5D0*0.06519D0
     scaling = ext_scaling ! CHANGE_AC
     !starting molecule number
     !init = 1
     init = ext_init ! CHANGE_AC
     !Sab constant factor
     !cbar = 0.5082
     cbar = ext_cbar ! CHANGE_AC
     !minao is an atomic overlap cutoff
     !a value of 1E-17 means all atoms are included
     !minao = 1.0D-17
     minao = ext_minao ! CHANGE_AC

!> dR is the finite difference to calculate the NACV
     dR = 0.0001D0*1.889725989D0


     !With surface-hopping, the program needs information about where in the simulation it is
     !i.e. if this is the first timestep, or if certain once-in-a-while actions need to be done
     !The file "inputget" contains this information
!     if (progtype == 'S') then
!          open(unit=39,file='./inputget',action='read')
!          read(39,*), input
!          close(39)
          !inrun is a logical which returns FALSE if this is the first timestep of the simulation
          !and TRUE otherwise
!          if (input == 'f') then
!               inrun = .FALSE.
!          else
               !inrun = .TRUE.

         !my_itimes = itimes 
!         print*, "here we are 001"
!         print*, "itimes", itimes 
!         print*, "here we are 002"
!         if (ASSOCIATED(itimes)) then
!            print*, "Associated!"
!            print*, 'itimes', itimes
!         else
!            print*, "Not associated"
!         endif
!         if (itimes.eq.1) then 
!           inrun = .FALSE.
!         else
!            inrun = .TRUE.
!         endif

!               inrun = aom_basket%inrun
!          end if
     !In classical path, it's more straightforward.
     !inrun is FALSE at t=2 (first timestep of calculation)
     !so we set it FALSE to begin with, and will set it TRUE as soon as we've done the first calculation
!     elseif (progtype == 'C') then
!          inrun = .FALSE.
!     end if

     !at present, input can also be 'c' which will cause the program to update neighbour lists,
     !or 'n', which will run the calculations but not propagate the electron any further.

     !for each orbital, there will be a pcoefficients file
     !allocate arrays to store filenames, then read in files
!     allocate(pcoefficients(orb))
!     if (.not. allocated(pcoefficients)) then
!          stop "Error: array 'pcoefficients' not allocated correctly."
!     end if
!     allocate(pfiles(orb))
!     if (.not. allocated(pfiles)) then
!          stop "Error: array 'pfiles' not allocated correctly."
!     end if
!     pcoefficients(1) = './coeffjhomop'
     !pcoefficients(2) = './coeffjlumop'
     !pcoefficients(3) = './coeffjlumo1p'
!     do d=1,orb
!          pfiles(d) = d+100
!OPEN          open(unit=pfiles(d),file=pcoefficients(d),action='read')
!     end do

!     open(unit=10,file=mdcrd,action='read')
!     if (progtype == 'S') then
!          open(unit=9,file=mdcrdm,action='read')
!          open(unit=11,file=mdcrdp,action='read')
!     end if
!OPEN     open(unit=17,file=scoefficients,action='read')
!     open(unit=23,file=bondlistfile,action='read')
     !relatively self-explanatory
     mat = orb*mol
     c = mol*atoms

    !if (itimes.le.1) then
    if (itimes.eq.0) then
       ALLOCATE(aom_basket%population(mat))
       ALLOCATE(aom_basket%ctin(mat))
       ALLOCATE(aom_basket%hamiltonian(mat,mat))
       allocate(aom_basket%connlist1(4,atoms))
       ALLOCATE(aom_basket%neighbourlist(1,mol,3,0:atoms))
       ALLOCATE(aom_basket%Stt(mat,mat))
       ALLOCATE(aom_basket%Stmt(mat,mat))
       ALLOCATE(aom_basket%Sttm(mat,mat))
       ALLOCATE(aom_basket%Stmtm(mat,mat))
       ALLOCATE(aom_basket%rm1(c,3))
       ALLOCATE(aom_basket%rm2(c,3))
       ALLOCATE(aom_basket%energies(mol))
       ALLOCATE(aom_basket%bneighbour(3,0:atoms))
       !ALLOCATE(aom_basket%nacv(3,mat,mat,2*mat*atoms))
       ALLOCATE(aom_basket%nacv(3,mat,mat,mat*atoms))
    endif

     !calculate timestep from femtoseconds in atomic units
     !molecular timestep is generally 0.5fs to improve energy conservation
!     deltatMD = (5.0D-16)/(2.4188843265D-17)
!     deltatMD = (ext_timestep*1.0D-15)/(2.4188843265D-17) !CHANGE_AC
     deltatMD = ext_timestep !CHANGE_AC
     !electronic timestep is a fraction of the molecular timestep
!     deltatE = deltatMD/5.0D0
     deltatE = deltatMD/ext_el_steps !CHANGE_AC
!print*, "deltatE", deltatE
!print*, "deltatMD", deltatMD
!print*, "ext_timestep", ext_timestep
!print*, "ext_el_steps", ext_el_steps

     !here are some useful conversion factors
     !1 Angstrom = bohrang Bohr
     bohrang = 1.889725989D0
     !mev = 1 Hartree in meV
     mev = 27211.0D0
     !kcal = 1 Hartree in kcal/mol
     kcal = 627.509D0
     !offset is a number taken away from the absolute value of the hamiltonian diagonal elements
     !the only thing important in the physics is the relative differences of the diagonals
     !but reducing the absolute values can make the RK4 method more stable
     !this option may be useful if you experience inexplicable discontinuities in your electron coefficients
     !offset = 2.433840754D0

     !deltaA is a diagonal matrix which allows arbitrary free energy difference between diabatic states
     allocate(deltaA(mat,mat))
     if (.not. allocated(deltaA)) then
          CPABORT("Error: array 'deltaA' not allocated correctly.")
          !stop "Error: array 'deltaA' not allocated correctly."
     end if
     do i=1,mat
          do j=1,mat
!               deltaA(i,j) = dcmplx(0.0D0, 0.0D0)
               deltaA(i,j) = cmplx(0.0D0, 0.0D0,dp)
          end do
     end do
     !change deltaA(k,k) to add a trap to site k
     !e.g. "deltaA(2,2) = dcmplx(0.0D0, 50.0D0/mev)" -> a 50meV 'trap' in orbital 2
     !e.g. "deltaA(2,2) = cmplx(0.0D0, 50.0D0/mev,dp)" -> a 50meV 'trap' in orbital 2

     !ct are the electron coefficients
     !ctin goes into the subroutine ZGEMV; ctout comes out of it (so essentially ctin = c(t), ctout = c(t+1))
     allocate(ctin(mat))
     if (.not. allocated(ctin)) then
          CPABORT("Error: array 'ctin' not allocated correctly.")
          !stop "Error: array 'ctin' not allocated correctly."
     end if
     allocate(ctout(mat))
     if (.not. allocated(ctout)) then
          CPABORT("Error: array 'ctout' not allocated correctly.")
          !stop "Error: array 'ctout' not allocated correctly."
     end if
     !ctotal is |c(t)|^2: this is the output for diabatic electronic state
     allocate(ctotal(mat))
     if (.not. allocated(ctotal)) then
          CPABORT("Error: array 'ctotal' not allocated correctly.")
          !stop "Error: array 'ctotal' not allocated correctly."
     end if

     !the stateno file contains a numeric representation of the adiabatic electronic state
     !open(unit=24,file='./stateno',action='read')
     !read(24,*) state
     !close(unit=24)
        if (itimes.gt.2) then
           old_state = aom_basket%new_state
        else
           old_state = 1
        endif

     !if the program has been previously called, a bunch of data will have been saved

     if (itimes.ge.3) then
!          !such as electronic coefficients from last time step
!          open(unit=12,file='./C.txt',position='append',action='write')
!          open(unit=18,file='./Cin_temp.txt',action='read')
!          read(18,*) ctin(1:mat)
          ctin(:) = aom_basket%ctin(:)
!          close(unit=18)
!          !and the overlaps we already calculated
!          open(unit=21,file='./stt',action='read')
!          open(unit=26,file='./sttm',action='read')
!          open(unit=27,file='./stmt',action='read')
!          open(unit=28,file='./stmtm',action='read')
     !if not, it may need to be initialised/calculated for the first time
     else if (itimes.eq.2) then
!OPEN          open(unit=12,file='./C.txt',action='write')
          !initialise ct such that ct(init) = 1.0, otherwise 0.0 (i.e. electron localised on molecule 'init')
          do i=1,mat
               ctin(i) = CMPLX(0.0D0, 0.0D0,dp)
               ctotal(i) = CMPLX(0.0D0, 0.0D0,dp)
          end do
          ctin(init*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
          ctotal(init*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
!WRITE          write(12,*) abs(ctotal(1:mat))
                aom_basket%population(1:mat) = abs(ctotal(1:mat)) 
                    CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,2.0D0)
                    CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,2.0D0)
                    CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,2.0D0)
     else
!OPEN          open(unit=12,file='./C.txt',action='write')
          !initialise ct such that ct(init) = 1.0, otherwise 0.0 (i.e. electron localised on molecule 'init')
          do i=1,mat
               ctin(i) = CMPLX(0.0D0, 0.0D0,dp)
               ctotal(i) = CMPLX(0.0D0, 0.0D0,dp)
          end do
          ctin(init*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
          ctotal(init*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
!WRITE          write(12,*) abs(ctotal(1:mat))
                aom_basket%population(1:mat) = abs(ctotal(1:mat)) 
     end if

!     open(unit=13,file='./hab',action='write')
!     open(unit=15,file='./habfull',position='append',action='write')
!     open(unit=18,file='./Cin_temp.txt',action='write')
!     open(unit=24,file='./stateno',action='write')
!     open(unit=30,file='sites',action='write',position='append')
!     open(unit=34,file='./forces',action='read')

     !this file is optional but may be useful if you want to get timing data on the program
     !you will need to insert timing statements into the program and have them write to this file
     !open(unit=19,file='./timings.txt',position='append',action='write')

     !c tells you how many atoms you will have per timestep
     !allocate(crds(3,c,times+1))
     allocate(crds(3,c,3))
     if (.not. allocated(crds)) then
          CPABORT("Error: array 'crds' not allocated correctly.")
          !stop "Error: array 'crds' not allocated correctly."
     end if

     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(Stt(mat,mat))
     if (.not. allocated(Stt)) then
          CPABORT("Error: array 'Stt' not allocated correctly.")
          !stop "Error: array 'Stt' not allocated correctly."
     end if
     allocate(Stpt(mat,mat))
     if (.not. allocated(Stpt)) then
          CPABORT("Error: array 'Stpt' not allocated correctly.")
          !stop "Error: array 'Stpt' not allocated correctly."
     end if
     allocate(Sttp(mat,mat))
     if (.not. allocated(Sttp)) then
          CPABORT("Error: array 'Stt' not allocated correctly.")
          !stop "Error: array 'Sttp' not allocated correctly."
     end if
     allocate(Stptp(mat,mat))
     if (.not. allocated(Stptp)) then
          CPABORT("Error: array 'Stptp' not allocated correctly.")
          !stop "Error: array 'Stptp' not allocated correctly."
     end if
     allocate(Stmt(mat,mat))
     if (.not. allocated(Stmt)) then
          CPABORT("Error: array 'Stmt' not allocated correctly.")
          !stop "Error: array 'Stmt' not allocated correctly."
     end if
     allocate(Sttm(mat,mat))
     if (.not. allocated(Sttm)) then
          CPABORT("Error: array 'Sttm' not allocated correctly.")
          !stop "Error: array 'Sttm' not allocated correctly."
     end if
     allocate(Stmtm(mat,mat))
     if (.not. allocated(Stmtm)) then
          CPABORT("Error: array 'Stmtm' not allocated correctly.")
          !stop "Error: array 'Stmtm' not allocated correctly."
     end if
     allocate(Sttdiff(mat,mat))
     if (.not. allocated(Sttdiff)) then
          CPABORT("Error: array 'Sttdiff' not allocated correctly.")
          !stop "Error: array 'Sttdiff' not allocated correctly."
     end if
     allocate(Stmtmdiff(mat,mat))
     if (.not. allocated(Stmtmdiff)) then
          CPABORT("Error: array 'Stmtmdiff' not allocated correctly.")
          !stop "Error: array 'Stmtmdiff' not allocated correctly."
     end if
     allocate(Stptpdiff(mat,mat))
     if (.not. allocated(Stptpdiff)) then
          CPABORT("Error: array 'Stptpdiff' not allocated correctly.")
          !stop "Error: array 'Stptpdiff' not allocated correctly."
     end if
     allocate(Sttpdiff(mat,mat))
     if (.not. allocated(Sttpdiff)) then
          CPABORT("Error: array 'Sttpdiff' not allocated correctly.")
          !stop "Error: array 'Sttpdiff' not allocated correctly."
     end if
     allocate(Sttmdiff(mat,mat))
     if (.not. allocated(Sttmdiff)) then
          CPABORT("Error: array 'Sttmdiff' not allocated correctly.")
          !stop "Error: array 'Sttmdiff' not allocated correctly."
     end if

     !then the complex ones i.e. H matrices; combinations of H and S

     !print *, "All done!",itimes,439
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H1 etc are always -i*H.
     !This is why they're complex variables

     allocate(H1(mat,mat))
     if (.not. allocated(H1)) then
          CPABORT("Error: array 'H1' not allocated correctly.")
          !stop "Error: array 'H1' not allocated correctly."
     end if
     allocate(H2save(mat,mat))
     allocate(H2(mat,mat))
     if (.not. allocated(H2)) then
          CPABORT("Error: array 'H2' not allocated correctly.")
          !stop "Error: array 'H2' not allocated correctly."
     end if
     allocate(H12(mat,mat))
     if (.not. allocated(H12)) then
          CPABORT("Error: array 'H12' not allocated correctly.")
          !stop "Error: array 'H12' not allocated correctly."
     end if
     allocate(Hdiff(mat,mat))
     if (.not. allocated(Hdiff)) then
          CPABORT("Error: array 'Hdiff' not allocated correctly.")
          !stop "Error: array 'Hdiff' not allocated correctly."
     end if
     allocate(Hsave(mat,mat))
     if (.not. allocated(Hsave)) then
          CPABORT("Error: array 'Hsave' not allocated correctly.")
          !stop "Error: array 'Hsave' not allocated correctly."
     end if
     allocate(X1(mat,mat))
     if (.not. allocated(X1)) then
          CPABORT("Error: array 'X1' not allocated correctly.")
          !stop "Error: array 'X1' not allocated correctly."
     end if
     allocate(X2(mat,mat))
     if (.not. allocated(X2)) then
          CPABORT("Error: array 'X2' not allocated correctly.")
          !stop "Error: array 'X2' not allocated correctly."
     end if
     allocate(X12(mat,mat))
     if (.not. allocated(X12)) then
          CPABORT("Error: array 'X12' not allocated correctly.")
          !stop "Error: array 'X12' not allocated correctly."
     end if

     !finally the K matrices, which are the terms in the Runge-Kutta integration
     !note K1 missing, because it's just a constant times X
     !Ktot is the final constant by which c(t) is multiplied
     allocate(Ktot(mat,mat))
     if (.not. allocated(Ktot)) then
          CPABORT("Error: array 'Ktot' not allocated correctly.")
          !stop "Error: array 'Ktot' not allocated correctly."
     end if
     allocate(K2(mat,mat))
     if (.not. allocated(K2)) then
          CPABORT("Error: array 'K2' not allocated correctly.")
          !stop "Error: array 'K2' not allocated correctly."
     end if
     allocate(K3(mat,mat))
     if (.not. allocated(K3)) then
          CPABORT("Error: array 'K3' not allocated correctly.")
          !stop "Error: array 'K3' not allocated correctly."
     end if
     allocate(K4(mat,mat))
     if (.not. allocated(K4)) then
          CPABORT("Error: array 'K4' not allocated correctly.")
          !stop "Error: array 'K4' not allocated correctly."
     end if

     !Imat is an identity matrix of the dimensions of the other square matrices
     allocate(Imat(mat,mat))
     if (.not. allocated(Imat)) then
          CPABORT("Error: array 'Imat' not allocated correctly.")
          !stop "Error: array 'Imat' not allocated correctly."
     end if
     !define once, then use whenever
     do j=1,mat
          do i=1,mat
               if (i == j) then
                    Imat(i,j) = 1.0D0
               else
                    Imat(i,j) = 0.0D0
               end if
          end do
     end do

     !print *, "All done!",itimes,515
     !bonds is an array that lists all the bonds between atoms in each molecule
     !it lists them as: atom 1, atom 2, bond type
!     allocate(bonds(3,bondsno))
!     if (.not. allocated(bonds)) then
!          stop "Error: array 'bonds' not allocated correctly."
!     end if
!     do i=1,bondsno
!
!          read(23,*) bonds(1,i), bonds(2,i), bonds(3,i)
!     end do

     !non-adiabatic coupling vectors
     !allocate(nacv(3,mat,mat,2*mat*atoms))
     allocate(nacv(3,mat,mat,mat*atoms))
     if (.not. allocated(nacv))then
          CPABORT("Error: array 'nacv' not allocated correctly.")
          !stop "Error: array 'nacv' not allocated correctly."
     end if
     nacv = 0.0D0 !CHANGE_AC 160420
     !diagonals of the Hamiltonian
     !allocate(siteenergies(times+1,mat))
     allocate(siteenergies(3,mat))
     if (.not. allocated(siteenergies)) then
          CPABORT("Error: array 'siteenergies' not allocated correctly.")
          !stop "Error: array 'siteenergies' not allocated correctly."
     end if

     !atomic coefficients for the overlap routine
     allocate(pcoeffs(atoms,orb))
     if (.not. allocated(pcoeffs)) then
          CPABORT("Error: array 'pcoeffs' not allocated correctly.")
          !stop "Error: array 'pcoeffs' not allocated correctly."
     end if
     !pcoeffsrun exists because the p-coefficients need to be renormalised for each structure
     allocate(pcoeffsrun(atoms,mat,3))
     if (.not. allocated(pcoeffsrun)) then
          CPABORT("Error: array 'pcoeffsrun' not allocated correctly.")
          !stop "Error: array 'pcoeffsrun' not allocated correctly."
     end if
     allocate(scoeffs(atoms))
     if (.not. allocated(scoeffs)) then
          CPABORT("Error: array 'scoeffs' not allocated correctly.")
          !stop "Error: array 'scoeffs' not allocated correctly."
     end if

!PRINT     print *, "Reading in coordinates..."
     !read in coordinates
     !this section of code probably needs to be updated depending on coordinate input format
     !some examples follow
!     if (mdcrdtype == 'SH') then
          !Surface hopping is assuming a very simple format from NAMD
          !each line has x, y, z coordinates of that atom
          !three files: R(t-dt), R(t), R(t+dt)
!          do i=1,c
!               read(9,*) crds(1,i,1), crds(2,i,1), crds(3,i,1)
!               read(10,*) crds(1,i,2), crds(2,i,2), crds(3,i,2)
!               read(11,*) crds(1,i,3), crds(2,i,3), crds(3,i,3)
!          end do
!     elseif (mdcrdtype == 'AMBER') then
          !AMBER .mdcrd format
          !This writes 10 coordinates to a line, such that 3 lines sequentially contains 10 atomic positions
!          read(10,*) filler
!          do t=1,times
!               do i=1,10,c
!                    read(10,*) crds(1,i,t), crds(2,i,t), crds(3,i,t), crds(1,i+1,t), crds(2,i+1,t), crds(3,i+1,t),&
!                    & crds(1,i+2,t), crds(2,i+2,t), crds(3,i+2,t), crds(1,i+3,t)
!                    read(10,*) crds(2,i+3,t), crds(3,i+3,t), crds(1,i+4,t), crds(2,i+4,t), crds(3,i+4,t), crds(1,i+5,t),&
!                    & crds(2,i+5,t), crds(3,i+5,t), crds(1,i+6,t), crds(2,i+6,t)
!                    read(10,*) crds(3,i+6,t), crds(1,i+7,t), crds(2,i+7,t), crds(3,i+7,t), crds(1,i+8,t), crds(2,i+8,t),&
!                    & crds(3,i+8,t), crds(1,i+9,t), crds(2,i+9,t), crds(3,i+9,t)
!               end do
!          end do
!     elseif (mdcrdtype == 'XYZ') then
          !standard .xyz format
!          do t=1,times
!               read(10,*) filler
!               read(10,*) filler
!               do i=1,c
!                    read(10,*) filler, crds(1,i,t), crds(2,i,t), crds(3,i,t)
!               end do
!          end do
!     else
!          print *, "This does not appear to be a valid coordinate file designation."
!          print *, "Please review the program's input parameters."
!          stop
!     end if

!CHANGE_AC
!READ THE COORDINATES
        do coord=1,3
              if (itimes.ge.2) then
!print*, "hey"
                 crds(coord,:,1) = PACK(aom_basket%rm2(:,coord), ext_pres)
                 crds(coord,:,2) = PACK(aom_basket%rm1(:,coord), ext_pres)
                 crds(coord,:,3) = PACK(positions_plus(:,coord), ext_pres)
              else if (itimes.eq.1) then 
                 crds(coord,:,2) = PACK(aom_basket%rm1(:,coord), ext_pres)
                 crds(coord,:,3) = PACK(positions_plus(:,coord), ext_pres)
              else if (itimes.eq.0) then 
                 crds(coord,:,3) = PACK(positions_plus(:,coord), ext_pres)
              endif
        enddo
!END CHANGE_AC
!print*, crds(:,1,3)
!print*, positions_plus(1,:)
!print*, ext_pres

 !    print *, "All done!",itimes,612
!CHANGE_AC
!READ THE COORDINATES
     !now we also read in the atomic coefficients
        ALLOCATE(pcoefftot(c))
        ALLOCATE(scoefftot(c))
        pcoefftot = ext_pcoeff
        scoefftot = ext_scoeff
!        pcoefftot = PACK(ext_pcoeff, ext_pres)
!        scoefftot = PACK(ext_scoeff, ext_pres)
!END CHANGE_AC
     do i=1,atoms
          do j=1,orb
!               read(pfiles(j),*) pcoeffs(i,j)
                       pcoeffs(i,j) = pcoefftot( i + (j-1)*orb)
          end do
!          read(17,*) scoeffs(i)
                scoeffs(i) = scoefftot(i)
     end do
     
!PRINT     print *, "Done! Calculating centres of mass..."

     !calculate CoMs at first timestep
     !coms(1 to 3,i) are the x, y and z coordinates of CoM of i
     allocate(coms(3,mol))
     do i=1,mol
          !initialise CoMs to 0
          do d=1,3
               coms(d,i) = 0.0D0
          end do
          do m=1,atoms
               do d=1,3
                    !coms(d,i) = coms(d,i) + (crds(d,(i-1)*atoms + m,1))/atoms
                    coms(d,i) = coms(d,i) + (crds(d,(i-1)*atoms + m,3))/atoms
               end do
          end do
     end do
!print*, "coms"
!print*, coms
     !construct neigbours list: 1 if close enough, 0 if not
!PRINT     print *, "Constructing nearest neighbours list..."
     allocate(NN(mat,mat))
     if (.not. allocated(NN)) then
          CPABORT("Error: array 'NN' not allocated correctly.")
          !stop "Error: neighbour list not allocated correctly."
     end if
     allocate(dist(mol,mol))
     if (.not. allocated(dist)) then
          CPABORT("Error: array 'dist' not allocated correctly.")
          !stop "Error: distance list not allocated correctly."
     end if
     do i=1,mol
          do j=1,mol
               dist(i,j) = ((coms(1,i) - coms(1,j))*(coms(1,i) - coms(1,j)) + (coms(2,i) - coms(2,j))*&
               &(coms(2,i) - coms(2,j)) + (coms(3,i) - coms(3,j))*(coms(3,i) - coms(3,j)))
               if (dist(i,j) < cutoff*cutoff) then
                    do d=1,orb
                         do m=1,orb
                              NN(orb*(i-1)+d,orb*(j-1)+m) = 1
                         end do
                    end do
               else
                    do d=1,orb
                         do m=1,orb
                              NN(orb*(i-1)+d,orb*(j-1)+m) = 0
                         end do
                    end do
               end if                    
          end do
     end do
!print*, "NN"
!print*, dist
!print*, NN
!PRINT     print *, "Nearest neighbours found."

     !allocate the arrays that will contain atomic coordinates for the overlap routine
     !first column is now atomic number of atom (e.g. 6 for carbon)
     !Naming convention: At = A(t), Atp = A(t+1), Atm = A(t-1), etc
     allocate(At(4,atoms))
     if (.not. allocated(At)) then
          CPABORT("Error: array 'At' not allocated correctly.")
          !stop "Error: array 'At' not allocated correctly."
     end if
     allocate(Atp(4,atoms))
     if (.not. allocated(Atp)) then
          CPABORT("Error: array 'Atp' not allocated correctly.")
          !stop "Error: array 'Atp' not allocated correctly."
     end if
     allocate(Atm(4,atoms))
     if (.not. allocated(Atm)) then
          CPABORT("Error: array 'Atm' not allocated correctly.")
          !stop "Error: array 'Atm' not allocated correctly."
     end if
     allocate(Bt(4,atoms))
     if (.not. allocated(Bt)) then
          CPABORT("Error: array 'Bt' not allocated correctly.")
          !stop "Error: array 'Bt' not allocated correctly."
     end if
     allocate(Btp(4,atoms))
     if (.not. allocated(Btp)) then
          CPABORT("Error: array 'Btp' not allocated correctly.")
          !stop "Error: array 'Btp' not allocated correctly."
     end if
     allocate(Btm(4,atoms))
     if (.not. allocated(Btm)) then
          CPABORT("Error: array 'Btm' not allocated correctly.")
          !stop "Error: array 'Btm' not allocated correctly."
     end if

     !array of diagonal atomic forces
!     allocate(F(3,mol*atoms,mol))
!     if (.not. allocated(F)) then
!          stop "Error: array 'F' not allocated correctly."
!     end if

     !print *, "All done!",itimes,715
     !elements contains the atomic number of each atom (this goes into At, Bt etc)
     !it is user-defined input
     allocate(elements(atoms))
     if (.not. allocated(elements)) then
          CPABORT("Error: array 'elements' not allocated correctly.")
          !stop "Error: array 'elements' not allocated correctly."
     end if
!OPEN     open(unit=44,file=elementfile,action='read')
!CHANGE_AC
     ALLOCATE(elementtot(c))
     elementtot = PACK(ext_elements, ext_pres)
!END CHANGE_AC
     do i=1,atoms
                elements(i) = elementtot(i)
!          read(44,*) elements(i)
     end do
!     close(44)

     !atomlist is a list of the atom numbers of the atoms of interest
     !i.e which atoms we will calculate the orbitals for
     !the information is contained in an input file called "atomlist"
!OPEN     open(unit=36,file='atomlist',action='read')
!     i = 0
!     stat = 0
!     read(36,*,iostat=stat) filler
!     do while(stat >= 0)
!          if (stat < 0) then
!               exit
!          else
!               i = i + 1
!          end if
!          read(36,*,iostat=stat) filler
!     end do
!     close(36)

!CHANGE_AC
        allocate(atomlist_mask(c))
        atomlist_mask =.NOT.((scoefftot.eq.0).AND.(pcoefftot.eq.0))
        allocate(atomlist(COUNT(atomlist_mask(1:atoms))))
        atomlist = PACK( (/ (i, i=1,atoms) /), atomlist_mask(1:atoms) )
    if (itimes.eq.0) then
       allocate(aom_basket%connlist2(4,size(atomlist)))
    endif
!END CHANGE_AC
        

!     allocate(atomlist(i)
!     if (.not. allocated(atomlist)) then
!          stop "Error: array 'atomlist' not allocated correctly."
!     end if
!OPEN     open(unit=36,file='atomlist',action='read')
!     do i=1,size(atomlist)
!          read(36,*) atomlist(i)
!     end do
!     close(36)

     !print *, "All done!",itimes,768
     !these arrays are data for the non-adiabatic coupling calculation
     !they help efficiency by telling the program which atoms are connected
     !and thus which pairs to carry out the calculation for
     allocate(connlist1(4,atoms))
     if (.not. allocated(connlist1)) then
          CPABORT("Error: array 'connlist1' not allocated correctly.")
          !stop "Error: array 'connlist1' not allocated correctly."
     end if
     allocate(connlist2(4,size(atomlist)))
     if (.not. allocated(connlist2)) then
          CPABORT("Error: array 'connlist2' not allocated correctly.")
          !stop "Error: array 'connlist2' not allocated correctly."
     end if
     !allocate(aneighbour(3,0:atoms))
     !allocate(aneighbour(4,0:atoms))
     allocate(aneighbour(0:atoms,0:atoms))
     if (.not. allocated(aneighbour)) then
          CPABORT("Error: array 'aneighbour' not allocated correctly.")
          !stop "Error: array 'aneighbour' not allocated correctly."
     end if
     !allocate(bneighbour(3,0:atoms))
     !allocate(bneighbour(4,0:atoms))
     allocate(bneighbour(0:atoms,0:atoms))
     if (.not. allocated(bneighbour)) then
          CPABORT("Error: array 'bneighbour' not allocated correctly.")
          !stop "Error: array 'bneighbour' not allocated correctly."
     end if
     !allocate(neighbourlist(mol,3,0:atoms))
     !allocate(neighbourlist(mol,4,0:atoms))
     allocate(neighbourlist(mol,0:atoms,0:atoms))
     if (.not. allocated(neighbourlist)) then
          CPABORT("Error: array 'neighbourlist' not allocated correctly.")
          !stop "Error: array 'neighbourlist' not allocated correctly."
     end if
     neighbourlist =0
     !we fill out the connlist arrays before calculating the nacv
     do i=1,atoms
          do j=1,3
               !At(j,i) = bohrang*crds(j,i,1)
               !At(j,i) = crds(j,i,1)
               At(j,i) = crds(j,i,3) !CHANGE_AC : 3 is the current timestep
          end do
          !print*, "a",At(:,i)
          !print*, "a",crds(:,i,3)
     end do
     !print *, "All done!",itimes,800
     !if (.not. inrun) then
!print*, "atomlist", atomlist

     if (itimes.eq.0) then
          !call connect_list2(At, atomlist, connlist1, connlist2)
          call connect_list2(At(1:3,:), atomlist, connlist1, connlist2)
!          open(unit=40,file='./connectlist1',action='write')
!          open(unit=41,file='./connectlist2',action='write')
!print*, connlist1(:,1), "connlist1 eq 0"
!print*, connlist1(:,2), "connlist1 eq 0"
!print*, connlist1(:,3), "connlist1 eq 0"
!print*, connlist1(:,4), "connlist1 eq 0"
!print*, connlist1(:,5), "connlist1 eq 0"
!print*, connlist1(:,6), "connlist1 eq 0"
!print*, connlist2(:,1), "connlist2 eq 0"
!print*, connlist2(:,2), "connlist2 eq 0"
          do i=1,atoms
!WRITE               write(40,*), connlist1(1:4,i)
                        aom_basket%connlist1(1:4,i) = connlist1(1:4,i)
!print*, "here"
!print*, "At", At(:,i)
!print*, connlist1(:,i)
          end do
          do i=1,size(atomlist)
!WRITE               write(41,*), connlist2(1:4,i)
                        aom_basket%connlist2(1:4,i) = connlist2(1:4,i)
          end do
!print*, aom_basket%connlist2, "here here"
!          close(unit=40)
!          close(unit=41)
     else if (itimes.gt.0) then
!          open(unit=40,file='./connectlist1',action='read')
!          open(unit=41,file='./connectlist2',action='read')
!          do i=1,atoms
!               read(40,*), connlist1(1:4,i)
                        connlist1(1:4,:) = aom_basket%connlist1(1:4,:)
!          end do
!print*, aom_basket%connlist2, "here here and and"
!          do i=1,size(atomlist)
!               read(41,*), connlist2(1:4,i)
                        connlist2(1:4,:) = aom_basket%connlist2(1:4,:)
!print*, atomlist
!print*, connlist1(:,1), "connlist1 gt 0"
!print*, connlist1(:,2), "connlist1 gt 0"
!print*, connlist1(:,3), "connlist1 gt 0"
!print*, connlist1(:,4), "connlist1 gt 0"
!print*, connlist1(:,5), "connlist1 gt 0"
!print*, connlist1(:,6), "connlist1 gt 0"
!print*, connlist2(:,1), "connlist2 gt 0"
!print*, connlist2(:,2), "connlist2 gt 0"
!          end do
!          close(unit=40)
!          close(unit=41)
     end if

     !used to give correct commands to ZGEMX subroutines later
     zgemm_type = 'N'

     !print *, "All done!",itimes,833
     !for each timestep: this loop is (probably) the main bulk of the code
!PRINT     print *, "Calling overlap routine for each pair along the trajectory (this may take a while)"
!     do t=2,times
t =2
           if (itimes.ge.2) then
!          if (progtype == 'C') then
!PRINT               print *, "Molecular timestep", t, "beginning."
!          end if
          !for each molecule
          do i=1,mol
               !create the array of all its atomic coordinates at t
               !numerical factor converts from Angstroms into Bohr radii
               do j=1,atoms
                    At(1,j) = dble(elements(j))
                    Atp(1,j) = dble(elements(j))
                    Atm(1,j) = dble(elements(j))
                    do d=2,4
                         At(d,j) = crds(d-1,(i-1)*atoms+j,t)
                         Atp(d,j) = crds(d-1,(i-1)*atoms+j,t+1)
                         Atm(d,j) = crds(d-1,(i-1)*atoms+j,t-1)
                         !At(d,j) = bohrang*crds(d-1,(i-1)*atoms+j,t)
                         !Atp(d,j) = bohrang*crds(d-1,(i-1)*atoms+j,t+1)
                         !Atm(d,j) = bohrang*crds(d-1,(i-1)*atoms+j,t-1)
                    end do
               end do
               !calculate the site energy for this molecule
               !at this timestep
!               call sitecalc(At, bonds, atoms, siteenergies(t,i))
                        siteenergies(t,i)   = aom_basket%energies(i) ! CHANGE_AC
                        siteenergies(t+1,i) = energies_plus(i) ! CHANGE_AC
!print*, energies_plus(i), "koko"
!print*, siteenergies(t,i), siteenergies(t+1,i), "koko"

               !then one timestep ahead: also calculate diagonal forces at this timestep
!               call sitecalc(Atp, bonds, atoms, siteenergies(t+1,i), F(:,:,i), mol, i)
!WRITE               write(30,*) i, siteenergies(t+1,i)
!                        aom_basket%hamiltonian(i,i) = siteenergies(t+1,i)

               !then for each other molecule
               do k=1,mol
                    !create the array of its atomic coordinates at both t and t+dt
                    do l=1,atoms
                         Bt(1,l) = dble(elements(l))
                         Btp(1,l) = dble(elements(l))
                         Btm(1,l) = dble(elements(l))
                         do d=2,4
                              Bt(d,l) =  crds(d-1,(k-1)*atoms+l,t)
                              Btp(d,l) = crds(d-1,(k-1)*atoms+l,t+1)
                              Btm(d,l) = crds(d-1,(k-1)*atoms+l,t-1)
                              !Bt(d,l) = bohrang*crds(d-1,(k-1)*atoms+l,t)
                              !Btp(d,l) = bohrang*crds(d-1,(k-1)*atoms+l,t+1)
                              !Btm(d,l) = bohrang*crds(d-1,(k-1)*atoms+l,t-1)
                         end do
                    end do
                    if (i == 1) then
!print*, "lololo", "k =",k
                         do d=1,orb
                              do l=1,atoms
                                   pcoeffsrun(l,orb*(k-1)+d,1) = pcoeffs(l,d)
                              end do
                              call calc_sab(Btm, connlist2, pcoeffsrun(:,orb*(k-1)+d,1))
                              do l=1,atoms
!print*, "pcoeffs l d",d,l,pcoeffs(l,d)
                                   pcoeffsrun(l,orb*(k-1)+d,2) = pcoeffs(l,d)
!print*, "pcoeffsrun l d",d,l,pcoeffsrun(l,orb*(k-1)+d,2)
                              end do
                              call calc_sab(Bt, connlist2, pcoeffsrun(:,orb*(k-1)+d,2))
                              do l=1,atoms
                                   pcoeffsrun(l,orb*(k-1)+d,3) = pcoeffs(l,d)
                              end do
                              call calc_sab(Btp, connlist2, pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
!print*, "pcoeffsrun"
!print*, pcoeffsrun(:,:,1)
!print*, pcoeffsrun(:,1,2)
!print*, pcoeffsrun(:,2,2)
!print*, "pcoeffsrun"
!                         if ((input == 'c') .or. (inrun .eqv. .FALSE.)) then
                         if (itimes.eq.2) then
                              do d=1,orb
                                   do m=1,orb
!print*, 947, d,m
!print*, At(1,:)
!print*, At(2,:)
!print*, At(3,:)
!print*, At(4,:)
                                        call sab_screen_new(At, atomlist, connlist1, pcoeffsrun(:,orb*(i-1)+d,2),&
                                        & atomlist, Bt, atomlist, connlist1, pcoeffsrun(:,orb*(k-1)+m,2),&
                                        & aneighbour(:,1:atoms), bneighbour(:,1:atoms), minao)
                                   end do
                              end do
                              neighbourlist(k,1:3,0) = [0, 0, 0]
                              neighbourlist(k,:,1:atoms) = bneighbour(:,1:atoms)
!print*, 954, "bneighbour", bneighbour(:,1:atoms)
!print*, "lololo itimes eq 2"
                              if (k == mol) then
!                                   open(unit=43,file='./bneighbour',action='write')
!print*, "lilili"
                                   do m=1,mol
                                        do j=1,atoms
!WRITE                                             write(43,*), neighbourlist(m,:,j)
!print*, neighbourlist(m,:,j), "lilili"
                                                                aom_basket%neighbourlist(1,m,:,j) = & 
                                                                        neighbourlist(m,:,j)
!print*, aom_basket%neighbourlist(m,:,j), "lilili"
                                        end do
                                   end do
                              !     close(43)
                              end if
                         else
!                              open(unit=43,file='./bneighbour',action='read')
                              do m=1,mol
                                   neighbourlist(m,1:3,0) = [0, 0, 0]
                                   do j=1,atoms
!print*, aom_basket%neighbourlist(m,:,j), "lalala"
                                                                neighbourlist(m,:,j)=&
                                                                aom_basket%neighbourlist(1,m,:,j)
!print*, neighbourlist(m,:,j), "lalala"
!                                        read(43,*), neighbourlist(m,:,j)
                                   end do
                              end do
!                              close(43)
                         end if !(inrun.eq.FALSE.)
                    end if  ! (i == 1)
                    !only call overlap if neighbours list says so
                    if (NN(i,k) == 1) then
                         !first the diagonals
                         if (i == k) then
                              do d=1,orb
                                   do m=1,orb
                                        !These diagonals are 1 by definition
                                        Stptp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
                                        !H1(orb*(i-1)+d,orb*(k-1)+d) = dcmplx(0.0D0, (-siteenergies(t,k)/kcal))
                                        H1(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(t,k)), dp)
                                        !H2(orb*(i-1)+d,orb*(k-1)+d) = dcmplx(0.0D0, (-siteenergies(t+1,k)/kcal))
                                        H2(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(t+1,k)),dp)
                                        call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Btp,&
                                        & connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                        Sttp(orb*(i-1)+d,orb*(k-1)+m) = S
!print*, "WOOHOO"
                                        call calc_sab(Atp, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3), S, Bt,&
                                        & connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                        Stpt(orb*(i-1)+d,orb*(k-1)+m) = S
!print*, "WAAHAA"
!print*, "yellow", pcoeffsrun(1:atoms,orb*(i-1)+d,2)
!print*, "i==k", i
                                        !call calc_dRSab(At, atomlist, connlist1, &
                                        call calc_dRSab(Atp, atomlist, connlist1, &
                                        &pcoeffsrun(1:atoms,orb*(i-1)+d,2), scoeffs, &
                                        &nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(k-1)*atoms+1:k*atoms), &
                                        &nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(i-1)*atoms+1:i*atoms), &
                                        !&nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,1:atoms), &
                                        !&nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,atoms+1:2*atoms), &
                                        !&Bt, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2), scoeffs, &
                                        &Btp, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2), scoeffs, &
                                        !&neighbourlist(i,:,:), neighbourlist(k,:,:))
                                        &neighbourlist(i,:,:), neighbourlist(k,:,:), dR)
!print*, "yellow", pcoeffsrun(1:atoms,orb*(i-1)+d,2)
                                   end do
                              end do

                              if (itimes.eq.2) then
                                   do d=1,orb
                                        do m=1,orb
!print*, At
!print*, connlist2
!print*, pcoeffsrun(1:atoms,orb*(i-1)+d,2)
!print*, Bt
!print*, pcoeffsrun(1:atoms,orb*(k-1)+m,2)
                                             call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2),&
                                             & S, Bt, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                             Stt(orb*(i-1)+d,orb*(k-1)+m) = S
!print*, "S", S
                                             call calc_sab(Atm, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,1),&
                                             & S, Bt, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                             Stmt(orb*(i-1)+d,orb*(k-1)+m) = S
                                             call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2),&
                                             & S, Btm, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                                             Sttm(orb*(i-1)+d,orb*(k-1)+m) = S
                                             call calc_sab(Atm, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,1),&
                                             & S, Btm, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                                             Stmtm(orb*(i-1)+d,orb*(k-1)+m) = S
                                        end do
                                   end do     
                              end if ! (.not.inrun)
                         !now the off-diagonals                         
                         else ! ( i == k)
                              !symmetric matrices
                              !no sense calling calc_sab twice as often as you need to
                              if (k < i) then
                                   do m=1,orb
                                        do d=1,orb
                                             H1(orb*(i-1)+d,orb*(k-1)+m) = H1(orb*(k-1)+m,orb*(i-1)+d)
                                             H2(orb*(i-1)+d,orb*(k-1)+m) = H2(orb*(k-1)+m,orb*(i-1)+d)
                                             Stt(orb*(i-1)+d,orb*(k-1)+m) = Stt(orb*(k-1)+m,orb*(i-1)+d)
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = Stptp(orb*(k-1)+m,orb*(i-1)+d)
                                             Stmtm(orb*(i-1)+d,orb*(k-1)+m) = Stmtm(orb*(k-1)+m,orb*(i-1)+d)
                                             !nacv(:,i,k,1:atoms) = nacv(:,k,i,atoms+1:2*atoms)
                                             !nacv(:,i,k,1:atoms) = 0.0D0
                                             !nacv(:,i,k,atoms+1:2*atoms) = nacv(:,k,i,1:atoms)
                                             !nacv(:,i,k,atoms+1:2*atoms) = 0.0D0
                                        end do
                                   end do
                              else
                                   do d=1,orb
                                        do m=1,orb
                                             call calc_sab(Atp, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                                             & S, Btp, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                             ! H2(orb*(i-1)+d,orb*(k-1)+m) = dcmplx(0.0D0, -scaling*S)
                                             H2(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)
!print*, "i,k",i,k
                                             !call calc_dRSab(At, atomlist, connlist1, &
                                             call calc_dRSab(Atp, atomlist, connlist1, &
                                             &pcoeffsrun(1:atoms,orb*(i-1)+d,2), scoeffs, &
                                             &nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(k-1)*atoms+1:k*atoms), &
                                             &nacv(1:3,orb*(k-1)+m,orb*(i-1)+d,(i-1)*atoms+1:i*atoms), &
                                             !&nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(k-1)*atoms+1:2*(k-1)*atoms), &
                                             !&nacv(1:3,orb*(k-1)+m,orb*(i-1)+d,(i-1)*atoms+1:2*(i-1)*atoms), &
                                             &Btp, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2),&
                                             !&Bt, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2),&
                                             !&scoeffs, neighbourlist(i,:,:), neighbourlist(k,:,:))
                                             &scoeffs, neighbourlist(i,:,:), neighbourlist(k,:,:), dR)
                                        end do
                                   end do

                                   !if (.not. inrun) then
                                   if (itimes.eq.2) then
                                        do d=1,orb
                                             do m=1,orb
                                                  call calc_sab(At, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Bt, connlist2,&
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                                  Stt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
!print*, cbar*S
                                                  ! H1(orb*(i-1)+d,orb*(k-1)+m) = dcmplx(0.0D0, -scaling*S)
                                                  H1(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)
                                                  call calc_sab(Atm, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, Btm, connlist2,&
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                                                  Stmtm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                                  call calc_sab(At, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Btm, connlist2,&
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                                                  Sttm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                                  call calc_sab(Atm, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, Bt, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                                  Stmt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                                  call calc_sab(Bt, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,2), S, Atm, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,1))
                                                  Sttm(orb*(k-1)+m,orb*(i-1)+d) = cbar*S
                                                  call calc_sab(Btm, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,1), S, At, connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,2))
                                                  Stmt(orb*(k-1)+m,orb*(i-1)+d) = cbar*S
                                             end do
                                        end do
                                   end if ! (.not.inrun)
                              end if ! ( k < i)

                              do d=1,orb
                                   do m=1,orb
                                        call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Btp,&
                                        & connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                        Sttp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                        call calc_sab(Atp, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3), S, Bt,&
                                        & connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                        Stpt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                        call calc_sab(Bt, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2), S, Atp,&
                                        & connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3))
                                        Sttp(orb*(k-1)+m,orb*(i-1)+d) = cbar*S
                                        call calc_sab(Btp, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3), S, At,&
                                        & connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2))
                                        Stpt(orb*(k-1)+m,orb*(i-1)+d) = cbar*S
                                   end do
                              end do
                         end if !(i==k)
                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   !H1(orb*(i-1)+d,orb*(k-1)+m) = dcmplx(0.0D0, 0.0D0)
                                   H1(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, 0.0D0,dp)
                                   !H2(orb*(i-1)+d,orb*(k-1)+m) = dcmplx(0.0D0, 0.0D0)
                                   H2(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, 0.0D0,dp)
                                   Stt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Sttp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Stptp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Stpt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Sttm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Stmt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                              end do
                         end do
                    end if !(NN(i,k) == 1
               end do !loop over k
          end do !loop over i
          !during a simulation, you can read in a lot of these overlaps from a previous timestep
          !if (inrun) then
          if (itimes.gt.2) then
               do i=1,mat
                    do j=1,mat
!                         read(21,*) Stt(j,i)
                                        Stt(j,i) = aom_basket%Stt(j,i)
                         if (i .ne. j) then
                              !H1(j,i) = dcmplx(0.0D0, -1.0D0*scaling*Stt(j,i)/cbar)
                              H1(j,i) = cmplx(0.0D0, -1.0D0*scaling*Stt(j,i)/cbar,dp)
                         end if
                                        Sttm(j,i)  = aom_basket%Sttm(j,i)
                                        Stmt(j,i)  = aom_basket%Stmt(j,i)
                                        Stmtm(j,i) = aom_basket%Stmtm(j,i)
!                         read(26,*) Sttm(j,i)
!                         read(27,*) Stmt(j,i)
!                         read(28,*) Stmtm(j,i)
                    end do
               end do
          end if !(inrun)
!print*, "Stptp", Stptp
!print*, "Stpt", Stpt
!print*, "Sttp", Sttp
!print*, "Stt", Stt
!print*, "Sttm", Sttm
!print*, "Stmt", Stmt
!print*, "Stmtm", Stmtm
!WRITE          write(13,*) -kcal*imagpart(H2(1,orb+1))
     !print *, "All done!",itimes,1106

!          close(unit=21)
!          close(unit=26)
!          close(unit=27)
!          close(unit=28)
!          open(unit=21,file='./stt',action='write')
!          open(unit=26,file='./sttm',action='write')
!          open(unit=27,file='./stmt',action='write')
!          open(unit=28,file='./stmtm',action='write')
         if (itimes.ge.2) then !CHANGE_AC
          do i=1,mat
               do j=1,mat
                                aom_basket%Stt(j,i)   = Stptp(j,i)
                                aom_basket%Sttm(j,i)  = Stpt(j,i)
                                aom_basket%Stmt(j,i)  = Sttp(j,i)
                                aom_basket%Stmtm(j,i) = Stt(j,i)
!                    write(21,*) Stptp(j,i)
!                    write(26,*) Stpt(j,i)
!                    write(27,*) Sttp(j,i)
!                    write(28,*) Stt(j,i)
               end do
          end do
         endif !CHANGE_AC
!               close(unit=21)
!          close(unit=26)
!          close(unit=27)
!          close(unit=28)

          nacv = cbar*nacv
!TO REMOVE AFTER

!PRINT          print *, "All pairs calculated."

          !It's sometimes necessary to just call this program to calculate forces, without propagating the electron.
!          if (input /= 'n') then
               !MD is done on a molecular timestep of femtoseconds: however electron dynamics happen on attoseconds
               !so for each molecular timestep, the program calculates for multiple electronic timesteps
               !during this process, the H and S matrices are found by linear interpolation between the original S1 and S2
               Hdiff = (deltatE/deltatMD)*(H2 - H1)
!print*, "Hdiff", Hdiff
               Sttdiff = (deltatE/deltatMD)*(Stptp - Stt)
               Sttpdiff = (deltatE/deltatMD)*(Sttp - Stmt)
               Sttmdiff = (deltatE/deltatMD)*(Stpt - Sttm)
               Stptpdiff = Sttdiff
               Stmtmdiff = (deltatE/deltatMD)*(Stt - Stmtm)
               Hsave = H1
               H2save=H2
!PRINT               print *, "Propagating the electron for one molecular timestep..."
               do t2=1,int(deltatMD/deltatE)

!print*, "Stt",Stt
                    call makeX(mat, H1, Stt, Sttp, Sttm, Stptp, Stmtm, X1, deltatMD)
!print*, "H1", H1
!print*, "X1",X1
                    !create timestepped equivalents
                    H12 = H1 + 0.5D0*Hdiff
                    Stt = Stt + 0.5D0*Sttdiff
                    Sttp = Sttp + 0.5D0*Sttpdiff
                    Sttm = Sttm + 0.5D0*Sttmdiff
                    Stptp = Stptp + 0.5D0*Stptpdiff
                    Stmtm = Stmtm + 0.5D0*Stmtmdiff
                    call makeX(mat, H12, Stt, Sttp, Sttm, Stptp, Stmtm, X12, deltatMD)
!print*, "H12", H12
!print*, "X12", X12
     
                    H2 = H12 + 0.5D0*Hdiff
                    Stt = Stt + 0.5D0*Sttdiff
                    Sttp = Sttp + 0.5D0*Sttpdiff
                    Sttm = Sttm + 0.5D0*Sttmdiff
                    Stptp = Stptp + 0.5D0*Stptpdiff
                    Stmtm = Stmtm + 0.5D0*Stmtmdiff
                    call makeX(mat, H2, Stt, Sttp, Sttm, Stptp, Stmtm, X2, deltatMD)
!print*, "H2", H2
!print*, "X2", X2
!print*, "vvvv"
!print*, mat
!print*, H2
!print*, Stt
!print*, Sttp
!print*, Sttm
!print*, Stptp
!print*, Stmtm
!print*, X2
!print*, deltatMD
!print*, "vvvv"
     
                    !K2 = dcmplx(0.0D0, 0.0D0)
                    K2 = cmplx(0.0D0, 0.0D0,dp)
                    !K3 = dcmplx(0.0D0, 0.0D0)
                    K3 = cmplx(0.0D0, 0.0D0,dp)
                    !K4 = dcmplx(0.0D0, 0.0D0)
                    K4 = cmplx(0.0D0, 0.0D0,dp)
          
                    !now that X has been found, use 4th order Runge-Kutta numerical integration method
                    !aiming to find operator such that c(t+1) = (1 + operator)*c
          
                    !there are four terms in X: each is an increasing multiple. so first term is proportional to X(t),
                    !second term is as X(t+deltat/2)*X(t), third term is as X(t+deltat/2)*X(t+deltat/2)*X(t),
                    !fourth term is as X(t+deltat)*X(t+deltat/2)*X(t+deltat/2)*X(t)
                    !this means each term is some X matrix times the previous term,
                    !so only three more calls to ZGEMM are needed
                    !the K matrices are the second, third and fourth Runge-Kutta terms
     
                    !alpha = dcmplx(deltatE/2.0D0, 0.0D0)
                    alpha = cmplx(deltatE/2.0D0, 0.0D0,dp)
                    !beta = dcmplx(0.0D0, 0.0D0)
                    beta = cmplx(0.0D0, 0.0D0,dp)
     
                    !second term
                    call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X12, mat, X1, &
                    &mat, beta, K2, mat)
                    K2 = K2 + X12
!print*, "K2", K2                    
                    !third term
                    !call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X12, mat, X2, &
                    call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X12, mat, K2, &
                    &mat, beta, K3, mat)
                    K3 = K3 + X12
!print*, "K3", K3 
                    !fourth term
                    !alpha = dcmplx(deltatE, 0.0D0)
                    alpha = cmplx(deltatE, 0.0D0,dp)
                    call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X2, mat, K3, &
                    &mat, beta, K4, mat)
                    K4 = K4 + X2
!print*, "K4", K4    
!print*, "deltatE", deltatE 
                    !put everything together
                    Ktot = Imat + (deltatE/6.0D0)*(X1 + 2*K2 + 2*K3 + K4)
!print*, "Ktot", Ktot
!print*, "wwww"
!print*, Imat
!print*, (deltatE/6.0D0)
!print*, X1
!print*, 2*K2
!print*, 2*K3
!print*, K4
!print*, Ktot
!print*, "wwww"
     
                    !timestep electron propagation and print out
                    !alpha = dcmplx(1.0D0, 0.0D0)
                    alpha = cmplx(1.0D0, 0.0D0,dp)
                    inc = 1
!print*, "ctin", ctin
                    call ZGEMV(zgemm_type, mat, mat, alpha, Ktot, mat, ctin, inc, beta, ctout, inc)
!print*, "ctin", ctin
!print*, "ctout", ctout
!print*, t2, zgemm_type, mat, mat, alpha, Ktot, mat, ctin, inc, beta, ctout, inc
!print*, itimes, ctout, "tttt"
!print*, itimes, ctout, "tttt"
     
                    ctmag = 0.0D0
                    !ctcomp = dcmplx(0.0D0, 0.0D0)
                    ctcomp = cmplx(0.0D0, 0.0D0,dp)
                    do i=1,mat
                         ctotal(i) = conjg(ctout(i))*ctout(i)
                         ctcomp = ctcomp + ctotal(i)
                    end do
                    ctmag = abs(ctcomp)
                    !rescale to conserve norm = 1
                    do i=1,mat
                         !ctout(i) = dcmplx((1.0D0/sqrt(ctmag))*realpart(ctout(i)), (1.0D0/sqrt(ctmag))*imagpart(ctout(i)))
                         ctout(i) = cmplx((1.0D0/sqrt(ctmag))*REAL(ctout(i)), (1.0D0/sqrt(ctmag))*AIMAG(ctout(i)),dp)
                    end do
                    do i=1,mat
                         ctotal(i) = conjg(ctout(i))*ctout(i)
                    end do
                    !if (t2 == int(deltatMD/deltatE)) then
!WRITE                         write(12,*) abs(ctotal(1:mat))
                                        aom_basket%population(1:mat) =&
                                            abs(ctotal(1:mat)) 
                    !end if
                    !update
                    H1 = H2
!print*, "ctout", ctout
                    ctin = ctout
!print*, "ctin", ctin
                    !ctout = dcmplx(0.0D0, 0.0D0)
                    ctout = cmplx(0.0D0, 0.0D0,dp)

                    times = DBLE(itimes) + t2*1/DBLE(ext_el_steps)   
                    CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,times)
                    CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,times)
                    CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,times)

               end do
!PRINT               print *, "Done."     
               H1 = Hsave
!          end if

!          if (progtype == 'C') then
!PRINT               print *, "Molecular timestep", t, "done."
!          else
!               print *, "Calculating surface hopping probabilities."

               !The adiabaticity subroutine diagonalises the diabatic Hamiltonian
               !It gives us back the eigenvalues and eigenvectors
!               call adiabaticity(-imagpart(H1), mat, C1, E1, .false.)
!               call adiabaticity(-imagpart(H2), mat, C2, E2, .true.)

               !The eigenvectors are needed to calculate the surface hopping probabilities
               !The surfacehop subroutine calculates these probabilites and decides on the
               !new (possibly unchanged) surface for simulation
!               if (input /= 'n') then
!                    call surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, ctin, C1, C2, deltatMD, state, decohere)
!               end if

               !The forces subroutine then calculates the adiabatic forces and prints them out
!               F = F/(kcal*bohrang)
!               call forces(state, -imagpart(H2), C2, nacv, NN, atoms, scaling, cbar, F, Stt)

               !Finally, collapse the wavefunction if called for
!               if (decohere) then
!                    open(unit=37,file='decohere',action='read')
!                    read(37,*), collapse
!                    close(37)
!                    if (collapse) then
!                         do i=1,mat
!                              ctin(i) = C2(state,i)
!                         end do
!                         
!                         write(37,*), 'f'
!                         close(37)
!                    end if
!               end if
!          end if
!WRITE          write(18,*) ctin(1:mat)
!print*, itimes, ctin, "tttt"
                aom_basket%ctin(1:mat) = ctin(1:mat)
!print*, itimes, aom_basket%ctin, "tttt"

!          close(unit=24)
!          close(unit=12)
!          close(unit=18)

!     end do !loop over t
      endif

!CHANGE_AC
!STORE THE COORDINATES
        if (itimes.ge.1) then
           do i=1,c
              do coord=1,3
                 aom_basket%rm2(i,coord) = aom_basket%rm1(i,coord)
                 aom_basket%rm1(i,coord) = positions_plus(i,coord)
              enddo
           enddo
           do i=1,mol
!print*, energies_plus(i), "keke"
              aom_basket%energies(i) = energies_plus(i)
           enddo
        else if (itimes.eq.0) then
           do i=1,c
              do coord=1,3
                 aom_basket%rm1(i,coord) = positions_plus(i,coord)
              enddo
           enddo
           do i=1,mol
!print*, energies_plus(i), "keke"
              aom_basket%energies(i) = energies_plus(i)
           enddo
        endif

!PREPARE THE PRINT
       !aom_basket%hamiltonian(i,i) = siteenergies(t+1,i)
!> H2 is in fact -i*Hamiltonian
!> TO_CHECK 16-03-19
       aom_basket%hamiltonian = - AIMAG(H2save)
       aom_basket%nacv        = nacv

!END CHANGE_AC
!PRINT     print *, "All done!",itimes
     !print *, "All done!",itimes

!     deallocate(pcoefficients, pfiles, pcoeffs, pcoeffsrun, scoeffs)
     deallocate(pcoeffs, pcoeffsrun, scoeffs)
     deallocate(ctin, ctout, ctotal)
     !deallocate(siteenergies, crds, coms, bonds, NN, At, Atm, Atp, Bt, Btm, Btp, atomlist, connlist1, connlist2)
     deallocate(siteenergies)
     deallocate(crds)
     deallocate(coms)
     deallocate(NN)
     deallocate(At)
     deallocate(Atm)
     deallocate(Atp)
     deallocate(Bt)
     deallocate(Btm)
     deallocate(Btp)
     deallocate(atomlist)
     deallocate(connlist1)
!print*, connlist2
     deallocate(connlist2)
     deallocate(H1, H2, H12, Hdiff, Hsave, X1, X2, X12, Ktot, K2, K3, K4)
     deallocate(Stt, Stptp, Stmtm, Sttp, Sttm, Stpt, Stmt)
     deallocate(Sttdiff, Stptpdiff, Stmtmdiff, Sttpdiff, Sttmdiff, nacv, dist)

    END SUBROUTINE construct_ham_X_SH_full_degen_classicalpath

END MODULE jacob_classicalpath
