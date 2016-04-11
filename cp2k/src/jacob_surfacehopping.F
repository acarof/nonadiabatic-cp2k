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

MODULE jacob_surfacehopping 

     USE kinds,                           ONLY: default_string_length,&
                                               dp


     USE jacob_overlapfinal,                 ONLY:  connect_list2, &
                                             calc_sab, &
                                             sab_screen_new, &
                                             calc_dRSab

     USE jacob_utils,                        ONLY: makeX, &
                                             adiabaticity, &
                                             surfacehop

     USE jacob_types,                     ONLY: aom_basket_type
     USE mixed_energy_types,              ONLY: mixed_force_type


#include "./base/base_uses.f90"


     IMPLICIT NONE

        PRIVATE
 
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='jacob_construct_ham_X_SH_full_degen_classicalpath'   

        PUBLIC :: construct_ham_X_SH

CONTAINS


        SUBROUTINE construct_ham_X_SH(nadiab, &
                   ext_atoms, norbitals, ext_cutoff, ext_scaling, &
                   ext_init, ext_cbar, ext_minao, ext_timestep, ext_el_steps, &
                   positions_plus, energies_plus, aom_basket, ext_elements, &
                   ext_pres, ext_scoeff, ext_pcoeff, itimes, collapse, is_hop, &
                   rescale_keyword, ext_dR, forces_adiab, global_forces)

        integer          :: nadiab,  norbitals, ext_init
        integer          :: ext_atoms, ext_el_steps
        REAL(KIND=dp) :: ext_cutoff, ext_scaling, ext_cbar, ext_minao, ext_timestep
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_plus
        REAL(kind=dp), DIMENSION(:),  POINTER    :: energies_plus
        TYPE(aom_basket_type)                    :: aom_basket
        integer, dimension(:),  INTENT(IN)        :: ext_elements
        LOGICAL, dimension(:), INTENT(IN)       :: ext_pres
        REAL(KIND=dp), dimension(:), ALLOCATABLE, INTENT(IN)  :: ext_scoeff,ext_pcoeff
        INTEGER, POINTER                                       :: itimes
        LOGICAL, INTENT(OUT)                         :: is_hop
        LOGICAL, INTENT(IN)                          :: collapse           
        CHARACTER(LEN=default_string_length)         :: rescale_keyword
        REAL(KIND=dp), INTENT(IN), OPTIONAL          :: ext_dR
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: forces_adiab
        TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces


        CHARACTER(len=*), PARAMETER :: routineN ='construct_ham_X_SH_full_degen_classicalpath', &
        routineP = moduleN//':'//routineN

     integer :: coord
     integer :: mol, i, j, k, l, m, d, atoms, c, t, t2, init, inc, state, mat, orb
     integer, dimension(:), allocatable :: atomlist, elements
     integer, dimension(:,:), allocatable :: NN, connlist1, connlist2, aneighbour, bneighbour
     integer, dimension(:,:,:), allocatable :: neighbourlist
     REAL(KIND=dp) :: cutoff, scaling, deltatMD, deltatE
     REAL(KIND=dp) :: S
     character(len=1) :: zgemm_type
     REAL(KIND=dp) :: ctmag, bohrang, mev, kcal, cbar, minao
     REAL(KIND=dp), dimension(:), allocatable :: scoeffs
     REAL(KIND=dp), dimension(:,:,:,:), allocatable :: nacv
     REAL(KIND=dp), dimension(:,:,:), allocatable :: crds, pcoeffsrun
     REAL(KIND=dp), dimension(:,:), allocatable :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm, Imat, pcoeffs
     REAL(KIND=dp), dimension(:,:), allocatable :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     REAL(KIND=dp), dimension(:,:), allocatable :: dist
     COMPLEX(KIND=dp) :: alpha, beta, ctcomp
     COMPLEX(KIND=dp), dimension(:), allocatable :: ctin, ctout, ctotal
     COMPLEX(KIND=dp), dimension(:,:), allocatable :: H1, H2, H12, Hdiff, Hsave, X1, X2, X12
     COMPLEX(KIND=dp), dimension(:,:), allocatable :: Ktot, K2, K3, K4
     REAL(KIND=dp), dimension(:,:), allocatable :: At, Atm, Atp, Bt, Btm, Btp, coms, siteenergies

     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pcoefftot, scoefftot
     INTEGER, DIMENSION(:), ALLOCATABLE :: elementtot
     LOGICAL, DIMENSION(:), ALLOCATABLE :: atomlist_mask


     REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: F
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: C1, C2
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: E1, E2
     !LOGICAL                                      :: collapse, decohere
     LOGICAL                                      :: decohere
     INTEGER                                      :: ext_state    
     REAL(KIND=dp)                                :: dR 
     INTEGER                                      :: imol, iatom, natom

     !We begin by reading in the input parameters
     !There are currently two modes this program can run in
     !In Classical Path Approximation mode, the program reads in a pre-computed trajectory
     !In Surface Hopping mode, the program runs alongside non-adiabatic MD
     !This variable is 'S' for Surface Hopping and 'C' for Classical Path
     !name of coordinate file(s)
     !technically only mdcrd matters for progtype == C
     !'type' of file so that program can interpret it
     !currently supported:
     ! - "AMBER" for AMBER .mdcrd file
     ! - "NAMD" for NAMD .dcd file
     ! - "XYZ" for .xyz file
     ! - "SH" is format printed out by NAMD surfhop code and should be default for progtype == S
!> Number of sites
     mol = nadiab ! CHANGE_AC
!> Number of atoms per site
        atoms = ext_atoms ! CHANGE_AC
     !for Classical Path Approximation, this is the number of timesteps in pre-computed trajectory
     !otherwise, 2
!> Number of orbitals per atom
     orb = norbitals ! CHANGE_AC
!> cutoff to define the neighborough list
!> In angstrom
     cutoff = ext_cutoff ! CHANGE_AC
     !scaling = 0.5D0*0.06519D0
!> Scaling factor Hab = scaling*Sab_bar
     scaling = ext_scaling ! CHANGE_AC
     !init = 1
!> Starting site number
     init = ext_init ! CHANGE_AC
     !cbar = 0.5082
!> cbar factor : Sab = cbar*Sab_bar
     cbar = ext_cbar ! CHANGE_AC
     !minao is an atomic overlap cutoff
     !a value of 1E-17 means all atoms are included
     !minao = 1.0D-17
     minao = ext_minao ! CHANGE_AC

!> dR is the finite difference to calculate the NACV
     if (PRESENT(ext_dR)) then 
        dR = ext_dR
     else
        dR = 0.0001D0*1.889725989D0
     endif

     !With surface-hopping, the program needs information about where in the simulation it is
     !i.e. if this is the first timestep, or if certain once-in-a-while actions need to be done
     !The file "inputget" contains this information

!> Size of the vectorial space where the wavefunctions are projected
     mat = orb*mol
!> Numbers of active atoms 
     c = mol*atoms
!> Numbers of atoms
     natom = SIZE(global_forces(1)%forces,2)

!> ALLOCATE aom_basket at timestep = 0
    if (itimes.eq.0) then
       ALLOCATE(aom_basket%population(mat))
       ALLOCATE(aom_basket%ctin(mat))
       ALLOCATE(aom_basket%hamiltonian(mat,mat))
       ALLOCATE(aom_basket%connlist1(4,atoms))
       ALLOCATE(aom_basket%neighbourlist(mol,3,0:atoms))
       ALLOCATE(aom_basket%Stt(mat,mat))
       ALLOCATE(aom_basket%Stmt(mat,mat))
       ALLOCATE(aom_basket%Sttm(mat,mat))
       ALLOCATE(aom_basket%Stmtm(mat,mat))
       ALLOCATE(aom_basket%rm1(c,3))
       ALLOCATE(aom_basket%rm2(c,3))
       ALLOCATE(aom_basket%energies(mol))
       ALLOCATE(aom_basket%bneighbour(3,0:atoms))
       ALLOCATE(aom_basket%nacv(3,mat,mat,mat*atoms))
    endif

     !calculate timestep from femtoseconds in atomic units
     !molecular timestep is generally 0.5fs to improve energy conservation
!     deltatMD = (5.0D-16)/(2.4188843265D-17)
!> Find deltaMD, nuclear timestep, in atomic units
     deltatMD = (ext_timestep*1.0D-15)/(2.4188843265D-17) !CHANGE_AC
     !electronic timestep is a fraction of the molecular timestep
!     deltatE = deltatMD/5.0D0
!> Find deltaE, electronic timestep, in atomic units
     deltatE = deltatMD/ext_el_steps !CHANGE_AC

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


     !ct are the electron coefficients
     !ctin goes into the subroutine ZGEMV; ctout comes out of it (so essentially ctin = c(t), ctout = c(t+1))
     allocate(ctin(mat))
     allocate(ctout(mat))
     !ctotal is |c(t)|^2: this is the output for diabatic electronic state
     allocate(ctotal(mat))

     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(Stt(mat,mat))
     allocate(Stpt(mat,mat))
     allocate(Sttp(mat,mat))
     allocate(Stptp(mat,mat))
     allocate(Stmt(mat,mat))
     allocate(Sttm(mat,mat))
     allocate(Stmtm(mat,mat))
     allocate(Sttdiff(mat,mat))
     allocate(Stmtmdiff(mat,mat))
     allocate(Stptpdiff(mat,mat))
     allocate(Sttpdiff(mat,mat))


     !c tells you how many atoms you will have per timestep
     !allocate(crds(3,c,times+1))
     allocate(crds(3,c,3))

     allocate(Sttpdiff(mat,mat))
     allocate(Sttmdiff(mat,mat))

     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H1 etc are always -i*H.
     !This is why they're complex variables

     allocate(H1(mat,mat))
     allocate(H2(mat,mat))
     allocate(H12(mat,mat))
     allocate(Hdiff(mat,mat))
     allocate(Hsave(mat,mat))
     allocate(X1(mat,mat))
     allocate(X2(mat,mat))
     allocate(X12(mat,mat))

     !finally the K matrices, which are the terms in the Runge-Kutta integration
     !note K1 missing, because it's just a constant times X
     !Ktot is the final constant by which c(t) is multiplied
     allocate(Ktot(mat,mat))
     allocate(K2(mat,mat))
     allocate(K3(mat,mat))
     allocate(K4(mat,mat))

     !Imat is an identity matrix of the dimensions of the other square matrices
     allocate(Imat(mat,mat))
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


     !non-adiabatic coupling vectors
     allocate(nacv(3,mat,mat,mat*atoms))
     !diagonals of the Hamiltonian
     allocate(siteenergies(3,mat))

     !atomic coefficients for the overlap routine
     allocate(pcoeffs(atoms,orb))
     !pcoeffsrun exists because the p-coefficients need to be renormalised for each structure
     allocate(pcoeffsrun(atoms,mat,3))
     allocate(scoeffs(atoms))

     allocate(pcoefftot(c))
     allocate(scoefftot(c))
     allocate(coms(3,mol))
     allocate(NN(mat,mat))
     allocate(dist(mol,mol))

     !allocate the arrays that will contain atomic coordinates for the overlap routine
     !first column is now atomic number of atom (e.g. 6 for carbon)
     !Naming convention: At = A(t), Atp = A(t+1), Atm = A(t-1), etc
     allocate(At(4,atoms))
     allocate(Atp(4,atoms))
     allocate(Atm(4,atoms))
     allocate(Bt(4,atoms))
     allocate(Btp(4,atoms))
     allocate(Btm(4,atoms))

     !array of diagonal atomic forces
     !allocate(F(3,mol*atoms,mol))
     allocate(F(3,natom,mol))

     !elements contains the atomic number of each atom (this goes into At, Bt etc)
     !it is user-defined input
     allocate(elements(atoms))
     allocate(elementtot(c))
     elementtot = PACK(ext_elements, ext_pres)
     allocate(atomlist_mask(c))
     allocate(atomlist(COUNT(atomlist_mask(1:atoms))))
    if (itimes.eq.0) then
       allocate(aom_basket%connlist2(4,size(atomlist)))
    endif

     !these arrays are data for the non-adiabatic coupling calculation
     !they help efficiency by telling the program which atoms are connected
     !and thus which pairs to carry out the calculation for
     allocate(connlist1(4,atoms))
     allocate(connlist2(4,size(atomlist)))
     allocate(aneighbour(3,0:atoms))
     allocate(bneighbour(3,0:atoms))
     allocate(neighbourlist(mol,3,0:atoms))


     !the stateno file contains a numeric representation of the adiabatic electronic state
        if (itimes.gt.1) then
           state = aom_basket%state
        else
           state = 1
        endif

     !if the program has been previously called, a bunch of data will have been saved

     if (itimes.ge.3) then
          ctin(:) = aom_basket%ctin(:)
     !if not, it may need to be initialised/calculated for the first time
     else if (itimes.eq.2) then
          !initialise ct such that ct(init) = 1.0, otherwise 0.0 (i.e. electron localised on molecule 'init')
          do i=1,mat
               ctin(i) = CMPLX(0.0D0, 0.0D0,dp)
               ctotal(i) = CMPLX(0.0D0, 0.0D0,dp)
          end do
          ctin(init*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
          ctotal(init*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
                aom_basket%population(1:mat) = abs(ctotal(1:mat)) 
     else
          !initialise ct such that ct(init) = 1.0, otherwise 0.0 (i.e. electron localised on molecule 'init')
          do i=1,mat
               ctin(i) = CMPLX(0.0D0, 0.0D0,dp)
               ctotal(i) = CMPLX(0.0D0, 0.0D0,dp)
          end do
          ctin(init*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
          ctotal(init*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
                aom_basket%population(1:mat) = abs(ctotal(1:mat)) 
     end if


!READ THE COORDINATES
        do coord=1,3
              if (itimes.ge.2) then
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
!END READ THE COORDINATES

!READ THE COEFFICIENTS
     !now we also read in the atomic coefficients
        pcoefftot = ext_pcoeff
        scoefftot = ext_scoeff
     do i=1,atoms
          do j=1,orb
                       pcoeffs(i,j) = pcoefftot( i + (j-1)*orb)
          end do
                scoeffs(i) = scoefftot(i)
     end do
     
!PRINT     print *, "Done! Calculating centres of mass..."

!CALCULATE THE CENTER OF MASS
     !calculate CoMs at first timestep
     !coms(1 to 3,i) are the x, y and z coordinates of CoM of i
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
!END CALCULATE THE CENTER OF MASS

     !construct neigbours list: 1 if close enough, 0 if not
!PRINT     print *, "Constructing nearest neighbours list..."
!CONSTRUCT NEIGHBOURS LIST
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
!END CONSTRUCT NEIGHBOURS LIST
!PRINT     print *, "Nearest neighbours found."


     do i=1,atoms
                elements(i) = elementtot(i)
     end do

     atomlist = PACK( (/ (i, i=1,atoms) /), atomlist_mask(1:atoms) )
     atomlist_mask =.NOT.((scoefftot.eq.0).AND.(pcoefftot.eq.0))
        

!CONSTRUCT THE CONNLIST
     !we fill out the connlist arrays before calculating the nacv
     do i=1,atoms
          do j=1,3
!CHANGE_AC: already in atomic units
               At(j,i) = crds(j,i,3) !CHANGE_AC : 3 is the current timestep
          end do
     end do
     !if (.not. inrun) then
     if (itimes.eq.0) then
          call connect_list2(At, atomlist, connlist1, connlist2)
          do i=1,atoms
                        aom_basket%connlist1(1:4,i) = connlist1(1:4,i)
          end do
          do i=1,size(atomlist)
                        aom_basket%connlist2(1:4,i) = connlist2(1:4,i)
          end do
     else if (itimes.gt.0) then
                        connlist1(1:4,:) = aom_basket%connlist1(1:4,:)
                        connlist2(1:4,:) = aom_basket%connlist2(1:4,:)
     end if
!END CONSTRUCT THE CONNLIST

     !used to give correct commands to ZGEMX subroutines later
     zgemm_type = 'N'

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
!> CHANGE_AC : already coordonates in atomic units
                         At(d,j) = crds(d-1,(i-1)*atoms+j,t)
                         Atp(d,j) = crds(d-1,(i-1)*atoms+j,t+1)
                         Atm(d,j) = crds(d-1,(i-1)*atoms+j,t-1)
                    end do
               end do
               !calculate the site energy for this molecule
               !at this timestep
                        siteenergies(t,i)   = aom_basket%energies(i) ! CHANGE_AC
                        siteenergies(t+1,i) = energies_plus(i) ! CHANGE_AC

               !then one timestep ahead: also calculate diagonal forces at this timestep
!               call sitecalc(Atp, bonds, atoms, siteenergies(t+1,i), F(:,:,i), mol, i)

               !then for each other molecule
               do k=1,mol
                    !create the array of its atomic coordinates at both t and t+dt
                    do l=1,atoms
                         Bt(1,l) = dble(elements(l))
                         Btp(1,l) = dble(elements(l))
                         Btm(1,l) = dble(elements(l))
                         do d=2,4
!> CHANGE_AC : already coordonates in atomic units
                              Bt(d,l) =  crds(d-1,(k-1)*atoms+l,t)
                              Btp(d,l) = crds(d-1,(k-1)*atoms+l,t+1)
                              Btm(d,l) = crds(d-1,(k-1)*atoms+l,t-1)
                         end do
                    end do
                    if (i == 1) then
                         do d=1,orb
                              do l=1,atoms
                                   pcoeffsrun(l,orb*(k-1)+d,1) = pcoeffs(l,d)
                              end do
                              call calc_sab(Btm, connlist2, pcoeffsrun(:,orb*(k-1)+d,1))
                              do l=1,atoms
                                   pcoeffsrun(l,orb*(k-1)+d,2) = pcoeffs(l,d)
                              end do
                              call calc_sab(Bt, connlist2, pcoeffsrun(:,orb*(k-1)+d,2))
                              do l=1,atoms
                                   pcoeffsrun(l,orb*(k-1)+d,3) = pcoeffs(l,d)
                              end do
                              call calc_sab(Btp, connlist2, pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
                         if (itimes.eq.2) then
                              do d=1,orb
                                   do m=1,orb
                                        call sab_screen_new(At, atomlist, connlist1, pcoeffsrun(:,orb*(i-1)+d,2),&
                                        & atomlist, Bt, atomlist, connlist1, pcoeffsrun(:,orb*(k-1)+m,2),&
                                        & aneighbour(:,1:atoms), bneighbour(:,1:atoms), minao)
                                   end do
                              end do
                              neighbourlist(k,1:3,0) = [0, 0, 0]
                              neighbourlist(k,:,1:atoms) = bneighbour(:,1:atoms)
                              if (k == mol) then
                                   do m=1,mol
                                        do j=1,atoms
                                                                aom_basket%neighbourlist(m,:,j) = & 
                                                                        neighbourlist(m,:,j)
                                        end do
                                   end do
                              end if
                         else
                              do m=1,mol
                                   neighbourlist(m,1:3,0) = [0, 0, 0]
                                   do j=1,atoms
                                                                neighbourlist(m,:,j)=&
                                                                aom_basket%neighbourlist(m,:,j)
                                   end do
                              end do
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
                                        H1(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(t,k)/kcal), dp)
                                        !H2(orb*(i-1)+d,orb*(k-1)+d) = dcmplx(0.0D0, (-siteenergies(t+1,k)/kcal))
                                        H2(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(t+1,k)/kcal),dp)
                                        call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Btp,&
                                        & connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                        Sttp(orb*(i-1)+d,orb*(k-1)+m) = S
                                        call calc_sab(Atp, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3), S, Bt,&
                                        & connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                        Stpt(orb*(i-1)+d,orb*(k-1)+m) = S
                                        call calc_dRSab(At, atomlist, connlist1, &
                                        &pcoeffsrun(1:atoms,orb*(i-1)+d,2), scoeffs, &
                                        &nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,1:atoms), &
                                        &nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,atoms+1:2*atoms), &
                                        &Bt, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2), scoeffs, &
                                        !&neighbourlist(i,:,:), neighbourlist(k,:,:))
                                        &neighbourlist(i,:,:), neighbourlist(k,:,:), dR)
                                   end do
                              end do

                              if (itimes.eq.2) then
                                   do d=1,orb
                                        do m=1,orb
                                             call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2),&
                                             & S, Bt, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                             Stt(orb*(i-1)+d,orb*(k-1)+m) = S
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
                                             H2(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)
                                             call calc_dRSab(At, atomlist, connlist1, &
                                             &pcoeffsrun(1:atoms,orb*(i-1)+d,2), scoeffs, &
                                             &nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(k-1)*atoms+1:2*(k-1)*atoms), &
                                             &nacv(1:3,orb*(k-1)+m,orb*(i-1)+d,(i-1)*atoms+1:2*(i-1)*atoms), &
                                             &Bt, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2),&
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
!GET INFORMATION FROM PREVIOUS TIMESTEP AND STORE SOME OF THEM
          if (itimes.gt.2) then
               do i=1,mat
                    do j=1,mat
                                        Stt(j,i) = aom_basket%Stt(j,i)
                         if (i .ne. j) then
                              H1(j,i) = cmplx(0.0D0, -1.0D0*scaling*Stt(j,i)/cbar,dp)
                         end if
                                        Sttm(j,i)  = aom_basket%Sttm(j,i)
                                        Stmt(j,i)  = aom_basket%Stmt(j,i)
                                        Stmtm(j,i) = aom_basket%Stmtm(j,i)
                    end do
               end do
         else if (itimes.ge.2) then !CHANGE_AC
               do i=1,mat
                  do j=1,mat
                                aom_basket%Stt(j,i)   = Stptp(j,i)
                                aom_basket%Sttm(j,i)  = Stpt(j,i)
                                aom_basket%Stmt(j,i)  = Sttp(j,i)
                                aom_basket%Stmtm(j,i) = Stt(j,i)
                  end do
               end do
         endif !CHANGE_AC
!END GET INFORMATION FROM PREVIOUS TIMESTEP AND STORE SOME OF THEM

          nacv = cbar*nacv


!PRINT          print *, "All pairs calculated."

!PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
          !It's sometimes necessary to just call this program to calculate forces, without propagating the electron.
!          if (input /= 'n') then
               !MD is done on a molecular timestep of femtoseconds: however electron dynamics happen on attoseconds
               !so for each molecular timestep, the program calculates for multiple electronic timesteps
               !during this process, the H and S matrices are found by linear interpolation between the original S1 and S2
               Hdiff = (deltatE/deltatMD)*(H2 - H1)
               Sttdiff = (deltatE/deltatMD)*(Stptp - Stt)
               Sttpdiff = (deltatE/deltatMD)*(Sttp - Stmt)
               Sttmdiff = (deltatE/deltatMD)*(Stpt - Sttm)
               Stptpdiff = Sttdiff
               Stmtmdiff = (deltatE/deltatMD)*(Stt - Stmtm)
               Hsave = H1
!PRINT               print *, "Propagating the electron for one molecular timestep..."
               do t2=1,int(deltatMD/deltatE)

                    call makeX(mat, H1, Stt, Sttp, Sttm, Stptp, Stmtm, X1, deltatMD)
                    !create timestepped equivalents
                    H12 = H1 + 0.5D0*Hdiff
                    Stt = Stt + 0.5D0*Sttdiff
                    Sttp = Sttp + 0.5D0*Sttpdiff
                    Sttm = Sttm + 0.5D0*Sttmdiff
                    Stptp = Stptp + 0.5D0*Stptpdiff
                    Stmtm = Stmtm + 0.5D0*Stmtmdiff
                    call makeX(mat, H12, Stt, Sttp, Sttm, Stptp, Stmtm, X12, deltatMD)
     
                    H2 = H12 + 0.5D0*Hdiff
                    Stt = Stt + 0.5D0*Sttdiff
                    Sttp = Sttp + 0.5D0*Sttpdiff
                    Sttm = Sttm + 0.5D0*Sttmdiff
                    Stptp = Stptp + 0.5D0*Stptpdiff
                    Stmtm = Stmtm + 0.5D0*Stmtmdiff
                    call makeX(mat, H2, Stt, Sttp, Sttm, Stptp, Stmtm, X2, deltatMD)
     
                    K2 = cmplx(0.0D0, 0.0D0,dp)
                    K3 = cmplx(0.0D0, 0.0D0,dp)
                    K4 = cmplx(0.0D0, 0.0D0,dp)
          
                    !now that X has been found, use 4th order Runge-Kutta numerical integration method
                    !aiming to find operator such that c(t+1) = (1 + operator)*c
          
                    !there are four terms in X: each is an increasing multiple. so first term is proportional to X(t),
                    !second term is as X(t+deltat/2)*X(t), third term is as X(t+deltat/2)*X(t+deltat/2)*X(t),
                    !fourth term is as X(t+deltat)*X(t+deltat/2)*X(t+deltat/2)*X(t)
                    !this means each term is some X matrix times the previous term,
                    !so only three more calls to ZGEMM are needed
                    !the K matrices are the second, third and fourth Runge-Kutta terms
     
                    alpha = cmplx(deltatE/2.0D0, 0.0D0,dp)
                    beta = cmplx(0.0D0, 0.0D0,dp)
     
                    !second term
                    call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X12, mat, X1, &
                    &mat, beta, K2, mat)
                    K2 = K2 + X12
                    
                    !third term
                    call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X12, mat, X2, &
                    &mat, beta, K3, mat)
                    K3 = K3 + X12
     
                    !fourth term
                    alpha = cmplx(deltatE, 0.0D0,dp)
                    call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X2, mat, K3, &
                    &mat, beta, K4, mat)
                    K4 = K4 + X2
     
                    !put everything together
                    Ktot = Imat + (deltatE/6.0D0)*(X1 + 2*K2 + 2*K3 + K4)
     
                    !timestep electron propagation and print out
                    !alpha = dcmplx(1.0D0, 0.0D0)
                    alpha = cmplx(1.0D0, 0.0D0,dp)
                    inc = 1
                    call ZGEMV(zgemm_type, mat, mat, alpha, Ktot, mat, ctin, inc, beta, ctout, inc)
     
                    ctmag = 0.0D0
                    ctcomp = cmplx(0.0D0, 0.0D0,dp)
                    do i=1,mat
                         ctotal(i) = conjg(ctout(i))*ctout(i)
                         ctcomp = ctcomp + ctotal(i)
                    end do
                    ctmag = abs(ctcomp)
                    !rescale to conserve norm = 1
                    do i=1,mat
                         ctout(i) = cmplx((1.0D0/sqrt(ctmag))*REAL(ctout(i)), (1.0D0/sqrt(ctmag))*AIMAG(ctout(i)),dp)
                    end do
                    do i=1,mat
                         ctotal(i) = conjg(ctout(i))*ctout(i)
                    end do
                    if (t2 == int(deltatMD/deltatE)) then
                                        aom_basket%population(1:mat) =&
                                            abs(ctotal(1:mat)) 
                    end if
                    !update
                    H1 = H2
                    ctin = ctout
                    ctout = cmplx(0.0D0, 0.0D0,dp)
               end do
!PRINT               print *, "Done."     
!END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
               H1 = Hsave

!DO THE HOPPING TEST AND CALCULATE THE FORCES
!          if (progtype == 'C') then
!PRINT               print *, "Molecular timestep", t, "done."
!          else
!PRINT               print *, "Calculating surface hopping probabilities."

               !The adiabaticity subroutine diagonalises the diabatic Hamiltonian
               !It gives us back the eigenvalues and eigenvectors
!               call adiabaticity(-imagpart(H1), mat, C1, E1, .false.)
!               call adiabaticity(-imagpart(H2), mat, C2, E2, .true.)
               call adiabaticity(-AIMAG(H1), mat, C1, E1, .false.)
               call adiabaticity(-AIMAG(H2), mat, C2, E2, .true.)

               !The eigenvectors are needed to calculate the surface hopping probabilities
               !The surfacehop subroutine calculates these probabilites and decides on the
               !new (possibly unchanged) surface for simulation
!               if (input /= 'n') then
                    !call surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, ctin, C1, C2, deltatMD, state, decohere)
                    call surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, ctin, C1, C2, deltatMD, state, decohere, &
                                    ext_state, is_hop)
                    aom_basket%state = ext_state
!               end if

               !The forces subroutine then calculates the adiabatic forces and prints them out
               !F = F/(kcal*bohrang)
!               call forces(state, -imagpart(H2), C2, nacv, NN, atoms, scaling, cbar, F, Stt)
               !call forces(state, -AIMAG(H2), C2, nacv, NN, atoms, scaling, cbar, F, Stt)
!CALCULATE THE FORCES IN THE NEW ADIABATIC STATE!
               !call forces(ext_state, -AIMAG(H2), C2, nacv, NN, atoms, scaling, cbar, F, Stt)
               DO imol=1,mol
                  DO iatom=1, natom
                     F(:,iatom,imol) =  global_forces(imol)%forces(:,iatom)
                  ENDDO
               ENDDO
               call forces(ext_state, -AIMAG(H2), C2, nacv, NN, atoms, scaling, cbar, F, Stt, forces_adiab, &
                           ext_pres)
               !call forces(state, -AIMAG(H2), C2, nacv, NN, atoms, scaling, cbar, F, Stt, ext_forces)

               !Finally, collapse the wavefunction if called for
               !if (decohere) then
               if (is_hop) aom_basket%futur_decohere = .TRUE.
!> Check if the collapse key word is TRUE in the input
               if (collapse) then
!                    open(unit=37,file='decohere',action='read')
!                    read(37,*), collapse
!                    close(37)

                    !if (collapse) then
!> Check if the surface_hop subroutine asks us to do the decoherence
                    if (aom_basket%futur_decohere) then 
                       if (decohere) then
                             do i=1,mat
                                 ctin(i) = C2(state,i)
                             end do
                             aom_basket%futur_decohere = .FALSE. 
!                            write(37,*), 'f'
!                            close(37)
                       end if
                    endif
               end if
!          end if
!END DO THE HOPPING TEST AND CALCULATE THE FORCES


!     end do !loop over t
      endif

!STORE THE COORDINATES
        if (itimes.ge.1) then
           do i=1,c
              do coord=1,3
                 aom_basket%rm2(i,coord) = aom_basket%rm1(i,coord)
                 aom_basket%rm1(i,coord) = positions_plus(i,coord)
              enddo
           enddo
           do i=1,mol
              aom_basket%energies(i) = energies_plus(i)
           enddo
        else if (itimes.eq.0) then
           do i=1,c
              do coord=1,3
                 aom_basket%rm1(i,coord) = positions_plus(i,coord)
              enddo
           enddo
           do i=1,mol
              aom_basket%energies(i) = energies_plus(i)
           enddo
        endif
!END STORE THE COORDINATES

!PREPARE THE PRINT
       aom_basket%ctin(1:mat) = ctin(1:mat)
       aom_basket%hamiltonian = AIMAG(H2)
       aom_basket%nacv        = nacv
!END PREPARE THE PRINT

!PRINT     print *, "All done!",itimes

!DEALLOCATE
     deallocate(pcoeffs, pcoeffsrun, scoeffs)
     deallocate(ctin, ctout, ctotal)
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
     deallocate(connlist2)
     deallocate(H1, H2, H12, Hdiff, Hsave, X1, X2, X12, Ktot, K2, K3, K4)
     deallocate(Stt, Stptp, Stmtm, Sttp, Sttm, Stpt, Stmt)
     deallocate(Sttdiff, Stptpdiff, Stmtmdiff, Sttpdiff, Sttmdiff, nacv, dist)
!END DEALLOCATE

    END SUBROUTINE construct_ham_X_SH

END MODULE jacob_surfacehopping
