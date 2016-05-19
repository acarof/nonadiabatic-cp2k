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

     USE jacob_utils,                        ONLY: makeX,  &
                                             adiabaticity, &
                                             surfacehop,   &
                                             calcforces,       &
                                             read_listcoeff

     USE jacob_types,                     ONLY: aom_basket_type
     USE mixed_energy_types,              ONLY: mixed_force_type
     USE input_section_types,             ONLY: section_vals_get,&
                                                section_vals_get_subs_vals,&
                                                section_vals_type,&
                                                section_vals_list_get,&
                                                section_vals_val_get
     USE cp_linked_list_val,              ONLY: &
         cp_sll_val_create, cp_sll_val_dealloc, cp_sll_val_get_el_at, &
         cp_sll_val_get_length, cp_sll_val_get_rest, cp_sll_val_insert_el_at, &
         cp_sll_val_next, cp_sll_val_p_type, cp_sll_val_rm_el_at, &
         cp_sll_val_set_el_at, cp_sll_val_type
     USE jacob_print,                     ONLY: adiab_print_coeff, &
                                                adiab_print_pop,   &
                                                adiab_print_matrix_density, &
                                                adiab_print_exact_forces



#include "./base/base_uses.f90"


     IMPLICIT NONE

        PRIVATE
 
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='jacob_construct_ham_X_SH_full_degen_classicalpath'   

        PUBLIC :: construct_ham_X_SH

CONTAINS

        SUBROUTINE construct_ham_X_SH(nadiab, aom_atoms, norbitals, MDtimestep, &
                                      positions_plus, energies_plus, kinetic_energies, &
                                      aom_basket, itimes, &
                                      is_hop, forces_adiab, global_forces, coupling_section, &
                                      adiab_section, my_bo_keyword, rescale_factor)


        INTEGER, INTENT(IN)                            :: nadiab,  norbitals 
        INTEGER, INTENT(IN)                            :: aom_atoms
        REAL(KIND=dp), INTENT(OUT)                     :: MDtimestep
        REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_plus
        REAL(kind=dp), DIMENSION(:),  POINTER          :: energies_plus
        REAL(kind=dp), DIMENSION(:),  POINTER          :: kinetic_energies
        TYPE(aom_basket_type), INTENT(INOUT)           :: aom_basket
        INTEGER, POINTER                               :: itimes
        LOGICAL, INTENT(OUT)                           :: is_hop
        REAL(kind=dp), DIMENSION(:,:),POINTER          :: forces_adiab
        TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces
        TYPE(section_vals_type), POINTER               :: coupling_section
        LOGICAL, INTENT(OUT)                           :: my_bo_keyword
        TYPE(section_vals_type), POINTER               :: adiab_section
        REAL(KIND=dp), INTENT(OUT)                     :: rescale_factor


        CHARACTER(len=*), PARAMETER :: routineN ='construct_ham_X_SH_full_degen_classicalpath', &
        routineP = moduleN//':'//routineN

     integer :: coord, t2
     integer :: mol, i, j, k, l, m, d, atoms, c,  first_diabat, inc, old_state, mat, orb
     integer, dimension(:), allocatable :: atomlist, elements
     integer, dimension(:,:), allocatable :: NN, connlist1, connlist2, aneighbour, bneighbour
     integer, dimension(:,:,:,:), allocatable :: neighbourlist
     REAL(KIND=dp) :: cutoff, scaling, deltatMD, deltatE
     REAL(KIND=dp) :: S
     character(len=1) :: zgemm_type
     REAL(KIND=dp) :: ctmag,  cbar, minao
     REAL(KIND=dp), dimension(:), allocatable :: scoeffs
     REAL(KIND=dp), dimension(:,:,:,:), allocatable :: nacv
     REAL(KIND=dp), dimension(:,:,:), allocatable :: crds, pcoeffsrun
     REAL(KIND=dp), dimension(:,:), allocatable :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm, Imat, pcoeffs
     REAL(KIND=dp), dimension(:,:), allocatable :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     REAL(KIND=dp), dimension(:,:), allocatable :: dist
     COMPLEX(KIND=dp) :: alpha, beta, ctcomp
     COMPLEX(KIND=dp), dimension(:), allocatable :: ctin, ctout, ctotal
     COMPLEX(KIND=dp), dimension(:,:), allocatable :: H1, H2, H12, Hdiff, Hsave, H2save, X1, X2, X12
     COMPLEX(KIND=dp), dimension(:,:), allocatable :: Ktot, K2, K3, K4
     REAL(KIND=dp), dimension(:,:), allocatable :: At, Atm, Atp, Bt, Btm, Btp, coms, siteenergies

     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pcoefftot, scoefftot
     LOGICAL, DIMENSION(:), ALLOCATABLE :: atomlist_mask


     REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: F
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: C1, C2
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: E1, E2
     LOGICAL                                      :: decohere, collapse, analytics
     INTEGER                                      :: new_state , ELsteps   
     REAL(KIND=dp)                                :: dR, times, delta_energy 
     INTEGER                                      :: imol, iatom, natom
     integer, dimension(:), ALLOCATABLE           :: ext_elements
     LOGICAL, dimension(:), ALLOCATABLE           :: ext_pres
     TYPE(cp_sll_val_type), POINTER               :: listcoeff
     INTEGER                                      :: first_adiabat
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: exact_forces
     CHARACTER(LEN=default_string_length)         :: rescale_keyword


!> cutoff to define the neighborough list
        CALL section_vals_val_get(coupling_section,"CUT_OFF",&
             r_val=cutoff)
!> Scaling factor Hab = scaling*Sab_bar
        CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
             r_val=scaling)
!> Starting site number
        CALL section_vals_val_get(coupling_section,"FIRST_DIABAT",&
             i_val=first_diabat)
!> cbar factor : Sab = cbar*Sab_bar
        CALL section_vals_val_get(coupling_section,"CBAR",&
             r_val=cbar)
     !minao is an atomic overlap cutoff
        CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
             r_val=minao)
        CALL section_vals_val_get(coupling_section,"ELECTRONIC_PARTIAL_STEP",&
             i_val=ELsteps)

     CALL section_vals_val_get(coupling_section,"NACV_INCREMENT",&
          r_val=dR)
     CALL section_vals_val_get(coupling_section,"COLLAPSE", l_val=collapse)
     CALL section_vals_val_get(coupling_section,"ANALYTICS", l_val=analytics)
     CALL section_vals_val_get(coupling_section,"FIRST_ADIABAT",&
          i_val=first_adiabat)


!> Number of sites
     mol = nadiab ! CHANGE_AC
!> Number of atoms per site
        atoms = aom_atoms ! CHANGE_AC
!> Number of orbitals per atom
     orb = norbitals ! CHANGE_AC
!> Size of the vectorial space where the wavefunctions are projected
     mat = orb*mol
!> Numbers of active atoms 
     c = mol*atoms
!> Numbers of atoms
     natom = SIZE(global_forces(1)%forces,2)
!> Find deltaMD, nuclear timestep, in atomic units
     deltatMD = MDtimestep !CHANGE_AC
!> electronic timestep is a fraction of the molecular timestep
     deltatE = deltatMD/ELsteps !CHANGE_AC

!> ALLOCATE aom_basket at timestep = 1
    if (itimes.eq.0) then
       ALLOCATE(aom_basket%population(mat))
       ALLOCATE(aom_basket%ctin(mat))
       ALLOCATE(aom_basket%hamiltonian(mat,mat))
       ALLOCATE(aom_basket%connlist1(4,atoms))
       ALLOCATE(aom_basket%neighbourlist(mol,mol,0:atoms,0:atoms))
       ALLOCATE(aom_basket%Stt(mat,mat))
       ALLOCATE(aom_basket%Stmt(mat,mat))
       ALLOCATE(aom_basket%Sttm(mat,mat))
       ALLOCATE(aom_basket%Stmtm(mat,mat))
       ALLOCATE(aom_basket%rm1(c,3))
       ALLOCATE(aom_basket%rm2(c,3))
       ALLOCATE(aom_basket%energies(mol))
       ALLOCATE(aom_basket%bneighbour(3,0:atoms))
       ALLOCATE(aom_basket%nacv(3,mat,mat,mat*atoms))
       aom_basket%futur_decohere = .FALSE.
    endif
     
     allocate(exact_forces(3,natom))
     allocate(ctin(mat))      !ctin goes into the subroutine ZGEMV
     IF (.NOT.ALLOCATED(ctin)) CPABORT("CTIN NOT ALLOCATED CORRECTLY.")
     allocate(ctout(mat))     !ctout comes out of it (so essentially ctin = c(t), ctout = c(t+1))
     allocate(ctotal(mat))    !ctotal is |c(t)|^2: this is the output for diabatic electronic state

     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(Stt(mat,mat))
     allocate(Stpt(mat,mat))
     allocate(Sttp(mat,mat))
     allocate(Stptp(mat,mat))
     allocate(Stmt(mat,mat))
     allocate(Sttm(mat,mat))
     allocate(Stmtm(mat,mat))
     allocate(Stptpdiff(mat,mat))
     allocate(Sttpdiff(mat,mat))
     allocate(Sttdiff(mat,mat))
     allocate(Sttmdiff(mat,mat))
     allocate(Stmtmdiff(mat,mat))

     allocate(crds(3,c,3))


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
     allocate(H2save(mat,mat))
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

     allocate(nacv(3,mat,mat,mat*atoms))      !non-adiabatic coupling vectors
     nacv = 0.0D0
     allocate(siteenergies(3,mat))      !diagonals of the Hamiltonian

     allocate(pcoeffs(atoms,orb)) !atomic coefficients for the overlap routine
     allocate(pcoeffsrun(atoms,mat,3)) !pcoeffsrun exists because the p-coefficients need to be renormalised for each structure
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

     allocate(F(3,natom,mol)) !array of diagonal atomic forces

     CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff)
     CALL read_listcoeff(listcoeff, ext_elements, ext_pres, &
                         scoefftot, pcoefftot, natom)
!     allocate(elements(atoms)) !elements contains the atomic number of each atom (this goes into At, Bt etc)
     allocate(elements(c))
     elements = PACK(ext_elements, ext_pres)
     allocate(atomlist_mask(c))
     atomlist_mask =.NOT.((scoefftot.eq.0).AND.(pcoefftot.eq.0))
     allocate(atomlist(COUNT(atomlist_mask(1:atoms))))
     atomlist = PACK( (/ (i, i=1,atoms) /), atomlist_mask(1:atoms) )
    if (itimes.eq.0) then
       allocate(aom_basket%connlist2(4,size(atomlist)))
    endif

     !these arrays are data for the non-adiabatic coupling calculation
     !they help efficiency by telling the program which atoms are connected
     !and thus which pairs to carry out the calculation for
     allocate(connlist1(4,atoms))
     allocate(connlist2(4,size(atomlist)))
     allocate(aneighbour(0:atoms,0:atoms))
     allocate(bneighbour(0:atoms,0:atoms))
     aneighbour = 0
     bneighbour = 0
     allocate(neighbourlist(mol,mol,0:atoms,0:atoms))
     neighbourlist = 0

!INITIALIZATION OR READING OF CTIN      
     if (itimes.eq.0) then
          !initialise ct such that ct(init) = 1.0, otherwise 0.0 (i.e. electron localised on molecule 'init')
          do i=1,mat
               ctin(i) = CMPLX(0.0D0, 0.0D0,dp)
               ctotal(i) = CMPLX(0.0D0, 0.0D0,dp)
          end do
          ctin(first_diabat*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
          ctotal(first_diabat*(orb-1)+1) = CMPLX(1.0D0, 0.0D0,dp)
          aom_basket%population(1:mat) = abs(ctotal(1:mat)) 
     else
          ctin(:) = aom_basket%ctin(:)
     end if
!END INITIALIZATION OR READING OF CTIN      


!READ THE COORDINATES, ENERGIES AND STATE
        do coord=1,3
              if (itimes.ge.2) then
                 crds(coord,:,1) = PACK(aom_basket%rm2(:,coord), ext_pres)
                 crds(coord,:,2) = PACK(aom_basket%rm1(:,coord), ext_pres)
                 crds(coord,:,3) = PACK(positions_plus(:,coord), ext_pres)
                 siteenergies(2,:) = aom_basket%energies(:) ! CHANGE_AC
                 siteenergies(3,:) = energies_plus(:) ! CHANGE_AC
                 old_state = aom_basket%new_state
              else if (itimes.eq.1) then 
                 crds(coord,:,1) = PACK(aom_basket%rm1(:,coord), ext_pres)
                 crds(coord,:,2) = PACK(aom_basket%rm1(:,coord), ext_pres)
                 crds(coord,:,3) = PACK(positions_plus(:,coord), ext_pres)
                 siteenergies(2,:) = aom_basket%energies(:) ! CHANGE_AC
                 siteenergies(3,:) = energies_plus(:) ! CHANGE_AC
                 old_state = aom_basket%new_state
              else if (itimes.eq.0) then 
                 crds(coord,:,1) = PACK(positions_plus(:,coord), ext_pres)
                 crds(coord,:,2) = PACK(positions_plus(:,coord), ext_pres)
                 crds(coord,:,3) = PACK(positions_plus(:,coord), ext_pres)
                 siteenergies(2,:) = 0.0D0 ! CHANGE_AC
                 siteenergies(3,:) = energies_plus(i) ! CHANGE_AC
                 old_state = 1
              endif
        enddo
!READ THE COORDINATES, ENERGIES AND STATE

!READ THE COEFFICIENTS
     do i=1,atoms
          do j=1,orb
                pcoeffs(i,j) = pcoefftot( i + (j-1)*orb)
          end do
          scoeffs(i) = scoefftot(i)
     end do
     
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
                    coms(d,i) = coms(d,i) + (crds(d,(i-1)*atoms + m,3))/atoms
               end do
          end do
     end do
!END CALCULATE THE CENTER OF MASS

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
        
!CONSTRUCT THE CONNLIST
     !we fill out the connlist arrays before calculating the nacv
     if (itimes.eq.0) then
         At(2:4,1:atoms) = crds(1:3,1:atoms,3)
         call connect_list2(At(2:4,:), atomlist, connlist1, connlist2)
         aom_basket%connlist1  = connlist1
         aom_basket%connlist2  = connlist2
     else
         connlist1 = aom_basket%connlist1
         connlist2 = aom_basket%connlist2
     endif
!END CONSTRUCT THE CONNLIST


     !used to give correct commands to ZGEMX subroutines later
     zgemm_type = 'N'

           if (itimes.ge.1) then
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
                         Atm(d,j) = crds(d-1,(i-1)*atoms+j,1)
                         At(d,j)  = crds(d-1,(i-1)*atoms+j,2)
                         Atp(d,j) = crds(d-1,(i-1)*atoms+j,3)
                    end do
               end do

               !then for each other molecule
               do k=1,mol
                    !create the array of its atomic coordinates at both t and t+dt
                    do l=1,atoms
                         Bt(1,l) = dble(elements(l))
                         Btp(1,l) = dble(elements(l))
                         Btm(1,l) = dble(elements(l))
                         do d=2,4
!> CHANGE_AC : already coordonates in atomic units
                              Btm(d,l) = crds(d-1,(k-1)*atoms+l,1)
                              Bt(d,l)  =  crds(d-1,(k-1)*atoms+l,2)
                              Btp(d,l) = crds(d-1,(k-1)*atoms+l,3)
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
                    end if
                         neighbourlist(:,:,:,:) = aom_basket%neighbourlist(:,:,:,:)
                         !if (itimes.eq.1) then
                         !     do d=1,orb
                         !          do m=1,orb
                         !               call sab_screen_new(At, atomlist, connlist1, pcoeffsrun(:,orb*(i-1)+d,2),&
                         !               & atomlist, Bt, atomlist, connlist1, pcoeffsrun(:,orb*(k-1)+m,2),&
                         !               & aneighbour(:,1:atoms), bneighbour(:,1:atoms), minao)
                         !          end do
                         !     end do
                              !neighbourlist(k,1:3,0) = [0, 0, 0]
                              !neighbourlist(k,:,1:atoms) = bneighbour(:,1:atoms)
                              !if (k == mol) then
                         !          do m=1,mol
                         !               do j=1,atoms
                         !                                       aom_basket%neighbourlist(m,:,j) = & 
                         !                                               bneighbour(m,:,j)
                         !               end do
                         !          end do
                              !end if
                         !else
                         !     do m=1,mol
                              !     neighbourlist(m,1:3,0) = [0, 0, 0]
                         !          do j=1,atoms
                         !                              bneighbourlist(:,j)= aom_basket%bneighbour(:,j)
                         !                              aneighbourlist(:,j)= aom_basket%aneighbour(:,j)
                         !          end do
                         !     end do
                         !end if !(inrun.eq.FALSE.)
                    !end if  ! (i == 1)
                    !only call overlap if neighbours list says so
                    if (NN(i,k) == 1) then
                         !first the diagonals
                         if (i == k) then
                              do d=1,orb
                                   do m=1,orb
                                        !These diagonals are 1 by definition
                                        Stptp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
                                        H1(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(2,k)), dp)
                                        H2(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(3,k)),dp)

                                        !call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Btp,&
                                        !& connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                        Sttp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0

                                        !call calc_sab(Atp, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3), S, Bt,&
                                        !& connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                        Stpt(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0

                                        !call calc_dRSab(Atp, atomlist, connlist1, &
                                        !&pcoeffsrun(1:atoms,orb*(i-1)+d,2), scoeffs, &
                                        !&nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(k-1)*atoms+1:k*atoms), &
                                        !&nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(i-1)*atoms+1:i*atoms), &
                                        !&Btp, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2), scoeffs, &
                                        !&neighbourlist(i,:,:), neighbourlist(k,:,:), dR)
                                   end do
                              end do

                              if (itimes.eq.1) then
                                   do d=1,orb
                                        do m=1,orb
                                             !call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2),&
                                             !& S, Bt, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                             Stt(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0

                                             !call calc_sab(Atm, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,1),&
                                             !& S, Bt, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                             Stmt(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0

                                             !call calc_sab(At, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2),&
                                             !& S, Btm, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                                             Sttm(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0

                                             !call calc_sab(Atm, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,1),&
                                             !& S, Btm, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                                             Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
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
                                             Sttp(orb*(i-1)+d,orb*(k-1)+m)  = Sttp(orb*(k-1)+m,orb*(i-1)+d)
                                             Stpt(orb*(i-1)+d,orb*(k-1)+m)  = Stpt(orb*(k-1)+m,orb*(i-1)+d)
                                        end do
                                   end do
                              else
                                   do d=1,orb
                                        do m=1,orb
                                             call calc_sab(Atp, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                                             & S, Btp, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                             H2(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)

                                             call calc_dRSab(Atp, atomlist, connlist1, &
                                             &pcoeffsrun(1:atoms,orb*(i-1)+d,2), scoeffs, &
                                             &nacv(1:3,orb*(k-1)+d,orb*(i-1)+m,(i-1)*atoms+1:i*atoms), &
                                             &nacv(1:3,orb*(i-1)+m,orb*(k-1)+d,(k-1)*atoms+1:k*atoms), &
                                             &Btp, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2),&
                                             !&scoeffs, neighbourlist(i,:,:), neighbourlist(k,:,:), dR)
                                             &scoeffs, neighbourlist(k,i,:,:), neighbourlist(i,k,:,:), dR)
                                        end do
                                   end do

                                   !if (.not. inrun) then
                                   if (itimes.eq.1) then
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

                                        !call calc_sab(Bt, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2), S, Atp,&
                                        !& connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3))
                                        !Sttp(orb*(k-1)+m,orb*(i-1)+d) = cbar*S

                                        !call calc_sab(Btp, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3), S, At,&
                                        !& connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2))
                                        !Stpt(orb*(k-1)+m,orb*(i-1)+d) = cbar*S
                                   end do
                              end do
                         end if !(i==k)
                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   H1(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, 0.0D0,dp)
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
          if (itimes.gt.1) then
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
         else if (itimes.eq.1) then !CHANGE_AC
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
               H2save = H2
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

                    times = DBLE(itimes) + t2*1/DBLE(ELsteps)
                    CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,times)
                    CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,times)
                    CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,times)

               end do
!PRINT               print *, "Done."     
!END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
               H1 = Hsave
               H2 = H2save

      ELSE IF (itimes.eq.0) THEN
          !for each molecule
          do i=1,mol
               !create the array of all its atomic coordinates at t
               !numerical factor converts from Angstroms into Bohr radii
               Atp(1,1:atoms) = dble(elements((i-1)*atoms+1:i*atoms))
               do j=1,atoms
                    Atp(2:4,j) = crds(1:3,(i-1)*atoms+j,3)
               end do

               !then for each other molecule
               do k=1,mol
                    !create the array of its atomic coordinates at both t and t+dt
                    Btp(1,1:atoms) = dble(elements((k-1)*atoms+1:k*atoms))
                    do l=1,atoms
                        Btp(2:4,l) = crds(1:3,(k-1)*atoms+l,3)
                    end do
                    if (i == 1) then
                         do d=1,orb
                              pcoeffsrun(1:atoms,orb*(k-1)+d,3) = pcoeffs(1:atoms,d)
                              call calc_sab(Btp, connlist2, pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
                    end if  ! (i == 1)
                         do d=1,orb
                                do m=1,orb
                                     call sab_screen_new(Atp, atomlist, connlist1, pcoeffsrun(:,orb*(i-1)+d,3),&
                                     & atomlist, Btp, atomlist, connlist1, pcoeffsrun(:,orb*(k-1)+m,3),&
                                     & aneighbour(:,1:atoms), bneighbour(:,1:atoms), minao)
                                end do
                         end do
                         !neighbourlist(k,1:3,0) = [0, 0, 0]
                         neighbourlist(i,k,:,1:atoms) = bneighbour(:,1:atoms)
                         neighbourlist(k,i,:,1:atoms) = aneighbour(:,1:atoms)
                    !else
                    !          do m=1,mol
                    !               neighbourlist(m,1:3,0) = [0, 0, 0]
                    !          end do
                    !only call overlap if neighbours list says so
                    if (NN(i,k) == 1) then
                         !first the diagonals
                         if (i == k) then
                              do d=1,orb
                                   do m=1,orb
                                        !These diagonals are 1 by definition
                                        H2(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(3,k)),dp)

                                        !call calc_dRSab(Atp, atomlist, connlist1, &
                                        !&pcoeffsrun(1:atoms,orb*(i-1)+d,3), scoeffs, &
                                        !&nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(k-1)*atoms+1:k*atoms), &
                                        !&nacv(1:3,orb*(i-1)+d,orb*(k-1)+m,(i-1)*atoms+1:i*atoms), &
                                        !&Btp, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,3), scoeffs, &
                                        !&neighbourlist(i,:,:), neighbourlist(k,:,:), dR)
                                   end do
                              end do

                         !now the off-diagonals                         
                         else ! ( i == k)
                              !symmetric matrices
                              !no sense calling calc_sab twice as often as you need to
                              if (k < i) then
                                   do m=1,orb
                                        do d=1,orb
                                             H2(orb*(i-1)+d,orb*(k-1)+m) = H2(orb*(k-1)+m,orb*(i-1)+d)
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = Stptp(orb*(k-1)+m,orb*(i-1)+d)
                                        end do
                                   end do
                              else
                                   do d=1,orb
                                        do m=1,orb
                                             call calc_sab(Atp, connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                                             & S, Btp, connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                             H2(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)

                                             call calc_dRSab(Atp, atomlist, connlist1, &
                                             &pcoeffsrun(1:atoms,orb*(i-1)+d,3), scoeffs, &
                                             &nacv(1:3,orb*(k-1)+d,orb*(i-1)+m,(i-1)*atoms+1:i*atoms), &
                                             &nacv(1:3,orb*(i-1)+m,orb*(k-1)+d,(k-1)*atoms+1:k*atoms), &
                                             &Btp, atomlist, connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,3),&
                                             !&scoeffs, neighbourlist(i,:,:), neighbourlist(k,:,:), dR)
                                             &scoeffs, neighbourlist(k,i,:,:), neighbourlist(i,k,:,:), dR)
                                        end do
                                   end do
                              endif 
                         end if !(i==k)
                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   H2(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, 0.0D0,dp)
                                   Stptp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                              end do
                         end do
                    end if !(NN(i,k) == 1
               end do !loop over k
         enddo !loop-over i
         DO t2=1,int(deltatMD/deltatE)
                    times = DBLE(itimes) + t2*1/DBLE(ELsteps)
                    CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,times)
                    CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,times)
                    CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,times)
         ENDDO
      endif !END OF THE BIG IF-LOOP itimes.ge.1



!DO THE HOPPING TEST 
               !The adiabaticity subroutine diagonalises the diabatic Hamiltonian
               !It gives us back the eigenvalues and eigenvectors
      IF (itimes.ge.1) call adiabaticity(-AIMAG(H1), mat, C1, E1, .false.)
      call adiabaticity(-AIMAG(H2), mat, C2, E2, .true.)

      !The eigenvectors are needed to calculate the surface hopping probabilities
      !The surfacehop subroutine calculates these probabilites and decides on the
      !new (possibly unchanged) surface for simulation
      IF (my_bo_keyword) THEN
          new_state = first_adiabat
          old_state = first_adiabat
          decohere = .FALSE.
          is_hop = .FALSE.
          rescale_factor = 1
       ELSE    
          IF (itimes.eq.0) THEN
              new_state = first_adiabat
              old_state = first_adiabat
              decohere = .FALSE.
              is_hop = .FALSE.
              rescale_factor = 1
          ELSE
              call surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, ctin, C1, C2, deltatMD, old_state, decohere,  &      
                              new_state, is_hop)
           !> Do the velocities rescaling if required
              IF (old_state.NE.new_state) THEN
                   CPASSERT(is_hop) !Be sure there is indeed a hop
                !> Check that there is enough kinetic energy to hop between old_state and new_state
                !> But first check kinetic energies are consistant between the different diabats
                   DO imol = 2, mol
                      IF (kinetic_energies(imol).NE.kinetic_energies(1)) THEN
                       CPABORT("The different diabatic FF don't share the same kinetic energies")
                      END IF
                   ENDDO
                   delta_energy = kinetic_energies(1)-(energies_plus(new_state)-energies_plus(old_state))
                !> If delta_energy < 0: frustrated hop, do nothing!
                   IF (delta_energy.lt.0) THEN
                      new_state = old_state
                      decohere  = .FALSE.
                      is_hop    = .FALSE.
                      rescale_factor = 1
                   ELSE
                      CALL section_vals_val_get(coupling_section,"METHOD_RESCALING",c_val=rescale_keyword)
                      SELECT CASE(rescale_keyword)
                      CASE("SIMPLE")
                   !> Factor to rescale the velocities (Eq. 49 of PAPER_1)
                           rescale_factor = SQRT( ( kinetic_energies(1) - (energies_plus(new_state) - &
                                            energies_plus(old_state) ) ) / kinetic_energies(1) )
                      CASE DEFAULT
                           CPABORT("There is no method to do the rescaling of the velocities.")
                      END SELECT
                   ENDIF
              ENDIF 
          ENDIF
       ENDIF
!END DO THE HOPPING TEST 

!CALCULATE THE FORCES IN THE NEW ADIABATIC STATE!
      IF (itimes.ge.1) THEN       
               DO imol=1,mol
                  DO iatom=1, natom
                     F(:,iatom,imol) =  global_forces(imol)%forces(:,iatom)
                  ENDDO
               ENDDO
               call calcforces(state=new_state,Hin=-AIMAG(H2), C=C2, nacv=nacv, NN=NN, atoms=atoms, & 
                           scaling=scaling, cbar=cbar, Fin=F, S=Stt, forces_adiab=forces_adiab, &
                           ext_pres=ext_pres, analytics=analytics, exact_forces=exact_forces,   &
                           ctin=ctin)
               IF (analytics) THEN
                   CALL adiab_print_exact_forces(adiab_section, exact_forces, natom, times)
               ENDIF
               !Finally, collapse the wavefunction if called for
               !if (decohere) then
               if (is_hop) aom_basket%futur_decohere = .TRUE.
!> Check if the collapse key word is TRUE in the input
               if (collapse) then
!> Check if the surface_hop subroutine asks us to do the decoherence
                    if (aom_basket%futur_decohere) then 
                       if (decohere) then
                             do i=1,mat
                                 ctin(i) = C2(old_state,i) 
                             end do
                             aom_basket%futur_decohere = .FALSE. 
                       end if
                    endif
               end if
       ENDIF
!END DO THE HOPPING TEST AND CALCULATE THE FORCES

!STORE THE COORDINATES
        if (itimes.ge.1) then
           aom_basket%rm2      = aom_basket%rm1
           aom_basket%rm1      = positions_plus
           aom_basket%energies = energies_plus
        else if (itimes.eq.0) then
           aom_basket%rm1      = positions_plus
           aom_basket%energies = energies_plus
           aom_basket%neighbourlist = neighbourlist
        endif
!END STORE THE COORDINATES

!PREPARE THE PRINT
       aom_basket%ctin        = ctin
       aom_basket%hamiltonian = AIMAG(H2)
       aom_basket%nacv        = nacv
       aom_basket%new_state = new_state
       aom_basket%old_state = old_state
!END PREPARE THE PRINT

!DEALLOCATE
     deallocate(ext_pres, ext_elements)
     deallocate(elements)
     deallocate(scoefftot, pcoefftot)
     deallocate(pcoeffs, pcoeffsrun, scoeffs)
     deallocate(ctin)
     deallocate(ctout)
     deallocate(ctotal)
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
     deallocate(atomlist_mask)
     deallocate(connlist1)
     deallocate(connlist2)
     deallocate(H1, H2, H12, Hdiff, Hsave, H2save,X1, X2, X12, Ktot, K2, K3, K4)
     deallocate(Stt, Stptp, Stmtm, Sttp, Sttm, Stpt, Stmt)
     deallocate(Sttdiff, Stptpdiff, Stmtmdiff, Sttpdiff, Sttmdiff, nacv, dist)
     deallocate(exact_forces)
     deallocate(F)
!END DEALLOCATE

    END SUBROUTINE construct_ham_X_SH

END MODULE jacob_surfacehopping
