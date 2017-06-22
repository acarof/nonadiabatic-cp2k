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

MODULE aom_main 

     USE kinds,                        ONLY: default_string_length,&
                                             dp


     USE aom_overlapfinal,             ONLY:  connect_list2, &
                                              calc_sab, &
                                              sab_screen_new, &
                                              calc_dRSab, &
                                              calc_rvecs, &
                                              calc_sab_pete

     USE aom_utils,                    ONLY: makeX,    &
                                             make_adiabX,      &
                                             adiabaticity,   &
                                             surf_hop_type,     &
                                             calcforces,     &
                                             read_listcoeff, &
                                             choose_firstadiabat, &
                                             simple_vel_rescaling, &
                                             translate_com, &
                                             calc_adiab_nace

     USE aom_types,                    ONLY: aom_env_type, &
                                                sh_print_env_type

     USE mixed_energy_types,           ONLY: mixed_force_type
     USE input_section_types,          ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_list_get,&
                                              section_vals_val_get

     USE cp_linked_list_val,           ONLY: &
         cp_sll_val_create, cp_sll_val_dealloc, cp_sll_val_get_el_at, &
         cp_sll_val_get_length, cp_sll_val_get_rest, cp_sll_val_insert_el_at, &
         cp_sll_val_next, cp_sll_val_p_type, cp_sll_val_rm_el_at, &
         cp_sll_val_set_el_at, cp_sll_val_type

     USE aom_print,                    ONLY: adiab_print_coeff,              &
                                             adiab_print_pop,                &
                                             adiab_print_adiabat_energies,   &
                                             adiab_print_matrix_density,     &
                                             adiab_print_exact_forces,       &
                                             adiab_print_active_state,       &
                                             adiab_print_sh_info,            &
                                             adiab_print_hamiltonian, &
                                             adiab_print_nacv, &
                                             adiab_print_pvecs



     USE atomic_kind_list_types,       ONLY: atomic_kind_list_type

     USE atomic_kind_types,            ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set

     USE distribution_1d_types,        ONLY: distribution_1d_type

     USE non_adiabatic_utils,            ONLY: decohere_wf,  &
                                               propagate_wf, &
                                               propagate_adiab_wf
                                          
     USE non_adiabatic_types,           ONLY: sh_env_type

     USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream

#include "./base/base_uses.f90"


     IMPLICIT NONE

        PRIVATE
 
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='aom_main'   

        PUBLIC :: construct_ham_X_SH

CONTAINS

    SUBROUTINE construct_ham_X_SH(nadiab, aom_atoms, norbitals, MDtimestep, &
                                  positions_tp,  &
                                  energies_tp,  &
                                  aom_env, itimes, &
                                  attempt_hop, global_forces, coupling_section, &
                                  adiab_section, my_bo_keyword, energy_adiabat, &
                                  atomic_kinds, local_particles,                &
                                  my_test_keyword)


     INTEGER, INTENT(IN)                            :: nadiab,  norbitals 
     INTEGER, INTENT(IN)                            :: aom_atoms
     REAL(KIND=dp), INTENT(OUT)                     :: MDtimestep
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     REAL(kind=dp), DIMENSION(:),  POINTER          :: energies_tp
     TYPE(aom_env_type), INTENT(INOUT)              :: aom_env
     INTEGER, POINTER                               :: itimes
     LOGICAL, INTENT(OUT)                           :: attempt_hop
     TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces
     TYPE(section_vals_type), POINTER               :: coupling_section
     LOGICAL, INTENT(IN)                            :: my_bo_keyword
     TYPE(section_vals_type), POINTER               :: adiab_section
     REAL(KIND=dp), INTENT(OUT)                     :: energy_adiabat
     TYPE(atomic_kind_list_type), POINTER           :: atomic_kinds
     TYPE(distribution_1d_type), POINTER            :: local_particles
     LOGICAL, INTENT(IN)                            :: my_test_keyword
     REAL(kind=dp)                                 :: K_tot
     INTEGER                                    :: nparticle_kind
     TYPE(atomic_kind_type), DIMENSION(:), POINTER                       :: atomic_kind_set



     CHARACTER(len=*), PARAMETER :: routineN ='construct_ham_X_SH_full_degen_classicalpath', &
     routineP = moduleN//':'//routineN

     TYPE(sh_env_type)                              :: sh_env
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: my_forces_adiab
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: my_forces_adiab_new
     integer :: coord
     integer :: mol, i, j, k, l, m, d, atoms, c,  active_state, mat, orb
     integer, dimension(:,:), allocatable :: NN,  aneighbour, bneighbour  
     integer, dimension(:,:,:,:), allocatable :: neighbourlist
     REAL(KIND=dp) :: cutoff_sites, cutoff_connect, scaling, deltatMD, deltatE
     REAL(KIND=dp) :: S
     REAL(KIND=dp)    :: cbar, minao
     REAL(KIND=dp), dimension(:,:,:,:), allocatable :: nacv
     REAL(KIND=dp), dimension(:,:,:), allocatable :: crds, pcoeffsrun
     REAL(KIND=dp), dimension(:,:), allocatable :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm
     REAL(KIND=dp), dimension(:,:), allocatable :: dist
     COMPLEX(KIND=dp), dimension(:), allocatable :: ctin
     COMPLEX(KIND=dp), dimension(:,:), allocatable :: H_t, H_tp
     REAL(KIND=dp), dimension(:,:), allocatable :: At, Atm, Atp, Bt, Btm, Btp, coms, siteenergies
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pcoefftot, scoefftot
     REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: F
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: C_t, C_tp
     COMPLEX(KIND=dp), dimension(:,:), allocatable :: X
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE_ad_t  !, NACE_ad_tp
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: E_t, E_tp
    ! LOGICAL                                      :: decohere, analytics
     LOGICAL                                      :: analytics
     INTEGER                                      :: new_state , ELsteps   
     REAL(KIND=dp)                                :: dR, times 
     REAL(KIND=dp)                                :: decohere_criterion
     INTEGER                                      :: imol, iatom, natom
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: exact_forces
     LOGICAL                                      :: keyword_translate
     TYPE(sh_print_env_type)                      :: sh_print_env
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_com_tp, velocities_com_tp
     CHARACTER(LEN=default_string_length)           :: deco_keyword, prop_wf_keyword 
      real(kind=dp), DIMENSION(:,:), allocatable :: PC_allrvecs
      real(kind=dp), allocatable :: PC_rvecs(:,:)
      integer :: o
      LOGICAL    :: do_random

     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_BETWEEN_SITES",&
             r_val=cutoff_sites)
     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_CONNECTIVITY",&
             r_val=cutoff_connect)
     !> Scaling factor Hab = scaling*Sab_bar
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
             r_val=scaling)
     !> cbar factor : Sab = cbar*Sab_bar
     CALL section_vals_val_get(coupling_section,"CBAR",&
             r_val=cbar)
     !minao is an atomic overlap cutoff
     CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
             r_val=minao)
     CALL section_vals_val_get(coupling_section,"ELECTRONIC_PARTIAL_STEP",&
             i_val=ELsteps)
     CALL section_vals_val_get(coupling_section,"DECOHERENCE_CRITERION",&
             r_val=decohere_criterion)
     CALL section_vals_val_get(coupling_section,"NACV_INCREMENT",&
             r_val=dR)
     CALL section_vals_val_get(coupling_section,"ANALYTICS", l_val=analytics)
     CALL section_vals_val_get(coupling_section,"RANDOM_DERIVATIVE", l_val=do_random)
     CALL section_vals_val_get(coupling_section,"CENTER_OF_MASS", l_val=keyword_translate)

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

     sh_print_env = aom_env%sh_print_env
     
     ALLOCATE(sh_env%E_tp(mat))
     ALLOCATE(sh_env%E_t(mat))
     ALLOCATE(sh_env%cttm(mat))
     ALLOCATE(sh_env%ctt(mat))
     ALLOCATE(sh_env%ctin(mat))

     ALLOCATE(E_t(mat))
     ALLOCATE(E_tp(mat))
     ALLOCATE(my_forces_adiab(3,natom))
     ALLOCATE(my_forces_adiab_new(3,natom))
     ALLOCATE(positions_com_tp(3, natom))
     ALLOCATE(velocities_com_tp(3, natom))
     allocate(exact_forces(3,natom))
     allocate(ctin(mat))      !ctin goes into the subroutine ZGEMV
     ALLOCATE(PC_allrvecs(5,mol*atoms))
     ALLOCATE(PC_rvecs(3,atoms))
     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(Stt(mat,mat))
     allocate(Stpt(mat,mat))
     allocate(Sttp(mat,mat))
     allocate(Stptp(mat,mat))
     allocate(Stmt(mat,mat))
     allocate(Sttm(mat,mat))
     allocate(Stmtm(mat,mat))
     allocate(crds(3,c,3))
     allocate(NACE_ad_t(mat,mat))
     !allocate(NACE_ad_tp(mat,mat))
     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H_t etc are always -i*H.
     !This is why they're complex variables
     allocate(H_t(mat,mat))
     allocate(H_tp(mat,mat))
     allocate(nacv(3,mat,mat,mat*atoms))      !non-adiabatic coupling vectors
     nacv = 0.0D0
     allocate(siteenergies(3,mat))      !diagonals of the Hamiltonian
     allocate(pcoeffsrun(atoms,mat,3)) !pcoeffsrun exists because the p-coefficients need to be renormalised for each structure
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
     !these arrays are data for the non-adiabatic coupling calculation
     !they help efficiency by telling the program which atoms are connected
     !and thus which pairs to carry out the calculation for
     allocate(aneighbour(0:atoms,0:atoms))
     allocate(bneighbour(0:atoms,0:atoms))
     aneighbour = 0
     bneighbour = 0
     allocate(neighbourlist(mol,mol,0:atoms,0:atoms))
     neighbourlist = 0


     !> OBTAIN COORDINATES, ENERGIES, STATES, etc. FROM PREVIOUS STEP
     CALL translate_com(positions_tp, keyword_translate, positions_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     CALL translate_com(aom_env%velocities_t, keyword_translate, velocities_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     do coord=1,3
        crds(coord,:,1) = PACK(aom_env%rm_tm(coord,:), aom_env%ext_pres)
        crds(coord,:,2) = PACK(aom_env%rm_t(coord,:), aom_env%ext_pres)
        crds(coord,:,3) = PACK(positions_com_tp(coord,:), aom_env%ext_pres)
        siteenergies(2,:) = aom_env%energies_t(:) ! CHANGE_AC
        siteenergies(3,:) = energies_tp(:) ! CHANGE_AC
     enddo
     if (itimes.gt.1) then   !!itimes strictly greater than 1
        Sttm  = aom_env%Sttm
        Stmt  = aom_env%Stmt
        Stmtm = aom_env%Stmtm
     end if
     Stt = aom_env%Stt
     active_state = aom_env%active_state
     H_t = aom_env%H_t
     ctin(:) = aom_env%ctin(:)
     aom_env%adiab_nace_t = aom_env%adiab_nace_tp


     !> CALCULATE THE CENTER OF MASS
     !calculate CoMs at first timestep
     !coms(1 to 3,i) are the x, y and z coordinates of CoM of i
     coms = 0
     do i=1,mol
          do m=1,atoms
               do d=1,3
                    coms(d,i) = coms(d,i) + (crds(d,(i-1)*atoms + m,3))/atoms
               end do
          end do
     end do


     !> CONSTRUCT NEIGHBOURS LIST
     do i=1,mol
          do j=1,mol
               dist(i,j) = ((coms(1,i) - coms(1,j))*(coms(1,i) - coms(1,j)) + (coms(2,i) - coms(2,j))*&
               &(coms(2,i) - coms(2,j)) + (coms(3,i) - coms(3,j))*(coms(3,i) - coms(3,j)))
               if (dist(i,j) < cutoff_sites*cutoff_sites) then
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


     !> CALCULATE THE HAMILTONIAN, THE STTP MATRIX AND THE NACV
     !for each molecule i
     do i=1,mol
               !create the array of all its atomic coordinates at t
               !numerical factor converts from Angstroms into Bohr radii
               do j=1,atoms
                    At(1,j)  = dble(aom_env%elements(j))
                    Atp(1,j) = dble(aom_env%elements(j))
                    Atm(1,j) = dble(aom_env%elements(j))
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
                         Bt(1,l)  = dble(aom_env%elements(l))
                         Btp(1,l) = dble(aom_env%elements(l))
                         Btm(1,l) = dble(aom_env%elements(l))
                         do d=2,4
                              !> CHANGE_AC : already coordonates in atomic units
                              Btm(d,l) = crds(d-1,(k-1)*atoms+l,1)
                              Bt(d,l)  = crds(d-1,(k-1)*atoms+l,2)
                              Btp(d,l) = crds(d-1,(k-1)*atoms+l,3)
                         end do
                    end do
                    if (i == 1) then
                         do d=1,orb
                              pcoeffsrun(:,orb*(k-1)+d,1) = aom_env%pcoeffs(:,d)
                              call calc_sab(Btm, aom_env%connlist2, pcoeffsrun(:,orb*(k-1)+d,1))

                              pcoeffsrun(:,orb*(k-1)+d,2) = aom_env%pcoeffs(:,d)
                              call calc_sab(Bt, aom_env%connlist2, pcoeffsrun(:,orb*(k-1)+d,2))

                              pcoeffsrun(:,orb*(k-1)+d,3) = aom_env%pcoeffs(:,d)
                              call calc_sab(Btp, aom_env%connlist2, pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
                    end if
                    neighbourlist(:,:,:,:) = aom_env%neighbourlist(:,:,:,:)
                    !only call overlap if neighbours list says so
                    if (NN(i,k) == 1) then
                         !first the diagonals
                         if (i == k) then
                              do d=1,orb
                                   do m=1,orb
                                        !These diagonals are 1 by definition
                                        !A bit of math to be convinced for Sttp,
                                        !Stpt: dSii/dt = 0
                                        Stptp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
                                        Sttp(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                                        Stpt(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                                        !H_t(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(2,k)), dp)
                                        H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(3,k)),dp)
                                   end do
                              end do
                              if (itimes.eq.1) then
                                   do d=1,orb
                                        do m=1,orb
                                        !These diagonals are 1 by definition
                                        !A bit of math to be convinced for Sttp,
                                        !Stpt: dSii/dt = 0
                                             Stt(orb*(i-1)+d,orb*(k-1)+m)   = 1.0D0
                                             Stmt(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                                             Sttm(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
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
                                             !H_t(orb*(i-1)+d,orb*(k-1)+m)   = H_t(orb*(k-1)+m,orb*(i-1)+d)
                                             H_tp(orb*(i-1)+d,orb*(k-1)+m)  = H_tp(orb*(k-1)+m,orb*(i-1)+d)
                                             Stt(orb*(i-1)+d,orb*(k-1)+m)   = Stt(orb*(k-1)+m,orb*(i-1)+d)
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = Stptp(orb*(k-1)+m,orb*(i-1)+d)
                                             Stmtm(orb*(i-1)+d,orb*(k-1)+m) = Stmtm(orb*(k-1)+m,orb*(i-1)+d)
                                             !A cross symmetry for Sttp, Stpt,
                                             !Sttm, Stmt
                                             Sttp(orb*(i-1)+d,orb*(k-1)+m)  = Stpt(orb*(k-1)+m,orb*(i-1)+d)
                                             Stpt(orb*(i-1)+d,orb*(k-1)+m)  = Sttp(orb*(k-1)+m,orb*(i-1)+d)
                                             IF (itimes.eq.1) THEN
                                                Sttm(orb*(i-1)+d,orb*(k-1)+m)  = Stmt(orb*(k-1)+m,orb*(i-1)+d)
                                                Stmt(orb*(i-1)+d,orb*(k-1)+m)  = Sttm(orb*(k-1)+m,orb*(i-1)+d)
                                             ENDIF
                                        end do
                                   end do
                              else
                                   do d=1,orb
                                        do m=1,orb
                                        !PC added print_rvecs =yes so p vectors and coords would be printed 
                                             call calc_sab(Atp, aom_env%connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                                             & S, Btp, aom_env%connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3),& 
                                             & print_rvecs='yes', molA_index=i, molB_index=k)
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                             H_tp(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)

                                             !Calcul of the NACV
                                             call calc_dRSab(Atp, aom_env%atomlist, aom_env%connlist1, &
                                             &pcoeffsrun(1:atoms,orb*(i-1)+d,2), aom_env%scoeffs, &
                                             &nacv(1:3,orb*(k-1)+d,orb*(i-1)+m,(i-1)*atoms+1:i*atoms), &
                                             &nacv(1:3,orb*(i-1)+m,orb*(k-1)+d,(k-1)*atoms+1:k*atoms), &
                                             &Btp, aom_env%atomlist, aom_env%connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,2),&
                                             &aom_env%scoeffs, neighbourlist(k,i,:,:), neighbourlist(i,k,:,:), dR,  &
                                              do_random, aom_env)
                                        end do
                                   end do
                                   do d=1,orb
                                        do m=1,orb
                                             call calc_sab(At, aom_env%connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Btp,&
                                             & aom_env%connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                             Sttp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S

                                             call calc_sab(Atp, aom_env%connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3), S, Bt,&
                                             & aom_env%connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                             Stpt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                        end do
                                   end do
                                   if (itimes.eq.1) then
                                        do d=1,orb
                                             do m=1,orb
                                                  !call calc_sab(At, aom_env%connlist2, &
                                                  !&pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Bt, aom_env%connlist2,&
                                                  !&pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                                  !Stt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                                  !H_t(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)

                                                  call calc_sab(Atm, aom_env%connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, Btm, aom_env%connlist2,&
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                                                  Stmtm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S


                                                  call calc_sab(At, aom_env%connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, Btm, aom_env%connlist2,&
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                                                  Sttm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S

                                                  call calc_sab(Atm, aom_env%connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, Bt, aom_env%connlist2, &
                                                  &pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                                                  Stmt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                             end do
                                        end do

                                   end if ! (.not.inrun)
                              end if ! ( k < i)

                         end if !(i==k)
                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   !H_t(orb*(i-1)+d,orb*(k-1)+m)   = cmplx(0.0D0, 0.0D0,dp)
                                   H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, 0.0D0,dp)
                                   Stt(orb*(i-1)+d,orb*(k-1)+m)   = 0.0D0
                                   Sttp(orb*(i-1)+d,orb*(k-1)+m)  = 0.0D0
                                   Stptp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Stpt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Sttm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   Stmt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                              end do
                         end do
                    end if !(NN(i,k) == 1
               end do !loop over k
               call calc_sab_pete(Atp, aom_env%connlist2, PC_rvecs)
               do o=1, size(PC_rvecs,2)
                  PC_allrvecs(1,((i-1)*atoms)+o) = i
                  PC_allrvecs(2,((i-1)*atoms)+o) = o
                  PC_allrvecs (3:5,((i-1)*atoms)+o) = PC_rvecs(:,o)
               end do
     end do !loop over i


     !> CONVERT THE NACV WITH THE GOOD UNITS
     nacv = cbar*nacv
     !The adiabaticity subroutine diagonalises the diabatic Hamiltonian
     !It gives us back the eigenvalues and eigenvectors
     !initialize Adiabatic energies that will be used for trivial hopping and
     !decoherence corrections
     CALL adiabaticity(-AIMAG(H_t), mat, C_t, E_t, .false.)
     sh_env%E_t = E_t
     CALL adiabaticity(-AIMAG(H_tp), mat, C_tp, E_tp, .true., sh_env)
     sh_env%E_tp = E_tp


     !> PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
     CALL calc_adiab_nace(Sttp, Stt, Stptp, Sttm, Stmtm, C_t, C_tp, deltatMD, aom_env)
     !NACE_ad_tp = aom_env%adiab_nace_tp
     NACE_ad_t = aom_env%adiab_nace_t
     
     !TEST
     CALL make_adiabX(mol, E_t, NACE_ad_t, X)
     !ENDTEST
     sh_env%cttm(:) = ctin(:)
     CALL section_vals_val_get(adiab_section,"REPRESENTATION", &
                                                  c_val=prop_wf_keyword)
     SELECT CASE(prop_wf_keyword)                                 
     CASE("DIABATIC_BASIS")
           CALL propagate_wf(H_t, H_tp,  Stt, Sttp, Stptp, Stpt, Stmt, & 
                             Sttm, Stmtm,    &
                             ctin, aom_env, nadiab,  norbitals, ELsteps, &
                             itimes, adiab_section, deltatMD, deltatE , times,    &
                             mat)
     CASE("ADIABATIC_BASIS")
           CALL propagate_adiab_wf(E_t, E_tp, ctin, aom_env, nadiab, norbitals, ELsteps, &
                                   itimes, adiab_section, deltatMD, deltatE , times,  & 
                                   mat)
     CASE("HS_CORRECTION")  
            print *, "TO BE ADDED"
     CASE DEFAULT                                           
           CPABORT("There is no method to propagate the wavefunction.") 
     END SELECT                                             
     sh_env%ctt(:) = ctin(:)


     !The eigenvectors are needed to calculate the surface hopping probabilities
     !The surfacehop subroutine calculates these probabilites and decides on the
     !new (possibly unchanged) surface for simulation
     !decohere = .FALSE.
     sh_env%decohere = .FALSE.
     attempt_hop = .FALSE.
     IF (my_test_keyword) THEN   
        attempt_hop = .TRUE.
        active_state = 1
        new_state = mod(active_state + 1, mat) 
        IF ((new_state).EQ.0) new_state = mat
     ELSE IF (.NOT.my_bo_keyword) THEN ! THEN IT IS A FSSH RUN
        call surf_hop_type(adiab_section, Stt, C_t, C_tp,   &
                                  deltatMD, active_state, sh_env,  &
                                  new_state, &
                                  attempt_hop, &
                                  decohere_criterion, &
                                  sh_print_env, aom_env)
     ENDIF


     !CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE RESULT OF THE SURFACE
     !HOP ROUTINE IS. 
     DO imol=1,mol
        DO iatom=1, natom
           F(:,iatom,imol) =  global_forces(imol)%forces(:,iatom)
        ENDDO
     ENDDO
     print*, "ACTIVE STATE = ", active_state
     call calcforces(state=active_state,Hin=-AIMAG(H_tp), C=C_tp, nacv=nacv, NN=NN, atoms=atoms, & 
                           scaling=scaling, cbar=cbar, Fin=F, forces_adiab=my_forces_adiab, &
                           ext_pres=aom_env%ext_pres, analytics=analytics, exact_forces=exact_forces)
     IF (attempt_hop) THEN
         call calcforces(state=new_state,Hin=-AIMAG(H_tp), C=C_tp, nacv=nacv, NN=NN, atoms=atoms, & 
                              scaling=scaling, cbar=cbar, Fin=F, forces_adiab=my_forces_adiab_new, &
                              ext_pres=aom_env%ext_pres, analytics=analytics, exact_forces=exact_forces)
     ENDIF


     !> APPLY DECOHERENCE CORRECTION
     !> initialize derived type variables
     sh_env%ctin=ctin
     sh_env%active_state=active_state
     CALL decohere_wf(adiab_section, sh_print_env, aom_env,sh_env, deco_keyword, itimes, &
                                    deltatMD, mol,velocities_com_tp, local_particles, atomic_kinds, nparticle_kind,    &
                                    atomic_kind_set, K_tot, aom_env%ext_pres)
     ctin=sh_env%ctin


     !> PRINT
     CALL adiab_print_hamiltonian(adiab_section, -AIMAG(H_tp), nadiab,norbitals,itimes, deltatMD)
     CALL adiab_print_nacv(adiab_section, nacv, nadiab, norbitals, aom_atoms,itimes, deltatMD)
     CALL adiab_print_adiabat_energies(adiab_section, E_tp, nadiab, deltatMD, times)
     CALL adiab_print_pvecs(adiab_section, PC_allrvecs, itimes, deltatMD, times, atoms, nadiab)
                         sh_print_env%times = times
                         sh_print_env%times = DBLE(itimes)
     IF (analytics) THEN
        times = DBLE(itimes)
        CALL adiab_print_exact_forces(adiab_section, exact_forces, natom, deltatMD, times)
     ENDIF


     !> STORE THE COORDINATES IF NOT SUCCESSFUL RUN
     aom_env%forces_adiab=my_forces_adiab
     aom_env%rm_tm      = aom_env%rm_t
     aom_env%rm_t       = positions_com_tp
     aom_env%energies_t = energies_tp
     aom_env%adiab_energies_t = E_t
     aom_env%ctin        = ctin
     aom_env%attempt_hop  = attempt_hop
     aom_env%H_t = H_tp
     aom_env%nacv        = nacv
     aom_env%new_state = new_state
     aom_env%active_state = active_state
     aom_env%Stt   = Stptp
     aom_env%Sttm  = Stpt
     aom_env%Stmt  = Sttp
     aom_env%Stmtm = Stt
     IF (attempt_hop) THEN
         aom_env%forces_adiab_new = my_forces_adiab_new
         aom_env%active_state = active_state
         aom_env%NN = NN
         aom_env%Stptp = Stptp
         aom_env%eigenvector_matrix = C_tp
         aom_env%diabatic_forces = F
         aom_env%scaling  = scaling
         aom_env%cbar     = cbar
         aom_env%energy_new = E_tp(aom_env%new_state)
     ENDIF
     energy_adiabat = E_tp(active_state)
     aom_env%sh_print_env = sh_print_env


     !> DEALLOCATE
     deallocate(scoefftot, pcoefftot)
     deallocate(pcoeffsrun)
     deallocate(ctin)
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
     deallocate(H_t, H_tp)
     deallocate(Stt, Stptp, Stmtm, Sttp, Sttm, Stpt, Stmt)
     deallocate(nacv, dist)
     deallocate(exact_forces)
     deallocate(F)
    END SUBROUTINE construct_ham_X_SH

END MODULE aom_main
