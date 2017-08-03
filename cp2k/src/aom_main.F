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
                                             calc_adiab_nace, &
                                             calc_hop_nacv, &
                                             cal_delta_energy_nacv

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
                                               propagate_adiab_wf, &
                                               prop_wf_HS_corr, & 
                                               calc_kinetic_tot, &
                                               classically_forbidden_hops
                                          
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
     !integer :: mol, i, j, k, l, m, d, atoms, c,  active_state, mat, orb
     integer :: mol, atoms, c, active_state, mat, orb
     integer, dimension(:,:), allocatable ::  aneighbour, bneighbour  
     !integer, dimension(:,:,:,:), allocatable :: neighbourlist
     REAL(KIND=dp) :: cutoff_sites, cutoff_connect, scaling, deltatMD, deltatE
     !REAL(KIND=dp) :: S
     REAL(KIND=dp)    :: cbar, minao
     !REAL(KIND=dp), dimension(:,:,:,:), allocatable :: nacv
     !REAL(KIND=dp), dimension(:,:,:), allocatable :: pcoeffsrun
!     REAL(kind=dp), DIMENSION(:,:), allocatable     :: hop_nacv
!     REAL(kind=dp)                               :: delta_energy

!    total_kin_energy, total_pot_energy
     REAL(KIND=dp), dimension(:), allocatable    :: energy_conservation
                                        !   total_kin_array, total_pot_array
     INTEGER, dimension(:), allocatable          :: class_forbidden_hop_tp
!     REAL(kind=dp)                               :: A, B
     COMPLEX(KIND=dp), dimension(:), allocatable :: ctin
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pcoefftot, scoefftot
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: C_tp
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: E_tp
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
      !real(kind=dp), allocatable :: PC_rvecs(:,:)
      !integer :: o
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

!COMMENT: if possible, transform all the mol, atoms, orb etc. in nadian,
!aom_atoms, norbitals, etc.
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
     

!COMMENT: one should check that all these arrays are used below (e.g. pcoefftot,
!scoefftot and they are de-allocated at the end of the subroutine.
     ALLOCATE(sh_env%E_tp(mat))
     ALLOCATE(sh_env%E_t(mat))
     ALLOCATE(sh_env%cttm(mat))
     ALLOCATE(sh_env%ctt(mat))
     ALLOCATE(sh_env%ctin(mat))
     !ALLOCATE(total_pot_array(mat)) 
     !ALLOCATE(total_kin_array(mat)) 
     !ALLOCATE(energy_conservation(mat)) 
     !ALLOCATE(class_forbidden_hop(mat)) 

     !ALLOCATE(E_t(mat))
     ALLOCATE(E_tp(mat))
     ALLOCATE(my_forces_adiab(3,natom))
     ALLOCATE(my_forces_adiab_new(3,natom))
     ALLOCATE(positions_com_tp(3, natom))
     ALLOCATE(velocities_com_tp(3, natom))
     allocate(exact_forces(3,natom))
     allocate(ctin(mat))      !ctin goes into the subroutine ZGEMV
     !ALLOCATE(PC_allrvecs(5,mol*atoms))
     !ALLOCATE(PC_rvecs(3,atoms))
     allocate(pcoefftot(c))
     allocate(scoefftot(c))
     !these arrays are data for the non-adiabatic coupling calculation
     !they help efficiency by telling the program which atoms are connected
     !and thus which pairs to carry out the calculation for
     allocate(aneighbour(0:atoms,0:atoms))
     allocate(bneighbour(0:atoms,0:atoms))
     aneighbour = 0
     bneighbour = 0
     !allocate(neighbourlist(mol,mol,0:atoms,0:atoms))
     !neighbourlist = 0


     !> OBTAIN COORDINATES, ENERGIES, STATES, etc. FROM PREVIOUS STEP
     CALL translate_com(positions_tp, keyword_translate, positions_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     CALL translate_com(aom_env%velocities_t, keyword_translate, velocities_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     do coord=1,3
        aom_env%crds(coord,:,1) = PACK(aom_env%rm_tm(coord,:), aom_env%ext_pres)
        aom_env%crds(coord,:,2) = PACK(aom_env%rm_t(coord,:), aom_env%ext_pres)
        aom_env%crds(coord,:,3) = PACK(positions_com_tp(coord,:), aom_env%ext_pres)
        aom_env%siteenergies(2,:) = aom_env%energies_t(:) ! CHANGE_AC
        aom_env%siteenergies(3,:) = energies_tp(:) ! CHANGE_AC
     enddo
     active_state = aom_env%active_state
     ctin(:) = aom_env%ctin(:)
     aom_env%adiab_nace_t = aom_env%adiab_nace_tp


     !> CALCULATE THE CENTER OF MASS
     !> CONSTRUCT NEIGHBOURS LIST
     CALL calc_coms_and_nlist(aom_env, atoms, mol, cutoff_sites, orb)

     !> CALCULATE THE HAMILTONIAN, THE STTP MATRIX AND THE NACV
     CALL calc_Xmatrix_elements(aom_env, itimes, atoms, mol, dR, do_random, cbar, &
                      scaling, PC_allrvecs, orb)


     !> CONVERT THE diabatic NACV WITH THE GOOD UNITS
     aom_env%nacv = cbar*aom_env%nacv

     !> HAMILTONIAN DIAGONALIZATION
     !The adiabaticity subroutine diagonalises the diabatic Hamiltonian
     !It gives us back the eigenvalues and eigenvectors
     !initialize Adiabatic energies that will be used for trivial hopping and
     !decoherence corrections
!COMMENT: in fact, C_t and E_t are already known at the previous timestep, maybe
!check aom_env%E_t, aom_env%C_t
!     CALL adiabaticity(-AIMAG(aom_env%H_t), mat, C_t, E_t, .false.)
!print *, "C_t", C_t
!print *, "E_t", E_t
     sh_env%E_t = aom_env%E_t
     CALL adiabaticity(-AIMAG(aom_env%H_tp), mat, C_tp, E_tp, .true., sh_env)
     sh_env%E_tp = E_tp

     !> ADIABATIC NACE CALCULATION
!COMMENT: do we really need to calculate the full NACE matrix? I guess one only
!need ADIAB_NACE(active_state, :)
     CALL calc_adiab_nace(aom_env%C_t, C_tp, deltatMD, aom_env)
     
     !TEST
     !CALL make_adiabX(mol, E_t, NACE_ad_t, X)
     !ENDTEST
      

     !> PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
     sh_env%cttm(:) = ctin(:)

     aom_env%HS_propagation = .FALSE.
      
     CALL section_vals_val_get(adiab_section,"RK_PROPAGATION", &
                                                  c_val=prop_wf_keyword)
     
     SELECT CASE(prop_wf_keyword)                                 
     CASE("DIABATIC_RK")
           CALL propagate_wf(ctin, aom_env, nadiab,  norbitals, ELsteps, &
                             itimes, adiab_section, deltatMD, deltatE , times,    &
                             mat)
     CASE("ADIABATIC_RK")
           CALL propagate_adiab_wf(aom_env%E_t, E_tp, ctin, aom_env, nadiab, norbitals, ELsteps, &
                                   itimes, adiab_section, deltatMD, deltatE , times,  & 
                                   mat)
     CASE("HS_CORRECTED_RK")  
           !> ADIABATIC NACV EVALUATION FOR HS CORRECTION
               
           
           CALL prop_wf_HS_corr(ctin, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times, &
                            mat, aom_env%class_forbidden_hop_t, active_state)
                                                     
           CALL classically_forbidden_hops(active_state, aom_env, mat, E_tp, &
                  velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                          atomic_kind_set, energy_conservation, class_forbidden_hop_tp)
!!TEST!
!if (active_state .eq. 1) then 
!    class_forbidden_hop_tp(1) = -1
!    class_forbidden_hop_tp(2) = 2
!else if (active_state .eq. 2) then
!    class_forbidden_hop_tp(1) = -1
!    class_forbidden_hop_tp(2) = -1
!end if  
!TEST!
!do i =1,mat
!class_forbidden_hop_tp(:) = -1
!end do
 
           aom_env%class_forbidden_hop_t = class_forbidden_hop_tp 
           aom_env%HS_propagation = .TRUE.
     CASE DEFAULT                                           
           CPABORT("There is no method to propagate the wavefunction.") 
     END SELECT                                             
     sh_env%ctt(:) = ctin(:)

     !> SURFACE HOPPING TEST
     !The eigenvectors are needed to calculate the surface hopping probabilities
     !The surfacehop subroutine calculates these probabilites and decides on the
     !new (possibly unchanged) surface for simulation
     sh_env%decohere = .FALSE.
     attempt_hop = .FALSE.
     IF (my_test_keyword) THEN   
        attempt_hop = .TRUE.
        active_state = 1
        new_state = mod(active_state + 1, mat) 
        IF ((new_state).EQ.0) new_state = mat
     ELSE IF (.NOT.my_bo_keyword) THEN ! THEN IT IS A FSSH RUN
            call surf_hop_type(adiab_section, aom_env%Stt, aom_env%C_t, C_tp,   &
                                 deltatMD, active_state, sh_env,  &
                                 new_state, &
                                 attempt_hop, &
                                 decohere_criterion, &
                                 sh_print_env, aom_env)

!       CALL section_vals_val_get(adiab_section,"RK_PROPAGATION", &
!                                                    c_val=prop_wf_keyword)
!       SELECT CASE(prop_wf_keyword)                                 
!       CASE("DIABATIC_RK")
!             call surf_hop_type(adiab_section, aom_env%Stt, C_t, C_tp,   &
!                                 deltatMD, active_state, sh_env,  &
!                                 new_state, &
!                                 attempt_hop, &
!                                 decohere_criterion, &
!                                 sh_print_env, aom_env)
!       CASE("ADIABATIC_RK")
!             call surf_hop_type(adiab_section, aom_env%Stt, C_t, C_tp,   &
!                                 deltatMD, active_state, sh_env,  &
!                                 new_state, &
!                                 attempt_hop, &
!                                 decohere_criterion, &
!                                 sh_print_env, aom_env)
!       CASE("HS_CORRECTED_RK")  
!             !> ADIABATIC NACV EVALUATION FOR HS CORRECTION
!                 
!             call surf_hop_type(adiab_section, aom_env%Stt, C_t, C_tp,   &
!                                 deltatMD, active_state, sh_env,  &
!                                 new_state, &
!                                 attempt_hop, &
!                                 decohere_criterion, &
!                                 sh_print_env, aom_env)
!             
!       CASE DEFAULT                                           
!            CPABORT("There is no method to propagate the wavefunction.") 
!       END SELECT                                             
     ENDIF


     !CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE RESULT OF THE SURFACE
     !HOP ROUTINE IS. 
     DO imol=1,mol
        DO iatom=1, natom
           aom_env%F(:,iatom,imol) =  global_forces(imol)%forces(:,iatom)
        ENDDO
     ENDDO
     print*, "ACTIVE STATE = ", active_state
     call calcforces(state=active_state,Hin=-AIMAG(aom_env%H_tp), C=C_tp, nacv=aom_env%nacv, &
                            NN=aom_env%NN, atoms=atoms, & 
                           scaling=scaling, cbar=cbar, Fin=aom_env%F, forces_adiab=my_forces_adiab, &
                           ext_pres=aom_env%ext_pres, analytics=analytics, exact_forces=exact_forces)
     IF (attempt_hop) THEN
         call calcforces(state=new_state,Hin=-AIMAG(aom_env%H_tp), C=C_tp, nacv=aom_env%nacv, &
                               NN=aom_env%NN, atoms=atoms, & 
                              scaling=scaling, cbar=cbar, Fin=aom_env%F, forces_adiab=my_forces_adiab_new, &
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
     CALL adiab_print_hamiltonian(adiab_section, -AIMAG(aom_env%H_tp), nadiab,norbitals,itimes, deltatMD)
     CALL adiab_print_nacv(adiab_section, aom_env%nacv, nadiab, norbitals, aom_atoms,itimes, deltatMD)
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
     aom_env%adiab_energies_t = aom_env%E_t
     aom_env%ctin        = ctin
     aom_env%attempt_hop  = attempt_hop
     aom_env%H_t = aom_env%H_tp
     aom_env%new_state = new_state
     aom_env%active_state = active_state
     aom_env%Stmtm = aom_env%Stt
     aom_env%Stt   = aom_env%Stptp
     aom_env%Sttm  = aom_env%Stpt
     aom_env%Stmt  = aom_env%Sttp
     aom_env%E_t = E_tp
     aom_env%C_t = C_tp

     !aom_env%Stmtm = Stt
     IF (attempt_hop) THEN
         aom_env%forces_adiab_new = my_forces_adiab_new
         aom_env%active_state = active_state
         !aom_env%NN = NN
         !aom_env%Stptp = Stptp
         aom_env%eigenvector_matrix = C_tp
         aom_env%diabatic_forces = aom_env%F
         aom_env%scaling  = scaling
         aom_env%cbar     = cbar
         aom_env%energy_new = E_tp(aom_env%new_state)
     ENDIF
     energy_adiabat = E_tp(active_state)
     aom_env%sh_print_env = sh_print_env
     

     !> DEALLOCATE
     !DEALLOCATE(energy_conservation, class_forbidden_hop)
     !DEALLOCATE(total_kin_array, total_pot_array)
     deallocate(scoefftot, pcoefftot)
     deallocate(ctin)
     deallocate(exact_forces)
    END SUBROUTINE construct_ham_X_SH


SUBROUTINE calc_coms_and_nlist(aom_env, atoms, mol, cutoff_sites, orb)
!COMMENT: one should follow the CP2K structure for subroutine: all the argument
!variable first, with a define INTENT, then the name of the subroutine and
!finally the local variable
     TYPE(aom_env_type), INTENT(INOUT)     :: aom_env
     REAL(KIND=dp), intent (IN)            :: cutoff_sites
     integer                               :: mol, i, j, m, d, atoms, orb

    
     !> CALCULATE THE CENTER OF MASS
     !calculate CoMs at first timestep
     !coms(1 to 3,i) are the x, y and z coordinates of CoM of i
     aom_env%coms = 0
     do i=1,mol
          do m=1,atoms
               do d=1,3
                    aom_env%coms(d,i) = aom_env%coms(d,i) + (aom_env%crds(d,(i-1)*atoms + m,3))/atoms
               end do
          end do
     end do


     !> CONSTRUCT NEIGHBOURS LIST
     do i=1,mol
          do j=1,mol
               aom_env%dist(i,j) = ((aom_env%coms(1,i) - aom_env%coms(1,j))*(aom_env%coms(1,i) - &
                                 aom_env%coms(1,j)) + (aom_env%coms(2,i) - aom_env%coms(2,j))*&
                      (aom_env%coms(2,i) - aom_env%coms(2,j)) + (aom_env%coms(3,i) - &
                          aom_env%coms(3,j))*(aom_env%coms(3,i) - aom_env%coms(3,j)))
               if (aom_env%dist(i,j) < cutoff_sites*cutoff_sites) then
                    do d=1,orb
                         do m=1,orb
                              aom_env%NN(orb*(i-1)+d,orb*(j-1)+m) = 1
                         end do
                    end do
               else
                    do d=1,orb
                         do m=1,orb
                              aom_env%NN(orb*(i-1)+d,orb*(j-1)+m) = 0
                         end do
                    end do
               end if                    
          end do
     end do
     

END SUBROUTINE calc_coms_and_nlist


SUBROUTINE calc_Xmatrix_elements(aom_env, itimes, atoms, mol, dR, do_random, cbar, &
                      scaling, PC_allrvecs, orb)


     TYPE(aom_env_type)                      :: aom_env
     INTEGER, POINTER                        :: itimes
     integer, intent(in)                     :: mol, orb
     integer, intent(IN)                     :: atoms
     REAL(KIND=dp), intent(in)                   :: dR  
     LOGICAL, intent(in)                         :: do_random
     REAL(KIND=dp), intent(in)                   :: cbar
     REAL(KIND=dp), intent(in) :: scaling
     !real(kind=dp), INTENT(IN) :: PC_rvecs
     !real(kind=dp), DIMENSION(:,:), intent(inout) :: PC_allrvecs
     !real(kind=dp), DIMENSION(:,:), intent(inout) :: PC_rvecs
     REAL(KIND=dp)                   :: S

     integer :: i, j, k, l, m, d, o
     real(kind=dp), DIMENSION(:,:), allocatable, intent(out) :: PC_allrvecs
     real(kind=dp), allocatable :: PC_rvecs(:,:)

     
     ALLOCATE(PC_allrvecs(5,mol*atoms))
     ALLOCATE(PC_rvecs(3,atoms))
!COMMENT: I guess Atm and Btm are only need for the first timestep (itimes = 1)
!one could duplicate this subroutine and create one only for the first timestep
     !for each molecule i
     do i=1,mol
               !create the array of all its atomic coordinates at t
               !numerical factor converts from Angstroms into Bohr radii
               do j=1,atoms
                    aom_env%At(1,j)  = dble(aom_env%elements(j))
                    aom_env%Atp(1,j) = dble(aom_env%elements(j))
                    aom_env%Atm(1,j) = dble(aom_env%elements(j))
                    do d=2,4
                         !> CHANGE_AC : already coordonates in atomic units
                         aom_env%Atm(d,j) = aom_env%crds(d-1,(i-1)*atoms+j,1)
                         aom_env%At(d,j)  = aom_env%crds(d-1,(i-1)*atoms+j,2)
                         aom_env%Atp(d,j) = aom_env%crds(d-1,(i-1)*atoms+j,3)
                    end do
               end do
               
               !then for each other molecule
               do k=1,mol
                    !create the array of its atomic coordinates at both t and t+dt
                    do l=1,atoms
                         aom_env%Bt(1,l)  = dble(aom_env%elements(l))
                         aom_env%Btp(1,l) = dble(aom_env%elements(l))
                         aom_env%Btm(1,l) = dble(aom_env%elements(l))
                         do d=2,4
                              !> CHANGE_AC : already coordonates in atomic units
                              aom_env%Btm(d,l) = aom_env%crds(d-1,(k-1)*atoms+l,1)
                              aom_env%Bt(d,l)  = aom_env%crds(d-1,(k-1)*atoms+l,2)
                              aom_env%Btp(d,l) = aom_env%crds(d-1,(k-1)*atoms+l,3)
                         end do
                    end do
                    if (i == 1) then
                         do d=1,orb
                              aom_env%pcoeffsrun(:,orb*(k-1)+d,1) = aom_env%pcoeffs(:,d)
                              call calc_sab(aom_env%Btm, aom_env%connlist2, &
                                                 aom_env%pcoeffsrun(:,orb*(k-1)+d,1))

                              aom_env%pcoeffsrun(:,orb*(k-1)+d,2) = aom_env%pcoeffs(:,d)
                              call calc_sab(aom_env%Bt, aom_env%connlist2, &
                                                 aom_env%pcoeffsrun(:,orb*(k-1)+d,2))

                              aom_env%pcoeffsrun(:,orb*(k-1)+d,3) = aom_env%pcoeffs(:,d)
                              call calc_sab(aom_env%Btp, aom_env%connlist2, & 
                                               aom_env%pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
                    end if
                    !neighbourlist(:,:,:,:) = aom_env%neighbourlist(:,:,:,:)
                    !only call overlap if neighbours list says so
                    if (aom_env%NN(i,k) == 1) then

                         !first the diagonals
                         if (i == k) then
                              CALL diagonal_elem(i,k,orb,itimes,aom_env)

                         !now the off-diagonals                         
                         else ! ( i == k)
                              CALL off_diagonal_elem(aom_env,i,k,orb,itimes, atoms, &
                                   dR, S, do_random, cbar, scaling)

                         end if !(i==k)

                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, 0.0D0,dp)
                                   aom_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = 0.0D0
                                   aom_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = 0.0D0
                                   aom_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   aom_env%Stpt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   aom_env%Sttm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   aom_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   aom_env%Stmt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                              end do
                         end do
                    end if !(NN(i,k) == 1


               end do !loop over k
               call calc_sab_pete(aom_env%Atp, aom_env%connlist2, PC_rvecs)
               do o=1, size(PC_rvecs,2)
                  PC_allrvecs(1,((i-1)*atoms)+o) = i
                  PC_allrvecs(2,((i-1)*atoms)+o) = o
                  PC_allrvecs (3:5,((i-1)*atoms)+o) = PC_rvecs(:,o)
               end do
     end do !loop over i

END SUBROUTINE calc_Xmatrix_elements

SUBROUTINE diagonal_elem(i,k,orb,itimes,aom_env)

     TYPE(aom_env_type)              :: aom_env
     INTEGER, POINTER                               :: itimes
     integer, intent(IN) :: i, k,orb

     integer ::  m, d


     do d=1,orb
          do m=1,orb

               !These diagonals are 1 by definition
               !A bit of math to be convinced for Sttp,
               !Stpt: dSii/dt = 0
               aom_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
               aom_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
               aom_env%Stpt(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
               aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, &
                                         (-aom_env%siteenergies(3,k)),dp)
          end do
     end do
     if (itimes.eq.1) then
          do d=1,orb
               do m=1,orb
               !These diagonals are 1 by definition
               !A bit of math to be convinced for Sttp,
               !Stpt: dSii/dt = 0
                    aom_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = 1.0D0
                    aom_env%Stmt(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                    aom_env%Sttm(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                    aom_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
               end do
          end do     
     end if ! (.not.inrun)

END SUBROUTINE diagonal_elem

SUBROUTINE off_diagonal_elem(aom_env,i,k,orb,itimes, atoms, &
                                   dR, S, do_random, cbar, scaling)

     TYPE(aom_env_type)                           :: aom_env
     INTEGER, POINTER                               :: itimes
     integer, intent(IN) :: i, k,orb, atoms
     !REAL(KIND=dp), dimension(:,:,:), intent(inout) :: pcoeffsrun
     REAL(KIND=dp), intent(in)                   :: dR  
     LOGICAL, intent(in)                         :: do_random
     REAL(KIND=dp), intent(in)                   :: cbar
     !integer, dimension(:,:,:,:), intent(in) :: neighbourlist
     REAL(KIND=dp), intent(in) :: scaling
     !REAL(KIND=dp), dimension(:,:,:,:), intent(out) :: nacv
     REAL(KIND=dp), intent(out)                  :: S

     integer ::  m, d
     !symmetric matrices
     !no sense calling calc_sab twice as often as you need to
     if (k < i) then
          do m=1,orb
               do d=1,orb
                    aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = aom_env%H_tp(orb*(k-1)+m,orb*(i-1)+d)
                    aom_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = &
                                             aom_env%Stt(orb*(k-1)+m,orb*(i-1)+d)
                    aom_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = &
                                            aom_env%Stptp(orb*(k-1)+m,orb*(i-1)+d)
                    aom_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = &
                                             aom_env%Stmtm(orb*(k-1)+m,orb*(i-1)+d)
                    !A cross symmetry for Sttp, Stpt,
                    !Sttm, Stmt
                    aom_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = &
                                             aom_env%Stpt(orb*(k-1)+m,orb*(i-1)+d)
                    aom_env%Stpt(orb*(i-1)+d,orb*(k-1)+m)  = & 
                                              aom_env%Sttp(orb*(k-1)+m,orb*(i-1)+d)
                    IF (itimes.eq.1) THEN
                       aom_env%Sttm(orb*(i-1)+d,orb*(k-1)+m)  = &
                                              aom_env%Stmt(orb*(k-1)+m,orb*(i-1)+d)
                       aom_env%Stmt(orb*(i-1)+d,orb*(k-1)+m)  = &
                                              aom_env%Sttm(orb*(k-1)+m,orb*(i-1)+d)
                    ENDIF
               end do
          end do
     else
          do d=1,orb
               do m=1,orb
               !PC added print_rvecs =yes so p vectors and coords would be printed 
                    call calc_sab(aom_env%Atp, aom_env%connlist2, &
                        aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                     S, aom_env%Btp, aom_env%connlist2, &
                               aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3),& 
                     print_rvecs='yes', molA_index=i, molB_index=k)
                    aom_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                    aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)

                    !Calcul of the NACV
                    call calc_dRSab(aom_env%Atp, aom_env%atomlist, aom_env%connlist1, &
                    aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,2), aom_env%scoeffs, &
                    aom_env%nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                    aom_env%nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                    aom_env%Btp, aom_env%atomlist, aom_env%connlist1, &
                    aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,2), &
                     aom_env%scoeffs, aom_env%neighbourlist(k,i,:,:), aom_env%neighbourlist(i,k,:,:), dR, &
                     do_random, aom_env)
               end do
          end do
          do d=1,orb
               do m=1,orb
                    call calc_sab(aom_env%At, aom_env%connlist2, &
                     aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, &
                     aom_env%Btp,&
                    & aom_env%connlist2, aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                    aom_env%Sttp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S

                    call calc_sab(aom_env%Atp, aom_env%connlist2, &
                         aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), S,& 
                     aom_env%Bt,&
                    & aom_env%connlist2, aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                    aom_env%Stpt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
               end do
          end do
          if (itimes.eq.1) then
               do d=1,orb
                    do m=1,orb

                         call calc_sab(aom_env%Atm, aom_env%connlist2, &
                         &aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, &
                             aom_env%Btm, aom_env%connlist2,&
                         &aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                         aom_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S


                         call calc_sab(aom_env%At, aom_env%connlist2, &
                         &aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, &
                                    aom_env%Btm, aom_env%connlist2,&
                         &aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                         aom_env%Sttm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S

                         call calc_sab(aom_env%Atm, aom_env%connlist2, &
                         &aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, &
                           aom_env%Bt, aom_env%connlist2, &
                         &aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                         aom_env%Stmt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                    end do
               end do

          end if ! (.not.inrun)
     end if ! ( k < i)


END SUBROUTINE off_diagonal_elem

END MODULE aom_main
