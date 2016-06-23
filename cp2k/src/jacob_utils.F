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

MODULE jacob_utils
   USE kinds,                           ONLY: default_string_length,&
                                               dp
   USE machine,                         ONLY: m_getpid
   USE jacob_types,                     ONLY: sh_print_basket_type

  USE cp_parser_methods,               ONLY: read_float_object

  USE cp_output_handling,              ONLY: cp_print_key_unit_nr,&
                                             cp_print_key_finished_output
  USE cp_log_handling,                 ONLY: cp_logger_type,&
                                             cp_logger_release,&
                                             cp_get_default_logger,&
                                             cp_to_string
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             get_results,&
                                             put_results,&
                                             test_for_result
  USE cp_result_types,                 ONLY: cp_result_p_type,&
                                             cp_result_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_list_get,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mixed_energy_types,              ONLY: mixed_force_type
  USE particle_list_types,             ONLY: particle_list_type
  USE virial_types,                    ONLY: virial_p_type,&
                                             virial_type,&
                                             zero_virial

  USE mixed_environment_utils,        ONLY: get_subsys_map_index
  USE mixed_environment_types,        ONLY: mixed_environment_type
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE cp_linked_list_val,              ONLY: &
       cp_sll_val_create, cp_sll_val_dealloc, cp_sll_val_get_el_at, &
       cp_sll_val_get_length, cp_sll_val_get_rest, cp_sll_val_insert_el_at, &
       cp_sll_val_next, cp_sll_val_p_type, cp_sll_val_rm_el_at, &
       cp_sll_val_set_el_at, cp_sll_val_type


  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set

  USE distribution_1d_types,           ONLY: distribution_1d_type


#include "./base/base_uses.f90"



IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'jacob_utils'

   !PUBLIC :: makeX, &
   PUBLIC :: makeX,            &
   !          sitecalc, &
             surfacehop,       &
   !          init_random_seed, &
   !          matinvsqrt, &
             adiabaticity,     &
             !forces, &
             calcforces,           &
             read_listcoeff,       & 
             choose_firstadiabat,  &
             simple_vel_rescaling, &
             calc_velocities_prime
   !          calcTDT

CONTAINS
  
!subroutine makeX(mol, H, Stt, Sttp, Sttm, Stptp, Stmtm, X, deltat)
subroutine makeX(mol, H, Stt, Sttp, Sttm, Stptp, Stmtm, X, deltat, NACE)
        !X is the name I give to the matrix in the solution to the Schroedinger equation:
        !dc/dt = X*c
        !makeX constructs an X matrix from inputs of S and H

        implicit none

!        interface
!                subroutine matinvsqrt(mol, invert, S)
!                        integer, intent(in) :: mol
!                        logical, intent(in) :: invert
!                        REAL(KIND=dp), dimension(:,:), intent(inout) :: S
!                end subroutine matinvsqrt
!!        end interface

        !integer :: stat, i, j
        integer, intent(in) :: mol
        REAL(KIND=dp), intent(in) :: deltat
        !REAL(KIND=dp), dimension(:), allocatable :: DGETR_IPIV, DGETRI_WORK, D
        REAL(KIND=dp), dimension(:,:), intent(in) :: Stt, Sttp, Sttm, Stptp, Stmtm
        REAL(KIND=dp), dimension(:,:), allocatable :: Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt
        REAL(KIND=dp), dimension(:,:), allocatable :: Sterm1, Sterm2, Splace1, Splace2
        COMPLEX(KIND=dp), dimension(:,:), intent(in) :: H
        COMPLEX(KIND=dp), dimension(:,:), intent(out) :: X
        REAL(KIND=dp), dimension(:,:), intent(out), optional    :: NACE   

        !allocate matrices
        allocate(Sttinvsqrt(mol,mol))
        if (.not. allocated(Sttinvsqrt)) then
                CPABORT("Error: array 'Sttinvsqrt' not allocated correctly.")
                !stop "Error: array 'Sttinvsqrt' not allocated correctly."
        end if
        allocate(Sttsqrt(mol,mol))
        if (.not. allocated(Sttsqrt)) then
                CPABORT("Error: array 'Sttsqrt' not allocated correctly.")
                !stop "Error: array 'Sttsqrt' not allocated correctly."
        end if
        allocate(Stptpinvsqrt(mol,mol))
        if (.not. allocated(Stptpinvsqrt)) then
                CPABORT("Error: array 'Stptpinvsqrt' not allocated correctly.")
                !stop "Error: array 'Stptpinvsqrt' not allocated correctly."
        end if
        allocate(Stmtminvsqrt(mol,mol))
        if (.not. allocated(Stmtminvsqrt)) then
                CPABORT("Error: array 'Stmtminvsqrt' not allocated correctly.")
                !stop "Error: array 'Stmtminvsqrt' not allocated correctly."
        end if
        allocate(Sterm1(mol,mol))
        if (.not. allocated(Sterm1)) then
                CPABORT("Error: array 'Sterm1' not allocated correctly.")
                !stop "Error: array 'Sterm1' not allocated correctly."
        end if
        allocate(Sterm2(mol,mol))
        if (.not. allocated(Sterm2)) then
                CPABORT("Error: array 'Sterm2' not allocated correctly.")
                !stop "Error: array 'Sterm2' not allocated correctly."
        end if
        allocate(Splace1(mol,mol))
        if (.not. allocated(Splace1)) then
                CPABORT("Error: array 'Splace1' not allocated correctly.")
                !stop "Error: array 'Splace1' not allocated correctly."
        end if
        allocate(Splace2(mol,mol))
        if (.not. allocated(Splace2)) then
                CPABORT("Error: array 'Splace2' not allocated correctly.")
                !stop "Error: array 'Splace2' not allocated correctly."
        end if

!print*, 'Sttm'
!print*, Sttm
!print*, 'Sttsqrt'
        Sttsqrt = Stt
        call matinvsqrt(mol, .FALSE., Sttsqrt)
!print*, Sttsqrt
!print*, "Sttinvsqrt"
        Sttinvsqrt = Stt
        call matinvsqrt(mol, .TRUE., Sttinvsqrt)
!print*, Sttinvsqrt
!print*, 'Stptpinvsqrt'
        Stptpinvsqrt = Stptp
        call matinvsqrt(mol, .TRUE., Stptpinvsqrt)
!print*, Stptpinvsqrt
!print*, 'Stmtminvsqrt'
        Stmtminvsqrt = Stmtm
        call matinvsqrt(mol, .TRUE., Stmtminvsqrt)
!print*, Stmtminvsqrt
        
        !subroutine DGEMM multiplies two matrices, from the BLAS
        !use to build the S terms in the X matrix

        !first term: Sttsqrt*(Stptpinvsqrt - Stmtmtinvsqrt)
        Splace1 = Stptpinvsqrt - Stmtminvsqrt
!print*, "Splace1"
!print*, Splace1
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Sttsqrt, mol, Splace1, mol, 0.0D0, Sterm1, mol)

        !second term: Sttinvsqrt*(Sttp - Sttm)*Sttinvsqrt
        Splace1 = Sttp - Sttm
!print*, "Splace1"
!print*, Splace1
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Splace1, mol, Sttinvsqrt, mol, 0.0D0, Splace2, mol)
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Sttinvsqrt, mol, Splace2, mol, 0.0D0, Sterm2, mol)
        X = H - (0.5D0/deltat)*(Sterm1 + Sterm2)

!> CHANGE_AC
        IF (PRESENT(NACE)) THEN
           NACE = (0.5D0/deltat)*(Sterm1 + Sterm2) 
!print*, "Sterm1"
!print*, Sterm1
!print*, "Sterm2"
!print*, Sterm2
!print*, "NACE"
!print*, NACE
        ENDIF

        deallocate(Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt, Sterm1, Sterm2, Splace1, Splace2)

end subroutine makeX

subroutine sitecalc(A, bondslist, atoms, site, F, mat, molno)

        !calculates the site energy of a molecule (in kcal/mol)
        
        !This subroutine makes some assumptions:
        !- the energy difference between charged and neutral molecule can be described in terms of
        !  bond length differences between the two
        !- all the molecules in the system are the same
        !- there are only two bond types per molecule
        !- both bond types have the same force constant k
        !This subroutine should be easily modified to remove the latter two assumptions but the
        !first two are currently integral to this method
        !
        !If the arguments F, mat and molno are present, this subroutine will also calculate
        !the diagonal forces on each atom.

        implicit none
        REAL(KIND=dp), dimension(:,:), intent(inout) :: A
        REAL(KIND=dp), dimension(:,:), optional, intent(out) :: F
        integer, optional, intent(in) :: mat, molno
        integer, dimension(:,:), intent(in) :: bondslist
        REAL(KIND=dp), intent(out) :: site
        !REAL(KIND=dp) :: k, r01n, r02n, r01c, r02c, e, dist, kcal, bohrang, F1
        REAL(KIND=dp) :: k, r01n, r02n, r01c, r02c, dist, kcal, bohrang, F1
        integer, intent(in) :: atoms
        integer :: i, j, l, a1, a2
        
        site = 0.0D0
        k = 1179.4D0
        r01n = 1.087D0
        r02n = 1.324D0
        r01c = 1.087D0
        r02c = 1.401D0
        kcal = 627.509D0
        bohrang = 1.889725989D0
        if(present(F)) then
                do i=1,mat*atoms
                        do j=1,3
                                F(j,i) = 0.0D0
                        end do
                end do
        end if

        !Coordinates come to us in Bohr, but since k, r01 etc are defined from MD (in Angstrom)
        !we need to convert, then later convert back
        A = A/bohrang

        do i=1,size(bondslist,2)
                a1 = bondslist(1,i)
                a2 = bondslist(2,i)
                dist = sqrt((A(2,a1) - A(2,a2))**2 + (A(3,a1) - A(3,a2))**2 + (A(4,a1) - A(4,a2))**2)
                if (bondslist(3,i) == 1) then
                        site = site + 0.5D0*k*((dist - r01c)*(dist - r01c) - (dist - r01n)*(dist - r01n))
                        if(present(F)) then
                                F1 = (k*((dist - r01c) - (dist - r01n)))
                                do j=1,3
                                        l = atoms*(molno-1)+a1
                                        F(j,l) = F(j,l) + (F1*(A(j+1,a2) - A(j+1,a1))/(dist))
                                        l = atoms*(molno-1)+a2
                                        F(j,l) = F(j,l) - (F1*(A(j+1,a2) - A(j+1,a1))/(dist))
                                end do
                        end if
                else if (bondslist(3,i) == 2) then
                        site = site + 0.5D0*k*((dist - r02c)*(dist - r02c) - (dist - r02n)*(dist - r02n))
                        if(present(F)) then
                                F1 = (k*((dist - r02c) - (dist - r02n)))
                                do j=1,3
                                        l = atoms*(molno-1)+a1
                                        F(j,l) = F(j,l) + (F1*(A(j+1,a2) - A(j+1,a1))/(dist))
                                        l = atoms*(molno-1)+a2
                                        F(j,l) = F(j,l) - (F1*(A(j+1,a2) - A(j+1,a1))/(dist))
                                end do
                        end if
                end if
        end do

        A = A*bohrang                

end subroutine sitecalc

!subroutine surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, c, C1, C2, deltat, state, decohere)
subroutine surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, c, C1, C2, deltat, old_state, & 
                      decohere, new_state, is_hop, decohere_criterion, sh_print_basket)

        !using the method outlined in John Tully's Electronic Transitions paper
        !calculate the probability of the molecular dynamics switching to evolve on a different potential surface

        implicit none

!        interface
!                subroutine matinvsqrt(mol, invert, S)
!                        integer, intent(in) :: mol
!                        logical, intent(in) :: invert
!                        REAL(KIND=dp), dimension(:,:), intent(inout) :: S
!                end subroutine matinvsqrt
!        end interface
        REAL(KIND=dp), INTENT(IN) :: decohere_criterion

        integer :: i, j, mol, k, l
        integer, intent(inout) :: old_state
        REAL(KIND=dp) :: hop, sumg
        REAL(KIND=dp), intent(in) :: deltat
        REAL(KIND=dp), dimension(:,:), allocatable :: dphi
        REAL(KIND=dp), dimension(:), allocatable :: b, g
        REAL(KIND=dp), dimension(:,:), intent(in) :: Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2
        REAL(KIND=dp), dimension(:,:), allocatable :: T1, T2, DT, TST, TDT, DS, DC
        COMPLEX(KIND=dp), dimension(:), intent(in) :: c
        COMPLEX(KIND=dp), dimension(:), allocatable :: ct
        COMPLEX(KIND=dp), dimension(:), allocatable :: ak
        COMPLEX(KIND=dp) :: cstate
        logical, intent(out) :: decohere
        INTEGER, INTENT(OUT)        :: new_state
        LOGICAL, INTENT(OUT)        :: is_hop
        TYPE(sh_print_basket_type), INTENT(OUT) :: sh_print_basket

        mol = size(Stt,1)

        !many allocations
        allocate(ak(mol))
        if (.not. allocated(ak)) then
                CPABORT("Error: array 'ak' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'ak' not allocated in subroutine 'surfacehop'"
        end if
        allocate(b(mol))
        if (.not. allocated(b)) then
                CPABORT("Error: array 'b' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'b' not allocated in subroutine 'surfacehop'"
        end if
        allocate(g(mol))
        if (.not. allocated(g)) then
                CPABORT("Error: array 'g' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'g' not allocated in subroutine 'surfacehop'"
        end if
        allocate(ct(mol))
        if (.not. allocated(ct)) then
                CPABORT("Error: array 'ct' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'ct' not allocated in subroutine 'surfacehop'"
        end if
        allocate(dphi(mol,mol))
        if (.not. allocated(dphi)) then
                CPABORT("Error: array 'dphi' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'dphi' not allocated in subroutine 'surfacehop'."
        end if
        allocate(T1(mol,mol))
        if (.not. allocated(T1)) then
                CPABORT("Error: array 'T1' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'T1' not allocated in subroutine 'surfacehop'."
        end if
        allocate(T2(mol,mol))
        if (.not. allocated(T2)) then
                CPABORT("Error: array 'T2' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'T2' not allocated in subroutine 'surfacehop'."
        end if
        allocate(DT(mol,mol))
        if (.not. allocated(DT)) then
                CPABORT("Error: array 'DT' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'DT' not allocated in subroutine 'surfacehop'."
        end if
        allocate(TST(mol,mol))
        if (.not. allocated(TST)) then
                CPABORT("Error: array 'TST' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'TST' not allocated in subroutine 'surfacehop'."
        end if
        allocate(TDT(mol,mol))
        if (.not. allocated(TDT)) then
                CPABORT("Error: array 'TDT' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'TDT' not allocated in subroutine 'surfacehop'."
        end if
        allocate(DS(mol,mol))
        if (.not. allocated(DS)) then
                !stop "Error: array 'DS' not allocated in subroutine 'surfacehop'."
                CPABORT("Error: array 'DS' not allocated in subroutine 'surfacehop'")
        end if
        allocate(DC(mol,mol))
        if (.not. allocated(DC)) then
                CPABORT("Error: array 'DC' not allocated in subroutine 'surfacehop'")
                !stop "Error: array 'DC' not allocated in subroutine 'surfacehop'."
        end if

        !use the matinvsqrt subroutine to find T = S^(-1/2)
        T1 = Stmtm
        T2 = Stptp
        call matinvsqrt(mol, .TRUE., T1)
        call matinvsqrt(mol, .TRUE., T2)

        !and do some straightforward first-order numerical derivatives
        DT = (0.5D0/deltat)*(T2 - T1)
        DS = (0.5D0/deltat)*(Sttp - Sttm)
        DC = (1.0D0/deltat)*(C2 - C1)
        
        !with these two loops, calculate the nacv in the time domain between adiabatic wavefunctions
        !by transforming from the diabatic overlaps
        do i=1,mol
                do j=1,mol
                        TST(i,j) = 0.0D0
                        TDT(i,j) = 0.0D0
                        do k=1,mol
                                do l=1,mol
                                        !TST(i,j) = TST(i,j) + T1(i,k)*Stt(k,l)*DT(j,l)
                                        !TDT(i,j) = TDT(i,j) + T1(i,k)*DS(k,l)*T1(j,l)
                                        TST(i,j) = TST(i,j) + T1(k,i)*Stt(k,l)*DT(l,j)
                                        TDT(i,j) = TDT(i,j) + T1(k,i)*DS(k,l)*T1(l,j)
                                end do
                        end do
                end do
        end do
        do i=1,mol
                do j=1,mol
                        dphi(i,j) = 0.0D0
                        do k=1,mol
                                !dphi(i,j) = dphi(i,j) + C1(i,k)*DC(j,k)
                                dphi(i,j) = dphi(i,j) + C2(k,i)*DC(k,j)
                                do l=1,mol
                                        !dphi(i,j) = dphi(i,j) + C1(i,k)*C1(j,l)*(TST(k,l) + TDT(k,l))
                                        dphi(i,j) = dphi(i,j) + C2(k,i)*C2(l,j)*(TST(k,l) + TDT(k,l))
                                end do
                        end do
                end do
        end do

        !now also calculate the adiabatic electronic coefficients by transforming from the diabatic
        !and fill out electronic density matrix elements a(i,j) for i = current adiabatic state
!        cstate = dcmplx(0.0D0, 0.0D0)
        cstate = cmplx(0.0D0, 0.0D0, dp)
        do i=1,mol
                !cstate = cstate + dcmplx(C1(state,i), 0.0D0)*c(i)
                cstate = cstate + cmplx(C2(i,old_state), 0.0D0,dp)*c(i)
        end do
        do i=1,mol
                !calculate a and b
                !ct(i) = dcmplx(0.0D0, 0.0D0)
                ct(i) = cmplx(0.0D0, 0.0D0,dp)
                do j=1,mol
                !        ct(i) = ct(i) + dcmplx(C1(i,j), 0.0D0)*c(j)
                        ct(i) = ct(i) + cmplx(C2(j,i), 0.0D0,dp)*c(j)
                end do
                ak(i) = ct(i)*conjg(cstate)
                b(i) = -2.0D0*REAL(conjg(ak(i))*dphi(i,old_state))
        end do

        sumg = 0.0D0
        !calculate g and sum(g)
        do i=1,mol
                g(i) = (deltat*b(i))/abs(ak(old_state))
                !as in Tully, set negative probabilities to 0
                if (g(i) < 0.0D0) then
                        g(i) = 0.0D0
                end if
                if (i .ne. old_state) then
                        sumg = sumg + g(i)
!                        print *, "Probability: ", i, g(i)
                end if
        end do
        !due to finite timestep, sometimes the probabilities may sum to greater than 1
        !if this happens, rescale them down to sum = 1
        if (sumg > 1.0D0) then
                do i=1,mol
                        g(i) = g(i)/sumg
                        g(old_state) = 0.0D0
                end do
!> CHANGE_AC
                sh_print_basket%renorm = .TRUE.
                sh_print_basket%sumprob   = sumg        
        else
!> CHANGE_AC
                sh_print_basket%renorm = .FALSE.
                sh_print_basket%sumprob   = sumg        
                g(old_state) = 1.0D0 - sumg
        end if


!> CHANGE_AC
        sh_print_basket%proba     = g

        !call a random number
        call init_random_seed()
        call random_number(hop)
        sh_print_basket%hop       = hop
        sh_print_basket%old_state = old_state

        !and see if this random number causes a surface hop
        sumg = 0.0D0
        do i=1,mol
                sumg = sumg + g(i)
                if (hop <= sumg) then
                        if (i .ne. old_state) then
!                                print *, "SURFACE HOP!"
                                new_state = i
                                !aom_basket%state = i
!                                write(24,*) i                        
!PRINT                                print *, "New state: ", i
                                decohere = .FALSE.
                                is_hop = .TRUE.
                                sh_print_basket%is_hop    = is_hop
                                sh_print_basket%new_state = new_state
                                return
                        elseif (i == old_state) then
!                                print *, "No surface hop."
                                new_state = old_state
                                is_hop = .FALSE.
                                sh_print_basket%is_hop    = is_hop
                                sh_print_basket%new_state = new_state
                                !aom_basket%state = state
!                                write(24,*) state
!                                print *, "State: ", state
                                
                                !if there has not been any surfacehop, we need to consider the possibility of the decoherence correction
                                !if the non-adiabatic coupling dphi to all other states is below a certain threshold
                                !then a wavefunction collapse will be enforced
                                do j=1,mol
                                        !if (abs(dphi(old_state,j)) > 1.0D-5 .and. old_state .ne. j) then
                                        if (abs(dphi(old_state,j)) > decohere_criterion .and. old_state .ne. j) then
                                                decohere = .FALSE.
                                                return
                                        end if
                                end do
                                decohere = .TRUE.
                                return
                        end if
                end if
        end do


 

        deallocate(ak)
        deallocate(b)
        deallocate(g)
        deallocate(ct)
        deallocate(dphi)
        deallocate(T1)
        deallocate(T2)
        deallocate(DT)
        deallocate(TST)
        deallocate(TDT)
        deallocate(DS)
        deallocate(DC)
end subroutine surfacehop

subroutine init_random_seed()

        !this subroutine just initialises a random seed for the random_number subroutine

        integer, allocatable :: seed(:)
        !integer :: i, n, un, istat, dt(8), pid, t(2), s
        integer :: i, n, dt(8), pid, t(2), s
        integer(8) :: count, tms
          
        call random_seed(size = n)
        allocate(seed(n))
        ! XOR:ing the current time and pid. The PID is
        ! useful in case one launches multiple instances of the same
        ! program in parallel.
        call system_clock(count)
        if (count /= 0) then
                t = transfer(count, t)
        else
                call date_and_time(values=dt)
                tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                        + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                        + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                        + dt(5) * 60 * 60 * 1000 &
                        + dt(6) * 60 * 1000 + dt(7) * 1000 &
                        + dt(8)
                t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
!        pid = getpid() + 1099279 ! Add a prime
        CALL m_getpid(pid)
        pid = pid + 1099279 ! Add a prime

        s = ieor(s, pid)
        if (n >= 3) then
                seed(1) = t(1) + 36269
                seed(2) = t(2) + 72551
                seed(3) = pid
                if (n > 3) then
                        seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                end if
        else
                seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
        call random_seed(put=seed)

end subroutine init_random_seed

subroutine matinvsqrt(mol, invert, S)

        !calculates either the square root or the inverse of the square root of a square matrix S
        !depending on whether invert is false or true respectively
        !to do this: diagonalise matrix, root eigenvalues, transform back

        implicit none
        
        integer, intent(in) :: mol
        integer :: i, j, stat
        logical, intent(in) :: invert
        REAL(KIND=dp), dimension(:,:), intent(inout) :: S
        REAL(KIND=dp), dimension(:,:), allocatable :: P, Pinv
        REAL(KIND=dp), dimension(:), allocatable :: D, DGETRI_WORK
        character(len=1) :: job, up

        allocate(D(mol))
        if (.not. allocated(D)) then
                CPABORT("Error: array 'D' not allocated correctly.")
                !stop "Error: array 'D' not allocated correctly."
        end if
        allocate(DGETRI_WORK(3*mol))
        if (.not. allocated(DGETRI_WORK)) then
                CPABORT("Error: array 'DGETRI_WORK' not allocated correctly.")
                !stop "Error: array 'DGETRI_WORK' not allocated correctly."
        end if
        allocate(P(mol,mol))
        if (.not. allocated(P)) then
                CPABORT("Error: array 'P' not allocated correctly.")
                !stop "Error: array 'P' not allocated correctly."
        end if
        allocate(Pinv(mol,mol))
        if (.not. allocated(Pinv)) then
                CPABORT("Error: array 'Pinv' not allocated correctly.")
                !stop "Error: array 'Pinv' not allocated correctly."
        end if

        !to calculate inverse square root of matrix, first diagonalise it
        job = 'V'
        up = 'U'
        !DSYEV is from the BLAS. Returns vector of eigenvalues D and overwrites S with matrix of eigenvectors
        call DSYEV(job, up, mol, S, mol, D, DGETRI_WORK, 3*mol, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
                !write(*,*) "Error with subroutine DSYEV: error status", stat
                !stop
        end if

        P = S
        !multiply P by diagonal matrix D^(-1/2)
        do i=1,mol
                if (invert) then
                        !D -> D^(-1/2)
                        D(i) = 1.0D0/sqrt(D(i))
                else
                        !D -> D^(1/2)
                        D(i) = sqrt(D(i))
                end if
                do j=1,mol
                        !Pinv = D*(P     )
                        Pinv(i,j) = D(i)*P(j,i)
                end do
        end do

        !call DGEMM from the BLAS to return S^(-1/2)
        !real symmetric matrix -> S = P*Pinv = P*D*(P     )
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, P, mol, Pinv, mol, 0.0D0, S, mol)

        deallocate(D, P, Pinv, DGETRI_WORK)

end subroutine matinvsqrt

subroutine adiabaticity(H, mol, C, E, print)

        !this subroutine will allow calculation of surface hopping probabilities
        !between N adiabatic surfaces

        !integer :: i, j, stat, LWORK
        integer :: stat, LWORK
        integer, intent(in) :: mol
        REAL(KIND=dp) :: kcal
        REAL(KIND=dp), dimension(:,:), intent(in) :: H
        REAL(KIND=dp), dimension(:,:), intent(out), allocatable :: C
        REAL(KIND=dp), dimension(:), intent(out), allocatable :: E
        character(len=1) :: job, up
        REAL(KIND=dp), dimension(:), allocatable :: DGETRI_WORK
        logical, intent(in) :: print

        !kcal = 1 Hartree in kcal/mol
        kcal = 627.509D0

        allocate(C(mol,mol))
        if (.not. allocated(C)) then
                CPABORT("Error: array 'C' not allocated correctly.")
                !stop "Error: array 'C' not allocated correctly."
        end if
        allocate(E(mol))
        if (.not. allocated(E)) then
                CPABORT("Error: array 'E' not allocated correctly.")
                !stop "Error: array 'E' not allocated correctly."
        end if
        allocate(DGETRI_WORK(3*size(C)-1))
        if (.not. allocated(DGETRI_WORK)) then
                CPABORT("Error: array 'DGETRI_WORK' not allocated correctly.")
                !stop "Error: array 'DGETRI_WORK' not allocated correctly."
        end if

        C = H

        job = 'V'
        up = 'U'
        !DSYEV is from the BLAS. Returns vector of eigenvalues D and overwrites S with matrix of eigenvectors
        !Original call establishes optimum size for speed purposes
        LWORK = -1
        call DSYEV(job, up, mol, C, mol, E, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
                !write(*,*) "Error with subroutine DSYEV: error status", stat
                !stop
        end if
        LWORK = max(3*size(C)-1, int(DGETRI_WORK(1)))

        deallocate(DGETRI_WORK)
        allocate(DGETRI_WORK(LWORK))
        if (.not. allocated(DGETRI_WORK)) then
                CPABORT("Error: array 'DGETRI_WORK' not allocated correctly.")
                !stop "Error: array 'DGETRI_WORK' not allocated correctly."
        end if

        call DSYEV(job, up, mol, C, mol, E, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
                !write(*,*) "Error with subroutine DSYEV: error status", stat
                !stop
        end if

!        if (print) then
!            open(unit=35,file="lambda",action='write')
!            do i=1,mol
!                    write(35,*) kcal*E(i)
!            end do
!            close(35)
!        end if

        deallocate(DGETRI_WORK)

end subroutine adiabaticity

!subroutine forces(state, Hin, C, nacv, NN, atoms, scaling, cbar, Fin, S)
!subroutine forces(state, Hin, C, nacv, NN, atoms, scaling, cbar, Fin, S, ext_forces)
subroutine calcforces(state, Hin, C, nacv, NN, atoms, scaling, cbar, Fin, S, forces_adiab, & 
                  ext_pres, analytics,  exact_forces,ctin)

        !per our full force expression derivation, this calculates the adiabatic forces
        !on each atom

!        interface
!                subroutine calcTDT(S, DS, T, DT, atoms, mol)
!                        REAL(KIND=dp), intent(in), dimension(:,:) :: S
!                        REAL(KIND=dp), intent(in), dimension(:,:,:,:) :: DS
!                        integer, intent(in) :: atoms, mol
!                        REAL(KIND=dp), intent(out), dimension(:,:) :: T
!                        REAL(KIND=dp), intent(out), dimension(:,:,:,:) :: DT
!                end subroutine calcTDT
!        end interface        

        COMPLEX(KIND=dp), DIMENSION(:), INTENT(INOUT), OPTIONAL   :: ctin
        integer, intent(in) :: state, atoms
        LOGICAL, INTENT(IN)  :: analytics
        REAL(KIND=dp), dimension(:,:), INTENT(OUT) :: exact_forces
        REAL(KIND=dp), dimension(:,:), POINTER                  :: forces_adiab
        REAL(KIND=dp), intent(in) :: scaling, cbar
        REAL(KIND=dp), intent(in), dimension(:,:) :: Hin, C, S
        REAL(KIND=dp), intent(in), dimension(:,:,:) :: Fin
        REAL(KIND=dp), intent(in), dimension(:,:,:,:) :: nacv
        integer, intent(in), dimension(:,:) :: NN

        integer :: i, j, k, mat, d
        !integer :: i, j, k, m, n, mat, d
!        integer :: i, j, k, m, mat, d
        REAL(KIND=dp) :: bohrang, kcal
        REAL(KIND=dp), dimension(:,:,:,:), allocatable :: DH, DT, my_DP
        !REAL(KIND=dp), dimension(:,:,:), allocatable :: F
        REAL(KIND=dp), dimension(:,:), allocatable :: F
        REAL(KIND=dp), dimension(:,:), allocatable :: H, T, nsave, DHtemp

        INTEGER                                                 :: natom
        LOGICAL, dimension(:), INTENT(IN)       :: ext_pres
        INTEGER                                                 :: aom_iatom, iatom, num_mol, nadiab, num_j

        natom = SIZE(Fin,2)
        nadiab= SIZE(Fin,3)
!print*, "NADIAB", nadiab
        mat = size(Hin,1)
        ctin = ctin

        !1 Angstrom = bohrang Bohr
        bohrang = 1.889725989D0
        !kcal = 1 Hartree in kcal/mol
        kcal = 627.509D0

        allocate(H(mat,mat))
        if (.not. allocated(H)) CPABORT("Error: array 'H' not allocated correctly.")
        allocate(T(mat,mat))
        if (.not. allocated(T)) CPABORT("Error: array 'T' not allocated correctly.")
        allocate(DHtemp(mat,mat))
        if (.not. allocated(DHtemp)) CPABORT("Error: array 'DHtemp' not allocated correctly.")
        !allocate(DH(mat,mat,mat*atoms,3))
        allocate(DH(mat,mat,natom,3))
        if (.not. allocated(DH)) CPABORT("Error: array 'DH' not allocated correctly.")
        allocate(DT(mat,mat,mat*atoms,3))
        if (.not. allocated(DT)) CPABORT("Error: array 'DT' not allocated correctly.")
        !allocate(my_DP(mat,mat,2*mat*atoms,3))
        !allocate(my_DP(mat,mat,mat*atoms,3)) CPABORT("Error: array 'my_DP' not allocated correctly.")
        allocate(my_DP(mat,mat,natom,3)) 
        if (.not. allocated(my_DP)) CPABORT("Error: array 'my_DP' not allocated correctly.")
        !allocate(F(3,mat,atoms))
        !allocate(F(3,mat,natom))
        allocate(F(3,natom))
!        allocate(exact_forces(3,natom))
        if (.not. allocated(F)) CPABORT("Error: array 'F' not allocated correctly.")
        !allocate(nsave(mat*atoms,3))
        allocate(nsave(natom,3))
        if (.not. allocated(nsave)) CPABORT("Error: array 'nsave' not allocated correctly.")

!        do d=1,3
!                my_DP(:,:,:,d) = nacv(d,:,:,:)
!        end do
        aom_iatom = 0
        do iatom=1, natom
           if (ext_pres(iatom)) then
              aom_iatom = aom_iatom +1
              do d=1,3
                  my_DP(:,:,iatom,d) = nacv(d,:,:,aom_iatom)
               enddo
           else 
              my_DP(:,:,:,iatom) = 0.0D0
           endif
        enddo
        CPASSERT(aom_iatom.EQ.SIZE(nacv,4))

        do i=1,mat
                do j=1,mat
                        if (j > i) then
                                my_DP(i,j,:,:) = my_DP(i,j,:,:) + my_DP(j,i,:,:)
                        else
                                my_DP(i,j,:,:) = my_DP(j,i,:,:)
                        end if
                end do
        end do

        !i is over each site
        do i=1,mat
                do j=1,mat
                        do k=1, natom
                                do d=1,3
                                        if (NN(i,j) == 1) then
                                                if (i == j) then
                                                        DH(i,j,k,d) = Fin(d,k,i)
                                                else
                                                       ! DH(i,j,k,d) = scaling*my_DP(i,j,k,d)/cbar
!> Correct the nacv
                                                        DH(i,j,k,d) = - scaling*my_DP(i,j,k,d)/cbar
                                                end if
                                        else
                                                DH(i,j,k,d) = 0.0D0
                                        end if
                                end do
                        end do
                end do
        end do

        !You can calculate the forces with a "shortcut formula"
        !This should take advantage of the sparsity of H for large systems

!if (shortcut) then
        IF (analytics) THEN
           IF (mat.NE.2) CPABORT("YOU CAN ONLY USE ANALYTICAL EXPRESSION FOR 2x2 SYSTEMS.")
           DO d=1,3
              DO j=1,natom
                 IF (state.EQ.1) THEN
                    exact_forces(d,j) = (DH(1,1,j,d)+DH(2,2,j,d))/2  -  &
                    ((Hin(2,2)-Hin(1,1))*(DH(2,2,j,d)-DH(1,1,j,d)) + 4*Hin(1,2)*DH(1,2,j,d) ) / &
                    (2*SQRT( (Hin(2,2)-Hin(1,1))**2 + 4*Hin(1,2)**2 ))
!print*, "ANALYTICS"
!print*, 4*Hin(1,2) / &
!                    (2*SQRT( (Hin(2,2)-Hin(1,1))**2 + 4*Hin(1,2)**2 ))
                 ELSE 
                    exact_forces(d,j) = (DH(1,1,j,d)+DH(2,2,j,d))/2.0  -  &
                    ((Hin(2,2)-Hin(1,1))*(DH(2,2,j,d)-DH(1,1,j,d)) + 4*Hin(1,2)*DH(1,2,j,d) ) / &
                    2.0*SQRT( (Hin(2,2)-Hin(1,1))**2 + 4*H(1,2)**2 )
                 ENDIF
              ENDDO 
           ENDDO
        ELSE
           exact_forces = 0.0D0
        ENDIF
        !begin shortcut formula 

        !conjugates and realpart commands have been omitted, as everything should be real
        !remember to change this if it ever becomes not the case!
        num_mol = 0
        do d=1,3
                num_j = -1
                do j=1,natom
!print*, "D",d,"J",j
                   F(d,j) = 0.0D0
                   IF (ext_pres(j)) THEN
                      num_j   = num_j + 1
!print*, "NUM_J", num_j
                      num_mol = int(num_j/atoms) + 1 
!                       num_mol = MODULO(num_mol,nadiab) + 1
!print*, "NUM_MOL", num_mol
                   ENDIF
                         do i=1,mat
                                F(d,j) = F(d,j) + (abs(C(i,state))**2)*DH(i,i,j,d)
                                IF (ext_pres(j)) THEN
                                   IF (i.ne.num_mol) THEN
                                      !F(d,j) = F(d,j) + 2.0D0*((C(m,state))*DH(i,m,j,d)*C(i,state))
                                      F(d,j) = F(d,j) + 2.0D0*((C(num_mol,state))*DH(i,num_mol,j,d)*C(i,state))
!print*, "NUM"
!print*, 2.0D0*((C(num_mol,state))*C(i,state))
!print*, 
                                   ENDIF
                                ENDIF
                         ENDDO
                ENDDO
        ENDDO
                !do j=1,atoms
                                !F(d,i,j) = 0.0D0
                                !do m=i,mat
                                        !F(d,i,j) = F(d,i,j) + (abs(C(m,state))**2)*DH(m,m,(i-1)*atoms+j,d)
                                !        if (i.eq. m) then
                                !                F(d,j) = F(d,j) + (abs(C(m,state))**2)*DH(m,m,j,d)
                                !        else
                                                !F(d,i,j) = F(d,i,j) + 2.0D0*((C(i,state))*DH(m,i,(i-1)*atoms+j,d)*C(m,state))
                                !                F(d,j) = F(d,j) + 2.0D0*((C(i,state))*DH(m,i,j,d)*C(m,state))
                                !        end if
                                                !These terms are here for 'historical interest'
                                                !We have shown that they should come out as exactly 0.
                                                !do m=1,mat        
                                                        !F(d,i,j) = F(d,i,j) - 2.0D0*((C(n,state))*H(n,i)*(C(m,state)*my_DP(m,i,j,d)) + &
                                                        !&(C(i,state)*C(n,state)*H(n,m)*my_DP(i,m,j,d)))
                                                !end do
!                                end do
!                        end do
!                end do
!        end do
!rint*, F(1,2)
!else
        !end shortcut formula

        !Alternatively you can just directly multiply the matrices, which works fine for small systems
        !For larger systems this probably becomes prohibitively expensive

        !begin matrix multiplication

!        do i=1,mat
!                do j=1,atoms
!                        do d=1,3
!                                do m=1,mat
!                                        do n=1,mat
!                                                DHtemp(m,n) = DH(m,n,(i-1)*atoms+j,d)
!                                        end do
!                                end do
!                                if (i == 1 .and. j == 1 .and. d == 1) then
!                                        print *, "DH", DHtemp
!                                end if
!                                DHtemp = matmul(DHtemp,C)
!                                DHtemp = matmul(transpose(C),DHtemp)
!                                F(d,i,j) = DHtemp(state,state)
!                                if (i == 1 .and. j == 1 .and. d == 1) then
!                                        print *, "DH", DHtemp
!                                end if
!                        end do
!                end do
!        end do
!endif        
        !end matrix multiplication        
        
        !Remember to re-add force conversion
!        open(unit=50,file='Fx',action='write')
!        open(unit=51,file='Fy',action='write')
!        open(unit=52,file='Fz',action='write')
        forces_adiab = F
        !forces_adiab = exact_forces
        !do i=1,mat
        !        do j=1,natom
!                        write(50,*) kcal*bohrang*F(1,i,j)
!                        write(51,*) kcal*bohrang*F(2,i,j)
!                        write(52,*) kcal*bohrang*F(3,i,j)
        !        end do
        !end do
!        close(50)
!        close(51)
!        close(52)

        deallocate(H, T, DHtemp, DH, DT, my_DP, F, nsave)
        !deallocat(ext_forces)
                
end subroutine calcforces

subroutine calcTDT(S, DS, T, DT, atoms, mol)

!        interface
!                subroutine matinvsqrt(mol, invert, S)
!                        integer, intent(in) :: mol
!                        logical, intent(in) :: invert
!                        REAL(KIND=dp), dimension(:,:), intent(inout) :: S
!                end subroutine matinvsqrt
!        end interface

        REAL(KIND=dp), dimension(:,:), intent(in) :: S
        integer, intent(in) :: atoms, mol
        REAL(KIND=dp), dimension(:,:), intent(out) :: T
        REAL(KIND=dp), dimension(:,:), allocatable :: L, Sinv, Wleft, Wright, Temp
        REAL(KIND=dp), dimension(:), allocatable :: Tdiag
        REAL(KIND=dp), dimension(:,:,:,:), allocatable :: W, Wtop
        REAL(KIND=dp), dimension(:,:,:,:), intent(in) :: DS
        REAL(KIND=dp), dimension(:,:,:,:), intent(out) :: DT
        !integer :: stat, LWORK, i, j, k, d
        integer :: stat, LWORK, i, j
        REAL(KIND=dp), dimension(:), allocatable :: DGETRI_WORK

        !print *, "S", S
        T = S
        call matinvsqrt(mol, .TRUE., T)
        !print *, "T", T

        allocate(W(3,mol,mol,mol*atoms))
        if (.not. allocated(W)) then
                CPABORT("Error: array 'W' not allocated correctly.")
                !stop "Error: array 'W' not allocated correctly."
        end if
        allocate(Wtop(3,mol,mol,mol*atoms))
        if (.not. allocated(Wtop)) then
                CPABORT("Error: array 'Wtop' not allocated correctly.")
                !stop "Error: array 'Wtop' not allocated correctly."
        end if
        allocate(L(mol,mol))
        if (.not. allocated(L)) then
                CPABORT("Error: array 'L' not allocated correctly.")
                !stop "Error: array 'L' not allocated correctly."
        end if
        allocate(Tdiag(mol))
        if (.not. allocated(Tdiag)) then
                CPABORT("Error: array 'Tdiag' not allocated correctly.")
                !stop "Error: array 'Tdiag' not allocated correctly."
        end if
        allocate(Sinv(mol,mol))
        if (.not. allocated(Sinv)) then
                CPABORT("Error: array 'Sinv' not allocated correctly.")
                !stop "Error: array 'Sinv' not allocated correctly."
        end if
        allocate(Wleft(mol,mol))
        if (.not. allocated(Wleft)) then
                CPABORT("Error: array 'Wleft' not allocated correctly.")
                !stop "Error: array 'Wleft' not allocated correctly."
        end if
        allocate(Wright(mol,mol))
        if (.not. allocated(Wright)) then
                CPABORT("Error: array 'Wright' not allocated correctly.")
                !stop "Error: array 'Wright' not allocated correctly."
        end if
        allocate(Temp(mol,mol))
        if (.not. allocated(Temp)) then
                CPABORT("Error: array  'Temp' not allocated correctly.")
                !stop "Error: array 'Temp' not allocated correctly."
        end if

        allocate(DGETRI_WORK(3*mol))
        if (.not. allocated(DGETRI_WORK)) then
                CPABORT("Error: array 'DGETRI_WORK' not allocated correctly.")
                !stop "Error: array 'DGETRI_WORK' not allocated correctly."
        end if

        Sinv = matmul(T,T)
        !print *, "Sinv", Sinv

        L = T        
        !DSYEV is from the BLAS. Returns vector of eigenvalues and overwrites input with matrix of eigenvectors
        LWORK = -1
        call DSYEV('V', 'U', mol, L, mol, Tdiag, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
                !write(*,*) "Error with subroutine DSYEV: error status", stat
                !print *, LWORK
                !stop
        end if
        LWORK = max(3*size(L)-1, int(DGETRI_WORK(1)))
        deallocate(DGETRI_WORK)
        allocate(DGETRI_WORK(LWORK))
        if (.not. allocated(DGETRI_WORK)) then
                CPABORT("Error: array 'DGETRI_WORK not allocated correctly.")
                !stop "Error: array 'DGETRI_WORK' not allocated correctly."
        end if
        call DSYEV('V', 'U', mol, L, mol, Tdiag, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
                !write(*,*) "Error with subroutine DSYEV: error status", stat
                !stop
        end if

        !print *, "T", T
        !print *, "Tdiag", Tdiag
        !print *, "L", L
        Wright = matmul(Sinv, L)
        Wleft = matmul(transpose(L), Sinv)
        !print *, "Wright", Wright
        !print *, "Wleft", Wleft

        Wtop = 0.0D0

        do i=1,3
                do j=1,mol*atoms
                        Wtop(i,:,:,j) = matmul(DS(i,:,:,j),Wright)
                        Wtop(i,:,:,j) = matmul(Wleft, Wtop(i,:,:,j))
                end do
        end do
        do i=1,mol
                do j=1,mol
                        W(:,i,j,:) = -Wtop(:,i,j,:)/(Tdiag(i) + Tdiag(j))
                end do
        end do
        do i=1,3
                do j=1,mol*atoms
                        DT(:,:,j,i) = matmul(W(i,:,:,j),transpose(L))
                        DT(:,:,j,i) = matmul(L,DT(:,:,j,i))
                end do
        end do

        deallocate(W, Wtop, L, Tdiag, Sinv, Wleft, Wright, Temp, DGETRI_WORK)

end subroutine calcTDT

   SUBROUTINE read_listcoeff(listcoeff, ext_elements, ext_pres, &
                            ext_scoeff, ext_pcoeff, natom)
    TYPE(cp_sll_val_type), POINTER                         :: listcoeff
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: ext_elements
    LOGICAL, DIMENSION(:), ALLOCATABLE,  INTENT(OUT)       :: ext_pres
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: ext_scoeff, ext_pcoeff
    INTEGER, INTENT(IN)                                    :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'read_listcoeff', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2*default_string_length)   :: line_att
    CHARACTER(len=default_string_length)     :: error_message
    !CHARACTER(len=default_string_length)     :: error_message, &
    !                                            my_default_index
                                                !unit_str
                                                !my_default_index, strtmp, &
                                                !unit_str
    INTEGER                                  :: default_id, end_c, &
    !INTEGER                                  :: default_id, end_c, handle, &
                                                !iatom, j, output_unit, &
                                                iatom, j,  &
                                                start_c, wrd
    LOGICAL                                  :: is_ok
    TYPE(val_type), POINTER                  :: val
    INTEGER                                  :: my_pres
    REAL(KIND=dp)                            :: my_ext_elements, my_pres2

    ALLOCATE(ext_elements(natom))
    ALLOCATE(ext_pres(natom))
    ALLOCATE(ext_scoeff(natom))
    ALLOCATE(ext_pcoeff(natom))

       DO iatom=1,natom
          ! we use only the first default_string_length characters of each line
          is_ok=cp_sll_val_next(listcoeff,val)
          CALL val_get(val,c_val=line_att)
          default_id = str2id(s2s(""))
          ! Read name and atomic coordinates
          start_c=1
          DO wrd=1,5
             DO j=start_c,LEN(line_att)
                IF (line_att(j:j)/=' ') THEN
                   start_c=j
                   EXIT
                END IF
             END DO
             end_c=LEN(line_att)+1
             DO j=start_c,LEN(line_att)
                IF (line_att(j:j)==' ') THEN
                   end_c=j
                   EXIT
                END IF
             END DO
             IF(LEN_TRIM(line_att(start_c:end_c-1)) == 0)&
                CALL cp_abort(__LOCATION__,&
                     "Incorrectly formatted input line for atom "//&
                     TRIM(ADJUSTL(cp_to_string(iatom)))//&
                     " found in COORD section. Input line: <"//&
                     TRIM(line_att)//"> ")
             SELECT CASE (wrd)
!             CASE (1)
!                 = str2id(s2s(line_att(start_c:end_c-1)))
             CASE (2)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       my_ext_elements,error_message)
                ext_elements(iatom) = INT(my_ext_elements)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (3)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       my_pres2,error_message)
                my_pres = INT(my_pres2)
                IF (my_pres.eq.0) THEN
                   ext_pres(iatom) = .FALSE.
                ELSE
                   ext_pres(iatom) = .TRUE.
                ENDIF
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (4)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       ext_scoeff(iatom),error_message)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (5)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       ext_pcoeff(iatom),error_message)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")

             END SELECT
             start_c = end_c
             IF (start_c > LEN_TRIM(line_att)) EXIT
          END DO
      END DO
   END SUBROUTINE read_listcoeff


   SUBROUTINE choose_firstadiabat(C, init, sh_print_basket, first_adiabat)
        REAL(KIND=dp), dimension(:,:), intent(in) :: C
        INTEGER, INTENT(IN)                       :: init
        TYPE(sh_print_basket_type), INTENT(OUT)   :: sh_print_basket
        INTEGER, INTENT(OUT)                      :: first_adiabat

        CHARACTER(len=*), PARAMETER :: routineN = 'choose_firstadiabat', &
        routineP = moduleN//':'//routineN
        
        INTEGER                                   :: i
        REAL( KIND =dp)                           :: random, sumc

        sh_print_basket%is_first = .TRUE. 
        CALL init_random_seed()
        CALL random_number(random)
        sh_print_basket%initrandom = random
        sumc = 0.0D0
        DO i=1, SIZE(C(init,:))
           sh_print_basket%pop_init(i) = C(init,i)*C(init,i)
        ENDDO
        DO i=1, SIZE(C(init,:))
           sumc = sumc + C(init,i)*C(init,i)
           IF (random.le.sumc) THEN
              first_adiabat = i
              sh_print_basket%first_adiabat = i
              return
           ENDIF
        ENDDO

   END SUBROUTINE choose_firstadiabat

   SUBROUTINE simple_vel_rescaling(velocities_prime, delta_energy, kinetic_energy_prime, &
                                   velocities_plus)
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_prime
        REAL(kind=dp), INTENT(IN)                :: delta_energy
        REAL(kind=dp), INTENT(IN)                :: kinetic_energy_prime
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_plus

        CHARACTER(len=*), PARAMETER :: routineN = 'simple_vel_rescaling', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: rescale_factor
!> Factor to rescale the velocities (Eq. 49 of
!PAPER_1)
        rescale_factor = SQRT(  delta_energy /kinetic_energy_prime )
print*, "rescale_factor 1", rescale_factor
        velocities_plus = rescale_factor*velocities_prime
print*, "velocities_plus"
        
print*, velocities_plus
   END SUBROUTINE simple_vel_rescaling

   SUBROUTINE calc_velocities_prime(velocities, forces_adiab, &
                                    atomic_kinds, deltatMD,   &
                                    local_particles,           &
                                    velocities_prime, kinetic_energy_prime)

        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: forces_adiab
        TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
        REAL(kind=dp), INTENT(IN)                :: deltatMD
        TYPE(distribution_1d_type), POINTER      :: local_particles
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_prime
        REAL(kind=dp), INTENT(OUT)                :: kinetic_energy_prime

        CHARACTER(len=*), PARAMETER :: routineN = 'simple_vel_rescaling', &
        routineP = moduleN//':'//routineN

        INTEGER                                   :: iparticle_kind, nparticle_kind
        REAL(kind=dp)                             :: dm,dt,mass
        INTEGER                                   :: nparticle_local, iparticle_local
        TYPE(atomic_kind_type), DIMENSION(:), &
          POINTER                                :: atomic_kind_set
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        INTEGER                                  :: iparticle
    
    dt = deltatMD
    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    kinetic_energy_prime = 0.0D0
    ! Update Verlet (second part)
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       dm = 0.5_dp*dt/mass
       nparticle_local = local_particles%n_el(iparticle_kind)
!print*, "dm = ", dm
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
!print*, "velocities:", iparticle
!print*, velocities( 1, iparticle)
!print*, velocities( 2, iparticle)
!print*, velocities( 3, iparticle)
             velocities_prime(1, iparticle) = velocities(1, iparticle) + &
                 dm*forces_adiab(1,iparticle)
             velocities_prime(2, iparticle) = velocities(2, iparticle) + &
                 dm*forces_adiab(2,iparticle)
             velocities_prime(3, iparticle) = velocities(3,iparticle) + &
                 dm*forces_adiab(3,iparticle)
!print*, "forces_adiab", iparticle
!print*, forces_adiab(1, iparticle)
!print*, forces_adiab(2, iparticle)
!print*, forces_adiab(3, iparticle)
!print*, "velocities_prime", iparticle
!print*, velocities_prime( 1, iparticle)
!print*, velocities_prime( 2, iparticle)
!print*, velocities_prime( 3, iparticle)
             kinetic_energy_prime = kinetic_energy_prime +   &
                 mass*( &
                 velocities_prime(1, iparticle)**2 + &
                 velocities_prime(2, iparticle)**2 + &
                 velocities_prime(3, iparticle)**2  )
       END DO
    END DO
    kinetic_energy_prime = 0.5*kinetic_energy_prime

   END SUBROUTINE calc_velocities_prime

END MODULE jacob_utils
