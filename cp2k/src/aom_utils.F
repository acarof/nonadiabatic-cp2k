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

MODULE aom_utils
   USE kinds,                          ONLY: default_string_length,&
                                               dp
   USE machine,                        ONLY: m_getpid
   USE aom_types,                      ONLY: sh_print_env_type

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

  USE mixed_environment_utils,         ONLY: get_subsys_map_index
  USE mixed_environment_types,         ONLY: mixed_environment_type
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

  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
!  USE integrator_utils,                ONLY: tmp_variables_type 

#include "./base/base_uses.f90"



IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'aom_utils'

   PUBLIC :: makeX,            &
             surfacehop,       &
             adiabaticity,     &
             calcforces,           &
             read_listcoeff,       & 
             choose_firstadiabat,  &
             simple_vel_rescaling, &
             check_hop_attempt

CONTAINS
  
subroutine makeX(mol, H, Stt, Sttp, Sttm, Stptp, Stmtm, X, deltat, NACE)
        !X is the name I give to the matrix in the solution to the Schroedinger equation:
        !dc/dt = X*c
        !makeX constructs an X matrix from inputs of S and H

        implicit none

        integer, intent(in) :: mol
        REAL(KIND=dp), intent(in) :: deltat
        REAL(KIND=dp), dimension(:,:), intent(in) :: Stt, Sttp, Sttm, Stptp, Stmtm
        REAL(KIND=dp), dimension(:,:), allocatable :: Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt
        REAL(KIND=dp), dimension(:,:), allocatable :: Sterm1, Sterm2, Splace1, Splace2
        COMPLEX(KIND=dp), dimension(:,:), intent(in) :: H
        COMPLEX(KIND=dp), dimension(:,:), intent(out) :: X
        REAL(KIND=dp), dimension(:,:), intent(out), optional    :: NACE   

        REAL(KIND=dp), dimension(:,:), allocatable     :: my_NACE   

        !allocate matrices
        allocate(Sttinvsqrt(mol,mol))
        allocate(Sttsqrt(mol,mol))
        allocate(Stptpinvsqrt(mol,mol))
        allocate(Stmtminvsqrt(mol,mol))
        allocate(Sterm1(mol,mol))
        allocate(Sterm2(mol,mol))
        allocate(Splace1(mol,mol))
        allocate(Splace2(mol,mol))
        allocate(my_NACE(mol,mol))

        Sttsqrt = Stt
        call matinvsqrt(mol, .FALSE., Sttsqrt)
        Sttinvsqrt = Stt
        call matinvsqrt(mol, .TRUE., Sttinvsqrt)
        Stptpinvsqrt = Stptp
        call matinvsqrt(mol, .TRUE., Stptpinvsqrt)
        Stmtminvsqrt = Stmtm
        call matinvsqrt(mol, .TRUE., Stmtminvsqrt)

        !subroutine DGEMM multiplies two matrices, from the BLAS
        !use to build the S terms in the X matrix

        !first term: Sttsqrt*(Stptpinvsqrt - Stmtmtinvsqrt)
        Splace1 = Stptpinvsqrt - Stmtminvsqrt
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Sttsqrt, mol, Splace1, mol, 0.0D0, Sterm1, mol)

        !second term: Sttinvsqrt*(Sttp - Sttm)*Sttinvsqrt
        Splace1 = Sttp - Sttm
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Splace1, mol, Sttinvsqrt, mol, 0.0D0, Splace2, mol)
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Sttinvsqrt, mol, Splace2, mol, 0.0D0, Sterm2, mol)
!ADD: SYMMETRIZE NACE
        my_NACE = (0.5D0/deltat)*(Sterm1 + Sterm2) 
        my_NACE = (my_NACE - TRANSPOSE(my_NACE)) / 2
        !X = H - (0.5D0/deltat)*(Sterm1 + Sterm2)
        X = H - my_NACE

!> CHANGE_AC
        IF (PRESENT(NACE)) THEN
           NACE = my_NACE 
        ENDIF

        deallocate(Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt, Sterm1, Sterm2, Splace1, Splace2)

end subroutine makeX


subroutine surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, c, C1, C2, deltat, old_state, & 
                      decohere, new_state, is_hop, decohere_criterion, sh_print_env)

        !using the method outlined in John Tully's Electronic Transitions paper
        !calculate the probability of the molecular dynamics switching to evolve on a different potential surface

        implicit none

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
        TYPE(sh_print_env_type), INTENT(OUT) :: sh_print_env

        mol = size(Stt,1)

        !many allocations
        allocate(ak(mol))
        allocate(b(mol))
        allocate(g(mol))
        allocate(ct(mol))
        allocate(dphi(mol,mol))
        allocate(T1(mol,mol))
        allocate(T2(mol,mol))
        allocate(DT(mol,mol))
        allocate(TST(mol,mol))
        allocate(TDT(mol,mol))
        allocate(DS(mol,mol))
        allocate(DC(mol,mol))

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
                                dphi(i,j) = dphi(i,j) + C2(k,i)*DC(k,j)
                                do l=1,mol
                                        dphi(i,j) = dphi(i,j) + C2(k,i)*C2(l,j)*(TST(k,l) + TDT(k,l))
                                end do
                        end do
                end do
        end do

        !now also calculate the adiabatic electronic coefficients by transforming from the diabatic
        !and fill out electronic density matrix elements a(i,j) for i = current adiabatic state
        cstate = cmplx(0.0D0, 0.0D0, dp)
        do i=1,mol
                cstate = cstate + cmplx(C2(i,old_state), 0.0D0,dp)*c(i)
        end do
        do i=1,mol
                !calculate a and b
                ct(i) = cmplx(0.0D0, 0.0D0,dp)
                do j=1,mol
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
                sh_print_env%renorm = .TRUE.
                sh_print_env%sumprob   = sumg        
        else
!> CHANGE_AC
                sh_print_env%renorm = .FALSE.
                sh_print_env%sumprob   = sumg        
                g(old_state) = 1.0D0 - sumg
        end if


!> CHANGE_AC
        sh_print_env%proba     = g

        !call a random number
        call init_random_seed()
        call random_number(hop)
        sh_print_env%hop       = hop
        sh_print_env%old_state = old_state

        !and see if this random number causes a surface hop
        sumg = 0.0D0
        do i=1,mol
                sumg = sumg + g(i)
                if (hop <= sumg) then
                        if (i .ne. old_state) then
                                new_state = i
!                                write(24,*) i                        
                                decohere = .FALSE.
                                is_hop = .TRUE.
                                sh_print_env%is_hop    = is_hop
                                sh_print_env%new_state = new_state
                                return
                        elseif (i == old_state) then
                                new_state = old_state
                                is_hop = .FALSE.
                                sh_print_env%is_hop    = is_hop
                                sh_print_env%new_state = new_state
!                                write(24,*) state
                                
                                !if there has not been any surfacehop, we need to consider the possibility of the decoherence correction
                                !if the non-adiabatic coupling dphi to all other states is below a certain threshold
                                !then a wavefunction collapse will be enforced
                                do j=1,mol
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
        allocate(DGETRI_WORK(3*mol))
        allocate(P(mol,mol))
        allocate(Pinv(mol,mol))

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

        integer :: stat, LWORK
        integer, intent(in) :: mol
        !REAL(KIND=dp) :: kcal
        REAL(KIND=dp), dimension(:,:), intent(in) :: H
        REAL(KIND=dp), dimension(:,:), intent(out), allocatable :: C
        REAL(KIND=dp), dimension(:), intent(out), allocatable :: E
        character(len=1) :: job, up
        REAL(KIND=dp), dimension(:), allocatable :: DGETRI_WORK
        logical, intent(in) :: print

        !kcal = 1 Hartree in kcal/mol
        !kcal = 627.509D0

        allocate(C(mol,mol))
        allocate(E(mol))
        allocate(DGETRI_WORK(3*size(C)-1))

        C = H

        job = 'V'
        up = 'U'
        !DSYEV is from the BLAS. Returns vector of eigenvalues D and overwrites S with matrix of eigenvectors
        !Original call establishes optimum size for speed purposes
        LWORK = -1
        call DSYEV(job, up, mol, C, mol, E, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
        end if
        LWORK = max(3*size(C)-1, int(DGETRI_WORK(1)))

        deallocate(DGETRI_WORK)
        allocate(DGETRI_WORK(LWORK))

        call DSYEV(job, up, mol, C, mol, E, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
        end if


        deallocate(DGETRI_WORK)

end subroutine adiabaticity

subroutine calcforces(state, Hin, C, nacv, NN, atoms, scaling, cbar, Fin, S, forces_adiab, & 
                  ext_pres, analytics,  exact_forces,ctin)

        !per our full force expression derivation, this calculates the adiabatic forces
        !on each atom


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
        !REAL(KIND=dp) :: bohrang, kcal
        REAL(KIND=dp), dimension(:,:,:,:), allocatable :: DH, DT, my_DP
        REAL(KIND=dp), dimension(:,:), allocatable :: F
        REAL(KIND=dp), dimension(:,:), allocatable :: H, T, nsave, DHtemp

        INTEGER                                                 :: natom
        LOGICAL, dimension(:), INTENT(IN)       :: ext_pres
        INTEGER                                                 :: aom_iatom, iatom, num_mol, nadiab, num_j

        natom = SIZE(Fin,2)
        nadiab= SIZE(Fin,3)
        mat = size(Hin,1)
        ctin = ctin

        !1 Angstrom = bohrang Bohr
        !bohrang = 1.889725989D0
        !kcal = 1 Hartree in kcal/mol
        !kcal = 627.509D0

        allocate(H(mat,mat))
        allocate(T(mat,mat))
        allocate(DHtemp(mat,mat))
        allocate(DH(mat,mat,natom,3))
        allocate(DT(mat,mat,mat*atoms,3))
        allocate(my_DP(mat,mat,natom,3)) 
        allocate(F(3,natom))
        allocate(nsave(natom,3))

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

        IF (analytics) THEN
           IF (mat.NE.2) CPABORT("YOU CAN ONLY USE ANALYTICAL EXPRESSION FOR 2x2 SYSTEMS.")
           DO d=1,3
              DO j=1,natom
                 IF (state.EQ.1) THEN
                    exact_forces(d,j) = (DH(1,1,j,d)+DH(2,2,j,d))/2  -  &
                    ((Hin(2,2)-Hin(1,1))*(DH(2,2,j,d)-DH(1,1,j,d)) + 4*Hin(1,2)*DH(1,2,j,d) ) / &
                    (2*SQRT( (Hin(2,2)-Hin(1,1))**2 + 4*Hin(1,2)**2 ))
                 ELSE 
                    exact_forces(d,j) = (DH(1,1,j,d)+DH(2,2,j,d))/2.0  +  &
                    ((Hin(2,2)-Hin(1,1))*(DH(2,2,j,d)-DH(1,1,j,d)) + 4*Hin(1,2)*DH(1,2,j,d) ) / &
                    (2*SQRT( (Hin(2,2)-Hin(1,1))**2 + 4*Hin(1,2)**2 ))
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
                   F(d,j) = 0.0D0
                   IF (ext_pres(j)) THEN
                      num_j   = num_j + 1
                      num_mol = int(num_j/atoms) + 1 
                   ENDIF
                         do i=1,mat
                                F(d,j) = F(d,j) + (abs(C(i,state))**2)*DH(i,i,j,d)
                                IF (ext_pres(j)) THEN
                                   IF (i.ne.num_mol) THEN
                                      F(d,j) = F(d,j) + 2.0D0*((C(num_mol,state))*DH(i,num_mol,j,d)*C(i,state))
                                   ENDIF
                                ENDIF
                         ENDDO
                ENDDO
        ENDDO
        forces_adiab = F

        deallocate(H, T, DHtemp, DH, DT, my_DP, F, nsave)
                
end subroutine calcforces

subroutine calcTDT(S, DS, T, DT, atoms, mol)


        REAL(KIND=dp), dimension(:,:), intent(in) :: S
        integer, intent(in) :: atoms, mol
        REAL(KIND=dp), dimension(:,:), intent(out) :: T
        REAL(KIND=dp), dimension(:,:), allocatable :: L, Sinv, Wleft, Wright, Temp
        REAL(KIND=dp), dimension(:), allocatable :: Tdiag
        REAL(KIND=dp), dimension(:,:,:,:), allocatable :: W, Wtop
        REAL(KIND=dp), dimension(:,:,:,:), intent(in) :: DS
        REAL(KIND=dp), dimension(:,:,:,:), intent(out) :: DT
        integer :: stat, LWORK, i, j
        REAL(KIND=dp), dimension(:), allocatable :: DGETRI_WORK

        T = S
        call matinvsqrt(mol, .TRUE., T)

        allocate(W(3,mol,mol,mol*atoms))
        allocate(Wtop(3,mol,mol,mol*atoms))
        allocate(L(mol,mol))
        allocate(Tdiag(mol))
        allocate(Sinv(mol,mol))
        allocate(Wleft(mol,mol))
        allocate(Wright(mol,mol))
        allocate(Temp(mol,mol))

        allocate(DGETRI_WORK(3*mol))

        Sinv = matmul(T,T)

        L = T        
        !DSYEV is from the BLAS. Returns vector of eigenvalues and overwrites input with matrix of eigenvectors
        LWORK = -1
        call DSYEV('V', 'U', mol, L, mol, Tdiag, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
        end if
        LWORK = max(3*size(L)-1, int(DGETRI_WORK(1)))
        deallocate(DGETRI_WORK)
        allocate(DGETRI_WORK(LWORK))
        if (.not. allocated(DGETRI_WORK)) then
                CPABORT("Error: array 'DGETRI_WORK not allocated correctly.")
        end if
        call DSYEV('V', 'U', mol, L, mol, Tdiag, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
        end if

        Wright = matmul(Sinv, L)
        Wleft = matmul(transpose(L), Sinv)

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
    INTEGER                                  :: default_id, end_c, &
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


   SUBROUTINE choose_firstadiabat(C, init, sh_print_env, first_adiabat)
        REAL(KIND=dp), dimension(:,:), intent(in) :: C
        INTEGER, INTENT(IN)                       :: init
        TYPE(sh_print_env_type), INTENT(OUT)   :: sh_print_env
        INTEGER, INTENT(OUT)                      :: first_adiabat

        CHARACTER(len=*), PARAMETER :: routineN = 'choose_firstadiabat', &
        routineP = moduleN//':'//routineN
        
        INTEGER                                   :: i
        REAL( KIND =dp)                           :: random, sumc

        sh_print_env%is_first = .TRUE. 
        CALL init_random_seed()
        CALL random_number(random)
        sh_print_env%initrandom = random
        sumc = 0.0D0
        DO i=1, SIZE(C(init,:))
           sh_print_env%pop_init(i) = C(init,i)*C(init,i)
        ENDDO
        DO i=1, SIZE(C(init,:))
           sumc = sumc + C(init,i)*C(init,i)
           IF (random.le.sumc) THEN
              first_adiabat = i
              sh_print_env%first_adiabat = i
              return
           ENDIF
        ENDDO

   END SUBROUTINE choose_firstadiabat

   SUBROUTINE simple_vel_rescaling(velocities, delta_energy, kinetic_energy_prime)
        !TYPE(tmp_variables_type), POINTER        :: tmp
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        REAL(kind=dp), INTENT(IN)                :: delta_energy
        REAL(kind=dp), INTENT(IN)                :: kinetic_energy_prime
        !REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_plus

        CHARACTER(len=*), PARAMETER :: routineN = 'simple_vel_rescaling', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: rescale_factor

!> Factor to rescale the velocities (Eq. 49 of PAPER_1)
        rescale_factor = SQRT(  delta_energy /kinetic_energy_prime )
        print*, "rescale_factor = ", rescale_factor
        velocities = rescale_factor*velocities

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
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             velocities_prime(1, iparticle) = velocities(1, iparticle) + &
                 dm*forces_adiab(1,iparticle)
             velocities_prime(2, iparticle) = velocities(2, iparticle) + &
                 dm*forces_adiab(2,iparticle)
             velocities_prime(3, iparticle) = velocities(3,iparticle) + &
                 dm*forces_adiab(3,iparticle)
             kinetic_energy_prime = kinetic_energy_prime +   &
                 mass*( &
                 velocities_prime(1, iparticle)**2 + &
                 velocities_prime(2, iparticle)**2 + &
                 velocities_prime(3, iparticle)**2  )
       END DO
    END DO
    kinetic_energy_prime = 0.5*kinetic_energy_prime

   END SUBROUTINE calc_velocities_prime


   SUBROUTINE check_hop_attempt(velocities, atomic_kind_set, nparticle_kind, local_particles, force_env)
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        !TYPE(tmp_variables_type), POINTER        :: tmp
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        INTEGER, INTENT(IN)                      :: nparticle_kind
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(force_env_type), POINTER            :: force_env


        CHARACTER(len=*), PARAMETER :: routineN = 'check_hop_attempt', &
        routineP = moduleN//':'//routineN

        INTEGER                                  :: iparticle_kind
        REAL(kind=dp)                            :: mass, kinetic_energy_prime
        REAL(kind=dp)                            :: delta_energy
        INTEGER                                  :: nparticle_local, iparticle_local
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        INTEGER                                  :: iparticle
        CHARACTER(LEN=default_string_length)     :: rescale_keyword
        INTEGER                                  :: old_state, new_state
        TYPE(section_vals_type), POINTER         :: mixed_section, adiab_section, coupling_section
        TYPE(section_vals_type), POINTER         :: force_env_section

        CALL force_env_get(force_env=force_env,&
                       force_env_section=force_env_section)
        mixed_section => section_vals_get_subs_vals(force_env_section,"MIXED")
        adiab_section => section_vals_get_subs_vals(mixed_section, "ADIABATIC")
        coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
        CALL section_vals_val_get(coupling_section,"METHOD_RESCALING",c_val=rescale_keyword)
        SELECT CASE(rescale_keyword)
        CASE("SIMPLE")
           kinetic_energy_prime = 0.0D0
           ! Update Verlet (second part)
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 kinetic_energy_prime = kinetic_energy_prime +   &
                    mass*( &
                    velocities(1, iparticle)**2 + &
                    velocities(2, iparticle)**2 + &
                    velocities(3, iparticle)**2  )
                    !tmp%vel(3, iparticle)**2  )
              END DO
           END DO
           kinetic_energy_prime = 0.5*kinetic_energy_prime
    
           new_state = force_env%mixed_env%aom_env%new_state    
           old_state = force_env%mixed_env%aom_env%active_state    
           delta_energy = kinetic_energy_prime-( &
                     force_env%mixed_env%aom_env%adiab_energies_t(new_state) - &
                     force_env%mixed_env%aom_env%adiab_energies_t(old_state))

print*, "kinetic_enrgy_prime = ", kinetic_energy_prime
print*, "energy_old_state = ", force_env%mixed_env%aom_env%adiab_energies_t(old_state)
print*, "energy_new_state = ", force_env%mixed_env%aom_env%adiab_energies_t(new_state)

           IF (delta_energy.lt.0) THEN
              force_env%mixed_env%aom_env%successful_hop = .FALSE.
           ELSE 
              force_env%mixed_env%aom_env%successful_hop = .TRUE.
              force_env%mixed_env%aom_env%active_state   = force_env%mixed_env%aom_env%new_state
              force_env%mixed_env%aom_env%new_state = -1
              CALL simple_vel_rescaling(velocities, delta_energy, kinetic_energy_prime)
           ENDIF
        CASE DEFAULT
                CPABORT("There is no method to do the rescaling of the velocities.")
        END SELECT


   END SUBROUTINE check_hop_attempt


END MODULE aom_utils
