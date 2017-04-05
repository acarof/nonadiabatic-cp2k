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
   USE aom_print,                      ONLY: adiab_print_sh_info,  &
                                             adiab_print_hop_nacv

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
!  USE kinds,                           ONLY: default_string_length,&
!                                            dp
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

  USE force_env_types,                 ONLY: force_env_get, &
                                             force_env_type

  USE non_adiabatic_types,             ONLY: sh_env_type

!  USE integrator_utils,                ONLY: tmp_variables_type 

#include "./base/base_uses.f90"



IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'aom_utils'

   PUBLIC :: makeX,            &
!             surfacehop,       &
             adiabaticity,     &
             calcforces,           &
             read_listcoeff,       & 
             choose_firstadiabat,  &
             simple_vel_rescaling, &
             check_hop_attempt,    &
             total_momentum,       &
!             translate_com
             translate_com,       &
             surf_hop_type

CONTAINS
  
subroutine makeX(mol, H, Stt, Sttp, Sttm, Stptp, Stmtm, X, deltat, NACE)
        !X is the name I give to the matrix in the solution to the Schroedinger equation:
        !dc/dt = X*sh_env%ctt
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
!           print *, 'LOOK NACE=', NACE 
        ENDIF

        deallocate(Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt, Sterm1, Sterm2, Splace1, Splace2)

end subroutine makeX


subroutine surf_hop_type(adiab_section, Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2, deltat, old_state, &
           sh_env, new_state, attempt_hop, decohere_criterion, sh_print_env)

!Choose whether doing surface hopping with trivial hopping correction or witgout
!it following the old code 
        
       implicit none
       TYPE(section_vals_type), POINTER               :: adiab_section
       REAL(KIND=dp), INTENT(IN) :: decohere_criterion
       integer, intent(in) :: old_state
       REAL(KIND=dp), intent(in) :: deltat
       REAL(KIND=dp), dimension(:,:), intent(in) :: Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2
       INTEGER, INTENT(OUT)        :: new_state
       LOGICAL, INTENT(OUT)        :: attempt_hop
       TYPE(sh_print_env_type)     :: sh_print_env
       TYPE(sh_env_type)           :: sh_env
       !COMPLEX(KIND=dp), dimension(:), intent(in) :: sh_env%ctt
       
       CHARACTER(len=*), PARAMETER :: routineN ='surf_hop_type', &
                               routineP = moduleN//':'//routineN

       CHARACTER(LEN=default_string_length)           :: surfhop_keyword
       REAL(KIND=dp) :: triv_threshold


       CALL section_vals_val_get(adiab_section,"SURF_HOP_CHOICE", c_val=surfhop_keyword)

       SELECT CASE(surfhop_keyword)
       CASE("TRIVIAL_HOP_CORRECT")

       CALL section_vals_val_get(adiab_section,"THRESHOLD_TRIVIAL_CROSS", r_val=triv_threshold)
       CALL modified_surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2, deltat, old_state, &
                        sh_env, new_state, attempt_hop, decohere_criterion, &
                       sh_print_env, triv_threshold)

PRINT *, 'TO BE CONTINUED'
       

       CASE("UNMODIFIED_SURF_HOP")
        
       CALL surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2, deltat, old_state, &
                        sh_env, new_state, attempt_hop, decohere_criterion, sh_print_env)
       
print *, "OKKKKKKKKKKK????????"

       CASE DEFAULT
         CPABORT("There is no method to add trivial crossing correction")
       
       END SELECT

END SUBROUTINE surf_hop_type

subroutine modified_surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2, deltat, old_state, &
              sh_env, new_state, attempt_hop, decohere_criterion, &
              sh_print_env, triv_threshold)

        !using the method outlined in John Tully's Electronic Transitions paper
        !calculate the probability of the molecular dynamics switching to evolve on a different potential surface

        implicit none

        REAL(KIND=dp), INTENT(IN) :: decohere_criterion
        integer, intent(in) :: old_state
        REAL(KIND=dp), intent(in) :: deltat
        REAL(KIND=dp), dimension(:,:), intent(in) :: Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2
        INTEGER, INTENT(OUT)        :: new_state
        LOGICAL, INTENT(OUT)        :: attempt_hop
        TYPE(sh_print_env_type), INTENT(OUT) :: sh_print_env
        TYPE(sh_env_type), INTENT(OUT) :: sh_env

        LOGICAL        :: trivial_correction, correct_for_negative
        integer :: i, j, mol, k, l
        REAL(KIND=dp) :: hop, sumg, ex_prob, min, g_exa, err
        REAL(KIND=dp) :: triv_threshold
        REAL(KIND=dp), dimension(:,:), allocatable :: dphi
        REAL(KIND=dp), dimension(:), allocatable :: b, g
        REAL(KIND=dp), dimension(:,:), allocatable :: T1, T2, DT, TST, TDT, DS, DC
       ! COMPLEX(KIND=dp), dimension(:), intent(in) :: sh_env%ctt
        COMPLEX(KIND=dp), dimension(:), allocatable :: ct, cttm_ad
        COMPLEX(KIND=dp), dimension(:), allocatable :: ak, akttm
        COMPLEX(KIND=dp) :: cstate, cstate_tm

        mol = size(Stt,1)

        !many allocations
        allocate(ak(mol))
        allocate(b(mol))
        allocate(g(mol))
        allocate(ct(mol))
        allocate(cttm_ad(mol))
        allocate(akttm(mol))        

        allocate(dphi(mol,mol))
        allocate(T1(mol,mol))
        allocate(T2(mol,mol))
        allocate(DT(mol,mol))
        allocate(TST(mol,mol))
        allocate(TDT(mol,mol))
        allocate(DS(mol,mol))
        allocate(DC(mol,mol))
        allocate(sh_env%enerdiff(mol)) 
        
        ALLOCATE(sh_print_env%proba(mol))

print *, "OKAAYYYYY"

!print *, "ADIAB ENERGY t-dt", sh_env%E_t
!print *, "A_EN t-dt OLD_STATE", sh_env%E_t(old_state)
print *, "ADIAB ENERGY t", sh_env%E_tp

        do i=1,mol
                if (i .ne. old_state) then 
                        sh_env%enerdiff(i) = abs(sh_env%E_tp(i) - sh_env%E_tp(old_state)) 
                else if (i .eq. old_state) then 
                        sh_env%enerdiff(old_state) = 0.0D0
                end if
        end do
print *, "sh_env%enerdiff = ", sh_env%enerdiff
print *, "old_state", old_state
       
        if (sh_env%enerdiff(1) == sh_env%enerdiff(old_state)) then
                min = sh_env%enerdiff(2)
                do i= 3,mol
                       if (sh_env%enerdiff(i) < min .and. i .ne. old_state) then
                               min = sh_env%enerdiff(i)
!                               sh_env%state_min = i
                        end if          
                end do
        else 
                min = sh_env%enerdiff(1)
                do i= 2,mol
                       if (sh_env%enerdiff(i) < min .and. i .ne. old_state) then
                               min = sh_env%enerdiff(i)
!                               sh_env%state_min = i
                        end if          
                end do
        end if


        do i =1,mol
                if (sh_env%enerdiff(i) == min) then
                       sh_env%state_min = i
                end if
        end do

print *, "MIN_ENER_DIFF", min
print *,"state_min", sh_env%state_min


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
        cstate_tm = cmplx(0.0D0, 0.0D0, dp)
        do i=1,mol
                cstate = cstate + cmplx(C2(i,old_state), 0.0D0,dp)*sh_env%ctt(i)
                cstate_tm = cstate_tm + cmplx(C1(i,old_state),0.0D0,dp)*sh_env%cttm(i)
        end do
print *, "cstate", cstate
print *, "cstate_tm", cstate_tm
        do i=1,mol
                !calculate a and b
                ct(i) = cmplx(0.0D0, 0.0D0,dp)
                cttm_ad(i) = cmplx(0.0D0, 0.0D0,dp)
                do j=1,mol
                        ct(i) = ct(i) + cmplx(C2(j,i), 0.0D0,dp)*sh_env%ctt(j)
                        cttm_ad(i) = cttm_ad(i) + cmplx(C1(j,i),0.0D0,dp)*sh_env%cttm(j)
                end do
!construct state density for t present and t-deltat from the active state to
!other state
print *, "ct(i)", ct(i)
print *, "cttm_ad(i)", cttm_ad(i)
                ak(i) = ct(i)*conjg(cstate)
                akttm(i) = cttm_ad(i)*conjg(cstate_tm)
print *, "ak(i)", ak(i)
print *, "akttm(i)", akttm(i)

                b(i) = -2.0D0*REAL(conjg(ak(i))*dphi(i,old_state))
!print*, "i b dphi", i, b(i), dphi(i, old_state)
        end do


!calc exact total surface hopping probability (also termed g_tot) according to Prezhdo and Wang 2014
        ex_prob = (abs(akttm(old_state))- abs(ak(old_state)))/abs(akttm(old_state))

!this probability can be negative when there is flux toward the current state
!meaning more probability to stay on the current state

        !if (ex_prob < 0.0D0) then
         !   ex_prob = 0.0D0
        !end if

print *, "ak(old_state)", (ak(old_state))
print *, "abs(akttm(old_state))", (akttm(old_state))
print *, "ABS(ak(old_state))", abs(ak(old_state))
print *, "ABS(akttm(old_state))", abs(akttm(old_state))
print *, "ex_prob =", ex_prob

! calculate sumg that can also be negative because each individual probability
! could be negative 

        sumg = 0.0D0
        !calculate g and sum(g)
        do i=1,mol
!print*, "deltat b ak", i, deltat, b(i), abs(ak(old_state))
!g(i) is the probability to hop from the current state (old_state) to any other
!state (i) 
                g(i) = (deltat*b(i))/abs(ak(old_state))
                
               ! if (g(i) < 0.0D0) then
               !         g(i) = 0.0D0
               ! end if

                if (i .ne. old_state) then
                        sumg = sumg + g(i)
                end if
        end do

print*, "sumg_tot", sumg
print *, "DIFF_SUMG_EXACT", ex_prob - sumg
print *, "Proba_gi_before_corrections: ", g

!EFFECTIVE TRIVIAL CROSSINGS CORRECTION

        trivial_correction = .FALSE.
        
        err = abs(sumg - ex_prob)
        if (err >= triv_threshold) then 
                trivial_correction = .TRUE.
                        if (trivial_correction) then  
                                g_exa = ex_prob - sumg + g(sh_env%state_min)
!                if (g_exa < 0.0D0) then
!                        g_exa = 0.0D0
!                end if
                                g(sh_env%state_min) = g_exa
!!                do i=1,mol
!!                        if (i == sh_env%state_min) then
!!                                g(i) = g_exa
!!                        end if
!!                end do 
                        end if
        end if
print *, "trivial_correction", trivial_correction
print *, "Prob_gi_CORRECT_for_trivial: ", g

!! CORRECTION FOR NEGATIVE NUMBER
       correct_for_negative = .FALSE.
       do i =1,mol
                if (g(i) < 0.0D0) then
                        g(i) = 0.0D0
                        correct_for_negative = .TRUE.
                end if
       end do
print *, "correct_for_negative", correct_for_negative
print *, "Prob_gi_CORRECT_for_NEGA: ", g

!! G RENORMALIZATION; BUILT A SUBROUTINE THAT DOES THIS
       sumg = 0.0D0
       do i=1,mol
               if (i .ne. old_state) then
                       sumg = sumg + g(i)
               end if
       end do
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





!        !due to finite timestep, sometimes the probabilities may sum to greater than 1
!        !if this happens, rescale them down to sum = 1
!        if (sumg > 1.0D0) then
!                do i=1,mol
!                        g(i) = g(i)/sumg
!                        g(old_state) = 0.0D0
!                end do
!!> CHANGE_AC
!                sh_print_env%renorm = .TRUE.
!                sh_print_env%sumprob   = sumg        
!        else
!!> CHANGE_AC
!                sh_print_env%renorm = .FALSE.
!                sh_print_env%sumprob   = sumg        
!                g(old_state) = 1.0D0 - sumg
!        end if

!print *, "sh_print_env%renorm", sh_print_env%renorm
!print *, "Prob_gi_RESCALING1: ", g  !this case assigns also g_old_state=
!1-g_i

!print*, "sum proba", sumg


!> CHANGE_AC
        sh_print_env%proba     = g

print *, "sh_print_env%renorm", sh_print_env%renorm
print *, "Prob_gi_after_RINORMALIZATION: ", g

!to move on in the trivial hopping correction we use rescaled g to avoid
!unphysical number larger than 1. 

        !call a random number
        call init_random_seed()
        call random_number(hop)
        sh_print_env%hop       = hop
print *, "RANDOM_NUMB", hop
        sh_print_env%old_state = old_state
        !and see if this random number causes a surface hop
        sumg = 0.0D0
        do i=1,mol
print *, "Prob_gi: ", g(i)
                sumg = sumg + g(i)
                if (hop <= sumg) then
                        if (i .ne. old_state) then
                                new_state = i
print *, "new_state", new_state
!                                write(24,*) i                        
                               ! decohere = .FALSE.
                                sh_env%decohere = .FALSE.
                                attempt_hop = .TRUE.
print *, "attempt_hop", attempt_hop
                                sh_print_env%attempt_hop    = attempt_hop
                                sh_print_env%new_state = new_state
                                EXIT
                        elseif (i == old_state) then
                                new_state = old_state
                                attempt_hop = .FALSE.
print *, "new_state", new_state
print *, "attempt_hop", attempt_hop
                                sh_print_env%attempt_hop    = attempt_hop
                                sh_print_env%new_state = new_state
!                                write(24,*) state
                                
                                !if there has not been any surfacehop, we need to consider the possibility of the decoherence correction
                                !if the non-adiabatic coupling dphi to all other states is below a certain threshold
                                !then a wavefunction collapse will be enforced
                                do j=1,mol
                                        if (abs(dphi(old_state,j)) > decohere_criterion .and. old_state .ne. j) then
                                               ! decohere = .FALSE.
                                                sh_env%decohere = .FALSE.
                                                EXIT
                                        end if
                                end do
                               ! decohere = .TRUE.
                                sh_env%decohere = .TRUE.
                                EXIT
                        end if
                end if
        end do


        CALL adiab_print_sh_info(sh_print_env, keyword = 'TULLY_TEST') 

        DEALLOCATE(sh_print_env%proba)

        deallocate(ak)
        deallocate(b)
        deallocate(G)
        deallocate(ct)
        deallocate(akttm)
        deallocate(cttm_ad)

        deallocate(dphi)
        deallocate(T1)
        deallocate(T2)
        deallocate(DT)
        deallocate(TST)
        deallocate(TDT)
        deallocate(DS)
        deallocate(DC)

end subroutine modified_surfacehop


subroutine surfacehop(Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2, deltat, old_state, &
              sh_env, new_state, attempt_hop, decohere_criterion, sh_print_env)

        !using the method outlined in John Tully's Electronic Transitions paper
        !calculate the probability of the molecular dynamics switching to evolve on a different potential surface

        implicit none

        REAL(KIND=dp), INTENT(IN) :: decohere_criterion
        integer, intent(in) :: old_state
        REAL(KIND=dp), intent(in) :: deltat
        REAL(KIND=dp), dimension(:,:), intent(in) :: Sttp, Stt, Stptp, Sttm, Stmtm, C1, C2
        INTEGER, INTENT(OUT)        :: new_state
        LOGICAL, INTENT(OUT)        :: attempt_hop
        TYPE(sh_print_env_type), INTENT(OUT) :: sh_print_env
        TYPE(sh_env_type), INTENT(OUT) :: sh_env

        integer :: i, j, mol, k, l
        REAL(KIND=dp) :: hop, sumg
        REAL(KIND=dp), dimension(:,:), allocatable :: dphi
        REAL(KIND=dp), dimension(:), allocatable :: b, g
        REAL(KIND=dp), dimension(:,:), allocatable :: T1, T2, DT, TST, TDT, DS, DC
        !COMPLEX(KIND=dp), dimension(:), intent(in) :: sh_env%ctt
        COMPLEX(KIND=dp), dimension(:), allocatable :: ct
        COMPLEX(KIND=dp), dimension(:), allocatable :: ak
        COMPLEX(KIND=dp) :: cstate

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

        ALLOCATE(sh_print_env%proba(mol))


        !use the matinvsqrt subroutine to find T = S^(-1/2)
        T1 = Stmtm
        T2 = Stptp
        call matinvsqrt(mol, .TRUE., T1)
        call matinvsqrt(mol, .TRUE., T2)

        !and do some straightforward first-order numerical derivatives
        DT = (0.5D0/deltat)*(T2 - T1)
        DS = (0.5D0/deltat)*(Sttp - Sttm)
        DC = (1.0D0/deltat)*(C2 - C1)
        
        print *, "NACE_FINTE_DIFF", DS

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
                cstate = cstate + cmplx(C2(i,old_state), 0.0D0,dp)*sh_env%ctt(i)
        end do
        do i=1,mol
                !calculate a and b
                ct(i) = cmplx(0.0D0, 0.0D0,dp)
                do j=1,mol
                        ct(i) = ct(i) + cmplx(C2(j,i), 0.0D0,dp)*sh_env%ctt(j)
                end do
                ak(i) = ct(i)*conjg(cstate)
                b(i) = -2.0D0*REAL(conjg(ak(i))*dphi(i,old_state))
!print*, "i b dphi", i, b(i), dphi(i, old_state)
        end do

        sumg = 0.0D0
        !calculate g and sum(g)
        do i=1,mol
!print*, "deltat b ak", i, deltat, b(i), abs(ak(old_state))
                g(i) = (deltat*b(i))/abs(ak(old_state))
print *, "Probability_INIT: ", i, g(i)
                !as in Tully, set negative probabilities to 0
                if (g(i) < 0.0D0) then
                        g(i) = 0.0D0
                end if
                if (i .ne. old_state) then
                        sumg = sumg + g(i)
!print*, "sum proba", sumg, g(i)
print *, "Probability_other_states: ", i, g(i)
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

print *, "Probability_RENORM: ", g

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
                               ! decohere = .FALSE.
                                sh_env%decohere = .FALSE.
                                attempt_hop = .TRUE.
                                sh_print_env%attempt_hop    = attempt_hop
                                sh_print_env%new_state = new_state
                                EXIT
                        elseif (i == old_state) then
                                new_state = old_state
                                attempt_hop = .FALSE.
                                sh_print_env%attempt_hop    = attempt_hop
                                sh_print_env%new_state = new_state
!                                write(24,*) state
                                
                                !if there has not been any surfacehop, we need to consider the possibility of the decoherence correction
                                !if the non-adiabatic coupling dphi to all other states is below a certain threshold
                                !then a wavefunction collapse will be enforced
                                do j=1,mol
                                        if (abs(dphi(old_state,j)) > decohere_criterion .and. old_state .ne. j) then
                                               ! decohere = .FALSE.
                                                sh_env%decohere = .FALSE.
                                                EXIT
                                        end if
                                end do
                               ! decohere = .TRUE.
                                sh_env%decohere = .TRUE.
                                EXIT
                        end if
                end if
        end do


        CALL adiab_print_sh_info(sh_print_env, keyword = 'TULLY_TEST') 

        DEALLOCATE(sh_print_env%proba)

        deallocate(ak)
        deallocate(b)
        deallocate(G)
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

subroutine adiabaticity(H, mol, C, E, print, sh_env)

        !this subroutine will allow calculation of surface hopping probabilities
        !between N adiabatic surfaces

        integer :: stat, LWORK
        integer, intent(in) :: mol
        !REAL(KIND=dp) :: kcal
        REAL(KIND=dp), dimension(:,:), intent(in) :: H
        REAL(KIND=dp), dimension(:,:), intent(out), allocatable :: C

        TYPE(sh_env_type), INTENT(INOUT), OPTIONAL               ::    sh_env
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
        
        if (print.AND.PRESENT(sh_env)) then 
        allocate(sh_env%C_tp(mol,mol))
                sh_env%C_tp=C
        end if 
        
         
        deallocate(DGETRI_WORK)

end subroutine adiabaticity

SUBROUTINE calc_hop_nacv(force_env, hop_nacv, ext_old_state, ext_new_state, symmetry, fast, &
                         exF_11, exF_1, exF_2, exF_3)
        TYPE(force_env_type), POINTER                            :: force_env
        REAL(kind=dp), DIMENSION(:,:), INTENT(OUT), ALLOCATABLE  :: hop_nacv
        INTEGER, OPTIONAL                                        :: ext_old_state, ext_new_state
        CHARACTER(len=*), OPTIONAL                               :: symmetry
        LOGICAL, OPTIONAL                                        :: fast
        REAL(kind=dp), DIMENSION(:,:), INTENT(OUT), ALLOCATABLE, OPTIONAL  :: exF_11, exF_1, exF_2, exF_3


        CHARACTER(len=*), PARAMETER :: routineN = 'calc_hop_nacv', &
        routineP = moduleN//':'//routineN

        INTEGER                                      :: aom_iatom, iatom, num_mol, nadiab, num_j, natom
        INTEGER                                      :: old_state, new_state
        REAL(KIND=dp)                                :: scaling, cbar, diff_energy
        REAL(KIND=dp)                                :: E_old_state, E_new_state
        REAL(KIND=dp), dimension(:,:,:), allocatable :: Fin
        REAL(KIND=dp), dimension(:,:,:,:), allocatable :: nacv, Wprime
        integer, dimension(:,:), allocatable :: NN
        integer :: i, j, k, mat, d, l, atoms
        REAL(KIND=dp), dimension(:,:,:,:), allocatable :: DH, my_DP
        REAL(KIND=dp), dimension(:,:), allocatable :: F_1, F_11, F_2, F_22, F_3
        REAL(KIND=dp), dimension(:), allocatable ::   Sdiag
        REAL(KIND=dp), dimension(:,:), allocatable :: S, C, T, X, LL, Y
        LOGICAL, dimension(:), allocatable         :: ext_pres
        REAL(KIND=dp)                              :: contribf1, contribf11, contribf2, contribf3
        REAL(KIND=dp)                              :: normf1, normf2, normf22, normf3, normf11, normnacv
        LOGICAL                                    :: my_fast
        REAL(KIND = dp)                            :: error_f1, error_f11, error_f2, error_f22, error_f3

        mat    = SIZE(force_env%mixed_env%aom_env%Stptp, 1)
        natom  = SIZE(force_env%mixed_env%aom_env%diabatic_forces, 2)
        nadiab = SIZE(force_env%mixed_env%aom_env%diabatic_forces, 3)
        atoms  = SIZE(force_env%mixed_env%aom_env%connlist1, 2)

        my_fast = .FALSE.
        IF (PRESENT(fast)) my_fast = fast

        allocate(T(mat,mat))
        allocate(NN(mat,mat))
        allocate(S(mat,mat))
        allocate(C(mat,mat))
        allocate(X(mat,mat))
        allocate(Y(mat,mat))
        allocate(LL(mat,mat))

        allocate(Sdiag(mat))
        allocate(ext_pres(natom))

        allocate(DH(mat,mat,natom,3))
        allocate(my_DP(mat,mat,natom,3)) 
        allocate(nacv(3,mat,mat, mat*atoms)) 
        allocate(Wprime(mat,mat,natom,3)) 
        allocate(F_1(3,natom))
        allocate(F_11(3,natom))
        allocate(F_2(3,natom))
        allocate(F_22(3,natom))
        allocate(F_3(3,natom))
        allocate(Fin(3,natom, mat))
        allocate(hop_nacv(3, natom)) 

        scaling = force_env%mixed_env%aom_env%scaling
        cbar = force_env%mixed_env%aom_env%cbar
        NN = force_env%mixed_env%aom_env%NN
        S = force_env%mixed_env%aom_env%Stptp
        C = force_env%mixed_env%aom_env%eigenvector_matrix
        
        nacv = force_env%mixed_env%aom_env%nacv
        Fin = force_env%mixed_env%aom_env%diabatic_forces
        ext_pres = force_env%mixed_env%aom_env%ext_pres       
        IF (PRESENT(ext_old_state)) THEN
           old_state = ext_old_state
        ELSE
           old_state = force_env%mixed_env%aom_env%active_state
        ENDIF
        IF (PRESENT(ext_new_state)) THEN
           new_state = ext_new_state 
        ELSE
           new_state = force_env%mixed_env%aom_env%new_state
        ENDIF
        E_old_state = force_env%mixed_env%aom_env%adiab_energies_t(old_state)
        E_new_state = force_env%mixed_env%aom_env%adiab_energies_t(new_state)

        IF (E_old_state.EQ.E_new_state) THEN
           diff_energy = 0.0D0
        ELSE
           diff_energy = 1 / (E_new_state - E_old_state)
        ENDIF  

    IF (.NOT.my_fast) THEN
        T = S
        CALL matinvsqrt(mat, .FALSE., T)
        X = matmul(T, C)
        CALL adiabaticity(T, mat, LL, Sdiag, .FALSE.)
        CALL adiabaticity(S, mat, LL, Sdiag, .FALSE.)

!CALCULATE Y = L^{dagger} U
        Y = matmul( transpose(LL), C)

!CALCULATE GRADIENT OF S: my_DP
        aom_iatom = 0
        do iatom=1, natom
           if (ext_pres(iatom)) then
              aom_iatom = aom_iatom +1
              do d=1,3
                  my_DP(:,:,iatom,d) = nacv(d,:,:,aom_iatom)
               enddo
           else 
              my_DP(:,:,iatom,:) = 0.0D0
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

!CALCULATE W'
        do d=1,3
           num_j = -1
           do j=1,natom
              IF (ext_pres(j)) THEN
                  num_j   = num_j + 1
                  num_mol = int(num_j/atoms) + 1
              ENDIF
              do k=1, mat
                 do l=1, mat
                    Wprime(k,l,j,d) = 0.0d0
                    do i=1,mat 
                       IF (ext_pres(j)) THEN
                          IF (i.ne.num_mol) THEN
                             Wprime(k,l,j,d) = Wprime(k,l,j,d) + my_DP(i,num_mol,j,d) &
                                               *( LL(num_mol,k)*LL(i,l)  + LL(num_mol,l) * LL(i, k))
                          ENDIF
                       ENDIF
                    ENDDO
                    Wprime(k, l, j, d) = - Wprime( k, l, j, d) / ( SQRT(Sdiag(k))*SQRT(Sdiag(l)) + Sdiag(l) )
                 ENDDO
              ENDDO
           ENDDO
        ENDDO          

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
     ELSE
        DH = 0.0D0
        do i=1,mat
                do j=1,mat
                        do k=1, natom
                                do d=1,3
                                        if (NN(i,j) == 1) then
                                                if (i == j) then
                                                        DH(i,j,k,d) = Fin(d,k,i)
                                                end if
                                        end if
                                end do
                        end do
                end do
        end do
     ENDIF


       IF (.NOT.my_fast) THEN
           num_mol = 0
           do d=1,3
              num_j = -1
              do j=1,natom
                 F_1(d,j) = 0.0D0
                 F_11(d,j) = 0.0D0
                 F_2(d,j) = 0.0D0
                 F_3(d,j) = 0.0D0
                 IF (ext_pres(j)) THEN
                    num_j   = num_j + 1
                    num_mol = int(num_j/atoms) + 1 
                 ENDIF
                 do i=1,mat
                    F_1(d,j) = F_1(d,j) + C(i,old_state)*DH(i,i,j,d)*C(i,new_state)*diff_energy
                    IF (ext_pres(j)) THEN
                       F_2(d, j) = F_2(d,j) + X(i, old_state) * nacv(d, i, num_mol, j) * X( num_mol, new_state) 
                       IF (i.ne.num_mol) THEN
                          F_11(d,j) = F_11(d,j) + diff_energy * DH(i,num_mol,j,d) * (              &
                                                 C(num_mol, new_state) * C(i, old_state) +   &
                                                 C(num_mol, old_state) * C(i,new_state) )
                       ENDIF
                       DO k = 1, mat
                          F_3(d,j) = F_3(d,j) + Y(i, old_state) * Wprime(i, k, j, d) * Y(k, new_state)
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO 


           hop_nacv = F_1 + F_11 + F_2 + F_3
           IF (PRESENT(symmetry)) THEN
              DO j=1, natom
                 DO d=1, 3
                    print*, "NACV", symmetry, hop_nacv(d,j) 
                 ENDDO
              ENDDO
           ENDIF
           IF (PRESENT(exF_1))  exF_1  = F_1
           IF (PRESENT(exF_11)) exF_11 = F_11
           IF (PRESENT(exF_2))  exF_2  = F_2
           IF (PRESENT(exF_3))  exF_3  = F_3
           IF (PRESENT(exF_1)) THEN
              F_22 = F_2 + F_3
              DO j=1, natom
                 normf1 = 0
                 normf11 = 0
                 normf22 = 0
                 normnacv = 0
                 contribf1 = 0
                 contribf11 = 0
                 contribf2 = 0
                 contribf3 = 0
                 DO d=1,3
                    IF ( mod(j,3).eq.1 ) THEN
                    print*, "NACVTOTCARBON = ", d, j, hop_nacv(d,j)
                    print*, "NACVF1CARBON   = ", d, j, F_1(d,j)
                    print*, "NACVFF1CARBON   = ", d, j, F_1(d,j) + F_11(d,j)
print*, "ERRORF1CARBON", SQRT( ( hop_nacv(d,j) - F_1(d,j) )**2 / hop_nacv(d,j)**2 )
print*, "ERRORFF1CARBON", SQRT( ( hop_nacv(d,j) - (F_1(d,j)+F_11(d,j) ))**2 / hop_nacv(d,j)**2 )
                    ELSE 
                    print*, "NACVTOTHYDROGEN = ", d, j, hop_nacv(d,j)
                    print*, "NACVF1HYDROGEN   = ", d, j, F_1(d,j)
                    print*, "NACVFF1HYDROGEN   = ", d, j, F_1(d,j) + F_11(d,j)
print*, "ERRORF1HYRDOGEN", SQRT( ( hop_nacv(d,j) - F_1(d,j) )**2 /  hop_nacv(d,j)**2 )
print*, "ERRORFF1HYDROGEN", SQRT( ( hop_nacv(d,j) - (F_1(d,j)+F_11(d,j) ))**2 / hop_nacv(d,j)**2 )
                    ENDIF 
                    !print*, "NACVF11  = ", d, j, F_11(d,j), symmetry
                    !print*, "NACVF2   = ", d, j, F_2(d,j), symmetry
                    !print*, "NACVF3   = ", d, j, F_3(d,j), symmetry
                    !print*, "NACVFtot = ", d, j, hop_nacv(d,j), symmetry
                    contribf1 = contribf1 + F_1(d,j) * hop_nacv(d,j)
                    normf1 = normf1 + F_1(d,j)*F_1(d,j)
                    contribf11 = contribf11 + F_11(d,j) * hop_nacv(d,j)
                    normf11 = normf11 + F_11(d,j)*F_11(d,j)
                    contribf2 = contribf2 + F_2(d,j) * hop_nacv(d,j)
                    normf2 = normf2 + F_2(d,j)*F_2(d,j)
                    normf22 = normf22 + F_22(d,j)*F_22(d,j)
                    contribf3 = contribf3 + F_3(d,j) * hop_nacv(d,j)
                    normf3 = normf3 + F_3(d,j)*F_3(d,j)
                    normnacv = normnacv + hop_nacv(d,j)*hop_nacv(d,j)
                 ENDDO
                 normnacv = SQRT(normnacv)
                 normf1 = SQRT(normf1)
                 normf11 = SQRT(normf11)
                 normf2 = SQRT(normf2)
                 normf22 = SQRT(normf22)
                 normf3 = SQRT(normf3)
                 error_f1 = SUM((hop_nacv(:,j)/normnacv - F_1(:,j)/normf1)**2) / 3
                 !print*, "ERRORF1 ATOM", j, error_f1
                 error_f1 = SUM((hop_nacv(:,j) - F_1(:,j))**2) / 3
                 !print*, "ERRORTOTF1 ATOM", j, error_f1
                 error_f1 = SUM((hop_nacv(:,j) - (F_1(:,j) + F_11(:,j)))**2) / 3
                 !print*, "ERRORTOTFF1 ATOM", j, error_f1
                 error_f11 = SUM((hop_nacv(:,j)/normnacv - F_11(:,j)/normf11)**2)  /3
                 error_f22 = SUM((hop_nacv(:,j)/normnacv - F_22(:,j)/normf22)**2)  /3
                 error_f2 = SUM((hop_nacv(:,j)/normnacv - F_2(:,j)/normf2)**2)     /3
                 error_f3 = SUM((hop_nacv(:,j)/normnacv - F_3(:,j)/normf3)**2)     /3
                 !print*, "ERRORF11 ATOM", j, error_f11
                 !print*, "ERRORF22 ATOM", j, error_f22
                 !print*, "ERRORF2 ATOM", j, error_f2
                 !print*, "ERRORF3 ATOM", j, error_f3
                 !print*, "CONTRIBF1 = ", contribf1 / normnacv,  "ATOM", j
                 !print*, "CONTRIBF11 = ", contribf11 / normnacv, "ATOM", j
                 !print*, "CONTRIBF2 = ", contribf2 / normnacv, "ATOM", j
                 !print*, "CONTRIBF3 = ", contribf3 / normnacv, "ATOM", j
              ENDDO
           ENDIF
        ELSE
           do d=1,3
                      do j=1,natom
                 F_1(d,j) = 0.0D0
                 do i=1,mat
                    F_1(d,j) = F_1(d,j) + C(i,old_state)*DH(i,i,j,d)*C(i,new_state)*diff_energy
                 ENDDO
              ENDDO
           ENDDO 

           hop_nacv = F_1 
        ENDIF
 
        deallocate(T)
        deallocate(DH)
        deallocate(my_DP)
        deallocate(F_1)
        deallocate(F_11)
        deallocate(F_2)
        deallocate(F_22)
        deallocate(F_3)
END SUBROUTINE calc_hop_nacv

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
              my_DP(:,:,iatom,:) = 0.0D0
           endif
        enddo
!       IF (aom_iatom.NE.SIZE(nacv,4)) THEN
!       ENDIF

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

        deallocate(T, DH, DT, my_DP, F, nsave)
                
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
        TYPE(sh_print_env_type), INTENT(INOUT)   :: sh_print_env
        INTEGER, INTENT(OUT)                      :: first_adiabat

        CHARACTER(len=*), PARAMETER :: routineN = 'choose_firstadiabat', &
        routineP = moduleN//':'//routineN
        
        INTEGER                                   :: i
        REAL( KIND =dp)                           :: random, sumc

        ALLOCATE(sh_print_env%pop_init(SIZE(C(init,:))))

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
              EXIT
           ENDIF
        ENDDO
        CALL adiab_print_sh_info(sh_print_env, keyword = 'START')
        DEALLOCATE(sh_print_env%pop_init)
   END SUBROUTINE choose_firstadiabat


   SUBROUTINE simple_vel_rescaling(velocities, delta_energy, kinetic_energy_prime, rescale_factor)
        !TYPE(tmp_variables_type), POINTER        :: tmp
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        REAL(kind=dp), INTENT(IN)                :: delta_energy
        REAL(kind=dp), INTENT(IN)                :: kinetic_energy_prime
         !REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_plus
        REAL(kind=dp), INTENT(OUT)               :: rescale_factor

        CHARACTER(len=*), PARAMETER :: routineN = 'simple_vel_rescaling', &
        routineP = moduleN//':'//routineN


 !> Factor to rescale the velocities (Eq. 49 of PAPER_1)
        rescale_factor = SQRT(  delta_energy /kinetic_energy_prime )
!print *, 'rescale_factor', rescale_factor

        velocities = rescale_factor*velocities
!print *, 'velocities', velocities 

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


   SUBROUTINE check_hop_attempt(velocities, atomic_kind_set, nparticle_kind, local_particles, &
 force_env)
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        !TYPE(tmp_variables_type), POINTER        :: tmp
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        INTEGER, INTENT(IN)                      :: nparticle_kind
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(force_env_type), POINTER            :: force_env

        CHARACTER(len=*), PARAMETER :: routineN = 'check_hop_attempt', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: delta_energy, rescale_factor
        REAL(kind=dp)                            :: kinetic_energy_prime, kinetic_energy_qsys, A, B
        REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE  :: hop_nacv, hop_nacv_oldold, hop_nacv_newnew
        REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE  :: hop_nacv_oldnew, hop_nacv_newold
        REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE  :: F_11, F_1, F_2, F_3
        CHARACTER(LEN=default_string_length)     :: rescale_keyword, nacv_keyword, reversal_keyword
        TYPE(section_vals_type), POINTER         :: mixed_section, adiab_section, coupling_section
        TYPE(section_vals_type), POINTER         :: force_env_section
        INTEGER                                  :: old_state, new_state
        TYPE(sh_print_env_type)                  :: sh_print_env
        LOGICAL                                  :: translate_keyword
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_com


        CALL force_env_get(force_env=force_env,&
                       force_env_section=force_env_section)
        mixed_section => section_vals_get_subs_vals(force_env_section,"MIXED")
        adiab_section => section_vals_get_subs_vals(mixed_section, "ADIABATIC")
        coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
        CALL section_vals_val_get(coupling_section,"METHOD_RESCALING",c_val=rescale_keyword)
        CALL section_vals_val_get(coupling_section,"CENTER_OF_MASS",l_val=translate_keyword)
        sh_print_env = force_env%mixed_env%aom_env%sh_print_env

        ALLOCATE(velocities_com(SIZE(velocities, 1), SIZE(velocities, 2)))
        CALL translate_com(velocities, translate_keyword, velocities_com, &
                           local_particles, nparticle_kind = nparticle_kind, atomic_kind_set = atomic_kind_set)
        SELECT CASE(rescale_keyword)
        CASE("SIMPLE")
           CALL cal_delta_energy(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                                 local_particles, delta_energy, kinetic_energy_prime)
!print *, 'kinetic_energy_prime is', kinetic_energy_prime          
           sh_print_env%old_state        = force_env%mixed_env%aom_env%active_state
           sh_print_env%new_state        = force_env%mixed_env%aom_env%new_state
           IF (delta_energy.lt.0) THEN
              force_env%mixed_env%aom_env%successful_hop = .FALSE.
           ELSE 
               !> Factor to rescale the velocities (Eq. 49 of PAPER_1)
              !rescale_factor = SQRT(  delta_energy /kinetic_energy_prime )
              !velocities = rescale_factor*velocities
              force_env%mixed_env%aom_env%successful_hop = .TRUE.
              force_env%mixed_env%aom_env%old_state      = force_env%mixed_env%aom_env%active_state
              force_env%mixed_env%aom_env%active_state   = force_env%mixed_env%aom_env%new_state
              force_env%mixed_env%aom_env%kinetic_energy = kinetic_energy_prime
              force_env%mixed_env%aom_env%delta_energy   = delta_energy
              CALL simple_vel_rescaling(velocities, delta_energy, kinetic_energy_prime, rescale_factor)
           ENDIF

           sh_print_env%adiab_energies_t = force_env%mixed_env%aom_env%adiab_energies_t
           sh_print_env%delta_energy     = delta_energy
           sh_print_env%successful_hop   = force_env%mixed_env%aom_env%successful_hop
           sh_print_env%kinetic_energy   = kinetic_energy_prime
           IF (sh_print_env%successful_hop) sh_print_env%rescale_factor = rescale_factor
           CALL adiab_print_sh_info(sh_print_env, keyword = 'SIMPLE_RESCALING')


        CASE("SIMPLE_QSYS")
           CALL cal_delta_ene_QSYS(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                         local_particles, delta_energy, kinetic_energy_qsys)
!print *, 'kinetic_energy_qsys is', kinetic_energy_qsys           
           sh_print_env%old_state        = force_env%mixed_env%aom_env%active_state
           sh_print_env%new_state        = force_env%mixed_env%aom_env%new_state
           IF (delta_energy.lt.0) THEN
              force_env%mixed_env%aom_env%successful_hop = .FALSE.
           ELSE 
               !> Factor to rescale the velocities (Eq. 49 of PAPER_1)
              !rescale_factor = SQRT(  delta_energy /kinetic_energy_prime )
              !velocities = rescale_factor*velocities
              force_env%mixed_env%aom_env%successful_hop = .TRUE.
              force_env%mixed_env%aom_env%old_state      = force_env%mixed_env%aom_env%active_state
              force_env%mixed_env%aom_env%active_state   = force_env%mixed_env%aom_env%new_state
              force_env%mixed_env%aom_env%kinetic_energy = kinetic_energy_qsys
              force_env%mixed_env%aom_env%delta_energy   = delta_energy
              CALL simple_vel_resc_qsys(atomic_kind_set, local_particles,  &
                       force_env, nparticle_kind, velocities, delta_energy, kinetic_energy_qsys, &
                         rescale_factor)

           ENDIF

           sh_print_env%adiab_energies_t = force_env%mixed_env%aom_env%adiab_energies_t
           sh_print_env%delta_energy     = delta_energy
           sh_print_env%successful_hop   = force_env%mixed_env%aom_env%successful_hop
           sh_print_env%kinetic_energy   = kinetic_energy_qsys
           IF (sh_print_env%successful_hop) sh_print_env%rescale_factor = rescale_factor
           CALL adiab_print_sh_info(sh_print_env, keyword = 'SIMPLE_RESCALING')


        CASE("NACV")
           CALL section_vals_val_get(coupling_section,"METHOD_ADIABATIC_NACV",c_val=nacv_keyword)
           CALL cal_delta_energy(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                                 local_particles, delta_energy, kinetic_energy_prime)
           sh_print_env%kinetic_energy   = kinetic_energy_prime
           sh_print_env%old_state        = force_env%mixed_env%aom_env%active_state
           sh_print_env%new_state        = force_env%mixed_env%aom_env%new_state
           IF ((delta_energy.lt.0).AND.(nacv_keyword.ne."TEST")) THEN
              force_env%mixed_env%aom_env%successful_hop = .FALSE.

              sh_print_env%adiab_energies_t = force_env%mixed_env%aom_env%adiab_energies_t
              sh_print_env%delta_energy     = delta_energy
              sh_print_env%successful_hop   = force_env%mixed_env%aom_env%successful_hop
              sh_print_env%old_state        = force_env%mixed_env%aom_env%active_state
              sh_print_env%new_state        = force_env%mixed_env%aom_env%new_state
              IF (sh_print_env%successful_hop) sh_print_env%rescale_factor = rescale_factor
              CALL adiab_print_sh_info(sh_print_env, keyword = 'SIMPLE_RESCALING')
           ELSE 
              SELECT CASE(nacv_keyword)
              CASE("TEST")
                 old_state = force_env%mixed_env%aom_env%active_state
                 new_state = force_env%mixed_env%aom_env%new_state
                 CALL calc_hop_nacv( force_env, hop_nacv_oldold, old_state, old_state, "OLDOLD")
                 CALL calc_hop_nacv( force_env, hop_nacv_newnew, new_state, new_state, "NEWNEW")
                 CALL calc_hop_nacv( force_env, hop_nacv_newold, new_state, old_state, "NEWOLD")
                 CALL calc_hop_nacv( force_env, hop_nacv_oldnew, old_state, new_state, "OLDNEW")
                 CALL adiab_print_hop_nacv(sh_print_env, hop_nacv_oldold, &
                                  hop_nacv_newnew, hop_nacv_newold, hop_nacv_oldnew) 
                 print*, "ERRORASYM = ", SUM( (hop_nacv_oldnew + hop_nacv_newold ) **2 )
                 hop_nacv = hop_nacv_oldnew
                 print*, "totalnacv =" , SUM(hop_nacv(1,:)),  SUM(hop_nacv(2,:)), SUM(hop_nacv(3,:))
                 print*, "totalnacvasym =" , SUM(hop_nacv_oldnew(1,:) - hop_nacv_newold(1,:)) / 2,  & 
                                             SUM(hop_nacv_oldnew(2,:) - hop_nacv_newold(2,:)) / 2,  &
                                             SUM(hop_nacv_oldnew(3,:) - hop_nacv_newold(3,:)) / 2
              CASE("CONTRIBUTION")
                 CALL calc_hop_nacv( force_env, hop_nacv, exF_11 = F_11, exF_1 = F_1, exF_2 = F_2, exF_3 =F_3)
                 print*, "totalnacv =" , SUM(hop_nacv(1,:)),  SUM(hop_nacv(2,:)), SUM(hop_nacv(3,:))
                 print*, "totalnacvF1 =" , SUM(F_1(1,:)),  SUM(F_1(2,:)), SUM(F_1(3,:))
                 print*, "totalnacvFF1 =" , SUM(F_1(1,:) + F_11(1,:)), &
                                            SUM(F_1(2,:) + F_11(2,:)), &
                                            SUM(F_1(3,:) + F_11(3,:))
                 CALL adiab_print_hop_nacv(sh_print_env, hop_nacv, F_11, F_1, F_2, F_3) 
                 CALL cal_delta_energy_nacv(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                                            local_particles, F_1 + F_11, delta_energy, A, B) 
print*, "DELTAFF1 = ", delta_energy
IF ((delta_energy).lt.0) THEN
     print*, "LOGFF1 = ", 0
ELSE
     print*, "LOGFF1 = ", 1
ENDIF
     
                 CALL cal_delta_energy_nacv(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                                            local_particles, F_1, delta_energy, A, B) 
print*, "DELTAF1 = ", delta_energy
IF ((delta_energy).lt.0) THEN
     print*, "LOGF1 = ", 0
ELSE
     print*, "LOGF1 = ", 1
ENDIF
     
              CASE("TOTAL")
                 CALL calc_hop_nacv( force_env, hop_nacv)
                 CALL adiab_print_hop_nacv(sh_print_env, hop_nacv) 
              CASE("FAST")
                 CALL calc_hop_nacv( force_env, hop_nacv, fast=.TRUE.)
                 CALL adiab_print_hop_nacv(sh_print_env, hop_nacv) 
              CASE DEFAULT
                CPABORT("There is no method to calculate the NACV.")
              END SELECT 
              CALL cal_delta_energy_nacv(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                                          local_particles, hop_nacv, delta_energy, A, B) 
IF ((delta_energy).lt.0) THEN
     print*, "LOGNACV = ", 0
ELSE
     print*, "LOGNACV = ", 1
ENDIF
              If (delta_energy.lt.0) THEN
                  force_env%mixed_env%aom_env%successful_hop = .FALSE.
                  CALL section_vals_val_get(coupling_section,"METHOD_REVERSAL",c_val=reversal_keyword)
                  CALL hop_nacv_vel_reversing(nparticle_kind, velocities,  & 
                                              force_env%mixed_env%aom_env%forces_adiab, &
                                              force_env%mixed_env%aom_env%forces_adiab_new, &
                                              hop_nacv, atomic_kind_set, local_particles, reversal_keyword, &
                                              sh_print_env)
              ELSE 
                  force_env%mixed_env%aom_env%successful_hop = .TRUE.
                  force_env%mixed_env%aom_env%old_state      = force_env%mixed_env%aom_env%active_state
                  force_env%mixed_env%aom_env%active_state   = force_env%mixed_env%aom_env%new_state
                  force_env%mixed_env%aom_env%kinetic_energy = kinetic_energy_prime
                  force_env%mixed_env%aom_env%delta_energy   = delta_energy
                  !force_env%mixed_env%aom_env%new_state = -1
                  CALL hop_nacv_vel_rescaling(nparticle_kind, velocities, delta_energy, A, B, hop_nacv, &
                                              atomic_kind_set, local_particles, rescale_factor)
              ENDIF

              sh_print_env%A = A
              sh_print_env%B = B
              sh_print_env%adiab_energies_t = force_env%mixed_env%aom_env%adiab_energies_t
              sh_print_env%delta_energy     = delta_energy
              sh_print_env%successful_hop   = force_env%mixed_env%aom_env%successful_hop
              IF (sh_print_env%successful_hop) THEN
                 sh_print_env%rescale_factor = rescale_factor
                 sh_print_env%momentum_shift(:) = rescale_factor *        & 
                          (/ SUM(hop_nacv(1,:)), SUM(hop_nacv(2,:)), SUM(hop_nacv(3,:)) /)
              ENDIF
              CALL adiab_print_sh_info(sh_print_env, keyword = 'NACV_RESCALING')
              !IF (.NOT.(sh_print_env%successful_hop)) THEN
              !   CALL adiab_print_sh_info(sh_print_env, keyword = 'NACV_REVERSING')
              !ENDIF 
           ENDIF 
        CASE DEFAULT
                CPABORT("There is no method to do the rescaling of the velocities.")
        END SELECT

   END SUBROUTINE check_hop_attempt

   SUBROUTINE hop_nacv_vel_reversing(nparticle_kind, velocities, forces_old, forces_new, hop_nacv, atomic_kind_set, &
                                     local_particles, reversal_keyword, sh_print_env)
        INTEGER, INTENT(IN)                      :: nparticle_kind
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)   :: forces_old, forces_new
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)   :: hop_nacv
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                             :: atomic_kind_set
        CHARACTER(LEN=default_string_length)     :: reversal_keyword
        TYPE(sh_print_env_type)                  :: sh_print_env

        CHARACTER(len=*), PARAMETER :: routineN = 'hop_nacv_vel_rescaling', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: product_pd, product_dfo, product_dfn
        REAL(kind=dp)                            :: scalar, norm
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        LOGICAL                                  :: reverse

!> TO EXPLAIN  

        IF ((reversal_keyword.ne."ALWAYS").OR.(reversal_keyword.ne."NEVER")) THEN
           product_pd = 0.0D0
           product_dfo = 0.0D0
           product_dfn = 0.0D0
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 product_pd = product_pd + mass * (                                   &
                                velocities(1, iparticle) * hop_nacv(1, iparticle) +   &            
                                velocities(2, iparticle) * hop_nacv(2, iparticle) +   &            
                                velocities(3, iparticle) * hop_nacv(3, iparticle) )               
                 product_dfo = product_dfo +  ( &
                                hop_nacv(1, iparticle) * forces_old(1, iparticle) +    &           
                                hop_nacv(2, iparticle) * forces_old(2, iparticle) +    &           
                                hop_nacv(3, iparticle) * forces_old(3, iparticle) )               
                 product_dfn = product_dfn +  ( &
                                hop_nacv(1, iparticle) * forces_new(1, iparticle) +    &           
                                hop_nacv(2, iparticle) * forces_new(2, iparticle) +    &           
                                hop_nacv(3, iparticle) * forces_new(3, iparticle) )               
              ENDDO
           END DO
        ENDIF

        reverse = .FALSE.
        SELECT CASE(reversal_keyword)
        CASE("NEVER")
           reverse = .FALSE.
        CASE("ALWAYS")
           reverse = .TRUE.
        CASE("TRHULAR")
           IF ((product_pd*product_dfn).lt.0) reverse = .TRUE.
        CASE("SUBOTNIK")
           IF (((product_pd*product_dfn).lt.0) .AND. ((product_dfn*product_dfo).lt.0)) reverse = .TRUE.
        CASE DEFAULT
           CPABORT("There is not method to reverse the velocities after a frustrated hop") 
        END SELECT

        IF (reverse) THEN
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 IF ((velocities(1, iparticle)**2 + velocities(2, iparticle)**2 + &
                      velocities(3, iparticle)**2 ).NE.0.0d0) THEN
                      scalar = SUM(velocities(:, iparticle) * hop_nacv(:, iparticle) )
                      norm = SUM( hop_nacv(:, iparticle)**2)
                      IF (norm.NE.0.0D0) THEN
                         velocities(:, iparticle) = velocities(:, iparticle) - 2 * (scalar/norm) * hop_nacv(:, iparticle)
                      ENDIF
                 ENDIF
              ENDDO
           END DO
        ENDIF
        
        sh_print_env%reversal_keyword = reversal_keyword
        sh_print_env%reverse = reverse
        sh_print_env%product_pd = product_pd
        sh_print_env%product_dfn = product_dfn
        sh_print_env%product_dfo = product_dfo
   END SUBROUTINE hop_nacv_vel_reversing



   SUBROUTINE hop_nacv_vel_rescaling(nparticle_kind, velocities, delta_energy, A, B, hop_nacv, atomic_kind_set, &
                                     local_particles, rescale_factor)
        INTEGER, INTENT(IN)                      :: nparticle_kind
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        REAL(kind=dp), INTENT(IN)                :: delta_energy, A, B
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)   :: hop_nacv
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                             :: atomic_kind_set
        REAL(kind=dp), INTENT(OUT)                :: rescale_factor

        CHARACTER(len=*), PARAMETER :: routineN = 'hop_nacv_vel_rescaling', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local 
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind

!> From Fabiano E, Keal TW, Thiel W. Implementation of surface hopping molecular dynamics using semiempirical methods. 
!>     Chem Phys 2008, 349:334–347.
!>  B = -B, compare to the paper

        IF (A == 0) THEN
           rescale_factor = 0
        ELSE
           IF (B < 0) THEN
              rescale_factor = (- B - SQRT(delta_energy))/(2 * A)
           ELSE
              rescale_factor = (- B + SQRT(delta_energy))/(2 * A)
           ENDIF
        ENDIF
        DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
              iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
              IF ((velocities(1, iparticle)**2 + velocities(2, iparticle)**2 + &
                   velocities(3, iparticle)**2 ).NE.0.0d0) THEN
                 velocities(:, iparticle ) = velocities(:, iparticle) + (rescale_factor/mass) * hop_nacv(:, iparticle)
              ENDIF
           ENDDO
        END DO
   END SUBROUTINE hop_nacv_vel_rescaling

   SUBROUTINE cal_delta_energy_nacv(nparticle_kind, velocities, force_env, atomic_kind_set, &
                                    local_particles, hop_nacv, delta_energy, A, B) 
        INTEGER, INTENT(IN)                      :: nparticle_kind
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(force_env_type), POINTER            :: force_env
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)   :: hop_nacv
        REAL(kind=dp), INTENT(OUT)                  :: delta_energy
        REAL(kind=dp), INTENT(OUT)                  :: A, B

        CHARACTER(len=*), PARAMETER :: routineN = 'cal_delta_energy_nacv', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind, new_state, old_state
        INTEGER                                  :: nparticle_local, coord
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind

!From Fabiano E, Keal TW, Thiel W. Implementation of surface hopping molecular dynamics using semiempirical methods. 
!     Chem Phys 2008, 349:334–347.
        
           A = 0
           B = 0
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 DO coord = 1, 3
                    IF ((velocities(1, iparticle)**2 + velocities(2, iparticle)**2 + &
                       velocities(3, iparticle)**2 ).NE.0.0d0) THEN
                       A = A + hop_nacv(coord, iparticle)**2 / mass
                    ELSE IF ( SUM(velocities).EQ.0) THEN
                       A = A + hop_nacv(coord, iparticle)**2 / mass
                    ENDIF
                    B = B + velocities(coord, iparticle) * hop_nacv(coord, iparticle)
                 ENDDO
              END DO
           END DO
           A = 0.5*A

           new_state = force_env%mixed_env%aom_env%new_state
           old_state = force_env%mixed_env%aom_env%active_state
           delta_energy = B**2 - 4 * A * (                                     &
                     force_env%mixed_env%aom_env%adiab_energies_t(new_state) - &
                     force_env%mixed_env%aom_env%adiab_energies_t(old_state))
   END SUBROUTINE cal_delta_energy_nacv
         
   SUBROUTINE cal_delta_energy(nparticle_kind, velocities, force_env, atomic_kind_set, &
                               local_particles, delta_energy, kinetic_energy_prime) 
        INTEGER, INTENT(IN)                      :: nparticle_kind
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(force_env_type), POINTER            :: force_env
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        REAL(kind=dp), INTENT(OUT)               :: delta_energy
        REAL(kind=dp), INTENT(OUT)               :: kinetic_energy_prime

        CHARACTER(len=*), PARAMETER :: routineN = 'delta_energy', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind, new_state, old_state
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind


           kinetic_energy_prime = 0.0D0
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
              END DO
           END DO
           kinetic_energy_prime = 0.5*kinetic_energy_prime

           new_state = force_env%mixed_env%aom_env%new_state
           old_state = force_env%mixed_env%aom_env%active_state
           delta_energy = kinetic_energy_prime-( &
                     force_env%mixed_env%aom_env%adiab_energies_t(new_state) - &
                     force_env%mixed_env%aom_env%adiab_energies_t(old_state))
   END SUBROUTINE cal_delta_energy


!---> CHANGE SG to calculate kinetic energy only for the quantum system         

   SUBROUTINE cal_delta_ene_QSYS(nparticle_kind, velocities, force_env, atomic_kind_set, &
                               local_particles, delta_energy, kinetic_energy_qsys) 
        INTEGER, INTENT(IN)                      :: nparticle_kind
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(force_env_type), POINTER            :: force_env
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        REAL(kind=dp), INTENT(OUT)               :: delta_energy
        REAL(kind=dp), INTENT(OUT)               :: kinetic_energy_qsys
!        LOGICAL, dimension(:), ALLOCATABLE, INTENT(IN)           :: ext_pres


        CHARACTER(len=*), PARAMETER :: routineN = 'cal_delta_ene_QSYS', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind, new_state, old_state
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind


           kinetic_energy_qsys = 0.0D0
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 if (force_env%mixed_env%aom_env%ext_pres(iparticle_local)) then 
                      iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                      kinetic_energy_qsys = kinetic_energy_qsys +   &
                         mass*( &
                         velocities(1, iparticle)**2 + &
                         velocities(2, iparticle)**2 + &
                         velocities(3, iparticle)**2  )
                 end if
              END DO
           END DO
          kinetic_energy_qsys = 0.5*kinetic_energy_qsys

           new_state = force_env%mixed_env%aom_env%new_state
           old_state = force_env%mixed_env%aom_env%active_state
           delta_energy = kinetic_energy_qsys-( &
                     force_env%mixed_env%aom_env%adiab_energies_t(new_state) - &
                     force_env%mixed_env%aom_env%adiab_energies_t(old_state))
   END SUBROUTINE cal_delta_ene_QSYS

   SUBROUTINE simple_vel_resc_qsys(atomic_kind_set,local_particles, force_env, nparticle_kind, & 
                                   velocities, delta_energy, kinetic_energy_qsys, rescale_factor)
        
        INTEGER, INTENT(IN)                      :: nparticle_kind 
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(force_env_type), POINTER            :: force_env
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        REAL(kind=dp), INTENT(IN)                :: delta_energy
        REAL(kind=dp), INTENT(IN)                :: kinetic_energy_qsys
        !REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_plus
        REAL(kind=dp), INTENT(OUT)               :: rescale_factor

        CHARACTER(len=*), PARAMETER :: routineN = 'simple_vel_resc_qys', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: iparticle_local, iparticle
        INTEGER                                  :: nparticle_local

!> Factor to rescale the velocities (Eq. 49 of PAPER_1)

        rescale_factor = SQRT(  delta_energy /kinetic_energy_qsys)

!print *, 'rescale_factor is', rescale_factor

           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 if (force_env%mixed_env%aom_env%ext_pres(iparticle_local)) then 
                      iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                         velocities(:, iparticle) =rescale_factor*velocities(:, iparticle)
                 end if
              END DO
           END DO
!print *, 'velocities are', velocities

   END SUBROUTINE simple_vel_resc_qsys



   SUBROUTINE total_momentum(velocities, atomic_kind_set, nparticle_kind, &
                               local_particles, force_env) 
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        INTEGER, INTENT(IN)                      :: nparticle_kind
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(force_env_type), POINTER            :: force_env
 
        CHARACTER(len=*), PARAMETER :: routineN = 'total_momentum', &
        routineP = moduleN//':'//routineN
 
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        REAL(kind=dp), DIMENSION(3)                       :: momentum


           momentum = 0.0D0 
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 momentum(:) = momentum(:) +   &
                    mass*velocities(:, iparticle)
              END DO
           END DO
           print*, "MOMENTUM = ", MOMENTUM
   END SUBROUTINE total_momentum
         
   SUBROUTINE translate_com(array, keyword_translate, translate_array, &
                            local_particles, atomic_kinds, nparticle_kind, atomic_kind_set)
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: array
        LOGICAL                                  :: keyword_translate
        REAL(KIND=dp), DIMENSION(:,:), POINTER   :: translate_array
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER, OPTIONAL                             :: atomic_kind_set

        CHARACTER(len=*), PARAMETER :: routineN = 'translate_com', &
        routineP = moduleN//':'//routineN

        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                             :: my_atomic_kind_set
        INTEGER                                  :: my_nparticle_kind
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        REAL(kind=dp), DIMENSION(3)              :: sum_array
        REAL(KIND=dp)                            :: sum_mass


        IF (keyword_translate) THEN
           IF (PRESENT(atomic_kinds)) THEN
              my_nparticle_kind = atomic_kinds%n_els
              my_atomic_kind_set => atomic_kinds%els
           ELSE
              IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
                 my_nparticle_kind = nparticle_kind
                 my_atomic_kind_set => atomic_kind_set
print*, "HERE"
              ENDIF
           ENDIF
           sum_array = 0.0D0
           sum_mass  = 0.0D0
           DO iparticle_kind=1, my_nparticle_kind
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 sum_array(:) = sum_array(:) +   &
                    mass*array(:, iparticle)
                 sum_mass = sum_mass + mass
              END DO
           END DO
           sum_array = sum_array / sum_mass

           DO iparticle_kind=1, my_nparticle_kind
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 translate_array(:, iparticle) = array(:, iparticle) - sum_array(:)
              END DO   
           END DO   
        ELSE
           translate_array = array
        ENDIF 

   END SUBROUTINE translate_com


END MODULE aom_utils
