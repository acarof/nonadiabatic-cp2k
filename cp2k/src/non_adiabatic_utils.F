!-----------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations           
!   !                                                                           
!   Copyright (C) 2000 - 2015  CP2K developers group                            
!   !                                                                           
!-----------------------------------------------------------------------------! 

!!!!!!!!!!!!!!!!!!!! CHANGE_SG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE non_adiabatic_utils
     
     USE kinds,                        ONLY: default_string_length, &
                                             dp                                    

     USE input_section_types,          ONLY:  section_vals_type, &
                                             section_vals_get, &
                                             section_vals_get_subs_vals, &
                                             section_vals_list_get, &
                                             section_vals_val_get

     USE aom_print,                    ONLY: adiab_print_coeff,              &
                                             adiab_print_pop,                &
                                             adiab_print_adiabat_energies,   &
                                             adiab_print_matrix_density,     &
                                             adiab_print_exact_forces,       &
                                             adiab_print_active_state,       &
                                             adiab_print_sh_info,            &
                                             adiab_print_hamiltonian, &
                      !                       adiab_print_nace,        &
                                             adiab_print_nacv        
     USE aom_types,                    ONLY: aom_env_type, &      
                                                sh_print_env_type 

     USE aom_utils,                    ONLY: makeX,  &
                                             adiabaticity,   &
                                            surfacehop,     &
                                             calcforces,     &
                                            read_listcoeff, &
                                            choose_firstadiabat, &
                                             simple_vel_rescaling, &
                                              translate_com


     USE non_adiabatic_types,           ONLY: sh_env_type

     USE atomic_kind_list_types,          ONLY: atomic_kind_list_type   

     USE atomic_kind_types,               ONLY: atomic_kind_type,&      
                                           get_atomic_kind,&       
                                           get_atomic_kind_set     
                                                                   
     USE distribution_1d_types,           ONLY: distribution_1d_type    



#include "./base/base_uses.f90" 

  IMPLICIT NONE
                                                                                        
        PRIVATE
                                                                                                       
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='non_adiabatic_utils'
                                                                                                       
        PUBLIC :: decohere_wf,      &  
                  propagate_wf

CONTAINS
 
!  SUBROUTINE decohere_wf(adiab_section, sh_print_env, aom_env,sh_env, deco_keyword, itimes, successful_hop, &
 !                         deltatMD, mol,velocities_t, local_particles, K_tot)
  SUBROUTINE decohere_wf(adiab_section, sh_print_env, aom_env,sh_env, deco_keyword, itimes, successful_hop, &
                         deltatMD, mol,velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, atomic_kind_set, K_tot)
  !SUBROUTINE decohere_wf(coupling_section, sh_print_env, aom_env,sh_env, deco_keyword)


        TYPE(section_vals_type), POINTER               :: adiab_section 
        TYPE(sh_print_env_type)                        :: sh_print_env
        TYPE(aom_env_type)                             :: aom_env
        TYPE(sh_env_type)                              :: sh_env
        REAL(KIND=dp), INTENT(IN)                      :: deltatMD
        REAL(KIND=dp)                                  :: C_param, E_0_param
       ! REAL(kind=dp), INTENT(IN)                      :: kinetic_energy_prime

        CHARACTER(len=*), PARAMETER :: routineN ='decohere_wf', &         
                                       routineP = moduleN//':'//routineN

        CHARACTER(LEN=default_string_length)           :: deco_keyword
        INTEGER, POINTER                               :: itimes
        LOGICAL                                        ::  successful_hop
        INTEGER                                        :: mol

!MODIFY
        REAL(kind=dp), DIMENSION(:,:), POINTER              :: velocities_com_tp
        TYPE(distribution_1d_type), POINTER                :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
                  POINTER, OPTIONAL                        :: atomic_kind_set

        REAL(kind=dp), INTENT(OUT)                         :: K_tot

                                                         
          CALL section_vals_val_get(adiab_section,"DECOHERENCE_CORRECTIONS",c_val=deco_keyword)
        
        SELECT CASE(deco_keyword)                                                          
        CASE("NO_DECO_CORR")

             PRINT *, 'NO DECOHERENCE CORRECTIONS'
    
        CASE("INSTANT_COLLAPSE")
          
          CALL wf_collapse(aom_env, sh_print_env, sh_env, itimes,successful_hop)
         ! CALL wf_collapse(itimes, sh_print_env, aom_env, collapse, successful_hop, decohere, ctin, C_tp, active_state)

        CASE("DAMPING")

          !CALL section_vals_val_get(adiab_section,"EDC_C",r_val=C_keyword)
          CALL section_vals_val_get(adiab_section,"EDC_C",r_val=C_param)
          !CALL section_vals_val_get(adiab_section,"EDC_E0",r_val=E_0_keyword)
          CALL section_vals_val_get(adiab_section,"EDC_E0",r_val=E_0_param)
!MODIFY     
          CALL calc_kinetic_tot(velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &                                               
                      atomic_kind_set, K_tot)
PRINT *, "HEYYYYYYYYYYYYYYYYYYYYYYY KINETIC IS", K_tot

          CALL wf_damping(C_param, E_0_param, sh_env, mol, deltatMD, K_tot)
 
        CASE DEFAULT                                           
          CPABORT("There is no method to add the decoherence.") 
        END SELECT                                             

  END SUBROUTINE decohere_wf

  SUBROUTINE wf_collapse(aom_env,sh_print_env, sh_env, itimes, successful_hop)
 ! SUBROUTINE wf_collapse(itimes, sh_print_env, aom_env, collapse, successful_hop, decohere, ctin, C_tp, active_state)
  
   ! COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)                 :: ctin
   ! REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)                  :: C_tp
   ! INTEGER, INTENT(IN)                                        :: active_state
    
    TYPE(sh_print_env_type)                        :: sh_print_env
    TYPE(aom_env_type)                             :: aom_env
    TYPE(sh_env_type)                              :: sh_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
    
    INTEGER, POINTER                               :: itimes
    LOGICAL                                        ::  successful_hop
   ! LOGICAL                                        :: decohere, collapse, successful_hop


 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (successful_hop) aom_env%futur_decohere = .TRUE.               
!> Check if the collapse key word is TRUE in the input                     
         !sh_print_env%collapse = sh_env%collapse
        ! if (sh_env%collapse .and. aom_env%futur_decohere .and. sh_env%decohere) then                                  
         if (aom_env%futur_decohere .and. sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
!print*, "WE ARE IN DECOHERE"                           
 !                 print*, sh_env%ctin(:)
 !                 print*, sh_env%active_state
!print*, "UNTIL HERE, I M FINE"
!stop
                  print*, sh_env%c_tp(:,:)
                     sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
                  aom_env%futur_decohere = .FALSE.                         
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF

!!!!sh_print_env%decohere = .FALSE.                                      
!!!!     IF (itimes.ge.1) THEN                                                
!!!!        !Finally, collapse the wavefunction if called for                 
!!!!        if (successful_hop) aom_env%futur_decohere = .TRUE.               
!!!!> Check if the collapse key word is TRUE in the input                     
!!!!        sh_print_env%collapse = collapse
!!!!        if (collapse .and. aom_env%futur_decohere .and. decohere) then                                  
!!!!> Check if the surface_hop subroutine asks us to do the decoherence       
!!!!                 sh_print_env%decohere = .TRUE.                           
!!!!                     ctin(:) = C_tp(: ,active_state)
!!!!                 aom_env%futur_decohere = .FALSE.                         
!!!!        end if                                                            
!!!!        CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
!!!!     ENDIF
!!!!
!!!!   print *, 'WF_COLLAPSE JUST CALLED, THIS WAS A TEST' 

  END SUBROUTINE wf_collapse                                                               
                                                 

  SUBROUTINE  wf_damping(C_param, E_0_param, sh_env, mol, deltatMD, K_tot)                    

    REAL(kind=dp), INTENT(IN)                    :: K_tot
    REAL(KIND=dp), INTENT(IN)                      :: deltatMD
    TYPE(sh_env_type)                              :: sh_env
!    TYPE(aom_env_type)                             :: aom_env
   ! REAL(kind=dp), INTENT(IN)                :: kinetic_energy_prime 

    CHARACTER(len=*), PARAMETER :: routineN ='wf_damping', &         
    routineP = moduleN//':'//routineN   
                                                        
    INTEGER  :: j,mol,i
    REAL(KIND=dp)   :: C_param, E_0_param
   ! REAL(KIND=dp), PARAMETER :: C, E_0
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tau
    REAL(KIND=dp)      :: sum!, sum1 
    COMPLEX(KIND=dp), dimension(:), allocatable :: ct
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop

allocate(pop(mol)) !TEST

allocate(ct(mol))  ! adiabatic states
allocate(tau(mol)) ! decoherence time (EDC method)

!set value C and E_O 
!C = 1.0D0
!E_0 = 0.1D0
!kinetic_energy_prime = 1.0D0

sum = 0.0D0

PRINT *, "CTIN diabatic BEFORE", sh_env%ctin
PRINT *, "ACTIVE_STATE IS", sh_env%active_state

do i=1,mol 
        if (i .ne. sh_env%active_state) then
            pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
            ct(i) = cmplx(0.0D0, 0.0D0,dp)                           
            tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))* &
                                     (C_param +(E_0_param/K_tot))
            PRINT *, "TAU is", tau(i) 
             do j=1,mol
                     ct(i) = ct(i) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*sh_env%ctin(j)
             end do

            PRINT *, "adiab_coeef_before damp is", ct(i)
            pop(i) = ct(i)*conjg(ct(i))

            PRINT *, "adiab_POP_before damp is", pop(i)
             ct(i) = ct(i)*exp(-(deltatMD/tau(i)))

             pop(i) = ct(i)*conjg(ct(i))

             PRINT *, "adiab_coeef_after damp is", ct(i)
             PRINT *, "adiab_POP_after damp is", pop(i)
             sum = sum + abs(ct(i))*abs(ct(i))
        end if
end do

pop(sh_env%active_state)= cmplx(0.0D0, 0.0D0, dp)
ct(sh_env%active_state) = cmplx(0.0D0, 0.0D0, dp)

do i=1,mol                                                       
       ! ct(sh_env%active_state) = ct(sh_env%active_state) +
       ! cmplx(C2(i,old_state), 0.0D0,dp)*c(i)  
        ct(sh_env%active_state) = ct(sh_env%active_state) + &
          cmplx(sh_env%C_tp(i,sh_env%active_state), 0.0D0,dp)*sh_env%ctin(i)
end do  
           PRINT *, "ACTIVE_adiab_coeef_before damp is", ct(sh_env%active_state)

pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))
           PRINT *, "ACTIVE_POP_adiab_coeef_before damp is", pop(sh_env%active_state)

        ct(sh_env%active_state)= ct(sh_env%active_state)* &
           sqrt((1-sum)/(abs(ct(sh_env%active_state))*abs(ct(sh_env%active_state))))

           PRINT *, "ACTIVE_adiab_coeef_AFTER damp is", ct(sh_env%active_state)

pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))
           PRINT *, "ACTIVE_POP_adiab_coeef_AFTER damp is", pop(sh_env%active_state)

do j=1,mol
       pop(j) = cmplx(0.0D0, 0.0D0,dp) !TEST
       do i=1,mol 
           sh_env%ctin(j) = sh_env%ctin(j) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*ct(i)
       end do
PRINT *, "CTIN_diabatic AFTER", sh_env%ctin(j)
      pop(j)=sh_env%ctin(j)*conjg(sh_env%ctin(j))
PRINT *, "POP_diabatic =", pop(j)
end do


!PRINT *, "CTIN_diabati-----", sh_env%ctin
!PRINT *, "POP_diabatic =", pop
!PRINT *, "ADIAB-------", ct

!do i=1,mol
 !  sum1 = sum1 + abs(ct(i))*abs(ct(i))
!end do 

!print *, "POPUL", sum1

deallocate(ct) 
deallocate(tau)
deallocate(pop)
    print *, 'I AM HEREEEEEEEEEEEEEEEEEEEEEEEE' 

!do i=1,mol
!        if (i .ne. active_state) then 
!            ct(i) = cmplx(0.0D0, 0.0D0,dp)    
!            tau(i)= (1/abs(E(i) - E(active_state)))*(C +(E_0/KE)) !NB solve zero 
!                                                                  !division
!             do j=1,mol
 !                    ct(i) = ct(i) + cmplx(C2(j,i), 0.0D0,dp)*c(j)
!             end do 
!             ct(i) = ct(i)*exp(-(deltatMD/tau(i)))        
!             sum = sum + ct(i)*ct(i)                                                        
!        end if
!end do 


!ct(active_state) = cmplx(0.0D0, 0.0D0, dp)
!do i=1,mol 
       ! ct(active_state) = ct(active_state) + cmplx(C2(i,old_state),
       ! 0.0D0,dp)*c(i)  
!        ct(active_state) = ct(active_state) + cmplx(C2(i,active_state),
!0.0D0,dp)*c(i)
!end do 
!        ct(active_state)=
!ct(active_state)*sqrt((1-sum)/(ct(active_state)*ct(active_state))
!
!do j=1,mol 
!       do i=1,mol                  
!           c(j) = c(j) + cmplx(C2(j,i), 0.0D0,dp)*ct(j)
!       end do
!end do

  END SUBROUTINE wf_damping 
 

  SUBROUTINE calc_kinetic_tot(array, local_particles, atomic_kinds, &
                                              nparticle_kind, atomic_kind_set, K_tot)
        REAL(kind=dp), DIMENSION(:,:), POINTER    :: array
        !LOGICAL                                  :: keyword_translate
        !REAL(KIND=dp), DIMENSION(:,:), POINTER   :: translate_array
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER, OPTIONAL                             :: atomic_kind_set

        REAL(kind=dp), INTENT(OUT)               :: K_tot

        CHARACTER(len=*), PARAMETER :: routineN = 'calc_kinetic_tot', &
        routineP = moduleN//':'//routineN

        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                             :: my_atomic_kind_set
        INTEGER                                  :: my_nparticle_kind
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
       ! REAL(kind=dp), DIMENSION(3), INTENT(OUT)              :: sum_array
       ! REAL(KIND=dp)                            :: sum_mass

PRINT *, "VELOCITIES ARRAY", array

           IF (PRESENT(atomic_kinds)) THEN
              my_nparticle_kind = atomic_kinds%n_els
              my_atomic_kind_set => atomic_kinds%els
           ELSE
              IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
                 my_nparticle_kind = nparticle_kind
                 my_atomic_kind_set => atomic_kind_set
              END IF
           END IF

PRINT *, "HEYYYY NPARTICLE IS =", nparticle_kind, "instead", my_nparticle_kind
          
          !K_tot = 0.0D0
          ! sum_array = 0.0D0
          ! sum_mass  = 0.0D0
           DO iparticle_kind= 1, my_nparticle_kind
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)

PRINT *, "nparticle_local", nparticle_local

              DO iparticle_local= 1, nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 K_tot = K_tot + mass*( &
                                   array(1, iparticle)**2 + &
                                   array(2, iparticle)**2 + &
                                   array(3, iparticle)**2  ) 
PRINT *, "K_tot is", K_tot         
                   
                 !sum_array(:) = sum_array(:) +   &
                   ! mass*array(:, iparticle)
                ! sum_mass = sum_mass + mass
              END DO
           END DO
           K_tot = 0.5*K_tot

! sum_array = sum_array / sum_mass

!           DO iparticle_kind=1, my_nparticle_kind
!              atomic_kind => my_atomic_kind_set(iparticle_kind)
!              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!              nparticle_local = local_particles%n_el(iparticle_kind)
!              DO iparticle_local=1,nparticle_local
!                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
!                 translate_array(:, iparticle) = array(:, iparticle) - sum_array(:)
!              END DO   
!           END DO   
!        ELSE
!           translate_array = array
!        ENDIF 

   END SUBROUTINE calc_kinetic_tot
 

  SUBROUTINE propagate_wf(H_t, H_tp,  Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm,     &
                          ctin, ctout, ctotal, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat)

     COMPLEX(KIND=dp), dimension(:,:), INTENT(IN)    :: H_t, H_tp
     REAL(KIND=dp), dimension(:,:), INTENT(INOUT)       :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm
     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)      :: ctin, ctout, ctotal
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 


     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf', &
     routineP = moduleN//':'//routineN

     REAL(KIND=dp), dimension(:,:), allocatable    :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_0, H_2, X1, X2, X12
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: Ktot, K2, K3, K4
     REAL(KIND=dp), dimension(:,:), allocatable      :: Imat 


     INTEGER                                         :: ELsteps   
     INTEGER, POINTER                                :: itimes
     TYPE(section_vals_type), POINTER                :: adiab_section
     TYPE(aom_env_type)                              :: aom_env
     REAL(KIND=dp)                                   :: deltatMD, deltatE
     REAL(KIND=dp)                                   :: ctmag
     REAL(KIND=dp)                                   :: times 
     COMPLEX(KIND=dp)                                :: alpha, beta, ctcomp
     integer                                         :: t2
     integer                                         :: i,j, inc, mat 
     character(len=1)                                :: zgemm_type


     !Various Allocation


     !Imat is an identity matrix of the dimensions of the other square matrices
     allocate(Imat(mat,mat))
     !define once, then use whenever
     Imat = 0.0D0
     FORALL(j = 1:mat) Imat(j,j) = 1     
     
     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(Stptpdiff(mat,mat))
     allocate(Sttpdiff(mat,mat))
     allocate(Sttdiff(mat,mat))
     allocate(Sttmdiff(mat,mat))
     allocate(Stmtmdiff(mat,mat))


     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H_t etc are always -i*H.
     !This is why they're complex variables
     allocate(NACE(mat,mat))
     allocate(Hdiff(mat,mat))
     allocate(H_0(mat,mat))
     allocate(H_1(mat,mat))
     allocate(H_2(mat,mat))
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

     !used to give correct commands to ZGEMX subroutines later
     zgemm_type = 'N'


             Hdiff = (deltatE/deltatMD)*(H_tp - H_t)
             Sttdiff = (deltatE/deltatMD)*(Stptp - Stt)
             Sttpdiff = (deltatE/deltatMD)*(Sttp - Stmt)
             Sttmdiff = (deltatE/deltatMD)*(Stpt - Sttm)
             Stptpdiff = Sttdiff
             Stmtmdiff = (deltatE/deltatMD)*(Stt - Stmtm)
             H_0  = H_t
             H_2 = H_tp

             do t2=1,int(deltatMD/deltatE)
                    call makeX(mat, H_t, Stt, Sttp, Sttm, Stptp, Stmtm, X1, deltatMD)


                    !create timestepped equivalents
                    H_1 = H_0 + 0.5D0*Hdiff
                    Stt = Stt + 0.5D0*Sttdiff
                    Sttp = Sttp + 0.5D0*Sttpdiff
                    Sttm = Sttm + 0.5D0*Sttmdiff
                    Stptp = Stptp + 0.5D0*Stptpdiff
                    Stmtm = Stmtm + 0.5D0*Stmtmdiff
                    call makeX(mat, H_1, Stt, Sttp, Sttm, Stptp, Stmtm, X12, deltatMD)
     
                    H_2 = H_1 + 0.5D0*Hdiff
                    Stt = Stt + 0.5D0*Sttdiff
                    Sttp = Sttp + 0.5D0*Sttpdiff
                    Sttm = Sttm + 0.5D0*Sttmdiff
                    Stptp = Stptp + 0.5D0*Stptpdiff
                    Stmtm = Stmtm + 0.5D0*Stmtmdiff
!                    call makeX(mat, H_2, Stt, Sttp, Sttm, Stptp, Stmtm, X2, deltatMD)
                    call makeX(mat, H_2, Stt, Sttp, Sttm, Stptp, Stmtm, X2, deltatMD, NACE)
     
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
                    call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X12, mat, K2, &
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
                                        aom_env%population(1:mat) =&
                                            abs(ctotal(1:mat)) 
                    end if

                    !update
                    H_0 = H_2

                    ctin = ctout
PRINT *, "diab_COEFF", ctin
                    ctout = cmplx(0.0D0, 0.0D0,dp)

                    times = DBLE(itimes) + t2*1/DBLE(ELsteps)
                    CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, &
                                           times)
                    CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,itimes, deltatMD, &
                                         times)
                    CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes,deltatMD, &
                                                    times)

             end do

!END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP

!DEALLOCATE
       DEALLOCATE(NACE)
       deallocate(H_0, Hdiff, H_1, H_2, X1, X2, X12, Ktot, K2, K3, K4)
       deallocate(Sttdiff, Stptpdiff, Stmtmdiff, Sttpdiff, Sttmdiff)
       deallocate(Imat)
!END DEALLOCATE

  END SUBROUTINE propagate_wf                                                                

END MODULE non_adiabatic_utils 
                       
