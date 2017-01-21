!-----------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations           
!   !                                                                           
!   Copyright (C) 2000 - 2015  CP2K developers group                            
!   !                                                                           
!-----------------------------------------------------------------------------! 

!!!!!!!!!!!!!!!!!!!! CHANGE_SG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE non_adiabatic_utils
     
     USE kinds,                        ONLY: default_string_length,&                                   
                                             dp                                                        

     USE input_section_types,          ONLY:  section_vals_type

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

     USE aom_utils,                    ONLY: makeX,    &
                                             adiabaticity,   &
                                             surfacehop,     &
                                             calcforces,     &
                                             read_listcoeff, &
                                             choose_firstadiabat, &
                                             simple_vel_rescaling, &
                                             translate_com
#include "./base/base_uses.f90" 

  IMPLICIT NONE
                                                                                        
        PRIVATE
                                                                                                       
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='non_adiabatic_utils'
                                                                                                       
        PUBLIC :: wf_collapse,      &
                  propagate_wf

CONTAINS

  SUBROUTINE wf_collapse(itimes, sh_print_env, aom_env, collapse, successful_hop, decohere, ctin, C_tp, active_state)
  
    COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)                 :: ctin
    REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)                  :: C_tp
    INTEGER, INTENT(IN)                                        :: active_state
    
    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
    
    INTEGER, POINTER                               :: itimes
    TYPE(sh_print_env_type)                        :: sh_print_env
    TYPE(aom_env_type)                             :: aom_env
    LOGICAL                                        :: decohere, collapse, successful_hop



 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (successful_hop) aom_env%futur_decohere = .TRUE.               
!> Check if the collapse key word is TRUE in the input                     
         sh_print_env%collapse = collapse
         if (collapse .and. aom_env%futur_decohere .and. decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.                           
                      ctin(:) = C_tp(: ,active_state)
                  aom_env%futur_decohere = .FALSE.                         
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF

!   print *, 'WF_COLLAPSE JUST CALLED, THIS WAS A TEST' 

  END SUBROUTINE wf_collapse                                                               
                                                 

!  SUBROUTINE test_sub()                           
!    implicit none                               
!    !this subroutine is just a test             
!                                               
!    print *, 'I AM HEREEEEEEEEEEEEEEEEEEEEEEEE' 

!  END SUBROUTINE test_sub 
  

  SUBROUTINE propagate_wf(H_t, H_tp,  Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm,     &
                          ctin, ctout, ctotal, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            active_state, mat)

     COMPLEX(KIND=dp), dimension(:,:), INTENT(IN)    :: H_t, H_tp
     REAL(KIND=dp), dimension(:,:), INTENT(INOUT)       :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm
     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)      :: ctin, ctout, ctotal
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 


     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf', &
     routineP = moduleN//':'//routineN

     REAL(KIND=dp), dimension(:,:), allocatable    :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
  !   REAL(KIND=dp), dimension(:,:), allocatable      :: NACE
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
     integer                                         :: i,j, inc, active_state, mat 
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
!     allocate(NACE(mat,mat))
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
                    call makeX(mat, H_2, Stt, Sttp, Sttm, Stptp, Stmtm, X2, deltatMD)
!                    call makeX(mat, H_2, Stt, Sttp, Sttm, Stptp, Stmtm, X2, deltatMD, NACE)
     
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
!       DEALLOCATE(NACE)
       deallocate(H_0, Hdiff, H_1, H_2, X1, X2, X12, Ktot, K2, K3, K4)
       deallocate(Sttdiff, Stptpdiff, Stmtmdiff, Sttpdiff, Sttmdiff)
       deallocate(Imat)
!END DEALLOCATE

  END SUBROUTINE propagate_wf                                                                

END MODULE non_adiabatic_utils 
                       
