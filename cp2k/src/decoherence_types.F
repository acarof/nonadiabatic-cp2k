!-----------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations           
!   !                                                                           
!   Copyright (C) 2000 - 2015  CP2K developers group                            
!   !                                                                           
!-----------------------------------------------------------------------------! 

MODULE decoherence_types
     
     USE kinds,                        ONLY: default_string_length,&                                   
                                             dp                                                        

     USE aom_print,                    ONLY: adiab_print_sh_info

     USE aom_types,                    ONLY: aom_env_type, &      
                                                sh_print_env_type 

#include "./base/base_uses.f90" 

  IMPLICIT NONE
                                                                                        
        PRIVATE
                                                                                                       
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='decoherence_types'
                                                                                                       
        PUBLIC :: wf_collapse

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
END MODULE decoherence_types                        
  

