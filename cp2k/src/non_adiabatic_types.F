
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   !
!   Copyright (C) 2000 - 2015  CP2K developers group
!   !
!-----------------------------------------------------------------------------!

!!!!!!!!!! CHANGE_SG !!!!!!!!!!!!!

MODULE non_adiabatic_types
   
  USE kinds,                           ONLY: default_string_length, &
                                             dp
  USE input_section_types,             ONLY: section_vals_type



IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'non_adiabatic_types'



   TYPE sh_env_type
!> Variables used for the decoherence methods
   
     LOGICAL                                     :: decohere

     !COMPLEX(KIND=dp), dimension(:), POINTER     :: ctin       
     COMPLEX(KIND=dp), dimension(:), POINTER     :: ctin, cttm, ctt       
     REAL(KIND=dp), DIMENSION(:,:), POINTER      :: C_tp          
     INTEGER                                     :: active_state, state_min  
     REAL(KIND=dp), dimension(:), POINTER        :: E_tp, E_t, enerdiff                                                                           


                    
   END TYPE sh_env_type
  
   PUBLIC :: sh_env_type


END MODULE non_adiabatic_types

