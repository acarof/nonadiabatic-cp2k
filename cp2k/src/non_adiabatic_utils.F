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
                                             adiab_print_nace,        &
                                             adiab_print_nacv        
     USE aom_types,                    ONLY: aom_env_type, &      
                                                sh_print_env_type 

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
                                             cal_delta_energy_nacv, &
                                             calc_diabatic_nace


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
                  propagate_wf,    &
                  propagate_adiab_wf, &
                  prop_wf_HS_corr, &
                  calc_kinetic_tot, &
                 classically_forbidden_hops

CONTAINS


SUBROUTINE classically_forbidden_hops(active_state, aom_env, mat, E_tp, &
                  velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                          atomic_kind_set, energy_conservation, class_forbidden_hop)

        !REAL(kind=dp), DIMENSION(:,:), POINTER              :: velocities_com_tp
     TYPE(distribution_1d_type), POINTER                :: local_particles
     TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
     INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
     TYPE(atomic_kind_type), DIMENSION(:), &
                   POINTER, OPTIONAL                        :: atomic_kind_set
 
    
     REAL(KIND=dp), DIMENSION(:), INTENT(IN)         ::  E_tp
     TYPE(aom_env_type)                                 :: aom_env
     INTEGER, INTENT(IN)                             :: active_state, mat
     REAL(kind=dp), DIMENSION(:,:), POINTER       :: velocities_com_tp


     REAL(KIND=dp), dimension(:), allocatable, INTENT(OUT)    :: energy_conservation
     INTEGER, dimension(:), allocatable, INTENT(OUT)          :: class_forbidden_hop

     INTEGER :: j
     REAL(kind=dp), DIMENSION(:,:), allocatable     :: hop_nacv
     REAL(kind=dp)                               :: delta_energy

!    total_kin_energy, total_pot_energy
                                        !   total_kin_array, total_pot_array
     REAL(kind=dp)                               :: A, B
     REAL(kind=dp)                               :: K_tot, E_cutoff  
     
     ALLOCATE(energy_conservation(mat)) 
     ALLOCATE(class_forbidden_hop(mat))


     CALL calc_kinetic_tot(velocities_com_tp, local_particles, atomic_kinds, &
                  nparticle_kind, atomic_kind_set, K_tot)
     E_cutoff = 0.0D0
     E_cutoff = E_tp(active_state) + K_tot  !THIS CUTOFF IS TO LARGE ALL THE
                                            !STATE ARE BELOW IT
!     E_cutoff = 4.6612978986598338E-002
print *, "KINETIC_ENERGY_QUNATUMM", K_tot  
print *, "ENERGY_ACTIVE_tp", E_tp(active_state)
print *,  "CUTOFF_ENERGY", E_cutoff
print *, "POTENTIAL_STATE9", E_tp(9)
     do j=1, mat
           if ((j .gt. active_state) .and. (E_tp(j) .le. E_cutoff)) then
               
              CALL calc_hop_nacv(aom_env, hop_nacv, active_state, j, fast=.TRUE.)
!print *, "ENERGY_ACTIVE_t", E_t(active_state)
!print *, "ENERGY_TO_J_t", E_t(j)
!print *, "ENERGY_ACTIVE_tp", E_tp(active_state)
print *, "ENERGY_TO_J_tp", E_tp(j)
print *, "velocities_com_tp", velocities_com_tp
              CALL evaluate_delta_energy_nacv(velocities_com_tp, local_particles, atomic_kinds, &
                                     nparticle_kind, atomic_kind_set, &
                                    aom_env, E_tp, hop_nacv, delta_energy, A, B, & 
                                              active_state, j)
              !total_kin_array(j) = total_kin_energy
              !total_pot_array(j) = total_pot_energy
              energy_conservation(j) = delta_energy
              if (energy_conservation(j) .lt. 0.0D0) then 
                   class_forbidden_hop(j) = j
              else 
                   class_forbidden_hop(j) = -1
              end if
           else if ((j .gt. active_state) .and. (E_tp(j) .gt. E_cutoff)) then
print *, "NOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"
              class_forbidden_hop(j) = j
              energy_conservation(j) = -2  
           else if (j .le. active_state) then
              !total_pot_array(j) = 0.0D0
              !total_kin_array(j) = 0.0D0 
              energy_conservation(j) = 0.0D0
              class_forbidden_hop(j) = -1
           end if  
print *, "energy_conservation(j)", energy_conservation(j)
print *, "ACTIVE_STATE_NAC", active_state
print *, "TO_STATE_j", j
print *, "HOP_NACV", hop_nacv
     end do
!print *, "KE_ARRAY", total_kin_array
!print *, "PE_ARRAY", total_pot_array
!print *, "DIFF", total_pot_array-total_kin_array
print *, "ENERGY_CONSERV_ARRAY", energy_conservation
print *, "CLASS_FORBIDD_ARRAY", class_forbidden_hop

END SUBROUTINE classically_forbidden_hops



   SUBROUTINE evaluate_delta_energy_nacv(velocities, local_particles, atomic_kinds, &
                                     nparticle_kind, atomic_kind_set, &
                                    aom_env,E_tp, hop_nacv, delta_energy, A, B, &
                                               i, j) 

        REAL(kind=dp), DIMENSION(:,:), POINTER             :: velocities
        TYPE(distribution_1d_type), POINTER                :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER, OPTIONAL                             :: atomic_kind_set
        TYPE(aom_env_type)                                 :: aom_env
        INTEGER, OPTIONAL, INTENT(IN)                      :: i, j
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)          :: hop_nacv
        REAL(kind=dp), INTENT(OUT)               :: delta_energy
!                                                   total_kin_energy, total_pot_energy 
        REAL(kind=dp), INTENT(OUT)                         :: A, B
        REAL(KIND=dp), DIMENSION(:), INTENT(IN)            ::  E_tp

        !LOGICAL, dimension(:), INTENT(IN)           :: ext_pres

        CHARACTER(len=*), PARAMETER :: routineN = 'evaluate_delta_energy_nacv', &
        routineP = moduleN//':'//routineN

        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                             :: my_atomic_kind_set
        !REAL(kind=dp)                            :: rescale_factor
        INTEGER                                  :: my_nparticle_kind
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind, new_state, old_state
        INTEGER                                  :: nparticle_local, coord
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind


!From Fabiano E, Keal TW, Thiel W. Implementation of surface hopping molecular dynamics using semiempirical methods. 
!     Chem Phys 2008, 349:334–347.


        IF (PRESENT(atomic_kinds)) THEN
           my_nparticle_kind = atomic_kinds%n_els
           my_atomic_kind_set => atomic_kinds%els
        ELSE
           IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
              my_nparticle_kind = nparticle_kind
              my_atomic_kind_set => atomic_kind_set
           END IF
        END IF


!print *, "HOP_NACV_INSIDE_DELT", hop_nacv
!print *, "VELOCITIES_INSIDE", velocities
           A = 0
           B = 0
           DO iparticle_kind=1, my_nparticle_kind
!print *, "nparticle_kind", my_nparticle_kind
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
!print *, "nparticle_local", nparticle_local
print *, "MASS", mass
!stop
              DO iparticle_local=1,nparticle_local
              ! if (ext_pres(iparticle_local)) then
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
!print *, "iparticle", iparticle
!stop
                 DO coord = 1, 3
                    IF ((velocities(1, iparticle)**2 + velocities(2, iparticle)**2 + &
                       velocities(3, iparticle)**2 ).NE.0.0d0) THEN
!print *, "VELOCITIES_INSIDE_1_2", velocities(1,1)
                       A = A + hop_nacv(coord, iparticle)**2 / mass
                    ELSE IF ( SUM(velocities).EQ.0) THEN
                       A = A + hop_nacv(coord, iparticle)**2 / mass
                    ENDIF
                    B = B + velocities(coord, iparticle) * hop_nacv(coord, iparticle)
                 ENDDO
              END DO
           END DO
           A = 0.5*A

           IF (PRESENT(j)) THEN
              new_state = j
           ELSE    
              new_state = aom_env%new_state
           ENDIF
           IF (PRESENT(i)) THEN
              old_state = i
           ELSE
              old_state = aom_env%active_state
           ENDIF
print *, "ENERGY_ACTIVE_IN_SUB", E_tp(old_state)
print *, "ENERGY_TO_j_IN_SUB", E_tp(new_state)
           delta_energy = B**2 - 4 * A * (  &
                     E_tp(new_state) - E_tp(old_state))

!>kinetic_energy calculation eq 23 Fabiano 2008
!> From Fabiano E, Keal TW, Thiel W. Implementation of surface hopping molecular
!dynamics using semiempirical methods. 
!>     Chem Phys 2008, 349:334–347.

!>  B = -B, compare to the paper
!print *, "A", A
!print *, "B", B
!print *, "delta_energy", delta_energy
!rescale_factor = (- B - SQRT(delta_energy))/(2 * A)
!print *, "RESCALE_BEFORE", rescale_factor 
!          IF (A == 0) THEN
!             rescale_factor = 0
!          ELSE
!             IF (B < 0) THEN
!                rescale_factor = (- B - SQRT(delta_energy))/(2 * A)
!             ELSE
!                rescale_factor = (- B + SQRT(delta_energy))/(2 * A)
!             ENDIF
!          ENDIF
!print *, "rescale_factor", rescale_factor
!
!          total_kin_energy = rescale_factor*rescale_factor*A + &
!                             rescale_factor*B
!          total_pot_energy = new_state - old_state
!        
!print *, "total_kin_energy", total_kin_energy
!print *, "total_pot_energy", total_pot_energy
!print *, "ENERGY_ACTIVE_IN_SUB", aom_env%adiab_energies_t(old_state)
!print *, "ENERGY_TO_j_IN_SUB", aom_env%adiab_energies_t(new_state)
         !  delta_energy = B**2 - 4 * A * (  &
          !           aom_env%adiab_energies_t(new_state) - &
           !          aom_env%adiab_energies_t(old_state))

   END SUBROUTINE evaluate_delta_energy_nacv




 
  SUBROUTINE decohere_wf(adiab_section, sh_print_env, aom_env,sh_env, deco_keyword, itimes, &
                         deltatMD, mol,velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                         atomic_kind_set, K_tot, ext_pres)


        TYPE(section_vals_type), POINTER               :: adiab_section 
        TYPE(sh_print_env_type)                        :: sh_print_env
        TYPE(aom_env_type)                             :: aom_env
        TYPE(sh_env_type)                              :: sh_env
        REAL(KIND=dp), INTENT(IN)                      :: deltatMD
        REAL(KIND=dp)                                  :: C_param, E_0_param
        LOGICAL, dimension(:), POINTER, INTENT(IN)     :: ext_pres

        CHARACTER(len=*), PARAMETER :: routineN ='decohere_wf', &         
                                       routineP = moduleN//':'//routineN

        CHARACTER(LEN=default_string_length)           :: deco_keyword
        CHARACTER(LEN=default_string_length)           :: wfprop_keyword
        INTEGER, POINTER                               :: itimes
        INTEGER                                        :: mol

!MODIFY
        REAL(kind=dp), DIMENSION(:,:), POINTER              :: velocities_com_tp
        TYPE(distribution_1d_type), POINTER                :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
                  POINTER, OPTIONAL                        :: atomic_kind_set

        REAL(kind=dp), INTENT(OUT)                         :: K_tot

        CALL section_vals_val_get(adiab_section,"REPRESENTATION", &
                                                     c_val=wfprop_keyword)
        CALL section_vals_val_get(adiab_section,"DECOHERENCE_CORRECTIONS",c_val=deco_keyword)
        
        SELECT CASE(wfprop_keyword)
        CASE("DIABATIC_BASIS")   
            SELECT CASE(deco_keyword)                                                          
            CASE("NO_DECO_CORR")
                 PRINT *, 'NO DECOHERENCE CORRECTIONS'
            CASE("INSTANT_COLLAPSE")
              CALL wf_collapse(aom_env, sh_print_env, sh_env, itimes)
            CASE("TRESH_ONLY_COLLAPSE")  
              CALL  wf_tresh_collaps(sh_print_env, sh_env, itimes)
            CASE("DAMPING")
              CALL section_vals_val_get(adiab_section,"EDC_C",r_val=C_param)
              CALL section_vals_val_get(adiab_section,"EDC_E0",r_val=E_0_param)
              CALL calc_kinetic_tot(velocities_com_tp, local_particles, atomic_kinds, &
                                    nparticle_kind, atomic_kind_set, K_tot, ext_pres)
              CALL wf_damping(C_param, E_0_param, sh_env, mol, deltatMD, K_tot, sh_print_env)
            CASE DEFAULT                                           
              CPABORT("There is no method to add the decoherence.") 
            END SELECT
        CASE("ADIABATIC_BASIS")
            SELECT CASE(deco_keyword)                                                          
            CASE("NO_DECO_CORR")
                 PRINT *, 'NO DECOHERENCE CORRECTIONS'
            CASE("INSTANT_COLLAPSE")
              !CALL wf_collapse(aom_env, sh_print_env, sh_env, itimes)
              CALL wf_collapse_adiab_prop(aom_env,sh_print_env, sh_env, itimes)
            CASE("TRESH_ONLY_COLLAPSE")  
              print *, "NOT IMPLEMENTED"
            CASE("DAMPING")
              CALL section_vals_val_get(adiab_section,"EDC_C",r_val=C_param)
              CALL section_vals_val_get(adiab_section,"EDC_E0",r_val=E_0_param)
              CALL calc_kinetic_tot(velocities_com_tp, local_particles, atomic_kinds, &
                                    nparticle_kind, atomic_kind_set, K_tot, ext_pres)
              !CALL wf_damping(C_param, E_0_param, sh_env, mol, deltatMD, K_tot, sh_print_env)
              CALL wf_damping_adiab_prop(C_param, E_0_param, sh_env, mol, deltatMD, & 
                                                              K_tot, sh_print_env)
            CASE DEFAULT                                           
              CPABORT("There is no method to add the decoherence.") 
            END SELECT
        CASE DEFAULT
           CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
        END SELECT

  END SUBROUTINE decohere_wf


  SUBROUTINE wf_collapse_adiab_prop(aom_env,sh_print_env, sh_env, itimes)
    TYPE(sh_print_env_type)                        :: sh_print_env
    TYPE(aom_env_type)                             :: aom_env
    TYPE(sh_env_type)                              :: sh_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
    
    INTEGER, POINTER                               :: itimes

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
!         print*, "future", aom_env%futur_decohere, sh_env%decohere
         if (aom_env%futur_decohere .and. sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
!                  print*, "deco!", sh_env%c_tp(:,:)
!The wafunction is collapsed in a pure state
                     sh_env%ctin(:) =  CMPLX(0.0D0, 0.0D0,dp)
                     sh_env%ctin(sh_env%active_state) = CMPLX(1.0D0, 0.0D0,dp)
                     !sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
                  aom_env%futur_decohere = .FALSE.                         
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF
  END SUBROUTINE wf_collapse_adiab_prop


  SUBROUTINE wf_collapse(aom_env,sh_print_env, sh_env, itimes)
    TYPE(sh_print_env_type)                        :: sh_print_env
    TYPE(aom_env_type)                             :: aom_env
    TYPE(sh_env_type)                              :: sh_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
    
    INTEGER, POINTER                               :: itimes

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
!         print*, "future", aom_env%futur_decohere, sh_env%decohere
         if (aom_env%futur_decohere .and. sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
!                  print*, "deco!", sh_env%c_tp(:,:)
                     sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
                  aom_env%futur_decohere = .FALSE.                         
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF
  END SUBROUTINE wf_collapse 


  SUBROUTINE wf_tresh_collaps(sh_print_env, sh_env, itimes)
    TYPE(sh_print_env_type)                        :: sh_print_env
    TYPE(sh_env_type)                              :: sh_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
    
    INTEGER, POINTER                               :: itimes

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
!                  print*, sh_env%c_tp(:,:)
                     sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF
  END SUBROUTINE wf_tresh_collaps

  SUBROUTINE  wf_damping_adiab_prop(C_param, E_0_param, sh_env, mol, deltatMD, K_tot, sh_print_env)                    
    REAL(KIND=dp), INTENT(IN)                  :: C_param, E_0_param
    TYPE(sh_env_type)                          :: sh_env
    INTEGER, INTENT(IN)                        :: mol
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    REAL(kind=dp)                              :: K_tot
    TYPE(sh_print_env_type)                    :: sh_print_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_damping', &         
    routineP = moduleN//':'//routineN   
                                                        
    !INTEGER  :: j,i
    INTEGER  :: i
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tau
    REAL(KIND=dp)      :: sum!, sum1 
    !COMPLEX(KIND=dp), dimension(:), allocatable :: ct
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop

    allocate(pop(mol)) !TEST

    !allocate(ct(mol))  ! adiabatic states
    allocate(tau(mol)) ! decoherence time (EDC method)
    allocate(sh_print_env%tau_deco(mol))

    sum = 0.0D0
!TEST WITHOUT KINETIC: useful if statment when E_0=0 and FROZEN HAMILTONIAN are
!used
    if (E_0_param == 0.0D0) then  
          K_tot = 1.0D0
    end if
!transfor the diabatic into the adiabatic coefficients and calculate tau
    do i=1,mol 
           if (i .ne. sh_env%active_state) then
              ! pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
              ! ct(i) = cmplx(0.0D0, 0.0D0,dp)                           
               tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))* &
                                        (C_param +(E_0_param/K_tot))
                print*, tau(i)
                !do j=1,mol
                !        ct(i) = ct(i) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*sh_env%ctin(j)
                !end do

!Damp the coefficients of the states different from the active state    
               !pop(i) = sh_env%ctin(i)*conjg(sh_env%ctin(i))
    
                !ct(i) = ct(i)
                sh_env%ctin(i) = sh_env%ctin(i)*exp(-(deltatMD/tau(i)))
    
                !pop(i) = sh_env%ctin(i)*conjg(sh_env%ctin(i))
    
                sum = sum + abs(sh_env%ctin(i))*abs(sh_env%ctin(i))
           else
                tau(i) = 0.0D0
           end if
    end do

    sh_print_env%tau_deco = tau
!transfor the active diabatic state into the active adiabatic state
    !pop(sh_env%active_state)= cmplx(0.0D0, 0.0D0, dp)
   ! ct(sh_env%active_state) = cmplx(0.0D0, 0.0D0, dp)
    !do i=1,mol                                                       
     !      ct(sh_env%active_state) = ct(sh_env%active_state) + &
      !       cmplx(sh_env%C_tp(i,sh_env%active_state), 0.0D0,dp)*sh_env%ctin(i)
    !end do  
    
   ! pop(sh_env%active_state) = sh_env%ctin(sh_env%active_state)*conjg(sh_env%ctin(sh_env%active_state))
    
!transfer population to the active state for keeping the detailed balance   
    if (sh_env%ctin(sh_env%active_state)== 0) then
           sh_env%ctin(sh_env%active_state) = sqrt(1-sum)
    else 
           sh_env%ctin(sh_env%active_state)= sh_env%ctin(sh_env%active_state)* &
              sqrt((1-sum)/(abs(sh_env%ctin(sh_env%active_state))*abs(sh_env%ctin(sh_env%active_state))))
    end if
    
    pop(sh_env%active_state) = sh_env%ctin(sh_env%active_state)*conjg(sh_env%ctin(sh_env%active_state))

!transfor back the adiabatic into the diabatic which, then, follow the
!electronic prpagation    
!   do j=1,mol
!         sh_env%ctin(j) = cmplx(0.0D0, 0.0D0,dp)                           
!         pop(j) = cmplx(0.0D0, 0.0D0,dp) !TEST
!         do i=1,mol 
!             sh_env%ctin(j) = sh_env%ctin(j) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*ct(i)
!         end do
!        pop(j)=sh_env%ctin(j)*conjg(sh_env%ctin(j))
!   end do
    
   CALL adiab_print_sh_info(sh_print_env, keyword = 'DAMPING')   
    
   ! deallocate(ct) 
    deallocate(tau)
    deallocate(sh_print_env%tau_deco)
    deallocate(pop)
  END SUBROUTINE wf_damping_adiab_prop 


  SUBROUTINE  wf_damping(C_param, E_0_param, sh_env, mol, deltatMD, K_tot, sh_print_env)                    
    REAL(KIND=dp), INTENT(IN)                  :: C_param, E_0_param
    TYPE(sh_env_type)                          :: sh_env
    INTEGER, INTENT(IN)                        :: mol
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    REAL(kind=dp)                              :: K_tot
    TYPE(sh_print_env_type)                    :: sh_print_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_damping', &         
    routineP = moduleN//':'//routineN   
                                                        
    INTEGER  :: j,i
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tau
    REAL(KIND=dp)      :: sum!, sum1 
    COMPLEX(KIND=dp), dimension(:), allocatable :: ct
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop

    allocate(pop(mol)) !TEST

    allocate(ct(mol))  ! adiabatic states
    allocate(tau(mol)) ! decoherence time (EDC method)
    allocate(sh_print_env%tau_deco(mol))

    sum = 0.0D0
!TEST WITHOUT KINETIC: useful if statment when E_0=0 and FROZEN HAMILTONIAN are
!used
    if (E_0_param == 0.0D0) then  
          K_tot = 1.0D0
    end if
!transfor the diabatic into the adiabatic coefficients and calculate tau
    do i=1,mol 
           if (i .ne. sh_env%active_state) then
               pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
               ct(i) = cmplx(0.0D0, 0.0D0,dp)                           
               tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))* &
                                        (C_param +(E_0_param/K_tot))
                print*, tau(i)
                do j=1,mol
                        ct(i) = ct(i) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*sh_env%ctin(j)
                end do

!damp the coefficients of the states different from the active state    
               pop(i) = ct(i)*conjg(ct(i))
    
                !ct(i) = ct(i)
                ct(i) = ct(i)*exp(-(deltatMD/tau(i)))
    
                pop(i) = ct(i)*conjg(ct(i))
    
                sum = sum + abs(ct(i))*abs(ct(i))
           else
                tau(i) = 0.0D0
           end if
    end do

    sh_print_env%tau_deco = tau
!transfor the active diabatic state into the active adiabatic state
    pop(sh_env%active_state)= cmplx(0.0D0, 0.0D0, dp)
    ct(sh_env%active_state) = cmplx(0.0D0, 0.0D0, dp)
    do i=1,mol                                                       
           ct(sh_env%active_state) = ct(sh_env%active_state) + &
             cmplx(sh_env%C_tp(i,sh_env%active_state), 0.0D0,dp)*sh_env%ctin(i)
    end do  
    
    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))
    
!transfer population to the active state for keeping the detailed balance   
    if (ct(sh_env%active_state)== 0) then
           ct(sh_env%active_state) = sqrt(1-sum)
    else 
           ct(sh_env%active_state)= ct(sh_env%active_state)* &
              sqrt((1-sum)/(abs(ct(sh_env%active_state))*abs(ct(sh_env%active_state))))
    end if
    
    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))

!transfor back the adiabatic into the diabatic which, then, follow the
!electronic prpagation    
    do j=1,mol
          sh_env%ctin(j) = cmplx(0.0D0, 0.0D0,dp)                           
          pop(j) = cmplx(0.0D0, 0.0D0,dp) !TEST
          do i=1,mol 
              sh_env%ctin(j) = sh_env%ctin(j) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*ct(i)
          end do
         pop(j)=sh_env%ctin(j)*conjg(sh_env%ctin(j))
    end do
    
   CALL adiab_print_sh_info(sh_print_env, keyword = 'DAMPING')   
    
    deallocate(ct) 
    deallocate(tau)
    deallocate(sh_print_env%tau_deco)
    deallocate(pop)
  END SUBROUTINE wf_damping 
 

  SUBROUTINE calc_kinetic_tot(array, local_particles, atomic_kinds, &
                                              nparticle_kind, atomic_kind_set, K_tot, ext_pres)
        REAL(kind=dp), DIMENSION(:,:), POINTER             :: array
        TYPE(distribution_1d_type), POINTER                :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER, OPTIONAL                             :: atomic_kind_set

        REAL(kind=dp), INTENT(OUT)                         :: K_tot
        LOGICAL, dimension(:), INTENT(IN), OPTIONAL           :: ext_pres

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

           IF (PRESENT(atomic_kinds)) THEN
              my_nparticle_kind = atomic_kinds%n_els
              my_atomic_kind_set => atomic_kinds%els
           ELSE
              IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
                 my_nparticle_kind = nparticle_kind
                 my_atomic_kind_set => atomic_kind_set
              END IF
           END IF
     
           K_tot = 0.0D0
           DO iparticle_kind= 1, my_nparticle_kind
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
print *, "MASS_INSIDE", mass
!calculate kinetic energy according to the COM since the array is velocities_COM
              DO iparticle_local= 1, nparticle_local
                 if ((PRESENT(ext_pres)) .and. (ext_pres(iparticle_local))) then
                        iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                        K_tot = K_tot + mass*( &
                                   array(1, iparticle)**2 + &
                                   array(2, iparticle)**2 + &
                                   array(3, iparticle)**2  ) 
                     ! end if
                 else 
                     iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                     K_tot = K_tot + mass*( &
                                array(1, iparticle)**2 + &
                                array(2, iparticle)**2 + &
                                array(3, iparticle)**2  )
                 end if 
              END DO
           END DO
           K_tot = 0.5*K_tot

print *, "VELOCITIES_COM", array
print *, 'K_TOT_INSIDE =', K_tot
   END SUBROUTINE calc_kinetic_tot
 
  !SUBROUTINE propagate_wf(H_t, H_tp,  Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm,     &
  !                       ctin, aom_env, nadiab,  norbitals, ELsteps, &
  !                         itimes, adiab_section, deltatMD, deltatE , times,      &
  !                          mat)
  SUBROUTINE propagate_wf(ctin, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat)

     !COMPLEX(KIND=dp), dimension(:,:), INTENT(IN)    :: H_t, H_tp
     !REAL(KIND=dp), dimension(:,:), INTENT(IN)       :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm
     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 


     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf', &
     routineP = moduleN//':'//routineN

     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctout, ctotal
     REAL(KIND=dp), dimension(:,:), allocatable      :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_0, H_2, X1, X2, X12
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: Ktot, K2, K3, K4
     REAL(KIND=dp), dimension(:,:), allocatable      :: Imat 
     REAL(KIND=dp), dimension(:,:), allocatable       :: Stt_int, Sttp_int, Stptp_int, &
                                                          Sttm_int, Stmtm_int


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
     allocate(Stt_int(mat,mat)) 
     !allocate(Stpt(mat,mat))
     allocate(Sttp_int(mat,mat))
     allocate(Stptp_int(mat,mat))
     !allocate(Stmt(mat,mat))
     allocate(Sttm_int(mat,mat))
     allocate(Stmtm_int(mat,mat))

     allocate(ctout(mat))
     allocate(ctotal(mat))

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


             Hdiff = (deltatE/deltatMD)*(aom_env%H_tp - aom_env%H_t)
             Sttdiff = (deltatE/deltatMD)*(aom_env%Stptp - aom_env%Stt)
             Sttpdiff = (deltatE/deltatMD)*(aom_env%Sttp - aom_env%Stmt)
             Sttmdiff = (deltatE/deltatMD)*(aom_env%Stpt - aom_env%Sttm)
             Stptpdiff = Sttdiff
             Stmtmdiff = (deltatE/deltatMD)*(aom_env%Stt - aom_env%Stmtm)
             H_0  = aom_env%H_t
             H_2 = aom_env%H_tp

!Assignment internal variabls used in the RK
             
             Stt_int = aom_env%Stt 
             Sttp_int = aom_env%Sttp
             Stptp_int = aom_env%Stptp
             !Stpt_int = Stpt!
             !Stmt_int = Stmt!
             Sttm_int = aom_env%Sttm
             Stmtm_int = aom_env%Stmtm

             do t2=1,int(deltatMD/deltatE)
                    call makeX(mat, aom_env%H_t, Stt_int, Sttp_int, Sttm_int, Stptp_int, &
                                             Stmtm_int, X1, deltatMD)



                    !create timestepped equivalents
                    H_1 = H_0 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff
                    call makeX(mat, H_1, Stt_int, Sttp_int, Sttm_int, Stptp_int, & 
                                                       Stmtm_int, X12, deltatMD)
     
                    H_2 = H_1 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff
                    call makeX(mat, H_2, Stt_int, Sttp_int, Sttm_int, Stptp_int, Stmtm_int, &
                                                 X2, deltatMD, NACE)
     
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

                    times = DBLE(itimes) + (t2-1)*1/DBLE(ELsteps)
                    CALL adiab_print_nace(adiab_section, NACE, nadiab, norbitals, itimes, deltatMD, times)
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
       deallocate(Stt_int, Sttp_int, Stptp_int, Sttm_int, Stmtm_int)
       deallocate(Imat)
!END DEALLOCATE

  END SUBROUTINE propagate_wf    


  SUBROUTINE prop_wf_HS_corr(ctin, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat, class_forbidden_hop)

     !COMPLEX(KIND=dp), dimension(:,:), INTENT(IN)    :: H_t, H_tp
     !REAL(KIND=dp), dimension(:,:), INTENT(IN)       :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm
     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
     INTEGER, dimension(:), INTENT(IN)               :: class_forbidden_hop

     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf', &
     routineP = moduleN//':'//routineN

     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctout, ctotal
     REAL(KIND=dp), dimension(:,:), allocatable      :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE
!     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE1
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_0, H_2, X1, X2, X12
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: Ktot, K2, K3, K4
     REAL(KIND=dp), dimension(:,:), allocatable      :: Imat 
     REAL(KIND=dp), dimension(:,:), allocatable       :: Stt_int, Sttp_int, Stptp_int, &
                                                          Sttm_int, Stmtm_int
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE       :: U_t, U_t1, U_t2
     COMPLEX(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: interm_hamil
     REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: interm_umatrix
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE         :: Eig_val_t, Eig_val_t1, Eig_val_t2


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
     allocate(interm_hamil(2,5,mat,mat))
     allocate(interm_umatrix(2,5,mat,mat))
     allocate(Stptpdiff(mat,mat))
     allocate(Sttpdiff(mat,mat))
     allocate(Sttdiff(mat,mat))
     allocate(Sttmdiff(mat,mat))
     allocate(Stmtmdiff(mat,mat))
     allocate(Stt_int(mat,mat)) 
     !allocate(Stpt(mat,mat))
     allocate(Sttp_int(mat,mat))
     allocate(Stptp_int(mat,mat))
     !allocate(Stmt(mat,mat))
     allocate(Sttm_int(mat,mat))
     allocate(Stmtm_int(mat,mat))

     allocate(ctout(mat))
     allocate(ctotal(mat))

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


             Hdiff = (deltatE/deltatMD)*(aom_env%H_tp - aom_env%H_t)
             Sttdiff = (deltatE/deltatMD)*(aom_env%Stptp - aom_env%Stt)
             Sttpdiff = (deltatE/deltatMD)*(aom_env%Sttp - aom_env%Stmt)
             Sttmdiff = (deltatE/deltatMD)*(aom_env%Stpt - aom_env%Sttm)
             Stptpdiff = Sttdiff
             Stmtmdiff = (deltatE/deltatMD)*(aom_env%Stt - aom_env%Stmtm)
             H_0  = aom_env%H_t
             H_2 = aom_env%H_tp

print *, "START_aom_env%H_t", aom_env%H_t
print *, "START_aom_env%H_tp", aom_env%H_tp
print *, "START_Hdiff", Hdiff

            !> LOOP TO BUILD INTERMEDIATE HAMILTONIANS AND CALCULATE 
            !> RESPECTIVE U MATRICES
           
            CALL adiabaticity(-AIMAG(aom_env%H_t), mat, U_t, Eig_val_t, .false.)
            do t2=1,int(deltatMD/deltatE)
                
                   !create 2 others intermediate hamiltonians usefull to build the X
                   !matrix in the following 
                   H_1 = H_0 + 0.5D0*Hdiff
                   interm_hamil(1,t2,:,:) = H_1(:,:)
                   CALL adiabaticity(-AIMAG(H_1), mat, U_t1, Eig_val_t1, .false.)
                   interm_umatrix(1,t2,:,:) = U_t1(:,:)
                   H_2 = H_1 + 0.5D0*Hdiff
                   interm_hamil(2,t2,:,:) = H_2(:,:)
                   CALL adiabaticity(-AIMAG(H_2), mat, U_t2, Eig_val_t2, .false.)
                   interm_umatrix(2,t2,:,:) = U_t2(:,:)

                   H_0 = H_2
            end do
print *, "interm_hamil", interm_hamil 
print *, "interm_umatrix", interm_umatrix 
                
             !H_0  = aom_env%H_t
             !H_2 = aom_env%H_tp

!Assignment internal variabls used in the RK
             
             Stt_int = aom_env%Stt 
             Sttp_int = aom_env%Sttp
             Stptp_int = aom_env%Stptp
             !Stpt_int = Stpt!
             !Stmt_int = Stmt!
             Sttm_int = aom_env%Sttm
             Stmtm_int = aom_env%Stmtm

             do t2=1,int(deltatMD/deltatE)
  
                    !CALL adiabaticity(-AIMAG(aom_env%H_t), mat, U_t, Eig_val_t, .false.)
                    CALL calc_diabatic_nace(mat, aom_env%H_t, Stt_int, Sttp_int, Sttm_int, &
                                                        Stptp_int, Stmtm_int, X1, deltatMD, NACE)
!print *, "LOOP_H_t", aom_env%H_t
!print *, "LOOP_H_0", H_0
print *, "LOOP_X_1", X1

!print *, "NACE_CALC_DIAB", NACE1
!                    call makeX(mat, aom_env%H_t, Stt_int, Sttp_int, Sttm_int, Stptp_int, &
!                                             Stmtm_int, X1, deltatMD, NACE)
!print *, "UMATRIX", U_t
!print *, "NACE_MAKEX_OLD", NACE
!print *, "class_forbidden_hop", class_forbidden_hop
     
                    !create timestepped equivalents
                    !H_1 = H_0 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff

!                    CALL adiabaticity(-AIMAG(H_1), mat, U_t1, Eig_val_t1, .false.)
                    call calc_diabatic_nace(mat, interm_hamil(1,t2,:,:), Stt_int, Sttp_int, Sttm_int, Stptp_int, & 
                                                       Stmtm_int, X12, deltatMD, NACE )

!print *, "LOOP_H_1", H_1
!print *, "LOOP_X_12", X12
!print *, "UMATRIX_1", U_t1

                   ! H_2 = H_1 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff

!                    CALL adiabaticity(-AIMAG(H_2), mat, U_t2, Eig_val_t2, .false.)
                    call calc_diabatic_nace(mat, interm_hamil(2,t2,:,:), Stt_int, Sttp_int, Sttm_int, Stptp_int, &
                                                     Stmtm_int, X2, deltatMD, NACE)
!print *, "LOOP_H_2", H_2     
!print *, "LOOP_X_2", X2
!print *, "UMATRIX_2", U_t2

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
                    !H_0 = H_2

                    ctin = ctout

                    ctout = cmplx(0.0D0, 0.0D0,dp)

                    times = DBLE(itimes) + (t2-1)*1/DBLE(ELsteps)
                    CALL adiab_print_nace(adiab_section, NACE, nadiab, norbitals, itimes, deltatMD, times)
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
       deallocate(Stt_int, Sttp_int, Stptp_int, Sttm_int, Stmtm_int)
       deallocate(Imat)
!END DEALLOCATE

  END SUBROUTINE prop_wf_HS_corr                                                                

                                                            

  SUBROUTINE propagate_adiab_wf(E_t, E_tp, ctin, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat)
     
     REAL(KIND=dp), DIMENSION(:), INTENT(IN)         :: E_t, E_tp 
     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
     TYPE(aom_env_type), intent(in)                  :: aom_env
 
 
     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf', &
     routineP = moduleN//':'//routineN
 
     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctout, ctotal
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACEdiff
     REAL(KIND=dp), dimension(:,:), allocatable     ::  NACEint
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: X1, X2, X12
     REAL(KIND=dp), DIMENSION(:), allocatable        :: E_1, Ediff, E_0, E_2
 
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: Ktot, K2, K3, K4
     REAL(KIND=dp), dimension(:,:), allocatable      :: Imat 
 
 
     INTEGER                                         :: ELsteps   
     INTEGER, POINTER                                :: itimes
     TYPE(section_vals_type), POINTER                :: adiab_section
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
     
     allocate(NACEdiff(mat,mat))
     allocate(NACEint(mat,mat)) 
 
     allocate(ctout(mat))
     allocate(ctotal(mat))
 
     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H_t etc are always -i*H.
     !This is why they're complex variables
     !allocate(NACE(mat,mat))
     allocate(Ediff(mat))
     allocate(E_0(mat))
     allocate(E_1(mat))
     allocate(E_2(mat))
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
 
 
             Ediff = (deltatE/deltatMD)*(E_tp - E_t)
             NACEdiff = (deltatE/deltatMD)*(aom_env%adiab_nace_tp - &
                                             aom_env%adiab_nace_t)
             
 
             E_0  = E_t
             E_2 = E_tp

             NACEint = aom_env%adiab_nace_t
 
             do t2=1,int(deltatMD/deltatE)
                    call make_adiabX(mat, E_t, NACEint, X1)
 
 
                    !create timestepped equivalents
                    E_1 = E_0 + 0.5D0*Ediff
                    NACEint = NACEint + 0.5D0*NACEdiff
                    
                    call make_adiabX(mat, E_1, NACEint, X12)

     
                    E_2 = E_1 + 0.5D0*Ediff
                    NACEint = NACEint + 0.5D0*NACEdiff

                    call make_adiabX(mat, E_2, NACEint, X2)
     
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
print *, "CT_BEFORE_RENORM", ctmag
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
                    E_0 = E_2
 
                    ctin = ctout
 
                    ctout = cmplx(0.0D0, 0.0D0,dp)
 
                    times = DBLE(itimes) + (t2-1)*1/DBLE(ELsteps)
                    CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, &
                                           times)
                    CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,itimes, deltatMD, &
                                         times)
                    CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes,deltatMD, &
                                                    times)
 
             end do
 
 !END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
 
 !DEALLOCATE
       !DEALLOCATE(NACE)
       deallocate(E_0, Ediff, E_1, E_2, X1, X2, X12, Ktot, K2, K3, K4)
       deallocate(Imat)
 !END DEALLOCATE
 
  END SUBROUTINE propagate_adiab_wf                                                                



END MODULE non_adiabatic_utils 
                       
