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
     USE parallel_rng_types,              ONLY: next_random_number    
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
                 classically_forbidden_hops, &
                 ida_collapse
                  

CONTAINS


SUBROUTINE classically_forbidden_hops(active_state, aom_env, mat, E_tp, &
                  velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                          atomic_kind_set, energy_conservation, class_forbidden_hop)

     INTEGER, INTENT(IN)                                :: active_state, mat
     TYPE(aom_env_type), INTENT(IN)                     :: aom_env
     REAL(KIND=dp), DIMENSION(:), INTENT(IN)         ::  E_tp
     REAL(kind=dp), DIMENSION(:,:), POINTER       :: velocities_com_tp
     TYPE(distribution_1d_type), POINTER                :: local_particles
     TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
     INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
     TYPE(atomic_kind_type), DIMENSION(:), &
                   POINTER, OPTIONAL                        :: atomic_kind_set
     REAL(KIND=dp), dimension(:), allocatable, INTENT(OUT)    :: energy_conservation
     INTEGER, dimension(:), allocatable, INTENT(OUT)          :: class_forbidden_hop


     INTEGER :: j
     REAL(kind=dp), DIMENSION(:,:), allocatable     :: hop_nacv
     REAL(kind=dp)                               :: delta_energy
     REAL(kind=dp)                               :: A, B
     REAL(kind=dp)                               :: K_tot, E_cutoff  
     
     ALLOCATE(energy_conservation(mat)) 
     ALLOCATE(class_forbidden_hop(mat))

     CALL calc_kinetic_tot(velocities_com_tp, local_particles, atomic_kinds, &
                  nparticle_kind, atomic_kind_set, K_tot)
     E_cutoff = 0.0D0
     E_cutoff = E_tp(active_state) + K_tot  !THIS CUTOFF IS TO LARGE ALL THE
                                            !STATE ARE BELOW IT
     do j=1, mat
           if ((j .gt. active_state) .and. (E_tp(j) .le. E_cutoff)) then
               
              CALL calc_hop_nacv(aom_env, hop_nacv, active_state, j, fast=.TRUE.)
              CALL evaluate_delta_energy_nacv(velocities_com_tp, local_particles, atomic_kinds, &
                                     nparticle_kind, atomic_kind_set, &
                                    aom_env, E_tp, hop_nacv, delta_energy, A, B, & 
                                              active_state, j)
              energy_conservation(j) = delta_energy
              if (energy_conservation(j) .lt. 0.0D0) then 
                   class_forbidden_hop(j) = j
              else 
                   class_forbidden_hop(j) = -1
              end if
           else if ((j .gt. active_state) .and. (E_tp(j) .gt. E_cutoff)) then
              class_forbidden_hop(j) = j
              energy_conservation(j) = -2  
           else if (j .le. active_state) then
              energy_conservation(j) = 0.0D0
              class_forbidden_hop(j) = -1
           end if  
     end do
END SUBROUTINE classically_forbidden_hops


!COMMENT: IS THIS SUBROUTINE REDUNDANT WITH cal_delta_energy IN aom_utils.F?
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
        REAL(kind=dp), INTENT(OUT)                         :: delta_energy
        REAL(kind=dp), INTENT(OUT)                         :: A, B
        REAL(KIND=dp), DIMENSION(:), INTENT(IN)            ::  E_tp


        CHARACTER(len=*), PARAMETER :: routineN = 'evaluate_delta_energy_nacv', &
        routineP = moduleN//':'//routineN

        TYPE(atomic_kind_type), ALLOCATABLE, TARGET, DIMENSION(:)     :: my_atomic_kind_set
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
           ALLOCATE(my_atomic_kind_set(SIZE(atomic_kinds%els)))
           my_atomic_kind_set = atomic_kinds%els
        ELSE
           IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
              my_nparticle_kind = nparticle_kind
              ALLOCATE(my_atomic_kind_set(SIZE(atomic_kind_set)))
              my_atomic_kind_set = atomic_kind_set
           END IF
        END IF

           A = 0
           B = 0
           DO iparticle_kind=1, my_nparticle_kind
              atomic_kind => my_atomic_kind_set(iparticle_kind)
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
           delta_energy = B**2 - 4 * A * (  &
                     E_tp(new_state) - E_tp(old_state))
   END SUBROUTINE evaluate_delta_energy_nacv




 
  SUBROUTINE decohere_wf(adiab_section, sh_print_env, aom_env,sh_env, deco_keyword, itimes, &
                         deltatMD, mol,velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                         atomic_kind_set, K_tot, ext_pres)

        TYPE(section_vals_type), POINTER               :: adiab_section 
        TYPE(sh_print_env_type)                        :: sh_print_env
        TYPE(aom_env_type)                             :: aom_env
        TYPE(sh_env_type)                              :: sh_env
        CHARACTER(LEN=default_string_length)           :: deco_keyword
        INTEGER, POINTER                               :: itimes
        REAL(KIND=dp), INTENT(IN)                      :: deltatMD
        INTEGER                                        :: mol
        REAL(kind=dp), DIMENSION(:,:), POINTER             :: velocities_com_tp
        TYPE(distribution_1d_type), POINTER                :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
                  POINTER, OPTIONAL                        :: atomic_kind_set

        REAL(kind=dp), INTENT(OUT)                         :: K_tot
        LOGICAL, dimension(:), POINTER, INTENT(IN)         :: ext_pres

        CHARACTER(len=*), PARAMETER :: routineN ='decohere_wf', &         
                                       routineP = moduleN//':'//routineN

!       logical                                            :: edc_deco_time
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: tau!, tau_forces
        CHARACTER(LEN=default_string_length)           :: wfprop_keyword
!        REAL(KIND=dp)                                  :: C_param, E_0_param

        CALL section_vals_val_get(adiab_section,"REPRESENTATION", &
                                                     c_val=wfprop_keyword)
        CALL section_vals_val_get(adiab_section,"DECOHERENCE_CORRECTIONS",c_val=deco_keyword)
        
        aom_env%ida_deco = .FALSE.
        
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
              CALL calc_kinetic_tot(velocities_com_tp, local_particles, atomic_kinds, &
                                    nparticle_kind, atomic_kind_set, K_tot, ext_pres)
              CALL eval_deco_tau(tau, adiab_section, sh_env, mol, K_tot, sh_print_env, &
                          local_particles, atomic_kinds,nparticle_kind, atomic_kind_set) 
                   
              CALL  wf_damping(tau, sh_env, mol, deltatMD, sh_print_env)  
            CASE("COLLAPSE_WITHIN_DECOTIME")
              CALL calc_kinetic_tot(velocities_com_tp, local_particles, atomic_kinds, &
                                    nparticle_kind, atomic_kind_set, K_tot, ext_pres)
              CALL eval_deco_tau(tau, adiab_section, sh_env, mol, K_tot, sh_print_env, &
                          local_particles, atomic_kinds,nparticle_kind, atomic_kind_set) 
                   
              CALL  collapse_within_deco_time(tau, mol, deltatMD, aom_env, sh_env)
            CASE("IDA")
               aom_env%ida_deco = .TRUE.
print *, "aom_env%ida_deco", aom_env%ida_deco
            CASE DEFAULT                                           
              CPABORT("There is no method to add the decoherence.") 
            END SELECT
        CASE("ADIABATIC_BASIS")
            SELECT CASE(deco_keyword)                                                          
            CASE("NO_DECO_CORR")
                 PRINT *, 'NO DECOHERENCE CORRECTIONS'
            CASE("INSTANT_COLLAPSE")
              CALL wf_collapse_adiab_prop(aom_env,sh_print_env, sh_env, itimes)
            CASE("TRESH_ONLY_COLLAPSE")  
              CPABORT("NOT IMPLEMENTED")
            CASE("DAMPING")
              CALL calc_kinetic_tot(velocities_com_tp, local_particles, atomic_kinds, &
                                    nparticle_kind, atomic_kind_set, K_tot, ext_pres)
              CALL eval_deco_tau(tau, adiab_section, sh_env, mol, K_tot, sh_print_env, &
                          local_particles, atomic_kinds,nparticle_kind, atomic_kind_set) 
              CALL wf_damping_adiab_prop(tau, sh_env, mol, deltatMD, sh_print_env)                    
            CASE DEFAULT                                           
              CPABORT("There is no method to add the decoherence.") 
            END SELECT
        CASE DEFAULT
           CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
        END SELECT
  END SUBROUTINE decohere_wf


  SUBROUTINE wf_collapse_adiab_prop(aom_env,sh_print_env, sh_env, itimes)
    TYPE(aom_env_type)                             :: aom_env
    TYPE(sh_print_env_type)                        :: sh_print_env
    TYPE(sh_env_type)                              :: sh_env
    INTEGER, POINTER                               :: itimes

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse_adiab_prop', &         
    routineP = moduleN//':'//routineN                                                                 

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (aom_env%futur_decohere .and. sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
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
    INTEGER, POINTER                               :: itimes

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
    

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (aom_env%futur_decohere .and. sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
                     sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
                  aom_env%futur_decohere = .FALSE.                         
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF
  END SUBROUTINE wf_collapse 


  SUBROUTINE ida_collapse(aom_env)
    TYPE(aom_env_type)                             :: aom_env
 
    CHARACTER(len=*), PARAMETER :: routineN ='ida_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
    
         if (aom_env%successful_hop) then
              aom_env%ctin(:) = aom_env%C_t(:, aom_env%new_state)
         else 
              aom_env%ctin(:) = aom_env%C_t(:, aom_env%active_state)
         end if                                                            
  END SUBROUTINE ida_collapse

 
  SUBROUTINE wf_tresh_collaps(sh_print_env, sh_env, itimes)
    TYPE(sh_print_env_type)                        :: sh_print_env
    TYPE(sh_env_type)                              :: sh_env
    INTEGER, POINTER                               :: itimes

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
                     sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF
  END SUBROUTINE wf_tresh_collaps


  SUBROUTINE  wf_damping_adiab_prop(tau, sh_env, mol, deltatMD, sh_print_env)                    
    REAL(KIND=dp), DIMENSION(:), INTENT(in)    :: tau!, tau_forces
    TYPE(sh_env_type)                          :: sh_env
    INTEGER, INTENT(IN)                        :: mol
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    TYPE(sh_print_env_type)                    :: sh_print_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_damping_adiab_prop', &         
    routineP = moduleN//':'//routineN   
                                                        
    INTEGER  :: i
    REAL(KIND=dp)      :: sum 
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop

    allocate(pop(mol)) !TEST

    sum = 0.0D0


!transfor the diabatic into the adiabatic coefficients and calculate tau
    do i=1,mol 
           if (i .ne. sh_env%active_state) then
!               tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))* &
!                                        (C_param +(E_0_param/K_tot))

!Damp the coefficients of the states different from the active state    
                sh_env%ctin(i) = sh_env%ctin(i)*exp(-(deltatMD/tau(i)))
                sum = sum + abs(sh_env%ctin(i))*abs(sh_env%ctin(i))
           end if
    end do

!transfer population to the active state for keeping the detailed balance   
    if (sh_env%ctin(sh_env%active_state)== 0) then
           sh_env%ctin(sh_env%active_state) = sqrt(1-sum)
    else 
           sh_env%ctin(sh_env%active_state)= sh_env%ctin(sh_env%active_state)* &
              sqrt((1-sum)/(abs(sh_env%ctin(sh_env%active_state))*abs(sh_env%ctin(sh_env%active_state))))
    end if
    
    pop(sh_env%active_state) = sh_env%ctin(sh_env%active_state)*conjg(sh_env%ctin(sh_env%active_state))
    deallocate(pop)
  END SUBROUTINE wf_damping_adiab_prop 

  
  SUBROUTINE  eval_deco_tau(tau, adiab_section, sh_env, mol, K_tot, sh_print_env, &
                          local_particles, atomic_kinds,nparticle_kind, atomic_kind_set)                    
!    LOGICAL, INTENT(IN)                        :: edc_deco_time
    TYPE(sh_env_type)                          :: sh_env
    INTEGER, INTENT(IN)                        :: mol
!    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    REAL(kind=dp)                              :: K_tot
    TYPE(sh_print_env_type)                    :: sh_print_env
    TYPE(distribution_1d_type), POINTER                :: local_particles
    TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
                POINTER, OPTIONAL                             :: atomic_kind_set
    INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: tau!, tau_forces
    TYPE(section_vals_type), POINTER               :: adiab_section



    CHARACTER(len=*), PARAMETER :: routineN ='eval_deco_tau', &         
    routineP = moduleN//':'//routineN   

    CHARACTER(LEN=default_string_length)           :: tau_param
    REAL(KIND=dp)                           :: C_param, E_0_param
    TYPE(atomic_kind_type), ALLOCATABLE, TARGET, DIMENSION(:)     :: my_atomic_kind_set
    INTEGER                                  :: my_nparticle_kind
    REAL(kind=dp)                            :: mass
    INTEGER                                  :: iparticle_kind
    INTEGER                                  :: nparticle_local
    REAL(kind=dp)                            :: mixing_threshold
    INTEGER                                  :: iparticle_local, iparticle
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    INTEGER  :: i, deg_of_freedom
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: adiab_nace_vector
    REAL(KIND=dp)      :: KbT, Temp, alpha, T_inp 
    REAL(kind=dp), DIMENSION(:,:),POINTER          :: my_forces_adiab
    LOGICAL                                  :: tau_schwartz


    ALLOCATE(my_forces_adiab(3,sh_env%natom))
    allocate(tau(mol)) ! decoherence time (EDC method)
    allocate(sh_print_env%tau_deco(mol))
    allocate(adiab_nace_vector(mol))

    adiab_nace_vector(:) = sh_env%adiab_nace(sh_env%active_state,:)
    !>number of degrees of freedom intensive EDC
    deg_of_freedom = (sh_env%active_atoms*3-6)


    IF (PRESENT(atomic_kinds)) THEN
       my_nparticle_kind = atomic_kinds%n_els
       ALLOCATE(my_atomic_kind_set(SIZE(atomic_kinds%els)))
       my_atomic_kind_set = atomic_kinds%els
    ELSE
       IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
          my_nparticle_kind = nparticle_kind
          ALLOCATE(my_atomic_kind_set(SIZE(atomic_kind_set)))
          my_atomic_kind_set = atomic_kind_set
       END IF
    END IF

!TEST TAU_SCHARTZ
tau_schwartz = .TRUE.
!END TEST
    
!   tau_schwartz = .FALSE.


!transfor the diabatic into the adiabatic coefficients and calculate tau
   
   CALL section_vals_val_get(adiab_section,"DECO_TIME",c_val=tau_param)
    do i=1,mol 
           if (i .ne. sh_env%active_state) then
              SELECT CASE(tau_param)
              CASE("EDC")
                    CALL section_vals_val_get(adiab_section,"EDC_C",r_val=C_param)
                    CALL section_vals_val_get(adiab_section,"EDC_E0",r_val=E_0_param)
                    !TEST WITHOUT KINETIC: useful if statment when E_0=0 and FROZEN HAMILTONIAN are
                    !used
                    if (E_0_param == 0.0D0) then  
                         K_tot = 1.0D0
                    end if
                    tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))* &
                                        (C_param+((E_0_param/K_tot)*deg_of_freedom))

              CASE("FORCES_BASED")

                    !Boltzmann factor to use in alpha parameter
                     KbT = (2*K_tot)/(3*sh_env%active_atoms)
                    !instead of using the exact boltzmann factor, at the moment the code uses Temp =
                    !298 K, maybe this point should be modified in future
                    
                    CALL section_vals_val_get(adiab_section,"THRESHOLD_TAU_FORCES",r_val=mixing_threshold)
                    CALL section_vals_val_get(adiab_section,"TEMPERATURE_FG_WIDTH",r_val=T_inp)
                    
                    Temp = 3.1693191275167785e-06*T_inp !conver in KbT
                      !Temp = 0.0009444571  !KbT at 298 K
                      ! mixing_threshold = 2.0D0
                    if ((abs(adiab_nace_vector(i))/(abs(sh_env%E_tp(i)-sh_env%E_tp(sh_env%active_state)))) & 
                                         .ge. mixing_threshold) then
                         !if there is mixing then we use force otherwise we use !pure dephasing
                         call calcforces(i ,-AIMAG(sh_env%H_tp), sh_env%C_tp, sh_env%nacv, &
                              sh_env%NN, sh_env%atoms, &
                              sh_env%scaling, sh_env%cbar, sh_env%F, my_forces_adiab, &
                             sh_env%ext_pres, sh_env%analytics, sh_env%exact_forces)
                         tau(i) = 0.0D0
                         DO iparticle_kind= 1, my_nparticle_kind
                                atomic_kind => my_atomic_kind_set(iparticle_kind)
                                CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
                                nparticle_local = local_particles%n_el(iparticle_kind)
                                DO iparticle_local= 1, nparticle_local
                                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                                   alpha = 6*mass*Temp                  
                                   if (tau_schwartz) then    
                                        tau(i) = tau(i) + &
                                        (dot_product((my_forces_adiab(:,iparticle)-sh_env%forces_active_state(:,iparticle)),&
                                        (my_forces_adiab(:,iparticle)-sh_env%forces_active_state(:,iparticle))))/(4*alpha)
                                   else
                                        tau(i) = tau(i) + &
                                        (sqrt(dot_product((my_forces_adiab(:,iparticle)-sh_env%forces_active_state(:,iparticle)),&
                                        (my_forces_adiab(:,iparticle)-sh_env%forces_active_state(:,iparticle)))))/(2*sqrt(alpha))
                                   end if
                                END DO
                         END DO
                         if (tau_schwartz) then
                             tau(i) = sqrt(1/tau(i))
                         else
                             tau(i) = 1/tau(i)
                         end if
                    else
                         tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))
                    end if
              CASE DEFAULT
                  CPABORT("There is no method to choose decoherence time.")
              END SELECT
           else
                tau(i) = 0.0D0
           end if
    end do

    sh_print_env%tau_deco = tau

print *, "TAU_NEW_EVAL", tau
    
    
    CALL adiab_print_sh_info(sh_print_env, keyword = 'DAMPING')   

!    deallocate(tau, adiab_nace_vector)
    deallocate(adiab_nace_vector)
    deallocate(sh_print_env%tau_deco)
    DEALLOCATE(my_forces_adiab)
  END SUBROUTINE eval_deco_tau


  SUBROUTINE collapse_within_deco_time(tau, mol, deltatMD, aom_env, sh_env)

    INTEGER, INTENT(IN)                        :: mol
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    TYPE(aom_env_type)                         :: aom_env
    TYPE(sh_env_type)                          :: sh_env
    REAL(KIND=dp), DIMENSION(:), INTENT(in)    :: tau!, tau_forces

    CHARACTER(len=*), PARAMETER :: routineN ='collapse_within_deco_time', &
    routineP = moduleN//':'//routineN

    INTEGER  :: i, j
    REAL(KIND=dp) :: random_num, collapse_probability
    COMPLEX(KIND=dp), dimension(:), allocatable :: ct
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop
    REAL(KIND=dp)      :: sum  


    allocate(pop(mol)) !TEST
    allocate(ct(mol))  ! adiabatic states

    sum = 0.0D0
    do i=1,mol
         if (i .ne. sh_env%active_state) then
             pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
             ct(i) = cmplx(0.0D0, 0.0D0,dp) 
             !> Tranform coefficients in the adiabatic basis 
             do j=1,mol
                    ct(i) = ct(i) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*sh_env%ctin(j)
             end do
             print *, "ct_SUBOTNIK_BEFRE", ct(i)                
             pop(i) = ct(i)*conjg(ct(i))
             print *, "POPU_BEFORE", pop(i)
        
             !>evaluate which states will be collpased
             collapse_probability = deltatMD/tau(i)
             random_num = next_random_number(aom_env%rng_stream)
             if (random_num .gt. collapse_probability) then 
                    ct(i) = ct(i)  !nothing happens
             else
                    ct(i) = 0.0D0
             end if
             sum = sum + abs(ct(i))*abs(ct(i))
         end if
    end do
 
    !>transfor the active diabatic state into the active adiabatic state
    pop(sh_env%active_state)= cmplx(0.0D0, 0.0D0, dp)
    ct(sh_env%active_state) = cmplx(0.0D0, 0.0D0, dp)
    do i=1,mol                                                       
           ct(sh_env%active_state) = ct(sh_env%active_state) + &
             cmplx(sh_env%C_tp(i,sh_env%active_state), 0.0D0,dp)*sh_env%ctin(i)
    end do  
    
    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))
    
    !>transfer population to the active state for keeping the detailed balance   
    if (ct(sh_env%active_state)== 0) then
           ct(sh_env%active_state) = sqrt(1-sum)
    else 
           ct(sh_env%active_state)= ct(sh_env%active_state)* &
              sqrt((1-sum)/(abs(ct(sh_env%active_state))*abs(ct(sh_env%active_state))))
    end if

    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))

    !>transfor back the adiabatic into the diabatic which, then, follow the
      !electronic prpagation    
    do j=1,mol
          sh_env%ctin(j) = cmplx(0.0D0, 0.0D0,dp)                           
          pop(j) = cmplx(0.0D0, 0.0D0,dp) !TEST
          do i=1,mol 
              sh_env%ctin(j) = sh_env%ctin(j) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*ct(i)
          end do
         pop(j)=sh_env%ctin(j)*conjg(sh_env%ctin(j))
    end do

    deallocate(pop) !TEST
    deallocate(ct)  ! adiabatic states

  END SUBROUTINE collapse_within_deco_time


  SUBROUTINE  wf_damping(tau, sh_env, mol, deltatMD, sh_print_env)                    
    REAL(KIND=dp), DIMENSION(:), INTENT(in)    :: tau!, tau_forces
    TYPE(sh_env_type)                          :: sh_env
    INTEGER, INTENT(IN)                        :: mol
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    TYPE(sh_print_env_type)                    :: sh_print_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_damping', &         
    routineP = moduleN//':'//routineN   

    INTEGER  :: j,i
    REAL(KIND=dp)      :: sum  
    COMPLEX(KIND=dp), dimension(:), allocatable :: ct
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop

    allocate(pop(mol)) !TEST
    allocate(ct(mol))  ! adiabatic states

    sum = 0.0D0
    do i=1,mol 
           if (i .ne. sh_env%active_state) then
               pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
               ct(i) = cmplx(0.0D0, 0.0D0,dp) 
               !> if statement to choose the decoherence time

               !> Tranform coefficients in the adiabatic basis 
                do j=1,mol
                        ct(i) = ct(i) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*sh_env%ctin(j)
                end do
                !>damp the coefficients of the states different from the active state    
                pop(i) = ct(i)*conjg(ct(i))
                !ct(i) = ct(i)
                ct(i) = ct(i)*exp(-(deltatMD/tau(i)))
                pop(i) = ct(i)*conjg(ct(i))
                sum = sum + abs(ct(i))*abs(ct(i))
           end if
    end do

    !>transfor the active diabatic state into the active adiabatic state
    pop(sh_env%active_state)= cmplx(0.0D0, 0.0D0, dp)
    ct(sh_env%active_state) = cmplx(0.0D0, 0.0D0, dp)
    do i=1,mol                                                       
           ct(sh_env%active_state) = ct(sh_env%active_state) + &
             cmplx(sh_env%C_tp(i,sh_env%active_state), 0.0D0,dp)*sh_env%ctin(i)
    end do  
    
    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))
    
    !>transfer population to the active state for keeping the detailed balance   
    if (ct(sh_env%active_state)== 0) then
           ct(sh_env%active_state) = sqrt(1-sum)
    else 
           ct(sh_env%active_state)= ct(sh_env%active_state)* &
              sqrt((1-sum)/(abs(ct(sh_env%active_state))*abs(ct(sh_env%active_state))))
    end if
    
    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))

    !>transfor back the adiabatic into the diabatic which, then, follow the
      !electronic prpagation    
    do j=1,mol
          sh_env%ctin(j) = cmplx(0.0D0, 0.0D0,dp)                           
          pop(j) = cmplx(0.0D0, 0.0D0,dp) !TEST
          do i=1,mol 
              sh_env%ctin(j) = sh_env%ctin(j) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*ct(i)
          end do
         pop(j)=sh_env%ctin(j)*conjg(sh_env%ctin(j))
    end do
    
!   CALL adiab_print_sh_info(sh_print_env, keyword = 'DAMPING')   
    
    deallocate(ct) 
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

        TYPE(atomic_kind_type), ALLOCATABLE, TARGET, DIMENSION(:)     :: my_atomic_kind_set
        INTEGER                                  :: my_nparticle_kind
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        logical                                  :: qsys_only

           IF (PRESENT(atomic_kinds)) THEN
              my_nparticle_kind = atomic_kinds%n_els
              ALLOCATE(my_atomic_kind_set(SIZE(atomic_kinds%els)))
              my_atomic_kind_set = atomic_kinds%els
           ELSE
              IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
                 my_nparticle_kind = nparticle_kind
                 ALLOCATE(my_atomic_kind_set(SIZE(atomic_kind_set)))
                 my_atomic_kind_set = atomic_kind_set
              END IF
           END IF

           qsys_only = .FALSE.
           if (PRESENT(ext_pres)) then 
               qsys_only = .TRUE.
           else 
               qsys_only = .FALSE.
           end if
     
           K_tot = 0.0D0
           DO iparticle_kind= 1, my_nparticle_kind
              NULLIFY(atomic_kind)
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
!calculate kinetic energy according to the COM since the array is velocities_COM
              DO iparticle_local= 1, nparticle_local
                 if (qsys_only) then 
                      if (ext_pres(iparticle_local)) then
!print *, "ext_pres", ext_pres
                        iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                        K_tot = K_tot + mass*( &
                                   array(1, iparticle)**2 + &
                                   array(2, iparticle)**2 + &
                                   array(3, iparticle)**2  ) 
                      end if
                 else 
!print *, "NOT_PRESENT"
                     iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                     K_tot = K_tot + mass*( &
                                array(1, iparticle)**2 + &
                                array(2, iparticle)**2 + &
                                array(3, iparticle)**2  )
                 end if 
              END DO
           END DO
           K_tot = 0.5*K_tot


   END SUBROUTINE calc_kinetic_tot
 

   SUBROUTINE propagate_wf(ctin, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat)

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

!COMMENT: INDENTATION
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
                    call makeX(mat, H_0, Stt_int, Sttp_int, Sttm_int, Stptp_int, &
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
!COMMENT: WHY CTOUT?
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
                            mat, class_forbidden_hop, active_state)

     !COMPLEX(KIND=dp), dimension(:,:), INTENT(IN)    :: H_t, H_tp
     !REAL(KIND=dp), dimension(:,:), INTENT(IN)       :: Stt, Sttp, Stptp, Stpt, Stmt, Sttm, Stmtm
     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
     INTEGER, dimension(:), INTENT(IN)               :: class_forbidden_hop
     integer, INTENT(IN)                             :: active_state, mat 

     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf', &
     routineP = moduleN//':'//routineN

     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctout, ctotal
     REAL(KIND=dp), dimension(:,:), allocatable      :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE , modified_diab_nace
!     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE1
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_0, H_2, X1, X2, X12!, Ytot
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: Ktot, K2, K3, K4
     REAL(KIND=dp), dimension(:,:), allocatable      :: Imat 
     REAL(KIND=dp), dimension(:,:), allocatable       :: Stt_int, Sttp_int, Stptp_int, &
                                                          Sttm_int, Stmtm_int
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE       :: U_t, U_t1, U_t2, U_tp, U_t0 
     COMPLEX(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: interm_hamil
     REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: interm_umatrix, udot_matrices 
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE         :: Eig_val_t, Eig_val_t1, Eig_val_t2, Eig_val_t0,& 
                                                          Eig_val_tp


     INTEGER                                         :: ELsteps   
     INTEGER, POINTER                                :: itimes
     TYPE(section_vals_type), POINTER                :: adiab_section
     TYPE(aom_env_type)                              :: aom_env
     REAL(KIND=dp)                                   :: deltatMD, deltatE
     REAL(KIND=dp)                                   :: ctmag
     REAL(KIND=dp)                                   :: times 
     COMPLEX(KIND=dp)                                :: alpha, beta, ctcomp
     integer                                         :: t2, dim_
     integer                                         :: i,j, inc 
     character(len=1)                                :: zgemm_type

!TEST
     COMPLEX(KIND=dp), dimension(:), allocatable     :: pop
     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctin_adiab
     REAL(KIND=dp), dimension(:,:), allocatable      :: Sterm1, UU, UDU, Sterm2, UDU_old_nace
!     integer :: k
!END TEST
     !Various Allocation


     !Imat is an identity matrix of the dimensions of the other square matrices
     allocate(Imat(mat,mat))
     !define once, then use whenever
     Imat = 0.0D0
     FORALL(j = 1:mat) Imat(j,j) = 1     
     
     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(interm_hamil(3,int(deltatMD/deltatE),mat,mat))
     allocate(interm_umatrix(3,int(deltatMD/deltatE),mat,mat))
     allocate(udot_matrices(3,int(deltatMD/deltatE),mat,mat))
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

!TEST
     allocate(ctin_adiab(mat))
     allocate(pop(mat))
    ! allocate(Ytot(mat,mat))
     allocate(Sterm1(mat,mat)) 
     allocate(UU(mat,mat)) 
     allocate(UDU(mat,mat)) 
     allocate(UDU_old_nace(mat,mat)) 
     allocate(Sterm2(mat,mat)) 
!END TEST
    
     !finally the K matrices, which are the terms in the Runge-Kutta integration
     !note K1 missing, because it's just a constant times X
     !Ktot is the final constant by which c(t) is multiplied
     allocate(Ktot(mat,mat))
     allocate(K2(mat,mat))
     allocate(K3(mat,mat))
     allocate(K4(mat,mat))

     !used to give correct commands to ZGEMX subroutines later
     zgemm_type = 'N'

!COMMENT: INDENTATION
             Hdiff = (deltatE/deltatMD)*(aom_env%H_tp - aom_env%H_t)
             Sttdiff = (deltatE/deltatMD)*(aom_env%Stptp - aom_env%Stt)
             Sttpdiff = (deltatE/deltatMD)*(aom_env%Sttp - aom_env%Stmt)
             Sttmdiff = (deltatE/deltatMD)*(aom_env%Stpt - aom_env%Sttm)
             Stptpdiff = Sttdiff
             Stmtmdiff = (deltatE/deltatMD)*(aom_env%Stt - aom_env%Stmtm)
             H_0  = aom_env%H_t
             H_2 = aom_env%H_tp


            !> LOOP TO BUILD INTERMEDIATE HAMILTONIANS AND CALCULATE 
            !> RESPECTIVE U MATRICES
!COMMENT: U_t, U_tp, Eig_val_t, Eig_val_tp are not already known?
!WITH FIVE ELECTRONIC TIMESTEP, ONE NEEDS THE FOLLOWING HAMILTONIANS AND U:
!H_t, H_t+deltaEL/2, H_t+deltaEL, H_t+1.5deltaEL, H_t+2deltaEL...
!AT THE END: NINE MATRICES, IN THE LOOP BELOW ONE CALCULATES 15 MATRICES. 
!INDEED: interm_hamil(3,1,:,:) = interm_hamil(1,2,:,:), etc.
!PROBABLY interm_hamil and intemr_umatrix could be only: 
!interm( 9, :, :)
            CALL adiabaticity(-AIMAG(aom_env%H_t), mat, U_t, Eig_val_t, .false.)
            CALL adiabaticity(-AIMAG(aom_env%H_tp), mat, U_tp, Eig_val_tp, .false.)
            do t2=1,int(deltatMD/deltatE)

                   !create 3 others intermediate hamiltonians usefull to build the X
                   !matrix in the following 
                   interm_hamil(1,t2,:,:) = H_0(:,:)
                   CALL adiabaticity(-AIMAG(H_0), mat, U_t0, Eig_val_t0, .false.)
                   interm_umatrix(1,t2,:,:) = U_t0(:,:)

                   H_1 = H_0 + 0.5D0*Hdiff
                   interm_hamil(2,t2,:,:) = H_1(:,:)
                   CALL adiabaticity(-AIMAG(H_1), mat, U_t1, Eig_val_t1, .false.)
                   interm_umatrix(2,t2,:,:) = U_t1(:,:)

                   H_2 = H_1 + 0.5D0*Hdiff
                   interm_hamil(3,t2,:,:) = H_2(:,:)
                   CALL adiabaticity(-AIMAG(H_2), mat, U_t2, Eig_val_t2, .false.)
                   interm_umatrix(3,t2,:,:) = U_t2(:,:)

                   H_0 = H_2
            end do

     
             !> CONSTRUCT Udot MATRICES
             
             !> 1st and 3rd arraies of the super-array:
             udot_matrices(1,1,:,:) = (interm_umatrix(2,1,:,:) - & 
                                    interm_umatrix(1,1,:,:))/(deltatE/2.0) !forwards
            ! udot_matrices(2,1,:,:) = (interm_umatrix(1,3,:,:) - &
             !                             interm_umatrix(1,1,:,:))/(deltatE) !centering
             do i =1,(int(deltatMD/deltatE)-1)
                   udot_matrices(3,i,:,:) =  (interm_umatrix(2,i+1,:,:) - & 
                                                     interm_umatrix(2,i,:,:))/(deltatE)
                   udot_matrices(1,i+1,:,:) = udot_matrices(3,i,:,:)
             end do 
                                                
             udot_matrices(3,int(deltatMD/deltatE),:,:) = (interm_umatrix(3,int(deltatMD/deltatE),:,:) - & 
                                       interm_umatrix(2,int(deltatMD/deltatE),:,:))/(deltatE/2.0) !backwards
             !> 2nd array of the super-array:
             do i =1,int(deltatMD/deltatE)
                   udot_matrices(2,i,:,:) = (interm_umatrix(3,i,:,:) - & 
                                         interm_umatrix(1,i,:,:))/(deltatE) !all centering
             end do



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
                 
                    dim_ = 1
                    !CALL adiabaticity(-AIMAG(aom_env%H_t), mat, U_t, Eig_val_t, .false.)
                    CALL calc_diabatic_nace(mat, interm_hamil(dim_,t2,:,:), Stt_int, Sttp_int, Sttm_int, &
                                                        Stptp_int, Stmtm_int, deltatMD, NACE)
                    !CALL calc_diabatic_nace(mat, aom_env%H_t, Stt_int, Sttp_int, Sttm_int, &
                     !                                   Stptp_int, Stmtm_int, X1, deltatMD, NACE)
                    CALL calc_modified_diab_nace(class_forbidden_hop, NACE, & 
                          udot_matrices, interm_umatrix, modified_diab_nace, mat, active_state, t2, dim_)
                    
                    CALL hs_corrected_makeX(interm_hamil, modified_diab_nace, t2, dim_, X1)


!                    call makeX(mat, aom_env%H_t, Stt_int, Sttp_int, Sttm_int, Stptp_int, &
!                                             Stmtm_int, X1, deltatMD, NACE)
     
                    !create timestepped equivalents
                    !H_1 = H_0 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff
                    
                    dim_ = 2 
!                    CALL adiabaticity(-AIMAG(H_1), mat, U_t1, Eig_val_t1, .false.)
                    CALL calc_diabatic_nace(mat, interm_hamil(dim_,t2,:,:), Stt_int, Sttp_int, Sttm_int, Stptp_int, & 
                                                       Stmtm_int, deltatMD, NACE)

                    CALL calc_modified_diab_nace(class_forbidden_hop, NACE, & 
                          udot_matrices, interm_umatrix, modified_diab_nace, mat, active_state, t2, dim_)
                    CALL hs_corrected_makeX(interm_hamil, modified_diab_nace, t2, dim_, X12)

                   ! H_2 = H_1 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff
                    
                    dim_ = 3  
!                    CALL adiabaticity(-AIMAG(H_2), mat, U_t2, Eig_val_t2, .false.)
                    call calc_diabatic_nace(mat, interm_hamil(dim_,t2,:,:), Stt_int, Sttp_int, Sttm_int, Stptp_int, &
                                                     Stmtm_int, deltatMD, NACE)
                    CALL calc_modified_diab_nace(class_forbidden_hop, NACE, & 
                          udot_matrices, interm_umatrix, modified_diab_nace, mat, active_state, t2, dim_)

                    CALL hs_corrected_makeX(interm_hamil, modified_diab_nace, t2, dim_, X2)

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


!!TEST> TO COMPARE OLD AND NEW NACE!

!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, interm_umatrix(dim_,t2,:,:), mat, modified_diab_nace, mat, 0.0D0, Sterm1, mat)
!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, Sterm1, mat, interm_umatrix(dim_,t2,:,:), mat, 0.0D0, UDU, mat)

!print *, "UDU", UDU
!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, interm_umatrix(dim_,t2,:,:), mat, NACE, mat, 0.0D0, Sterm2, mat)
!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, Sterm2, mat, interm_umatrix(dim_,t2,:,:), mat, 0.0D0, UDU_old_nace, mat)
!print *, "UDU_old_nace", UDU_old_nace

!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, interm_umatrix(dim_,t2,:,:), mat, & 
!                                         udot_matrices(dim_,t2,:,:), mat, 0.0D0, UU, mat)

!print *, "UU", UU
!print *, "NEW_NACE_IN_ADIAB", UDU+UU
!print *, "OLD_NACE_IN_ADIAB", UDU_old_nace+UU

!!do i=1,mat
!!  do j=1,mat
!!    Sterm2(i,j)=0
!!    do k=1,mat
!!    Sterm2(i,j) = Sterm2(i,j) + interm_umatrix(dim_,t2,i,k)*modified_diab_nace(k,j)
!!    end do
!!  end do
!!end do

!!print *, "Sterm2", Sterm2

!do i = 1,mat
!   pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
!   ctin_adiab(i) = cmplx(0.0D0, 0.0D0,dp)
!   do j = 1, mat

!       ctin_adiab(i) = ctin_adiab(i) + cmplx(interm_umatrix(dim_,t2,j,i), 0.0D0,dp)*ctin(j)
!   end do
!   pop(i) = ctin_adiab(i)*conjg(ctin_adiab(i))
!end do

!print *, "ct_ADIAB", ctin_adiab
!print *, "ct_SECOND_STATE", ctin_adiab(2)
!print *, "POP_SECOND", pop(2)
!print *, "POP_SECOND_REAL", REAL(pop(2))

!END_TEST

                    CALL adiab_print_nace(adiab_section, NACE, nadiab, norbitals, itimes, deltatMD, times)
                    CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, &
                                           times)
                    CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,itimes, deltatMD, &
                                         times)
                    CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes,deltatMD, &
                                                    times)

             end do

!stop

!END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP

!DEALLOCATE

       DEALLOCATE(NACE)
       deallocate(interm_hamil,interm_umatrix, udot_matrices)
       deallocate(H_0, Hdiff, H_1, H_2, X1, X2, X12, Ktot, K2, K3, K4)
       deallocate(Sttdiff, Stptpdiff, Stmtmdiff, Sttpdiff, Sttmdiff)
       deallocate(Stt_int, Sttp_int, Stptp_int, Sttm_int, Stmtm_int)
       deallocate(Imat)
!END DEALLOCATE

  END SUBROUTINE prop_wf_HS_corr                                                                

!COMMENT: IF THIS SUBROUTINE IS ONLY TWO LINES, WHY NOT WRITE IT DIRECTLY IN THE
!MAIN ONE
  SUBROUTINE hs_corrected_makeX(interm_hamil, modified_diab_nace, t2, dim_, Ytot)
      
     COMPLEX(KIND=dp), DIMENSION(:,:,:,:), INTENT(IN)     :: interm_hamil
     integer, INTENT(IN)                                  :: t2,dim_ 
     !COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_0, H_2, X1, X2, X12
     COMPLEX(KIND=dp), dimension(:,:), intent(out) :: Ytot
     
     REAL(KIND=dp), dimension(:,:)                 :: modified_diab_nace

     modified_diab_nace = (modified_diab_nace - TRANSPOSE(modified_diab_nace)) / 2

!print *, "modified_diab_nace_MOD", modified_diab_nace

     Ytot = interm_hamil(dim_,t2,:,:) - modified_diab_nace


  END SUBROUTINE hs_corrected_makeX
                                                            
  SUBROUTINE calc_modified_diab_nace(class_forbidden_hop, NACE, udot_matrices, interm_umatrix, & 
                                            modified_diab_nace, mat, active_state, t2, dim_)


     REAL(KIND=dp), DIMENSION(:,:,:,:), INTENT(IN)     :: interm_umatrix, udot_matrices 
     INTEGER, dimension(:), INTENT(IN)          :: class_forbidden_hop
     REAL(KIND=dp), dimension(:,:), INTENT(IN)      :: NACE
     REAL(KIND=dp), dimension(:,:), allocatable, INTENT(OUT)     :: modified_diab_nace
     integer, INTENT(IN)                                         :: active_state,mat,t2,dim_ 
     

     integer                                         :: i,j, k, l
     REAL(KIND=dp), dimension(:,:), allocatable      :: modified_ad_nace
     REAL(KIND=dp), dimension(:,:), allocatable      :: rotated_ad_nace
     
     

     allocate(modified_ad_nace(mat,mat))
     allocate(rotated_ad_nace(mat,mat))
     allocate(modified_diab_nace(mat,mat))


     !> CALCULATE MODIFIED ADIABATIC NACE 
     modified_ad_nace = 0.0D0
!COMMENT: THE DO LOOP IS UNNECESSARY, BECAUSE ONLY i=active_state IS POSSIBLE
!     do i=1,mat
!         if (i .eq. active_state) then
             do j=1,mat
                  if((j .ge. active_state) .and. (j .eq. class_forbidden_hop(j))) then
                      modified_ad_nace(active_state,j) = 0.0D0
                         do k=1,mat
                             modified_ad_nace(active_state,j) = & 
                                modified_ad_nace(active_state,j) + &
                                           interm_umatrix(dim_,t2,k,active_state)*udot_matrices(dim_,t2,k,j)
                             do l=1,mat
                                modified_ad_nace(active_state,j) = & 
                                   modified_ad_nace(active_state,j) + interm_umatrix(dim_,t2,k,active_state)* &
                                        NACE(k,l)*interm_umatrix(dim_,t2,l,j)
!print  *, "NACE(k,l)", NACE(k,l)
 
                             end do
                         end do
                  !else 
                      !modified_ad_nace(i,j) = -modified_ad_nace(j,i)
                  end if 
             end do
         !else 
          !   modified_ad_nace(i,j) = 0.0D0
!         end if
!     end do

     !> SIMMETRIZATION OF modified_ad_nace
     do i=1,mat
         do j=1,mat
              if(j .lt. i) then 
                 modified_ad_nace(i,j) = -modified_ad_nace(j,i)
              end if
         end do
     end do    
  
!print *, "modified_ad_nace", modified_ad_nace
      
     !rotated_ad_nace = modified_ad_nace

     !> ROTATE BACK modified_ad_nace INTO DIABATIC BASIS
!COMMENT: AS modified_ad_nace(k,l) = 0 IF k != active_state, A SHORTCUT FORMULA
!CAN BE USED. THE LOOP OVER k IS PROBABLY NOT USEFUL. 
   do i=1,mat
       do j=1,mat
            rotated_ad_nace(i,j) = 0.0D0
               do k=1,mat
                    do l=1,mat
                         rotated_ad_nace(i,j) =  rotated_ad_nace(i,j) + &
                   interm_umatrix(dim_,t2,i,k)*modified_ad_nace(k,l)*interm_umatrix(dim_,t2,j,l)
!print *, "rotated_ad_nace", rotated_ad_nace
                    end do
               end do
       end do
   end do

!print *, "ROTATED_ARRAY", rotated_ad_nace
     !> CALCULATION OF THE NEW MODIFIED DIABATIC NACV
!NB:THERE WAS A PROBLEM WITH THE SIGN IN THE CORRECTION I JUST CHANGE THIS SIGN
!HERE IN THE modified_diab_nace, BUT I STILL HAVE TO FIND WHERE THE EFFECTIVE
!ERROR IN SIGN IS
!     modified_diab_nace = rotated_ad_nace + NACE 
     modified_diab_nace = NACE - rotated_ad_nace

!print *, "NACE", NACE
!print *, "modified_diab_nace_NEW", modified_diab_nace
   
 
 
     deallocate(modified_ad_nace, rotated_ad_nace)

  END SUBROUTINE calc_modified_diab_nace

  SUBROUTINE propagate_adiab_wf(E_t, E_tp, ctin, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat)
     
     REAL(KIND=dp), DIMENSION(:), INTENT(IN)         :: E_t, E_tp 
     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
     TYPE(aom_env_type), intent(in)                  :: aom_env
 
 
     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_adiab_wf', &
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
                    call make_adiabX(mat, E_0, NACEint, X1)
 
 
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
                       
