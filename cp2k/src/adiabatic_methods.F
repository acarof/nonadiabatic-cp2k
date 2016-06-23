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
!>      carof, 2015.12.11: creation, put only ffsh_aom subroutines
! *****************************************************************************


MODULE adiabatic_methods

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

  USE jacob_types,                     ONLY: aom_basket_type 
!  USE timings,                         ONLY: timeset,&
!                                             timestop
!  USE jacob_classicalpath, ONLY: &
!                 construct_ham_X_SH_full_degen_classicalpath

  USE jacob_surfacehopping,            ONLY:  &
                 construct_ham_X_SH
  USE cp_units,                        ONLY: cp_unit_from_cp2k

  USE jacob_print,                     ONLY: adiab_print_hamiltonian, &
                                             adiab_print_nace,        &
                                             adiab_print_nacv,        &
                                             adiab_print_stateno
  USE distribution_1d_types,           ONLY: distribution_1d_type

  USE jacob_utils,                     ONLY: read_listcoeff

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set


#include "./base/base_uses.f90"

IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'adiabatic_methods'
   
   PUBLIC :: adiab_fssh_routines
   !PUBLIC :: adiab_fssh_routines,          &
   !          adiab_classicalpath_routines

CONTAINS

! *****************************************************************************
!> \brief Call the Fewest Surface Hopping Algorithm with the Analaytic Overlap
!>  \brief Method
!sections/environments
!> \param particles_mix ...
!> \param virial_mix ...
!> \param results_mix ...
!> \param global_forces ...
!> \param virials ...
!> \param results ...
!> \param factor ...
!> \param iforce_eval ...
!> \param nforce_eval ...
!> \param map_index ...
!> \param mapping_section ...
!> \param overwrite ...
!> \author Antoine Carof - 12.2015
! *****************************************************************************

!    SUBROUTINE  adiab_classicalpath_routines(particles_mix, virial_mix, results_mix,   &
!                global_forces, energies_plus, virials, results, nforce_eval,&
!                map_index, mapping_section, adiab_section, mixed_env, my_itimes, &
!                frozen_hamil)

!    TYPE(particle_list_type), POINTER        :: particles_mix
!    TYPE(virial_type), POINTER               :: virial_mix
!    TYPE(cp_result_type), POINTER            :: results_mix
!    TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces
!    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies_plus
!    TYPE(virial_p_type), DIMENSION(:), &
!      POINTER                                :: virials
!    TYPE(cp_result_p_type), DIMENSION(:), &
!      POINTER                                :: results
!    INTEGER, INTENT(IN)                      :: nforce_eval
!    INTEGER, DIMENSION(:), POINTER           :: map_index
!    TYPE(section_vals_type), POINTER         :: mapping_section
!    TYPE(section_vals_type), POINTER         :: adiab_section
!    TYPE(mixed_environment_type), OPTIONAL, &
!      POINTER                                :: mixed_env
!    INTEGER, POINTER                         :: my_itimes
!    LOGICAL, OPTIONAL                        :: frozen_hamil

!    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_classicalpath_routines', &
!      routineP = moduleN//':'//routineN

!    INTEGER                                  :: iparticle, jparticle, natom,iforce_eval
!    INTEGER                                  :: aom_natom, nadiab
!    INTEGER                                  :: norbitals, ext_init, ext_atoms
!    INTEGER                                  :: ext_el_steps, init
!    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies
!    REAL(kind=dp)                            :: ext_cutoff, ext_scaling
!    REAL(KIND=dp)                            :: ext_timestep
!    REAL(kind=dp)                            :: ext_cbar, ext_minao
!    REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions, velocities
    !REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_plus, positions_minus
!    REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_plus
!    REAL(kind=dp)                            :: factor
!    TYPE(section_vals_type), POINTER         :: coupling_section
!    CHARACTER(LEN=default_string_length)     :: coupling_keyword
!    TYPE(cp_sll_val_type), POINTER           :: listcoeff
!    TYPE(aom_basket_type)                    :: aom_basket
!    integer, dimension(:), ALLOCATABLE       :: ext_elements
!    LOGICAL, dimension(:), ALLOCATABLE       :: ext_pres
!    REAL(KIND=dp), dimension(:), ALLOCATABLE :: ext_scoeff, ext_pcoeff
    !COMPLEX(KIND=dp), DIMENSION(:),  POINTER :: coefficients
!    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    !REAL(KIND=dp), DIMENSION(:),POINTER      :: populations
!    REAL(KIND=dp), DIMENSION(:,:,:,:),POINTER :: nacv
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!    INTEGER, POINTER                         :: itimes
!    LOGICAL                                  :: my_frozen_hamil

!    my_frozen_hamil = .FALSE.
!    IF (PRESENT(frozen_hamil)) my_frozen_hamil = frozen_hamil

!    print*, "la", my_itimes
 !   itimes => my_itimes
 !       if (ASSOCIATED(itimes)) then
 !   print*, "ici", itimes
!    natom = SIZE(global_forces(1)%forces,2)
!print*, energies_plus, "pppp"
!    ALLOCATE(energies(nforce_eval))
!    ALLOCATE(energies_plus(nforce_eval))
!print*, energies_plus, "pppp"
!    ALLOCATE(positions(natom,3))
!    ALLOCATE(positions_plus(natom,3))
!    ALLOCATE(velocities(natom,3))
!    factor = 1

!> Find the method to calculate the off-diagonal elements
!    CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
!    CALL get_subsys_map_index(mapping_section, natom, 1, nforce_eval, map_index)
!> Choose the method to calculate the off-diagonal elements
!    SELECT CASE(coupling_keyword)
!    CASE("AOM")
!> Get the information for the previous timesteps
!        aom_basket = mixed_env%aom_basket
!        DO iparticle = 1, natom
!           jparticle = map_index(iparticle)
!           positions_plus(jparticle,:)     = particles_mix%els(jparticle)%r(:)
!           positions(jparticle,:)  = particles_mix%els(jparticle)%rm1(:)
!           positions_m2(jparticle,:)  = particles_mix%els(jparticle)%rm2(:)
!        END DO
!        energies_plus(:) = particles_mix%els(1)%energies(:)
!> Get the coupling section and all the information 
!        coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
!        CALL section_vals_val_get(coupling_section,"_DEFAULT_KEYWORD_",&
!             n_rep_val=aom_natom)
!        CALL section_vals_val_get(coupling_section,"NUMBER_DIABATIC_STATES",&
!             i_val=nadiab)
!        CALL section_vals_val_get(coupling_section,"NUMBER_ATOMS_PER_SITE",&
!             i_val=ext_atoms)
!        CALL section_vals_val_get(coupling_section,"NUMBER_ORBITALS",&
!             i_val=norbitals)
!        CALL section_vals_val_get(coupling_section,"CUT_OFF",&
!             r_val=ext_cutoff)
!        CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
!             r_val=ext_scaling)
!        CALL section_vals_val_get(coupling_section,"FIRST_DIABAT",&
!             i_val=ext_init)
!        CALL section_vals_val_get(coupling_section,"CBAR",&
!             r_val=ext_cbar)
!        CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
!             r_val=ext_minao)
!        CALL section_vals_val_get(coupling_section,"NUCLEAR_TIMESTEP",&
!             r_val=ext_timestep)
!        CALL section_vals_val_get(coupling_section,"ELECTRONIC_PARTIAL_STEP",&
!             i_val=ext_el_steps)
!        CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff)
!        CALL read_listcoeff(listcoeff, ext_elements, ext_pres, & 
!                            ext_scoeff, ext_pcoeff, natom)
!    itimes => my_itimes
        !if (ASSOCIATED(itimes)) then
          ! print*, "ASSOCIATED 1"
          ! print*, "itimes", itimes
        !endif
!> Main CALL: calculate the Hamiltonian, propagate the electronic dynamics
!        CALL construct_ham_X_SH_full_degen_classicalpath(nadiab=nadiab, &
!             ext_atoms=ext_atoms, norbitals=norbitals, &
!             ext_cutoff = ext_cutoff, ext_scaling = ext_scaling, &
!             ext_init = ext_init, ext_cbar = ext_cbar, &
!             ext_minao = ext_minao, ext_timestep = ext_timestep, &
!             ext_el_steps = ext_el_steps, &
!             positions_plus=positions_plus, &
!             energies_plus = energies_plus, &
!             aom_basket = aom_basket, ext_elements=ext_elements, &
!             ext_pres = ext_pres, ext_scoeff=ext_scoeff, &
!             ext_pcoeff=ext_pcoeff, itimes=itimes, adiab_section=adiab_section)
!print*, 'here we are 02'
            
        !PRINT COEFFICIENTS, POPULATIONS AND HAMILTONIAN
!        ALLOCATE(coefficients(norbitals*nadiab))
!        ALLOCATE(populations(norbitals*nadiab))
!        ALLOCATE(hamiltonian(norbitals*nadiab, norbitals*nadiab))
!        ALLOCATE(nacv(3, norbitals*nadiab, norbitals*nadiab,norbitals*nadiab*ext_atoms))

!> Get the renewed information from aom_basket
!        coefficients = aom_basket%ctin
!print*, coefficients, "mmmm"
!        populations  = aom_basket%population
!        hamiltonian  = aom_basket%hamiltonian
!        nacv         = aom_basket%nacv
        
!print*, hamiltonian, "mmmm"
!> Print the coefficients, the populations, the Hamiltonian, the density matrix
!> and the NACV
!        CALL adiab_print_coeff(adiab_section, coefficients,nadiab,norbitals,itimes)
!        CALL adiab_print_pop(adiab_section, populations, nadiab,norbitals,itimes)
!        if (itimes.ge.2) then
!           CALL adiab_print_hamiltonian(adiab_section, hamiltonian, nadiab,norbitals,itimes)
!        CALL adiab_print_matrix_density(adiab_section, coefficients, nadiab,norbitals,itimes)
!           CALL adiab_print_nacv(adiab_section, nacv, nadiab, norbitals, ext_atoms,itimes)        
!        endif
 
!> IMPORTANT: store all aom_basket for the future timesteps
!        mixed_env%aom_basket=aom_basket
!    CASE DEFAULT
!        CPABORT("There is no method to calculate the coupling")
!    END SELECT

!> IMPORTANT : GIVE BACK THE FORCES TO FIST
!    DO iparticle = 1, natom
!      jparticle = map_index(iparticle)
!> DEPEND IF WE WORK IN THE FROZEN HAMILTONIAN CASE OR NOT
!      IF (my_frozen_hamil) THEN
!          particles_mix%els(jparticle)%f(:)= 0.0D0
!          particles_mix%els(jparticle)%v(:)= 0.0D0
!      ELSE 
!          CALL section_vals_val_get(coupling_section,"FIRST_DIABAT",&
!               i_val=init)
!          particles_mix%els(jparticle)%f(:)= global_forces(init)%forces(:,iparticle)
!print*,  "part forces", particles_mix%els(jparticle)%f(:) 
!      ENDIF
!      positions(jparticle,:)  = particles_mix%els(jparticle)%r(:)
!    END DO


!    iforce_eval = nforce_eval
    ! Mixing Virial
    !Copy-Paste of mixed_map_forces
    !Probably useless
    ! But overwrite is useless
!    IF (virial_mix%pv_availability) THEN
!       IF (overwrite) CALL zero_virial(virial_mix,reset=.FALSE.)
!       virial_mix%pv_total      = virial_mix%pv_total + factor*virials(iforce_eval)%virial%pv_total
!       virial_mix%pv_kinetic    = virial_mix%pv_kinetic + factor*virials(iforce_eval)%virial%pv_kinetic
!       virial_mix%pv_virial     = virial_mix%pv_virial + factor*virials(iforce_eval)%virial%pv_virial
!       virial_mix%pv_xc         = virial_mix%pv_xc + factor*virials(iforce_eval)%virial%pv_xc
!       virial_mix%pv_fock_4c    = virial_mix%pv_fock_4c + factor*virials(iforce_eval)%virial%pv_fock_4c
!       virial_mix%pv_constraint = virial_mix%pv_constraint + factor*virials(iforce_eval)%virial%pv_constraint
!    END IF

    ! Deallocate map_index array
!    IF (ASSOCIATED(map_index)) THEN
!       DEALLOCATE(map_index)
!    END IF
    !DEALLOCATE(coefficients, populations, hamiltonian, nacv)
!    DEALLOCATE(hamiltonian, nacv)

!   END SUBROUTINE adiab_classicalpath_routines




    SUBROUTINE  adiab_fssh_routines(particles_mix, virial_mix, results_mix,      &
                global_forces, energies_plus, virials, results, nforce_eval,     &
                map_index, mapping_section, adiab_section, mixed_env, my_itimes, &
                kinetic_energy, method_keyword, final_energy, atomic_kinds,      &
                local_particles)

    TYPE(particle_list_type), POINTER        :: particles_mix
    TYPE(virial_type), POINTER               :: virial_mix
    TYPE(cp_result_type), POINTER            :: results_mix
    TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces
    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies_plus
    REAL(KIND=dp), INTENT(OUT)               :: kinetic_energy
    TYPE(virial_p_type), DIMENSION(:), &
      POINTER                                :: virials
    TYPE(cp_result_p_type), DIMENSION(:), &
      POINTER                                :: results
    INTEGER, INTENT(IN)                      :: nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: map_index
    TYPE(section_vals_type), POINTER         :: mapping_section
    TYPE(section_vals_type), POINTER         :: adiab_section
    TYPE(mixed_environment_type), OPTIONAL, &
      POINTER                                :: mixed_env
    INTEGER, POINTER                         :: my_itimes
    CHARACTER(LEN=default_string_length)     :: method_keyword
    REAL(KIND=dp), INTENT(OUT)               :: final_energy
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(distribution_1d_type), POINTER      :: local_particles


    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_classicalpath_routines', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, jparticle, natom,iforce_eval
    INTEGER                                  :: nadiab
    INTEGER                                  :: norbitals, aom_atoms, init
    REAL(KIND=dp)                            :: MDtimestep
    REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_plus
!    REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
!    REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_plus
    REAL(kind=dp), DIMENSION(:,:),POINTER    :: forces_adiab
    REAL(kind=dp)                            :: factor
    TYPE(section_vals_type), POINTER         :: coupling_section
    CHARACTER(LEN=default_string_length)     :: coupling_keyword
    TYPE(aom_basket_type)                    :: aom_basket
    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: nace
    REAL(KIND=dp), DIMENSION(:,:,:,:),POINTER :: nacv
    INTEGER, POINTER                         :: itimes
    INTEGER                                  :: old_state, new_state
    REAL(KIND=dp)                            :: rescale_factor, energy_adiabat
    LOGICAL                                  :: is_hop
    LOGICAL                                  :: my_bo_keyword

    natom = SIZE(global_forces(1)%forces,2)
    ALLOCATE(positions_plus(natom,3))
!    ALLOCATE(velocities(natom,3))
!    ALLOCATE(velocities_plus(natom,3))
    ALLOCATE(forces_adiab(3,natom))
    factor = 1

    IF (method_keyword.eq."FSSH") THEN
       my_bo_keyword = .FALSE.
    ELSE 
       my_bo_keyword = .TRUE.
    ENDIF
    !IF (PRESENT(bo_keyword)) my_bo_keyword = bo_keyword


!> Find the method to calculate the off-diagonal elements
    CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
    CALL get_subsys_map_index(mapping_section, natom, 1, nforce_eval, map_index)
!> Choose the method to calculate the off-diagonal elements
    SELECT CASE(coupling_keyword)
    CASE("AOM")
!> Get the information for the previous timesteps
        aom_basket = mixed_env%aom_basket
        IF (my_itimes.eq.0) THEN
           ALLOCATE(aom_basket%velocities(3, natom))
        ENDIF
        DO iparticle = 1, natom
           jparticle = map_index(iparticle)
           positions_plus(jparticle,:)        = particles_mix%els(jparticle)%r(:)
           IF (my_itimes.eq.0) THEN 
              aom_basket%velocities(:, jparticle) = particles_mix%els(jparticle)%v(:)
            ENDIF
!           velocities(jparticle,:)           = particles_mix%els(jparticle)%v(:)
        END DO
!> Get the coupling section and all the information 
        coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
        CALL section_vals_val_get(coupling_section,"NUMBER_DIABATIC_STATES",&
             i_val=nadiab)
        CALL section_vals_val_get(coupling_section,"NUMBER_ATOMS_PER_SITE",&
             i_val=aom_atoms)
        CALL section_vals_val_get(coupling_section,"NUMBER_ORBITALS",&
             i_val=norbitals)
        CALL section_vals_val_get(coupling_section,"NUCLEAR_TIMESTEP",&
             r_val=MDtimestep)
        itimes => my_itimes
!> Main CALL: calculate the Hamiltonian, propagate the electronic dynamics
        CALL construct_ham_X_SH(nadiab, aom_atoms, norbitals, MDtimestep, &
             positions_plus, energies_plus, kinetic_energy, aom_basket, itimes, &
             is_hop, global_forces, coupling_section, adiab_section,& 
             my_bo_keyword, rescale_factor, energy_adiabat, atomic_kinds, &
             local_particles)
            
        !PRINT COEFFICIENTS, POPULATIONS AND HAMILTONIAN
        ALLOCATE(hamiltonian(norbitals*nadiab, norbitals*nadiab))
        ALLOCATE(nace(norbitals*nadiab, norbitals*nadiab))
        ALLOCATE(nacv(3, norbitals*nadiab, norbitals*nadiab,norbitals*nadiab*aom_atoms))

!> Get the renewed information from aom_basket
        hamiltonian  = aom_basket%hamiltonian
        nace         = aom_basket%nace
        nacv         = aom_basket%nacv
        new_state    = aom_basket%new_state
        old_state    = aom_basket%old_state
        forces_adiab = aom_basket%forces_adiab
!> Print the coefficients, the populations, the Hamiltonian, the density matrix
!> and the NACV
        CALL adiab_print_hamiltonian(adiab_section, hamiltonian, nadiab,norbitals,itimes, MDtimestep)
        CALL adiab_print_nace(adiab_section, nace, nadiab,norbitals,itimes, MDtimestep)
        CALL adiab_print_nacv(adiab_section, nacv, nadiab, norbitals, aom_atoms,itimes, MDtimestep)        
        CALL adiab_print_stateno(adiab_section, new_state, aom_atoms,itimes, MDtimestep)


!> IMPORTANT: store all aom_basket for the future timesteps
        mixed_env%aom_basket=aom_basket
!> IMPORTANT

    CASE DEFAULT
        CPABORT("There is no method to calculate the coupling")
    END SELECT

!> IMPORTANT : GIVE BACK THE FORCES TO FIST
    DO iparticle = 1, natom
      jparticle = map_index(iparticle)
      SELECT CASE(method_keyword)
      CASE("CLASSICAL_PATH")
          CALL section_vals_val_get(coupling_section,"FIRST_DIABAT",&
               i_val=init)
          particles_mix%els(jparticle)%f(:)= global_forces(init)%forces(:,iparticle)
          final_energy = energies_plus(init)
      CASE("FROZEN_HAMILTONIAN")
          particles_mix%els(jparticle)%f(:)= 0.0D0
          particles_mix%els(jparticle)%v(:)= 0.0D0
          final_energy = 0.0D0
      CASE("BORN_OPPENHEIMER")
        particles_mix%els(jparticle)%f(:)= forces_adiab(:,iparticle)
        final_energy = energy_adiabat
      CASE("FSSH")
        particles_mix%els(jparticle)%f(:)= forces_adiab(:,iparticle)
        final_energy = energy_adiabat
!> If hop, rescale the velocities (frustrated hop present in the
!> construct_ham_X_SH subroutine)
!        IF (is_hop) THEN
!              particles_mix%els(jparticle)%v(:)   = velocities_plus(iparticle,:) 
!print*, "HER"
!print*, velocities_plus
!        ENDIF
      END SELECT
!      IF (itimes.eq.0) THEN
!> Starting site number
!        CALL section_vals_val_get(coupling_section,"FIRST_DIABAT",&
!             i_val=init)
!       particles_mix%els(jparticle)%f(:)= global_forces(init)%forces(:,iparticle) 
!      ELSE     
!        particles_mix%els(jparticle)%f(:)= forces_adiab(:,iparticle)
!      ENDIF
    END DO

    iforce_eval = nforce_eval
    ! Mixing Virial
    !Copy-Paste of mixed_map_forces
    !Probably useless
    ! But overwrite is useless
    IF (virial_mix%pv_availability) THEN
!       IF (overwrite) CALL zero_virial(virial_mix,reset=.FALSE.)
       virial_mix%pv_total      = virial_mix%pv_total + factor*virials(iforce_eval)%virial%pv_total
       virial_mix%pv_kinetic    = virial_mix%pv_kinetic + factor*virials(iforce_eval)%virial%pv_kinetic
       virial_mix%pv_virial     = virial_mix%pv_virial + factor*virials(iforce_eval)%virial%pv_virial
       virial_mix%pv_xc         = virial_mix%pv_xc + factor*virials(iforce_eval)%virial%pv_xc
       virial_mix%pv_fock_4c    = virial_mix%pv_fock_4c + factor*virials(iforce_eval)%virial%pv_fock_4c
       virial_mix%pv_constraint = virial_mix%pv_constraint + factor*virials(iforce_eval)%virial%pv_constraint
    END IF

    ! Deallocate map_index array
    IF (ASSOCIATED(map_index)) THEN
       DEALLOCATE(map_index)
    END IF
    
    DEALLOCATE(hamiltonian, nacv, nace)
    DEALLOCATE(positions_plus, forces_adiab)

   END SUBROUTINE adiab_fssh_routines

END MODULE adiabatic_methods

