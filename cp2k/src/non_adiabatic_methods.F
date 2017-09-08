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


MODULE non_adiabatic_methods

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

  USE aom_types,                       ONLY: aom_env_type 

  USE aom_main,                        ONLY: construct_ham_X_SH

  USE distribution_1d_types,           ONLY: distribution_1d_type

  USE aom_utils,                       ONLY: read_listcoeff

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE aom_initialization,              ONLY: create_aom_env,   &
                                             propagate_aom_namd_first_step

#include "./base/base_uses.f90"

IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'non_adiabatic_methods'
   
   PUBLIC :: adiab_fssh_routines

CONTAINS

    SUBROUTINE  adiab_fssh_routines(particles_mix, virial_mix, results_mix,      &
                global_forces, energies_tp, virials, results, nforce_eval,     &
                map_index, mapping_section, adiab_section, mixed_env, my_itimes, &
                method_keyword, final_energy, atomic_kinds,      &
                local_particles)

    TYPE(particle_list_type), POINTER        :: particles_mix
    TYPE(virial_type), POINTER               :: virial_mix
    TYPE(cp_result_type), POINTER            :: results_mix
    TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces
    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies_tp
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
    REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_tp
    REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_tp
    REAL(kind=dp), DIMENSION(:,:),POINTER    :: forces_adiab
    REAL(kind=dp)                            :: factor
    TYPE(section_vals_type), POINTER         :: coupling_section
    CHARACTER(LEN=default_string_length)     :: coupling_keyword
    TYPE(aom_env_type)                       :: aom_env
    INTEGER, POINTER                         :: itimes
    REAL(KIND=dp)                            :: energy_adiabat
    LOGICAL                                  :: is_hop
    LOGICAL                                  :: my_bo_keyword
    LOGICAL                                  :: my_test_keyword

    natom = SIZE(global_forces(1)%forces,2)
    ALLOCATE(positions_tp(3,natom))
    ALLOCATE(velocities_tp(3,natom))
    ALLOCATE(forces_adiab(3,natom))
    factor = 1


    my_bo_keyword = .TRUE.
    my_test_keyword = .FALSE.
    IF (method_keyword.eq."FSSH") THEN
       my_bo_keyword = .FALSE.
    ELSE IF (method_keyword.eq."TEST_HOP") THEN
       my_bo_keyword = .FALSE.
       my_test_keyword = .TRUE. 
    ENDIF


    !> Find the method to calculate the off-diagonal elements
    CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
    CALL get_subsys_map_index(mapping_section, natom, 1, nforce_eval, map_index)
    !> Choose the method to calculate the off-diagonal elements
    SELECT CASE(coupling_keyword)
    CASE("AOM")
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
        DO iparticle = 1, natom
           jparticle = map_index(iparticle)
           positions_tp(:,jparticle)   = particles_mix%els(jparticle)%r(:)
           IF (my_itimes.eq.0) THEN
                velocities_tp(:, jparticle) = particles_mix%els(jparticle)%v(:)
           ENDIF
        END DO
        IF (my_itimes.eq.0) THEN
           CALL create_aom_env(nadiab, aom_atoms, norbitals, natom, &
                               positions_tp, coupling_section, adiab_section, & 
                               aom_env)
           CALL propagate_aom_namd_first_step(nadiab, aom_atoms, norbitals, MDtimestep, &
                     positions_tp, velocities_tp, energies_tp, aom_env, itimes, &
                     is_hop, global_forces, coupling_section, adiab_section,& 
                     my_bo_keyword, energy_adiabat, atomic_kinds, &
                     local_particles, my_test_keyword, method_keyword)
        ELSE
           !> Get the information for the previous timesteps
           aom_env = mixed_env%aom_env
           !> Velocities are taken directly from the integrator, note from
           !> velocities_tp
           !aom_env%velocities_t = velocities_tp
           IF (aom_env%successful_hop) THEN
              aom_env%forces_adiab = aom_env%forces_adiab_new
              energy_adiabat = aom_env%energy_new
              aom_env%futur_decohere = .TRUE.
              aom_env%successful_hop = .FALSE.
              aom_env%forces_adiab_new = 0.0D0
              aom_env%new_state = -1
              aom_env%old_state = -1
           ELSE
              !> Main CALL: calculate the Hamiltonian, propagate the electronic dynamics
              CALL construct_ham_X_SH(nadiab, aom_atoms, norbitals, MDtimestep, &
                            positions_tp, energies_tp, aom_env, itimes, &
                            is_hop, global_forces, coupling_section, adiab_section,& 
                            my_bo_keyword, energy_adiabat, atomic_kinds, &
                            local_particles, my_test_keyword)
           ENDIF
        ENDIF
        !> IMPORTANT: store all aom_env for the future timesteps
        mixed_env%aom_env=aom_env
    CASE DEFAULT
        CPABORT("There is no method to calculate the coupling")
    END SELECT


    !> IMPORTANT : GIVE BACK THE FORCES TO FIST
    forces_adiab = aom_env%forces_adiab
    DO iparticle = 1, natom
      jparticle = map_index(iparticle)
      SELECT CASE(method_keyword)
      CASE("CLASSICAL_PATH")
          CALL section_vals_val_get(adiab_section,"FIRST_DIABAT",&
               i_val=init)
          particles_mix%els(jparticle)%f(:)= global_forces(init)%forces(:,iparticle)
          final_energy = energies_tp(init)
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
      CASE("TEST_HOP")
        particles_mix%els(jparticle)%f(:)= 0.0D0
        final_energy = energy_adiabat
      CASE("GALILEAN")
        particles_mix%els(jparticle)%f(:)= 0.0D0
        final_energy = energy_adiabat
      END SELECT
    END DO

    iforce_eval = nforce_eval
    ! Mixing Virial
    !Copy-Paste of mixed_map_forces
    !Probably useless
    ! But overwrite is useless
    IF (virial_mix%pv_availability) THEN
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
    DEALLOCATE(positions_tp, velocities_tp, forces_adiab)
   END SUBROUTINE adiab_fssh_routines

END MODULE non_adiabatic_methods

