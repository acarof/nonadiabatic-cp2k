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
  USE jacob_classicalpath, ONLY: &
                 construct_ham_X_SH_full_degen_classicalpath

  USE jacob_surfacehopping,            ONLY:  &
                 construct_ham_X_SH

#include "./base/base_uses.f90"

IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'adiabatic_methods'
   
   PUBLIC :: adiab_fssh_routines,&
             adiab_classicalpath_routines 

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

    SUBROUTINE  adiab_classicalpath_routines(particles_mix, virial_mix, results_mix,   &
                global_forces, energies_plus, virials, results, nforce_eval,&
                map_index, mapping_section, adiab_section, mixed_env, my_itimes)

    TYPE(particle_list_type), POINTER        :: particles_mix
    TYPE(virial_type), POINTER               :: virial_mix
    TYPE(cp_result_type), POINTER            :: results_mix
    TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces
    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies_plus
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

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_classicalpath_routines', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, jparticle, natom,iforce_eval
    INTEGER                                  :: aom_natom, nadiab
    INTEGER                                  :: norbitals, ext_init, ext_atoms
    INTEGER                                  :: ext_el_steps
!    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies
    REAL(kind=dp)                            :: ext_cutoff, ext_scaling
    REAL(KIND=dp)                            :: ext_timestep
    REAL(kind=dp)                            :: ext_cbar, ext_minao
!    REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions, velocities
    !REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_plus, positions_minus
    REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_plus
    REAL(kind=dp)                            :: factor
    TYPE(section_vals_type), POINTER         :: coupling_section
    CHARACTER(LEN=default_string_length)     :: coupling_keyword
    TYPE(cp_sll_val_type), POINTER           :: listcoeff
    TYPE(aom_basket_type)                    :: aom_basket
    integer, dimension(:), ALLOCATABLE       :: ext_elements
    LOGICAL, dimension(:), ALLOCATABLE       :: ext_pres
    REAL(KIND=dp), dimension(:), ALLOCATABLE :: ext_scoeff, ext_pcoeff
    COMPLEX(KIND=dp), DIMENSION(:),  POINTER :: coefficients
    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    REAL(KIND=dp), DIMENSION(:),POINTER      :: populations
    REAL(KIND=dp), DIMENSION(:,:,:,:),POINTER :: nacv
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, POINTER                         :: itimes

 !   print*, "la", my_itimes
 !   itimes => my_itimes
 !       if (ASSOCIATED(itimes)) then
 !   print*, "ici", itimes
    natom = SIZE(global_forces(1)%forces,2)
!print*, energies_plus, "pppp"
!    ALLOCATE(energies(nforce_eval))
!    ALLOCATE(energies_plus(nforce_eval))
!print*, energies_plus, "pppp"
!    ALLOCATE(positions(natom,3))
    ALLOCATE(positions_plus(natom,3))
!    ALLOCATE(velocities(natom,3))
    factor = 1

!> Find the method to calculate the off-diagonal elements
    CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
    CALL get_subsys_map_index(mapping_section, natom, 1, nforce_eval, map_index)
!> Choose the method to calculate the off-diagonal elements
    SELECT CASE(coupling_keyword)
    CASE("AOM")
!> Get the information for the previous timesteps
        aom_basket = mixed_env%aom_basket
        DO iparticle = 1, natom
           jparticle = map_index(iparticle)
           positions_plus(jparticle,:)     = particles_mix%els(jparticle)%r(:)
!           positions(jparticle,:)  = particles_mix%els(jparticle)%rm1(:)
!           positions_m2(jparticle,:)  = particles_mix%els(jparticle)%rm2(:)
        END DO
!        energies_plus(:) = particles_mix%els(1)%energies(:)
!> Get the coupling section and all the information 
        coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
        CALL section_vals_val_get(coupling_section,"_DEFAULT_KEYWORD_",&
             n_rep_val=aom_natom)
        CALL section_vals_val_get(coupling_section,"NUMBER_DIABATIC_STATES",&
             i_val=nadiab)
        CALL section_vals_val_get(coupling_section,"NUMBER_ATOMS_PER_SITE",&
             i_val=ext_atoms)
        CALL section_vals_val_get(coupling_section,"NUMBER_ORBITALS",&
             i_val=norbitals)
        CALL section_vals_val_get(coupling_section,"CUT_OFF",&
             r_val=ext_cutoff)
        CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
             r_val=ext_scaling)
        CALL section_vals_val_get(coupling_section,"INIT_SITES",&
             i_val=ext_init)
        CALL section_vals_val_get(coupling_section,"CBAR",&
             r_val=ext_cbar)
        CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
             r_val=ext_minao)
        CALL section_vals_val_get(coupling_section,"NUCLEAR_TIMESTEP",&
             r_val=ext_timestep)
        CALL section_vals_val_get(coupling_section,"ELECTRONIC_PARTIAL_STEP",&
             i_val=ext_el_steps)
        CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff)
        CALL read_listcoeff(listcoeff, ext_elements, ext_pres, & 
                            ext_scoeff, ext_pcoeff, natom)
    itimes => my_itimes
        !if (ASSOCIATED(itimes)) then
          ! print*, "ASSOCIATED 1"
          ! print*, "itimes", itimes
        !endif
!> Main CALL: calculate the Hamiltonian, propagate the electronic dynamics
        CALL construct_ham_X_SH_full_degen_classicalpath(nadiab=nadiab, &
             ext_atoms=ext_atoms, norbitals=norbitals, &
             ext_cutoff = ext_cutoff, ext_scaling = ext_scaling, &
             ext_init = ext_init, ext_cbar = ext_cbar, & 
             ext_minao = ext_minao, ext_timestep = ext_timestep, &
             ext_el_steps = ext_el_steps, &
             positions_plus=positions_plus, &
             energies_plus = energies_plus, &
             aom_basket = aom_basket, ext_elements=ext_elements, &
             ext_pres = ext_pres, ext_scoeff=ext_scoeff, &
             ext_pcoeff=ext_pcoeff, itimes=itimes)
!print*, 'here we are 02'
            
        !PRINT COEFFICIENTS, POPULATIONS AND HAMILTONIAN
        ALLOCATE(coefficients(norbitals*nadiab))
        ALLOCATE(populations(norbitals*nadiab))
        ALLOCATE(hamiltonian(norbitals*nadiab, norbitals*nadiab))
        ALLOCATE(nacv(3, norbitals*nadiab, norbitals*nadiab,2*norbitals*nadiab*ext_atoms))

!> Get the renewed information from aom_basket
        coefficients = aom_basket%ctin
print*, coefficients, "mmmm"
        populations  = aom_basket%population
        hamiltonian  = aom_basket%hamiltonian
        nacv         = aom_basket%nacv
print*, hamiltonian, "mmmm"
!> Print the coefficients, the populations, the Hamiltonian, the density matrix
!> and the NACV
        CALL adiab_print_coeff(adiab_section, coefficients, nadiab,norbitals)
        CALL adiab_print_pop(adiab_section, populations, nadiab,norbitals)
        CALL adiab_print_hamiltonian(adiab_section, hamiltonian, nadiab,norbitals)
        CALL adiab_print_matrix_density(adiab_section, coefficients, nadiab,norbitals)
        CALL adiab_print_nacv(adiab_section, nacv, nadiab, norbitals, ext_atoms)        

!> IMPORTANT: store all aom_basket for the future timesteps
        mixed_env%aom_basket=aom_basket
    CASE DEFAULT
        CPABORT("There is no method to calculate the coupling")
    END SELECT

!> IMPORTANT : GIVE BACK THE FORCES TO FIST
    DO iparticle = 1, natom
      jparticle = map_index(iparticle)
      particles_mix%els(jparticle)%f(:)= global_forces(1)%forces(:,iparticle)
!      positions(jparticle,:)  = particles_mix%els(jparticle)%r(:)
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

   END SUBROUTINE adiab_classicalpath_routines




    SUBROUTINE  adiab_fssh_routines(particles_mix, virial_mix, results_mix,      &
                global_forces, energies_plus, virials, results, nforce_eval,     &
                map_index, mapping_section, adiab_section, mixed_env, my_itimes, &
                kinetic_energies)

    TYPE(particle_list_type), POINTER        :: particles_mix
    TYPE(virial_type), POINTER               :: virial_mix
    TYPE(cp_result_type), POINTER            :: results_mix
    TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces
    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies_plus
    REAL(KIND=dp), DIMENSION(:), POINTER     :: kinetic_energies
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

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_classicalpath_routines', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, jparticle, natom,iforce_eval
    INTEGER                                  :: aom_natom, nadiab
    INTEGER                                  :: norbitals, ext_init, ext_atoms
    INTEGER                                  :: ext_el_steps
    REAL(kind=dp)                            :: ext_cutoff, ext_scaling
    REAL(KIND=dp)                            :: ext_timestep
    REAL(kind=dp)                            :: ext_cbar, ext_minao
    REAL(kind=dp), DIMENSION(:,:),POINTER    :: positions_plus
    REAL(kind=dp)                            :: factor
    TYPE(section_vals_type), POINTER         :: coupling_section
    CHARACTER(LEN=default_string_length)     :: coupling_keyword
    TYPE(cp_sll_val_type), POINTER           :: listcoeff
    TYPE(aom_basket_type)                    :: aom_basket
    integer, dimension(:), ALLOCATABLE       :: ext_elements
    LOGICAL, dimension(:), ALLOCATABLE       :: ext_pres
    REAL(KIND=dp), dimension(:), ALLOCATABLE :: ext_scoeff, ext_pcoeff
    COMPLEX(KIND=dp), DIMENSION(:),  POINTER :: coefficients
    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    REAL(KIND=dp), DIMENSION(:),POINTER      :: populations
    REAL(KIND=dp), DIMENSION(:,:,:,:),POINTER :: nacv
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: collapse
    INTEGER                                  :: old_state, new_state
    REAL(KIND=dp)                            :: scalefactor
    LOGICAL                                  :: is_hop

    natom = SIZE(global_forces(1)%forces,2)
    ALLOCATE(positions_plus(natom,3))
    factor = 1
!> Save the previous state
    old_state = aom_basket%state


!> Find the method to calculate the off-diagonal elements
    CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
    CALL get_subsys_map_index(mapping_section, natom, 1, nforce_eval, map_index)
!> Choose the method to calculate the off-diagonal elements
    SELECT CASE(coupling_keyword)
    CASE("AOM")
!> Get the information for the previous timesteps
        aom_basket = mixed_env%aom_basket
        DO iparticle = 1, natom
           jparticle = map_index(iparticle)
           positions_plus(jparticle,:)     = particles_mix%els(jparticle)%r(:)
        END DO
!> Get the coupling section and all the information 
        coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
        CALL section_vals_val_get(coupling_section,"_DEFAULT_KEYWORD_",&
             n_rep_val=aom_natom)
        CALL section_vals_val_get(coupling_section,"NUMBER_DIABATIC_STATES",&
             i_val=nadiab)
        CALL section_vals_val_get(coupling_section,"NUMBER_ATOMS_PER_SITE",&
             i_val=ext_atoms)
        CALL section_vals_val_get(coupling_section,"NUMBER_ORBITALS",&
             i_val=norbitals)
        CALL section_vals_val_get(coupling_section,"CUT_OFF",&
             r_val=ext_cutoff)
        CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
             r_val=ext_scaling)
        CALL section_vals_val_get(coupling_section,"INIT_SITES",&
             i_val=ext_init)
        CALL section_vals_val_get(coupling_section,"CBAR",&
             r_val=ext_cbar)
        CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
             r_val=ext_minao)
        CALL section_vals_val_get(coupling_section,"NUCLEAR_TIMESTEP",&
             r_val=ext_timestep)
        CALL section_vals_val_get(coupling_section,"ELECTRONIC_PARTIAL_STEP",&
             i_val=ext_el_steps)
        CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff)
        CALL read_listcoeff(listcoeff, ext_elements, ext_pres, & 
                            ext_scoeff, ext_pcoeff, natom)
        CALL section_vals_val_get(coupling_section,"COLLAPSE", l_val=collapse)
    itimes => my_itimes
!> Main CALL: calculate the Hamiltonian, propagate the electronic dynamics
        CALL construct_ham_X_SH(nadiab=nadiab, &
             ext_atoms=ext_atoms, norbitals=norbitals, &
             ext_cutoff = ext_cutoff, ext_scaling = ext_scaling, &
             ext_init = ext_init, ext_cbar = ext_cbar, & 
             ext_minao = ext_minao, ext_timestep = ext_timestep, &
             ext_el_steps = ext_el_steps, &
             positions_plus=positions_plus, &
             energies_plus = energies_plus, &
             aom_basket = aom_basket, ext_elements=ext_elements, &
             ext_pres = ext_pres, ext_scoeff=ext_scoeff, &
             ext_pcoeff=ext_pcoeff, itimes=itimes, &
             collapse = collapse, is_hop=is_hop)
            
        !PRINT COEFFICIENTS, POPULATIONS AND HAMILTONIAN
        ALLOCATE(coefficients(norbitals*nadiab))
        ALLOCATE(populations(norbitals*nadiab))
        ALLOCATE(hamiltonian(norbitals*nadiab, norbitals*nadiab))
        ALLOCATE(nacv(3, norbitals*nadiab, norbitals*nadiab,2*norbitals*nadiab*ext_atoms))

!> Get the renewed information from aom_basket
        coefficients = aom_basket%ctin
        populations  = aom_basket%population
        hamiltonian  = aom_basket%hamiltonian
        nacv         = aom_basket%nacv
        new_state        = aom_basket%state

!> Print the coefficients, the populations, the Hamiltonian, the density matrix
!> and the NACV
        CALL adiab_print_coeff(adiab_section, coefficients, nadiab,norbitals)
        CALL adiab_print_pop(adiab_section, populations, nadiab,norbitals)
        CALL adiab_print_hamiltonian(adiab_section, hamiltonian, nadiab,norbitals)
        CALL adiab_print_matrix_density(adiab_section, coefficients, nadiab,norbitals)
        CALL adiab_print_nacv(adiab_section, nacv, nadiab, norbitals, ext_atoms)        
        CALL adiab_print_stateno(adiab_section, new_state, ext_atoms)

!> Do the velocities rescaling if required
        IF (old_state.NE.new_state) THEN
           CPASSERT(is_hop) !Be sure there is indeed a hop
!> Check that all the kinetic energies are the same for the different force
!> fields
          DO iforce_eval = 2, nforce_eval
             IF (kinetic_energies(iforce_eval).NE.kinetic_energies(1)) THEN
                CPABORT("The different diabatic FF don't share the same kinetic energies")
             END IF
          ENDDO
!Factor to rescale the velocities (Eq. 49 of PAPER_1)
          scalefactor = SQRT( ( kinetic_energies(1) - (energies_plus(new_state) - &
                        energies_plus(old_state) ) ) / kinetic_energies(1) ) 
           DO iparticle = 1, natom
              jparticle = map_index(iparticle)
              particles_mix%els(jparticle)%v(:)   = &
              scalefactor*particles_mix%els(jparticle)%v(:)
           END DO
        ENDIF

!> IMPORTANT: store all aom_basket for the future timesteps
        mixed_env%aom_basket=aom_basket
    CASE DEFAULT
        CPABORT("There is no method to calculate the coupling")
    END SELECT

!> IMPORTANT : GIVE BACK THE FORCES TO FIST
    DO iparticle = 1, natom
      jparticle = map_index(iparticle)
      particles_mix%els(jparticle)%f(:)= global_forces(1)%forces(:,iparticle)
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

   END SUBROUTINE adiab_fssh_routines



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
    !CHARACTER(len=default_string_length)     :: error_message, &
    !                                            my_default_index
                                                !unit_str
                                                !my_default_index, strtmp, &
                                                !unit_str
    INTEGER                                  :: default_id, end_c, &
    !INTEGER                                  :: default_id, end_c, handle, &
                                                !iatom, j, output_unit, &
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
!print*, "b", my_pres2, my_pres
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


   SUBROUTINE adiab_print_coeff(adiab_section, coefficients, nadiab, norbitals)
    TYPE(section_vals_type), POINTER         :: adiab_section
    COMPLEX(KIND=dp), DIMENSION(:),  POINTER :: coefficients
    INTEGER                                  :: nadiab, norbitals
!    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
!    INTEGER                                  :: i
    INTEGER                                  :: iadiab, iorbitals
!    TYPE(particle_type), POINTER             :: particle_set(:)
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    CHARACTER(LEN=default_string_length)     :: my_format
    TYPE(cp_logger_type), POINTER            :: logger


    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%COEFFICIENTS",&
               extension=".xyz", middle_name="coeff", log_filename=.FALSE.)
   
    IF (unit_nr>0) THEN

!       ncoeff = SIZE(coefficients)

       WRITE ( unit_nr, "(2I8)" ) nadiab, norbitals
       WRITE ( unit_nr, '(A)' ) 'Diab. State, Orb., Re(C), I(C) '

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
!          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
!                  name=atm_name)
!          my_format="(A4,"
         my_format="(I4,1X,I4,"
!          name = TRIM(atm_name)
         WRITE(unit_nr,TRIM(my_format)//"1X,2F20.10)") iadiab, iorbitals, &
                   REAL(coefficients((iadiab-1)*norbitals + iorbitals)),   &
                   AIMAG(coefficients((iadiab-1)*norbitals + iorbitals))
         ENDDO
       ENDDO

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%COEFFICIENTS")
   END SUBROUTINE adiab_print_coeff


   SUBROUTINE adiab_print_pop(adiab_section, populations, nadiab, norbitals)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:),POINTER :: populations
    INTEGER                                  :: nadiab, norbitals
!    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_pop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
!    INTEGER                                  :: i
    INTEGER                                  :: iadiab, iorbitals
!    TYPE(particle_type), POINTER             :: particle_set(:)
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    CHARACTER(LEN=4)                         :: name
    CHARACTER(LEN=default_string_length)     :: my_format
    TYPE(cp_logger_type), POINTER            :: logger


    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%POPULATIONS",&
               extension=".xyz", middle_name="pop",log_filename=.FALSE.)

    IF (unit_nr>0) THEN

!       ncoeff = SIZE(coefficients)

       WRITE ( unit_nr, "(2I8)" ) nadiab, norbitals
       WRITE ( unit_nr, '(A)' ) 'Diab. State, Orb., Pop. '

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
!          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
!                  name=atm_name)
!          my_format="(A4,"
         my_format="(I4,1X,I4,"
!          name = TRIM(atm_name)
         WRITE(unit_nr,TRIM(my_format)//"1X,F20.10)") iadiab, iorbitals, &
                   populations((iadiab-1)*norbitals + iorbitals)
         ENDDO
       ENDDO

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%POPULATIONS")
   END SUBROUTINE adiab_print_pop


   SUBROUTINE adiab_print_hamiltonian(adiab_section, hamiltonian, nadiab, &
              norbitals)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    INTEGER                                  :: nadiab, norbitals
!    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_hamiltonian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
!    INTEGER                                  :: i
    INTEGER                                  :: iadiab, iorbitals
 !   INTEGER                                  :: jorbitals
!    TYPE(particle_type), POINTER             :: particle_set(:)
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    CHARACTER(LEN=default_string_length)     ::  my_format, my_format2
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%HAMILTONIAN",&
               extension=".xyz",middle_name="hamilt", log_filename=.FALSE.,&
               is_new_file=new_file)


    IF (unit_nr>0) THEN

!       ncoeff = SIZE(coefficients)

       WRITE ( unit_nr, "(2I8)" ) nadiab, norbitals
       WRITE ( unit_nr, '(A)' ) 'Hamiltonian: (nadiab*norbitals)^2'

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
!            DO jadiab=1, nadiab
!               DO jorbitals=1,norbitals 
!          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
!                  name=atm_name)
!          my_format="(A4,"
         my_format="(I4,1X,I4,"
         write(my_format2,"(I8)") nadiab*norbitals
!          name = TRIM(atm_name)
         WRITE(unit_nr,TRIM(my_format)//"1X,"//TRIM(my_format2)//"F20.10)") iadiab, iorbitals, &
                           hamiltonian( (iadiab-1)*norbitals + iorbitals,:) 
print*, hamiltonian( (iadiab-1)*norbitals + iorbitals,:), "bbbb"
!              ENDDO
!            ENDDO
          ENDDO
       ENDDO

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%HAMILTONIAN")
   END SUBROUTINE adiab_print_hamiltonian


   SUBROUTINE adiab_print_matrix_density(adiab_section, coefficients, nadiab, &
              norbitals)
    TYPE(section_vals_type), POINTER         :: adiab_section
    COMPLEX(KIND=dp), DIMENSION(:),  POINTER :: coefficients
    !REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    INTEGER                                  :: nadiab, norbitals
!    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_matrix_density', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:,:),  POINTER :: matrix
    INTEGER                                  :: unit_nr
!    INTEGER                                  :: i
    INTEGER                                  :: iadiab, iorbitals
    INTEGER                                  :: jadiab, jorbitals
!    TYPE(particle_type), POINTER             :: particle_set(:)
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    CHARACTER(LEN=default_string_length)     ::  my_format, my_format2
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file

    ALLOCATE(matrix(nadiab*norbitals, nadiab*norbitals))

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%DENSITY_MATRIX",&
               extension=".xyz",middle_name="mat", log_filename=.FALSE.,&
               is_new_file=new_file)


    IF (unit_nr>0) THEN

!       ncoeff = SIZE(coefficients)

       WRITE ( unit_nr, "(2I8)" ) nadiab, norbitals
       WRITE ( unit_nr, '(A)' ) 'Matrix Density: (nadiab*norbitals)^2'

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
             DO jadiab=1, nadiab
                DO jorbitals=1, norbitals
                   matrix( (iadiab-1)*norbitals + iorbitals, (jadiab-1)*norbitals + jorbitals ) = &
                   coefficients( (iadiab-1)*norbitals + iorbitals)* &
                   CONJG(coefficients( (jadiab-1)*norbitals + jorbitals)) 
 print*, "hwre", coefficients( (iadiab-1)*norbitals + iorbitals)
 print*, "here", matrix( (iadiab-1)*norbitals + iorbitals, (jadiab-1)*norbitals + jorbitals ) 
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
!            DO jadiab=1, nadiab
!               DO jorbitals=1,norbitals 
!          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
!                  name=atm_name)
!          my_format="(A4,"
         my_format="(I4,1X,I4,"
         write(my_format2,"(I8)") 2*nadiab*norbitals
!          name = TRIM(atm_name)
         WRITE(unit_nr,TRIM(my_format)//"1X,"//TRIM(my_format2)//"F20.10)") iadiab, iorbitals, &
                           matrix( (iadiab-1)*norbitals + iorbitals,:)
!print*, hamiltonian( (iadiab-1)*norbitals + iorbitals,:), "bbbb"
!              ENDDO
!            ENDDO
          ENDDO
       ENDDO

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%HAMILTONIAN")
   END SUBROUTINE adiab_print_matrix_density


   SUBROUTINE adiab_print_nacv(adiab_section, nacv, nadiab, &
              norbitals,ext_atoms)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:,:,:),POINTER :: nacv
    INTEGER                                  :: nadiab, norbitals, ext_atoms

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_nacv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    INTEGER                                  :: iatom
    INTEGER                                  :: iadiab, iorbitals
    INTEGER                                  :: jadiab, jorbitals
    CHARACTER(LEN=default_string_length)     ::  my_format
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file


    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%NACV",&
               extension=".xyz",middle_name="nacv", log_filename=.FALSE.,&
               is_new_file=new_file)


    IF (unit_nr>0) THEN

       WRITE ( unit_nr, "(2I8)" ) ext_atoms, nadiab, norbitals
       WRITE ( unit_nr, '(A)' ) 'Atoms, Adiab., Orb.,NACV'

       DO iatom=1,ext_atoms
          DO iadiab=1, nadiab
             DO iorbitals=1,norbitals
                DO jadiab=1, nadiab
                   DO jorbitals=1, norbitals
                      my_format="(I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,3F20.10)"
                      WRITE(unit_nr,TRIM(my_format)) &
                             iatom, iadiab, iorbitals, jadiab, jorbitals,   &
                             nacv(:, (iadiab-1)*norbitals+iorbitals,        & 
                             (jadiab-1)*norbitals+jorbitals, iatom)
                   ENDDO
                ENDDO 
             ENDDO
          ENDDO
       ENDDO

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%HAMILTONIAN")
   END SUBROUTINE adiab_print_nacv


   SUBROUTINE adiab_print_stateno(adiab_section, state, ext_atoms)
    TYPE(section_vals_type), POINTER         :: adiab_section
    INTEGER                                  :: state
    INTEGER                                  :: ext_atoms

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_stateno', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file


    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%STATE",&
               extension=".xyz",middle_name="state", log_filename=.FALSE.,&
               is_new_file=new_file)

    IF (unit_nr>0) THEN

       WRITE ( unit_nr, "(2I8)" ) ext_atoms
       WRITE ( unit_nr, '(A)' ) 'State'
       WRITE ( unit_nr, "(I8)"  )  state

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%HAMILTONIAN")
   END SUBROUTINE adiab_print_stateno


END MODULE adiabatic_methods

