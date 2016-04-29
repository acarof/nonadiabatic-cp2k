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


MODULE jacob_print

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
  USE cp_units,                        ONLY: cp_unit_from_cp2k


#include "./base/base_uses.f90"

IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'jacob_print'
   
   PUBLIC :: adiab_print_coeff,          &
             adiab_print_pop,            &
             adiab_print_matrix_density, &
             adiab_print_hamiltonian,    &
             adiab_print_nacv,           &
             adiab_print_stateno,        &
             adiab_print_exact_forces

CONTAINS


   SUBROUTINE adiab_print_coeff(adiab_section, coefficients, nadiab, norbitals,times)
    TYPE(section_vals_type), POINTER         :: adiab_section
    COMPLEX(KIND=dp), DIMENSION(:)           :: coefficients
    INTEGER                                  :: nadiab, norbitals
    REAL(KIND=dp)                            :: times
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

       WRITE ( unit_nr, "( A,A,I8,A,I8)" ) "Psi", "; nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,F10.1)' )   " i = ",times

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


   SUBROUTINE adiab_print_exact_forces(adiab_section, exact_forces,natom,times)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:),INTENT(OUT) :: exact_forces
    INTEGER, INTENT(OUT)                      :: natom
    REAL(KIND=dp), INTENT(OUT)                :: times
!    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_exact_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr, iatom
!    INTEGER                                  :: i
!    INTEGER                                  :: iadiab, iorbitals
!    TYPE(particle_type), POINTER             :: particle_set(:)
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    CHARACTER(LEN=4)                         :: name
!    CHARACTER(LEN=default_string_length)     :: my_format
    TYPE(cp_logger_type), POINTER            :: logger


    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%EXACT_FORCES",&
               extension=".xyz", middle_name="exfrc",log_filename=.FALSE.)

    IF (unit_nr>0) THEN

!       ncoeff = SIZE(coefficients)

       WRITE ( unit_nr, "(I8)" )          natom
       WRITE ( unit_nr, '(A,F10.1)' )   " i = ",times
       DO iatom=1,natom 
          WRITE(unit_nr,"(3F20.10)")  exact_forces(:,iatom)
       ENDDO

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%EXACT_FORCES")
   END SUBROUTINE adiab_print_exact_forces



   SUBROUTINE adiab_print_pop(adiab_section, populations, nadiab, norbitals,times)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:)              :: populations
    INTEGER                                  :: nadiab, norbitals
    REAL(KIND=dp)                            :: times
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

       WRITE ( unit_nr, "( A,A,I8,A,I8)" ) "Pop.", "; nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,F10.1)' )   " i = ",times

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
              norbitals, itimes)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    INTEGER                                  :: nadiab, norbitals
    INTEGER, POINTER                         :: itimes
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
    REAL(KIND=dp)                            :: unit_conv
    CHARACTER(LEN=default_string_length)     :: unit_str 

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%HAMILTONIAN",&
               extension=".xyz",middle_name="hamilt", log_filename=.FALSE.,&
               is_new_file=new_file)


    IF (unit_nr>0) THEN

!       ncoeff = SIZE(coefficients)
       CALL section_vals_val_get(adiab_section,"PRINT%HAMILTONIAN%UNIT",&
            c_val=unit_str)
       unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))

       WRITE ( unit_nr, "( A,A,A,I8,A,I8)" ) "H (",TRIM(unit_str),");  nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,I8)' )   " i = ",itimes

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
                           unit_conv*hamiltonian( (iadiab-1)*norbitals + iorbitals,:) 
!print*, hamiltonian( (iadiab-1)*norbitals + iorbitals,:), "bbbb"
!              ENDDO
!            ENDDO
          ENDDO
       ENDDO

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%HAMILTONIAN")
   END SUBROUTINE adiab_print_hamiltonian


   SUBROUTINE adiab_print_matrix_density(adiab_section, coefficients, nadiab, &
              norbitals, times)
    TYPE(section_vals_type), POINTER         :: adiab_section
    COMPLEX(KIND=dp), DIMENSION(:)           :: coefficients
    !REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    INTEGER                                  :: nadiab, norbitals
    REAL(KIND=dp)                            :: times
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

       WRITE ( unit_nr, "( A,A,I8,A,I8)" ) "Mat. Dens.", "; nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,F10.1)' )   " i = ",times

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
             DO jadiab=1, nadiab
                DO jorbitals=1, norbitals
                   matrix( (iadiab-1)*norbitals + iorbitals, (jadiab-1)*norbitals + jorbitals ) = &
                   coefficients( (iadiab-1)*norbitals + iorbitals)* &
                   CONJG(coefficients( (jadiab-1)*norbitals + jorbitals)) 
! print*, "hwre", coefficients( (iadiab-1)*norbitals + iorbitals)
! print*, "here", matrix( (iadiab-1)*norbitals + iorbitals, (jadiab-1)*norbitals + jorbitals ) 
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
              norbitals,ext_atoms, itimes)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:,:,:),POINTER :: nacv
    INTEGER                                  :: nadiab, norbitals, ext_atoms
    INTEGER, POINTER                         :: itimes

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_nacv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    !INTEGER                                  :: iatom
    !INTEGER                                  :: iadiab, iorbitals
    !INTEGER                                  :: jadiab, jorbitals
    !CHARACTER(LEN=default_string_length)     ::  my_format
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file
    INTEGER                                  :: i,j,k
    REAL(KIND=dp)                            :: unit_conv
    CHARACTER(LEN=default_string_length)     :: unit_str


    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%NACV",&
               extension=".xyz",middle_name="nacv", log_filename=.FALSE.,&
               is_new_file=new_file)


    IF (unit_nr>0) THEN

       CALL section_vals_val_get(adiab_section,"PRINT%NACV%UNIT",&
            c_val=unit_str)
       unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))

        WRITE ( unit_nr, "( A,A,A,I8,A,I8)" ) "NACV (",TRIM(unit_str),"); nadiab =", nadiab,"; norbitals =",&
                                           norbitals
        WRITE ( unit_nr, '(A,I8)' )   " i = ",itimes
        DO i=1,nadiab*norbitals*ext_atoms
           DO j=1,3
              WRITE( unit_nr, *) "Atom = ",i,"Dimension = ", j
              DO k=1,norbitals*nadiab
                 WRITE( unit_nr, *) nacv(j, k, 1:(norbitals*nadiab), i)
              ENDDO
           ENDDO
        ENDDO
           
!       WRITE ( unit_nr, "(2I8)" ) ext_atoms, nadiab, norbitals
!       WRITE ( unit_nr, '(A)' ) 'Atoms, Adiab., Orb.,NACV'

!       DO iatom=1,ext_atoms
!          DO iadiab=1, nadiab
!             DO iorbitals=1,norbitals
!                DO jadiab=1, nadiab
!                   DO jorbitals=1, norbitals
!                      my_format="(I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,3F20.10)"
!                      WRITE(unit_nr,TRIM(my_format)) &
!                             iatom, iadiab, iorbitals, jadiab, jorbitals,   &
!                             nacv(:, (iadiab-1)*norbitals+iorbitals,        & 
!                             (jadiab-1)*norbitals+jorbitals, iatom)
!                   ENDDO
!                ENDDO 
!             ENDDO
!          ENDDO
!       ENDDO

    END IF
!    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%HAMILTONIAN")
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%NACV")
   END SUBROUTINE adiab_print_nacv


   SUBROUTINE adiab_print_stateno(adiab_section, state, ext_atoms,itimes)
    TYPE(section_vals_type), POINTER         :: adiab_section
    INTEGER                                  :: state
    INTEGER                                  :: ext_atoms
    INTEGER, POINTER                         :: itimes

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

       WRITE ( unit_nr, '(A)' ) 'Present Adiabatic State'
       WRITE ( unit_nr, '(A,I8)' )   " i = ",itimes
       WRITE ( unit_nr, '(A,I8)' )   " State = ", state

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%STATE")
   END SUBROUTINE adiab_print_stateno


END MODULE jacob_print

