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

  USE jacob_types,                     ONLY: aom_basket_type, &
                                             sh_print_basket_type
!  USE timings,                         ONLY: timeset,&
!                                             timestop
  USE cp_units,                        ONLY: cp_unit_from_cp2k


#include "./base/base_uses.f90"

IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'jacob_print'
   
   PUBLIC :: adiab_print_coeff,          &
             adiab_print_pop,            &
             adiab_print_adiabat_energies,&
             adiab_print_matrix_density, &
             adiab_print_hamiltonian,    &
             adiab_print_nace,           &
             adiab_print_nacv,           &
             adiab_print_stateno,        &
             adiab_print_exact_forces,   &
             adiab_print_active_state,   &
             adiab_print_sh_info

CONTAINS


   SUBROUTINE adiab_print_coeff(adiab_section, coefficients, nadiab, norbitals,itimes, deltatMD, times)
    TYPE(section_vals_type), POINTER           :: adiab_section
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN) :: coefficients
    INTEGER, INTENT(IN)                        :: nadiab, norbitals, itimes
    REAL(KIND=dp), INTENT(IN)                  :: times, deltatMD

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: unit_conv_fs, time
    INTEGER                                  :: iadiab, iorbitals
    CHARACTER(LEN=default_string_length)     :: my_format
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%COEFFICIENTS",&
               extension=".xyz", middle_name="coeff", log_filename=.FALSE.)
    unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")
   
    IF (unit_nr>0) THEN

       time = times*deltatMD*unit_conv_fs
       WRITE ( unit_nr, "( A,A,I8,A,I8)" ) "Psi", "; nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",itimes,", time = ",time, " (fs)"

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
         my_format="(I4,1X,I4,"
         WRITE(unit_nr,TRIM(my_format)//"1X,2F20.10)") iadiab, iorbitals, &
                   REAL(coefficients((iadiab-1)*norbitals + iorbitals)),   &
                   AIMAG(coefficients((iadiab-1)*norbitals + iorbitals))
         ENDDO
       ENDDO

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%COEFFICIENTS")
   END SUBROUTINE adiab_print_coeff




   SUBROUTINE adiab_print_exact_forces(adiab_section, exact_forces,natom, deltatMD, times)
    TYPE(section_vals_type), POINTER          :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:),INTENT(IN ) :: exact_forces
    INTEGER, INTENT(IN)                       :: natom
    REAL(KIND=dp), INTENT(IN)                 :: deltatMD, times

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_exact_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr, iatom
    REAL(KIND=dp)                            :: unit_conv_fs, time
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%EXACT_FORCES",&
               extension=".xyz", middle_name="exfrc",log_filename=.FALSE.)
    unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")

    IF (unit_nr>0) THEN

       time = times*deltatMD*unit_conv_fs
       WRITE ( unit_nr, "(I8)" )          natom
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",INT(times),", time = ",time, " (fs)" 
       DO iatom=1,natom 
          WRITE(unit_nr,"(3F20.10)")  exact_forces(:,iatom)
       ENDDO

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%EXACT_FORCES")
   END SUBROUTINE adiab_print_exact_forces



   SUBROUTINE adiab_print_adiabat_energies(adiab_section, energies, nadiab, deltatMD, times)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: energies
    INTEGER, INTENT(IN)                      :: nadiab
    REAL(KIND=dp), INTENT(IN)                :: times
    REAL(KIND=dp), INTENT(IN)                :: deltatMD

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_pop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: unit_conv_fs, unit_conv, time
    INTEGER                                  :: iadiab
    CHARACTER(LEN=default_string_length)     :: my_format
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(LEN=default_string_length)     :: unit_str 

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%ADIABAT_ENERGIES",&
               extension=".xyz", middle_name="adiab",log_filename=.FALSE.)
    unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")

    IF (unit_nr>0) THEN
       CALL section_vals_val_get(adiab_section,"PRINT%ADIABAT_ENERGIES%UNIT",&
            c_val=unit_str)
       unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))
       unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")
       time = times*deltatMD*unit_conv_fs

       WRITE ( unit_nr, "( A,A,A,I8,A,I8)" ) "Adiab. Energies (",TRIM(unit_str),");  nadiab =", nadiab
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",INT(times),", time = ",time, " (fs)"

       DO iadiab=1, nadiab
         my_format="(I4,"
         WRITE(unit_nr,TRIM(my_format)//"1X,F20.10)") iadiab, energies(iadiab) 
       ENDDO

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%ADIABAT_ENERGIES")
   END SUBROUTINE adiab_print_adiabat_energies




   SUBROUTINE adiab_print_pop(adiab_section, populations, nadiab, norbitals,itimes, deltatMD, times)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: populations
    INTEGER, INTENT(IN)                      :: nadiab, norbitals, itimes
    REAL(KIND=dp),INTENT(IN)                 :: deltatMD, times

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_pop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: unit_conv_fs, time
    INTEGER                                  :: iadiab, iorbitals
    CHARACTER(LEN=default_string_length)     :: my_format
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%POPULATIONS",&
               extension=".xyz", middle_name="pop",log_filename=.FALSE.)
    unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")

    IF (unit_nr>0) THEN

       time = unit_conv_fs*times*deltatMD
       WRITE ( unit_nr, "( A,A,I8,A,I8)" ) "Pop.", "; nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",itimes,", time = ",time, " (fs)" 

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
         my_format="(I4,1X,I4,"
         WRITE(unit_nr,TRIM(my_format)//"1X,F20.10)") iadiab, iorbitals, &
                   populations((iadiab-1)*norbitals + iorbitals)
         ENDDO
       ENDDO

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%POPULATIONS")
   END SUBROUTINE adiab_print_pop



   SUBROUTINE adiab_print_active_state(adiab_section, state, energy, vector, & 
              nadiab, norbitals, itimes, deltatMD)
    TYPE(section_vals_type), POINTER         :: adiab_section
    INTEGER, INTENT(IN)                      :: state
    REAL(KIND=dp), INTENT(IN)                :: energy
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector
    INTEGER, INTENT(IN)                      :: nadiab, norbitals
    INTEGER, POINTER                         :: itimes
    REAL(KIND=dp),INTENT(IN)                 :: deltatMD

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_active_state', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    INTEGER                                  :: iadiab, iorbitals
    CHARACTER(LEN=default_string_length)     :: my_format
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file
    REAL(KIND=dp)                            :: unit_conv, unit_conv_fs, time
    CHARACTER(LEN=default_string_length)     :: unit_str

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%ACTIVE_STATE_INFO",&
               extension=".xyz",middle_name="active", log_filename=.FALSE.,&
               is_new_file=new_file)

    IF (unit_nr>0) THEN

       CALL section_vals_val_get(adiab_section,"PRINT%ACTIVE_STATE_INFO%UNIT",&
            c_val=unit_str)
       unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))
       unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")
       time = unit_conv_fs*itimes*deltatMD

       WRITE ( unit_nr, "(A,I8,A,I8)" ) "Active State Info, nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",itimes,", time = ",time, " (fs)"
       WRITE ( unit_nr, '(A,I8)' )   " State = ", state
       WRITE ( unit_nr, '(A,A,A,F20.10)' )   " Energy (",TRIM(unit_str),") = ", energy
       WRITE ( unit_nr, '(A)' )   " Eigenvector = "
       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
         my_format="(I4,1X,I4,"
         WRITE(unit_nr,TRIM(my_format)//"1X,F20.10)") iadiab, iorbitals, &
                   vector((iadiab-1)*norbitals + iorbitals)
         ENDDO
       ENDDO
       WRITE ( unit_nr, '(A)' )   "  "
       WRITE ( unit_nr, '(A)' )   "  "

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%ACTIVE_STATE_INFO")
   END SUBROUTINE adiab_print_active_state




   SUBROUTINE adiab_print_hamiltonian(adiab_section, hamiltonian, nadiab, &
              norbitals, itimes, deltatMD)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: hamiltonian
    INTEGER, INTENT(IN)                      :: nadiab, norbitals
    INTEGER, POINTER                         :: itimes
    REAL(KIND=dp),INTENT(IN)                 :: deltatMD

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_hamiltonian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    INTEGER                                  :: iadiab, iorbitals
    CHARACTER(LEN=default_string_length)     ::  my_format, my_format2
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file
    REAL(KIND=dp)                            :: unit_conv, unit_conv_fs, time
    CHARACTER(LEN=default_string_length)     :: unit_str 

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%HAMILTONIAN",&
               extension=".xyz",middle_name="hamilt", log_filename=.FALSE.,&
               is_new_file=new_file)

    IF (unit_nr>0) THEN

       CALL section_vals_val_get(adiab_section,"PRINT%HAMILTONIAN%UNIT",&
            c_val=unit_str)
       unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))
       unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")
       time = unit_conv_fs*itimes*deltatMD

       WRITE ( unit_nr, "( A,A,A,I8,A,I8)" ) "H (",TRIM(unit_str),");  nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",itimes,", time = ",time, " (fs)" 

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
         my_format="(I4,1X,I4,"
         write(my_format2,"(I8)") nadiab*norbitals
         WRITE(unit_nr,TRIM(my_format)//"1X,"//TRIM(my_format2)//"F20.10)") iadiab, iorbitals, &
                           unit_conv*hamiltonian( (iadiab-1)*norbitals + iorbitals,:) 
          ENDDO
       ENDDO

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%HAMILTONIAN")
   END SUBROUTINE adiab_print_hamiltonian



   SUBROUTINE adiab_print_nace(adiab_section, nace, nadiab, &
              norbitals, itimes, deltatMD)
    TYPE(section_vals_type), POINTER         :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:),POINTER    :: nace
    INTEGER, INTENT(IN)                      :: nadiab, norbitals
    INTEGER, POINTER                         :: itimes
    REAL(KIND=dp),INTENT(IN)                 :: deltatMD

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_nace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    INTEGER                                  :: iadiab, iorbitals
    CHARACTER(LEN=default_string_length)     ::  my_format, my_format2
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file
    REAL(KIND=dp)                            :: unit_conv, unit_conv_fs, time
    CHARACTER(LEN=default_string_length)     :: unit_str

    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%NACE",&
               extension=".xyz",middle_name="nace", log_filename=.FALSE.,&
               is_new_file=new_file)

    IF (unit_nr>0) THEN

       CALL section_vals_val_get(adiab_section,"PRINT%NACE%UNIT",&
            c_val=unit_str)
       unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))
       unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")
       time = unit_conv_fs*itimes*deltatMD

       WRITE ( unit_nr, "( A,A,A,I8,A,I8)" ) "NACE (",TRIM(unit_str),");  nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",itimes,", time = ",time, " (fs)"

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
         my_format="(I4,1X,I4,"
         write(my_format2,"(I8)") nadiab*norbitals
         WRITE(unit_nr,TRIM(my_format)//"1X,"//TRIM(my_format2)//"F20.10)") iadiab, iorbitals, &
                           unit_conv*nace( (iadiab-1)*norbitals + iorbitals,:)
          ENDDO
       ENDDO

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%NACE")
   END SUBROUTINE adiab_print_nace




   SUBROUTINE adiab_print_matrix_density(adiab_section, coefficients, nadiab, &
              norbitals, itimes, deltatMD, times)
    TYPE(section_vals_type), POINTER           :: adiab_section
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN) :: coefficients
    INTEGER, INTENT(IN)                        :: nadiab, norbitals, itimes
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD, times

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_matrix_density', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:,:),  POINTER :: matrix
    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: unit_conv_fs, time
    INTEGER                                  :: iadiab, iorbitals
    INTEGER                                  :: jadiab, jorbitals
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

       unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")
       time = unit_conv_fs*times*deltatMD
       WRITE ( unit_nr, "( A,A,I8,A,I8)" ) "Mat. Dens.", "; nadiab =", nadiab, "; norbitals =",&
                                   norbitals
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",itimes,", time = ",time, " (fs)"

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
             DO jadiab=1, nadiab
                DO jorbitals=1, norbitals
                   matrix( (iadiab-1)*norbitals + iorbitals, (jadiab-1)*norbitals + jorbitals ) = &
                   coefficients( (iadiab-1)*norbitals + iorbitals)* &
                   CONJG(coefficients( (jadiab-1)*norbitals + jorbitals)) 
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DO iadiab=1, nadiab
          DO iorbitals=1,norbitals
         my_format="(I4,1X,I4,"
         write(my_format2,"(I8)") 2*nadiab*norbitals
         WRITE(unit_nr,TRIM(my_format)//"1X,"//TRIM(my_format2)//"F20.10)") iadiab, iorbitals, &
                           matrix( (iadiab-1)*norbitals + iorbitals,:)
          ENDDO
       ENDDO

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%HAMILTONIAN")
   END SUBROUTINE adiab_print_matrix_density




   SUBROUTINE adiab_print_nacv(adiab_section, nacv, nadiab, &
              norbitals,ext_atoms, itimes, deltatMD)
    TYPE(section_vals_type), POINTER          :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:,:,:),POINTER :: nacv
    INTEGER, INTENT(IN)                       :: nadiab, norbitals, ext_atoms
    INTEGER, POINTER                          :: itimes
    REAL(KIND=dp), INTENT(IN)                 :: deltatMD

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_nacv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file
    INTEGER                                  :: i,j,k
    REAL(KIND=dp)                            :: unit_conv, time, unit_conv_fs
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
       unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")
       time = unit_conv_fs*itimes*deltatMD

        WRITE ( unit_nr, "( A,A,A,I8,A,I8)" ) "NACV (",TRIM(unit_str),"); nadiab =", nadiab,"; norbitals =",&
                                           norbitals
        WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",itimes,", time = ",time, " (fs)"
        DO i=1,nadiab*norbitals*ext_atoms
           DO j=1,3
              WRITE( unit_nr, *) "Atom = ",i,"Dimension = ", j
              DO k=1,norbitals*nadiab
                 WRITE( unit_nr, *) nacv(j, k, 1:(norbitals*nadiab), i)
              ENDDO
           ENDDO
        ENDDO
           
    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%NACV")
   END SUBROUTINE adiab_print_nacv




   SUBROUTINE adiab_print_stateno(adiab_section, state, ext_atoms,itimes, deltatMD)
    TYPE(section_vals_type), POINTER         :: adiab_section
    INTEGER, INTENT(IN)                      :: state
    INTEGER, INTENT(IN)                      :: ext_atoms
    INTEGER, POINTER                         :: itimes
    REAL(KIND=dp), INTENT(IN)                :: deltatMD

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_stateno', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: unit_conv_fs, time
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: new_file


    NULLIFY(logger)
    logger => cp_get_default_logger()

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%STATE",&
               extension=".xyz",middle_name="state", log_filename=.FALSE.,&
               is_new_file=new_file)
    unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")

    IF (unit_nr>0) THEN

       time = itimes*deltatMD*unit_conv_fs
       WRITE ( unit_nr, '(A)' ) 'Present Adiabatic State'
       WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",itimes,", time = ",time, " (fs)" 
       WRITE ( unit_nr, '(A,I8)' )   " State = ", state

    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%STATE")
   END SUBROUTINE adiab_print_stateno


   SUBROUTINE adiab_print_sh_info(adiab_section, sh_print_basket, nadiab, norbitals, deltatMD, times)
    TYPE(section_vals_type), POINTER         :: adiab_section
    TYPE(sh_print_basket_type), INTENT(IN)   :: sh_print_basket
    INTEGER, INTENT(IN)                      :: nadiab, norbitals
    REAL(KIND=dp), INTENT(IN)                :: times
    REAL(KIND=dp), INTENT(IN)                :: deltatMD

    CHARACTER(len=*), PARAMETER :: routineN = 'adiab_print_pop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    !REAL(KIND=dp)                            :: unit_conv_fs, unit_conv, time
    REAL(KIND=dp)                            :: unit_conv_fs, time
    INTEGER                                  :: iadiab
    CHARACTER(LEN=default_string_length)     :: my_format
    TYPE(cp_logger_type), POINTER            :: logger

     REAL(KIND = dp), DIMENSION(:), POINTER      :: proba, pop_init
     REAL(KIND = dp)                             :: hop
     LOGICAL                                     :: renorm
     REAL(KIND = dp)                             :: sumg
     REAL(KIND = dp)                             :: factor
     INTEGER                                     :: old_state
     INTEGER                                     :: new_state, first_adiabat
     REAL(KIND = dp)                             :: delta_energy, initrandom
     LOGICAL                                     :: is_hop
     LOGICAL                                     :: decohere, is_first
   
    NULLIFY(logger)
    logger => cp_get_default_logger()

    ALLOCATE(proba(nadiab))
    ALLOCATE(pop_init(nadiab))

    unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%SH_INFO",&
               extension=".log", middle_name="sh",log_filename=.FALSE.)
    unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")

    IF (unit_nr>0) THEN
       !CALL section_vals_val_get(adiab_section,"PRINT%SH_INFO%UNIT",&
       !     c_val=unit_str)
       !unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))
       unit_conv_fs = cp_unit_from_cp2k(1.0_dp,"fs")
       time = times*deltatMD*unit_conv_fs

       sumg         = sh_print_basket%sumprob
       renorm       = sh_print_basket%renorm
       proba        = sh_print_basket%proba
       hop          = sh_print_basket%hop
       is_hop       = sh_print_basket%is_hop
       old_state    = sh_print_basket%old_state
       new_state    = sh_print_basket%new_state
       delta_energy = sh_print_basket%delta_energy
       factor       = sh_print_basket%factor
       decohere     = sh_print_basket%decohere
       is_first     = sh_print_basket%is_first
       initrandom   = sh_print_basket%initrandom
       pop_init     = sh_print_basket%pop_init
       first_adiabat= sh_print_basket%first_adiabat    

       IF (is_first) THEN
          WRITE( unit_nr, "(A)")   "LET'S START. "
          WRITE( unit_nr, "(A)")   "Adiabatic Initial Population = "
          DO iadiab=1, nadiab
             my_format="(I4,"
             WRITE(unit_nr,TRIM(my_format)//"1X,F20.10)") iadiab, pop_init(iadiab)
          ENDDO
          WRITE( unit_nr, "(A,F20.10)" ) "Initial Random Number = ", initrandom
          WRITE( unit_nr, "(A,I8)" ) "First Adiabatic State = ", first_adiabat
          WRITE( unit_nr, "(A)") "  "
       ELSE        
          WRITE ( unit_nr, "( A,I8,A,I8)" ) "Surface Hopping Info,  nadiab =", nadiab, &
                                                " , norbitals = ", norbitals
          WRITE ( unit_nr, '(A,I0,A,F20.3,A,F0.10)' )   " i = ",INT(times),", time = ",time, " (fs)"
          WRITE( unit_nr, "(A,F20.10)" ) "Sum of probability = ", sumg

          IF (renorm)  WRITE( unit_nr, "(A)" ) "RE-NORMALIZATION!"

          WRITE( unit_nr, "(A)" ) "Probability = "
          DO iadiab=1, nadiab
             my_format="(I4,"
             WRITE(unit_nr,TRIM(my_format)//"1X,F20.10)") iadiab, proba(iadiab)
          ENDDO
          WRITE( unit_nr, "(A,F20.10)" ) "Random number = ", hop
          IF (is_hop) THEN
             WRITE( unit_nr, "(A,I8,A,I8)" ) "ATTEMPT between old_state = ", old_state, " and new_state = ", new_state
             WRITE( unit_nr, "(A,F20.10)" ) "Energy Differences = ", delta_energy
             IF (delta_energy.ge.0) THEN
                WRITE( unit_nr, "(A)") "SUCCESSFUL HOP!"
                WRITE( unit_nr, "(A,F20.10)" ) "Rescaling factor = ", factor
             ELSE
                WRITE( unit_nr, "(A)") "FRUSTRATED HOP: NO HOP!"
             ENDIF
          ELSE 
             WRITE( unit_nr, "(A)") "NO HOP!"
          ENDIF   
          IF (decohere) WRITE( unit_nr, "(A)") "DECOHERENCE"       
          WRITE( unit_nr, "(A)")     "  "
          WRITE ( unit_nr, '(A)' )   "  "
       ENDIF   
    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,adiab_section,"PRINT%SH_INFO")
    DEALLOCATE(proba)
   END SUBROUTINE adiab_print_sh_info


END MODULE jacob_print

