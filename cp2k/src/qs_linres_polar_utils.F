!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Polarizability calculation by dfpt
!>      Initialization of the polar_env, 
!>      Perturbation Hamiltonian by application of the Berry phase operator to psi0
!>      Write output
!>      Deallocate everything
!> periodic Raman SL February 2013
!> \note
! *****************************************************************************
MODULE qs_linres_polar_utils
  USE bibliography,                    ONLY: Luber2014,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm,&
                                             cp_fm_trace
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_get_default_io_unit,&
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             put_results
  USE cp_result_types,                 ONLY: cp_result_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: twopi
  USE physcon,                         ONLY: angstrom
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_methods,               ONLY: linres_solver
  USE qs_linres_types,                 ONLY: get_polar_env,&
                                             linres_control_type,&
                                             polar_env_create,&
                                             polar_env_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: polar_env_cleanup, polar_env_init, polar_response, polar_polar, polar_print

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_polar_utils'

CONTAINS

! *****************************************************************************
!> \brief Initialize the polar environment
!> \param polar_env ...
!> \param p_env ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE polar_response(polar_env,p_env,qs_env)
   
    TYPE(polar_env_type)                     :: polar_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_response', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, ispin, nao, &
                                                nmo, nspins, output_unit
    LOGICAL                                  :: do_raman, should_stop
    REAL(dp)                                 :: chk, fro
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: h1_psi0, psi0_order, psi1
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: dBerry_psi0, psi1_dBerry
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(section_vals_type), POINTER         :: lr_section, polar_section

    CALL timeset(routineN,handle)
    !
    NULLIFY(dft_control, linres_control, lr_section, polar_section)
    NULLIFY(logger, mpools, psi1, h1_psi0, mo_coeff, para_env)
    NULLIFY(tmp_fm_struct, psi1_dBerry,dBerry_psi0)

    logger => cp_get_default_logger()
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES")
    polar_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%POLAR")

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog")
    IF(output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "*** Self consistent optimization of the response wavefunctions ***"
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    mpools=mpools,&
                    linres_control=linres_control,&
                    mos=mos,&
                    para_env=para_env)

    nspins = dft_control%nspins

    CALL get_polar_env(polar_env=polar_env, &
                      do_raman=do_raman)
    !
    ! allocate the vectors
    ALLOCATE(psi0_order(nspins))
    ALLOCATE(psi1(nspins),h1_psi0(nspins))
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       psi0_order(ispin)%matrix => mo_coeff
       CALL cp_fm_get_info(mo_coeff,ncol_global=nmo,nrow_global=nao)
       NULLIFY(tmp_fm_struct,psi1(ispin)%matrix,h1_psi0(ispin)%matrix)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                               ncol_global=nmo,&
                               context=mo_coeff%matrix_struct%context)
       CALL cp_fm_create(psi1(ispin)%matrix,tmp_fm_struct)
       CALL cp_fm_create(h1_psi0(ispin)%matrix,tmp_fm_struct)
       CALL cp_fm_struct_release(tmp_fm_struct)
    ENDDO
    chk = 0.0_dp
    !
    IF(do_raman) THEN
       CALL get_polar_env(polar_env=polar_env, &
                      psi1_dBerry=psi1_dBerry, &
                      dBerry_psi0=dBerry_psi0)
       DO idir = 1,3
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") &
               "Response to the perturbation operator Berry phase_"//ACHAR(idir+119)
          ENDIF
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1_dBerry(idir,ispin)%matrix,0.0_dp)
          ENDDO
          !
          !Initial guess for psi1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(dBerry_psi0(idir,ispin)%matrix,h1_psi0(ispin)%matrix)
          ENDDO
          !
          !
          linres_control%lr_triplet = .FALSE. ! we do singlet response
          linres_control%do_kernel  = .TRUE.  ! we do coupled response
          linres_control%converged  = .FALSE.
          CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order, should_stop)
          !
          !
          ! copy the response
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_dBerry(idir,ispin)%matrix)
             CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro)
             chk = chk + fro
          ENDDO
          !
         IF(output_unit>0) THEN
           WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
         ENDIF
          !
          ! Write the result in the restart file
       ENDDO ! idir
    ENDIF  !raman

    !
    ! print the checksum
 !   IF(output_unit>0) THEN
 !      WRITE(output_unit,'(T2,A,E23.16)') 'polar| response: CheckSum =',chk
 !   ENDIF
    !
    !
    ! clean up
    DO ispin = 1,nspins
       CALL cp_fm_release(psi1(ispin)%matrix)
       CALL cp_fm_release(h1_psi0(ispin)%matrix)
    ENDDO
    DEALLOCATE(psi1,h1_psi0,psi0_order)
    !
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
                                     "PRINT%PROGRAM_RUN_INFO")
    !
    CALL timestop(handle)
    !
  END SUBROUTINE polar_response


! *****************************************************************************
!> \brief ...
!> \param polar_env ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE polar_polar(polar_env,qs_env)

    TYPE(polar_env_type)                     :: polar_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_polar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, nspins, &
                                                output_unit, z
    LOGICAL                                  :: do_raman
    REAL(dp), DIMENSION(:, :), POINTER       :: polar, polar_tmp
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: dBerry_psi0, mo_derivs, &
                                                psi1_dBerry
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL timeset(routineN,handle)

    NULLIFY(cell,dft_control,polar,psi1_dBerry,logger)
    NULLIFY(mos, dBerry_psi0, mo_derivs)
    logger => cp_get_default_logger()
    output_unit= cp_logger_get_default_io_unit(logger)

    CALL get_qs_env(qs_env=qs_env,&
                   cell=cell,&
                   dft_control=dft_control,&
                   mos=mos)

    nspins = dft_control%nspins

    CALL get_polar_env(polar_env=polar_env,&
                      do_raman=do_raman)
!
    IF(do_raman) THEN

        CALL cite_reference(Luber2014)
    
        CALL get_polar_env(polar_env=polar_env,&
                      psi1_dBerry=psi1_dBerry, &
                      dBerry_psi0=dBerry_psi0, &
                      polar=polar, &
                      mo_derivs=mo_derivs)
    !
    !
!    polar_section => section_vals_get_subs_vals(qs_env%input, &
!         & "PROPERTIES%LINRES%POLAR")
    !
    ! Initialize
    ALLOCATE ( polar_tmp( 3, 3))
    polar_tmp(:,:)=0.0_dp

    DO i = 1,3 ! directions of electric field
     DO z=1,3 !dipole directions
       DO ispin=1,dft_control%nspins
!SL compute trace 
            CALL cp_fm_trace(psi1_dBerry(i,ispin)%matrix,mo_derivs(z,ispin)%matrix,polar_tmp(i,z))
         polar_tmp(i,z)=polar_tmp(i,z) + polar_tmp(i,z)
       END DO
          polar_tmp(i,z) = polar_tmp(i,z)/(twopi * twopi)
      END DO
     END DO !spin
   
     polar(:,:) = MATMUL (MATMUL (cell%hmat,polar_tmp),TRANSPOSE(cell%hmat)) 
!SL evtl maxocc instead?
      IF (dft_control%nspins .EQ. 1) THEN
        polar(:,:) = 2.0_dp * polar(:,:)
      END IF

    !
       IF(ASSOCIATED(polar_tmp)) THEN
         DEALLOCATE(polar_tmp)
       END IF
       
    ENDIF  !raman
    !
    ! Finalize
    CALL timestop(handle)
    !
  END SUBROUTINE polar_polar

! *****************************************************************************
!> \brief ...
!> \param polar_env ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE polar_print(polar_env,qs_env)
    TYPE(polar_env_type)                     :: polar_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'polar_print', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description
    INTEGER                                  :: output_unit, unit_p
    LOGICAL                                  :: do_raman
    REAL(dp), DIMENSION(:, :), POINTER       :: polar
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_result_type), POINTER            :: results
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: polar_section

    NULLIFY(logger,dft_control,para_env,results)
    CALL get_qs_env(qs_env=qs_env, results=results, dft_control=dft_control,&
                    para_env=para_env)

    logger => cp_get_default_logger()
    output_unit= cp_logger_get_default_io_unit(logger)

    polar_section => section_vals_get_subs_vals(qs_env%input, &
         "PROPERTIES%LINRES%POLAR")

    CALL get_polar_env(polar_env=polar_env, polar=polar, do_raman=do_raman)

    description = "[POLAR]"
    CALL cp_results_erase(results,description=description)
    CALL put_results(results,description=description,values=polar(:,:))

    IF(do_raman)  THEN
      IF(BTEST(cp_print_key_should_output(logger%iter_info,polar_section,&
             "PRINT%POLAR_MATRIX"),cp_p_file)) THEN

        unit_p=cp_print_key_unit_nr(logger,polar_section,"PRINT%POLAR_MATRIX",&
                  extension=".data",middle_name="raman",log_filename=.FALSE.)
        IF(unit_p>0) THEN 
          IF(unit_p/=output_unit) THEN 
              WRITE(unit_p,*)
              WRITE (unit_p,'(T10,A)') 'POLARIZABILITY TENSOR (atomic units):'
              WRITE (unit_p,'(T10,A,3F15.5)')"xx,yy,zz", polar(1,1),polar(2,2),polar(3,3)
              WRITE (unit_p,'(T10,A,3F15.5)')"xy,xz,yz",polar(1,2),polar(1,3),polar(2,3)
              WRITE (unit_p,'(T10,A,3F15.5)')"yx,zx,zy",polar(2,1),polar(3,1),polar(3,2)
              WRITE (unit_p,'(T10,A)') 'POLARIZABILITY TENSOR (Angstrom^3):'
              WRITE (unit_p,'(T10,A,3F15.5)')"xx,yy,zz", polar(1,1)*angstrom**3,&
                      polar(2,2)*angstrom**3,polar(3,3)*angstrom**3
              WRITE (unit_p,'(T10,A,3F15.5)')"xy,xz,yz",polar(1,2)*angstrom**3,& 
                      polar(1,3)*angstrom**3,polar(2,3)*angstrom**3
              WRITE (unit_p,'(T10,A,3F15.5)')"yx,zx,zy",polar(2,1)*angstrom**3,& 
                      polar(3,1)*angstrom**3,polar(3,2)*angstrom**3
            CALL cp_print_key_finished_output(unit_p,logger,polar_section, & 
                "PRINT%POLAR_MATRIX")
          ENDIF
        ENDIF
      ENDIF
      IF(output_unit >0) THEN 
        WRITE (output_unit,*)
        WRITE (output_unit,'(T10,A)') 'POLARIZABILITY TENSOR (atomic units):'
        WRITE (output_unit,'(T10,A,3F15.5)')"xx,yy,zz", polar(1,1),polar(2,2),polar(3,3)
        WRITE (output_unit,'(T10,A,3F15.5)')"xy,xz,yz",polar(1,2),polar(1,3),polar(2,3)
        WRITE (output_unit,'(T10,A,3F15.5)')"yx,zx,zy",polar(2,1),polar(3,1),polar(3,2)
        WRITE (output_unit,'(T10,A)') 'POLARIZABILITY TENSOR (Angstrom^3):'
        WRITE (output_unit,'(T10,A,3F15.5)')"xx,yy,zz", polar(1,1)*angstrom**3 ,&
                polar(2,2)*angstrom**3 ,polar(3,3)*angstrom**3 
        WRITE (output_unit,'(T10,A,3F15.5)')"xy,xz,yz",polar(1,2)*angstrom**3 ,&
                polar(1,3)*angstrom**3 ,polar(2,3 )*angstrom**3
        WRITE (output_unit,'(T10,A,3F15.5)')"yx,zx,zy",polar(2,1)*angstrom**3 ,&
                polar(3,1)*angstrom**3 ,polar(3,2 )*angstrom**3
      END IF
    END IF
    !
    !
  END SUBROUTINE polar_print


! *****************************************************************************
!> \brief Initialize the polar environment
!> \param polar_env ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE polar_env_init(polar_env,qs_env)
    !
    TYPE(polar_env_type)                     :: polar_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, ispin, m, nao, &
                                                nmo, nspins, output_unit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(section_vals_type), POINTER         :: lr_section, polar_section

!

    CALL timeset(routineN,handle)


    NULLIFY(linres_control)
    NULLIFY(logger,polar_section)

    logger => cp_get_default_logger()
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES")


    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                                     extension=".linresLog")

   IF(polar_env%ref_count /= 0) THEN
      CALL polar_env_cleanup(polar_env)
   ENDIF

   IF(output_unit>0) THEN
      WRITE(output_unit,"(/,T20,A,/)") "*** Start polarizability calculation ***"
      WRITE(output_unit,"(T10,A,/)") "Initialization of the polar environment"
   ENDIF

   CALL polar_env_create(polar_env)
  !
   polar_section => section_vals_get_subs_vals(qs_env%input, &
                  "PROPERTIES%LINRES%POLAR")

   CALL get_qs_env(qs_env=qs_env,&
                  dft_control=dft_control,&
                  matrix_s=matrix_s, &
                  linres_control=linres_control,&
                  mos=mos)
   !
   !
    nspins = dft_control%nspins
   !
    CALL section_vals_val_get(polar_section,"DO_RAMAN",l_val=polar_env%do_raman)
   !
   !
   !
   ALLOCATE(polar_env%polar(3,3))
    polar_env%polar(:,:) = 0.0_dp
  !
      ! allocation
   ALLOCATE(polar_env%dBerry_psi0(3,nspins), &
           polar_env%psi1_dBerry(3,nspins), polar_env%mo_derivs(3,nspins))
   DO ispin = 1,nspins
      CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
      CALL cp_fm_get_info(mo_coeff,ncol_global=m,nrow_global=nao)

      NULLIFY(tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                              ncol_global=m,&
                              context=mo_coeff%matrix_struct%context)
      DO idir = 1,3
         NULLIFY(polar_env%psi1_dBerry(idir,ispin)%matrix,polar_env%dBerry_psi0(idir,ispin)%matrix,&
                 polar_env%mo_derivs(idir,ispin)%matrix)
         CALL cp_fm_create(polar_env%psi1_dBerry(idir,ispin)%matrix,tmp_fm_struct)
         CALL cp_fm_create(polar_env%dBerry_psi0(idir,ispin)%matrix,tmp_fm_struct)
         CALL cp_fm_create(polar_env%mo_derivs(idir,ispin)%matrix,tmp_fm_struct)
      ENDDO
      CALL cp_fm_struct_release(tmp_fm_struct)

   !
   END DO

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
                                     "PRINT%PROGRAM_RUN_INFO")

    CALL timestop(handle)

  END SUBROUTINE polar_env_init

! *****************************************************************************
!> \brief Deallocate the polar environment
!> \param polar_env ...
!> \par History
! *****************************************************************************
  SUBROUTINE polar_env_cleanup(polar_env)

    TYPE(polar_env_type)                     :: polar_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idir, ispin

    polar_env%ref_count = polar_env%ref_count - 1
    IF(polar_env%ref_count == 0 ) THEN
       IF(ASSOCIATED(polar_env%polar)) THEN
          DEALLOCATE(polar_env%polar)
       ENDIF
       IF(ASSOCIATED(polar_env%dBerry_psi0)) THEN
          DO idir = 1,SIZE(polar_env%dBerry_psi0,1)
             DO ispin = 1,SIZE(polar_env%dBerry_psi0,2)
                CALL cp_fm_release(polar_env%dBerry_psi0(idir,ispin)%matrix)
             ENDDO
          ENDDO
          DEALLOCATE(polar_env%dBerry_psi0)
       ENDIF
       IF(ASSOCIATED(polar_env%mo_derivs)) THEN
          DO idir = 1,SIZE(polar_env%mo_derivs,1)
             DO ispin = 1,SIZE(polar_env%mo_derivs,2)
                CALL cp_fm_release(polar_env%mo_derivs(idir,ispin)%matrix)
             ENDDO
          ENDDO
          DEALLOCATE(polar_env%mo_derivs)
       ENDIF
       IF(ASSOCIATED(polar_env%psi1_dBerry)) THEN
          DO idir = 1,SIZE(polar_env%psi1_dBerry,1)
             DO ispin = 1,SIZE(polar_env%psi1_dBerry,2)
                CALL cp_fm_release(polar_env%psi1_dBerry(idir,ispin)%matrix)
             ENDDO
          ENDDO
          DEALLOCATE(polar_env%psi1_dBerry)
       ENDIF
       !



    ENDIF ! ref count

  END SUBROUTINE polar_env_cleanup

END MODULE qs_linres_polar_utils
