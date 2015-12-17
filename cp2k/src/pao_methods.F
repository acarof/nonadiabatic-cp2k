!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Methods used by pao_main.F
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_methods
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_binary_read, cp_dbcsr_checksum, &
       cp_dbcsr_col_block_sizes, cp_dbcsr_copy_into_existing, &
       cp_dbcsr_create, cp_dbcsr_desymmetrize, cp_dbcsr_distribution, &
       cp_dbcsr_frobenius_norm, cp_dbcsr_get_block_p, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_reserve_diag_blocks, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_trace, cp_dbcsr_type
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type,&
                                             cp_to_string
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dm_ls_scf_methods,               ONLY: density_matrix_trs4,&
                                             ls_scf_init_matrix_S
  USE dm_ls_scf_qs,                    ONLY: matrix_qs_to_ls
  USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                             ls_scf_env_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE machine,                         ONLY: m_walltime
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_max
  USE pao_input,                       ONLY: pao_exp_param
  USE pao_param,                       ONLY: pao_update_AB
  USE pao_types,                       ONLY: pao_env_type
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                             qs_ks_env_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_methods'

  PUBLIC :: pao_print_atom_info, pao_build_orthogonalizer, pao_build_selector
  PUBLIC :: pao_test_convergence, pao_calc_CG_beta
  PUBLIC :: pao_calc_trace_ph, pao_check_trace_ps
  PUBLIC :: pao_read_preopt_dm
  PUBLIC :: pao_dm_trs4, pao_calc_outer_grad_lnv
  PUBLIC :: pao_check_grad_full

CONTAINS

! *****************************************************************************
!> \brief Prints a one line summary for each atom.
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_print_atom_info(pao)
    TYPE(pao_env_type), POINTER              :: pao

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_print_atom_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, natoms
    INTEGER, DIMENSION(:), POINTER           :: pao_basis, param_cols, &
                                                param_rows, pri_basis

    pri_basis => cp_dbcsr_row_block_sizes(pao%matrix_Y)
    pao_basis => cp_dbcsr_col_block_sizes(pao%matrix_Y)

    CPASSERT(SIZE(pao_basis) == SIZE(pri_basis))
    natoms = SIZE(pao_basis)

    param_rows => cp_dbcsr_row_block_sizes(pao%matrix_X)
    param_cols => cp_dbcsr_col_block_sizes(pao%matrix_X)
    CPASSERT(SIZE(param_rows)==natoms.AND.SIZE(param_cols)==natoms)

    DO iatom=1, natoms
      IF(pao%iw>0) WRITE (pao%iw,"(A,I7,T20,A,I3,T45,A,I3,T65,A,I3)")&
                         " PAO| atom: ", iatom,&
                         " prim_basis: ", pri_basis(iatom),&
                         " pao_basis: ", pao_basis(iatom),&
                         " pao_params: ", (param_cols(iatom)*param_rows(iatom))
    ENDDO
  END SUBROUTINE pao_print_atom_info


! *****************************************************************************
!> \brief Constructs matrix_N and its inverse.
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE pao_build_orthogonalizer(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_build_orthogonalizer', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, iatom, &
                                                j, k, N, natoms
    LOGICAL                                  :: found
    REAL(dp)                                 :: v, w
    REAL(dp), DIMENSION(:), POINTER          :: evals
    REAL(dp), DIMENSION(:, :), POINTER       :: A, block_N, block_N_inv, &
                                                block_S
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    matrix_s=matrix_s,&
                    natom=natoms)

    CALL cp_dbcsr_init(pao%matrix_N)
    CALL cp_dbcsr_create(pao%matrix_N, template=matrix_s(1)%matrix,  name="PAO matrix_N")
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_N)

    CALL cp_dbcsr_init(pao%matrix_N_inv)
    CALL cp_dbcsr_create(pao%matrix_N_inv, template=matrix_s(1)%matrix, name="PAO matrix_N_inv")
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_N_inv)

    CALL cp_dbcsr_iterator_start(iter, pao%matrix_N)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_N)
       IF(arow /= acol) CPABORT("got non-diagonal block")
       iatom = arow

       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_N_inv, row=arow, col=acol, block=block_N_inv, found=found)
       CPASSERT(ASSOCIATED(block_N_inv))

       CALL cp_dbcsr_get_block_p(matrix=matrix_s(1)%matrix, row=arow, col=acol, block=block_S, found=found)
       CPASSERT(ASSOCIATED(block_S))

       IF(SIZE(block_S,1)/=SIZE(block_S,2)) CPABORT("Expected symmetric block")
       N = SIZE(block_S,1)
       ALLOCATE(A(N,N), evals(N))

       ! take square root of atomic overlap matrix
       A(:,:) = block_S(:,:)
       CALL diamat_all(A, evals) !afterwards A contains the eigenvectors
       DO k=1, N
          ! NOTE: To maintain a consistent notation with the Berghold paper,
          ! the "_inv" is swapped: N^{-1}=sqrt(S); N=sqrt(S)^{-1}
          w = 1.0_dp / SQRT(evals(k))
          v = SQRT(evals(k))
          DO i=1, N
             DO j=1, N
                block_N(i,j)     = block_N(i,j)     + w*A(i,k)*A(j,k)
                block_N_inv(i,j) = block_N_inv(i,j) + v*A(i,k)*A(j,k)
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE(A, evals)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE pao_build_orthogonalizer


! *****************************************************************************
!> \brief Build rectangular matrix to converert between primary and PAO basis.
!> \param qs_env ...
!> \param matrix_Y ...
! *****************************************************************************
  SUBROUTINE pao_build_selector(qs_env, matrix_Y)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_type)                      :: matrix_Y

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_build_selector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, ikind, &
                                                nblocks, pao_basis_size
    INTEGER, DIMENSION(:), POINTER           :: blk_sizes_aux, blk_sizes_pri
    REAL(dp), DIMENSION(:, :), POINTER       :: block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    NULLIFY(ks_env, matrix_s)
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    matrix_s=matrix_s,&
                    particle_set=particle_set,&
                    qs_kind_set=qs_kind_set)

    blk_sizes_pri => cp_dbcsr_col_block_sizes(matrix_s(1)%matrix)
    nblocks = SIZE(blk_sizes_pri)
    IF(nblocks /= SIZE(particle_set, 1)) CPABORT("nblocks wrong")

    ALLOCATE(blk_sizes_aux(nblocks))
    DO i=1, nblocks
      CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), pao_basis_size=pao_basis_size)
      IF(pao_basis_size < 1) CPABORT("PAO_BASIS_SIZE too small.")
      blk_sizes_aux(i) = pao_basis_size
    ENDDO

    CALL cp_dbcsr_init(matrix_Y)
    CALL cp_dbcsr_create(matrix_Y,&
                         template=matrix_s(1)%matrix,&
                         matrix_type="N",&
                         row_blk_size=blk_sizes_pri,&
                         col_blk_size=blk_sizes_aux,&
                         name="PAO matrix_Y")
    DEALLOCATE(blk_sizes_aux)

    CALL cp_dbcsr_reserve_diag_blocks(matrix_Y)

    CALL cp_dbcsr_iterator_start(iter, matrix_Y)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block)
       IF(SIZE(block,2)>SIZE(block,1)) CPABORT("Unexpected block sizes")
       block = 0.0_dp
       DO i=1, SIZE(block, 2)
          block(i,i) = 1.0_dp
       ENDDO
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE pao_build_selector


! *****************************************************************************
!> \brief Test whether the PAO optimization has reached convergence
!> \param pao ...
!> \param ls_scf_env ...
!> \param new_energy ...
!> \param icycle ...
!> \param is_converged ...
! *****************************************************************************
  SUBROUTINE pao_test_convergence(pao, ls_scf_env, new_energy, icycle, is_converged)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    REAL(KIND=dp), INTENT(IN)                :: new_energy
    INTEGER, INTENT(IN)                      :: icycle
    LOGICAL, INTENT(OUT)                     :: is_converged

    REAL(KIND=dp)                            :: energy_diff, loop_eps, &
                                                norm_G, now, time_diff

    is_converged = .FALSE.

    ! calculate progress
    energy_diff = new_energy - pao%energy_prev
    pao%energy_prev = new_energy
    now = m_walltime()
    time_diff = now - pao%step_start_time
    pao%step_start_time = now

    norm_G = cp_dbcsr_frobenius_norm(pao%matrix_G)
    IF(pao%iw>0) WRITE(pao%iw,*) "PAO| energy:", new_energy
    IF(pao%iw>0) WRITE(pao%iw,*) "PAO| gradient norm:", norm_G

    IF(pao%istep>1) THEN
       IF(pao%iw>0) WRITE(pao%iw,*) "PAO| energy improvement:", energy_diff
       ! IF(energy_diff>0.0_dp) CPWARN("PAO| energy increased")

       loop_eps = norm_G / ls_scf_env%nelectron_total
       IF(loop_eps < pao%eps_pao) &
         is_converged = .TRUE.

       ! safeguard against noisy gradient
       IF(icycle > 1 .AND. pao%linesearch%step_size < pao%min_step) &
         is_converged = .TRUE.

       ! print one-liner
       IF(pao%iw>0) WRITE (pao%iw,'(A,I6,11X,F20.9,1X,E10.3,1X,E10.3,1X,F6.3)')&
                    " PAO| step ",&
                    pao%istep,&
                    new_energy,&
                    loop_eps,&
                    pao%linesearch%step_size,& !prev step, which let to the current energy
                    time_diff
    ENDIF
  END SUBROUTINE pao_test_convergence


! *****************************************************************************
!> \brief Conjugate gradient algorithm
!> \param pao ...
!> \param beta ...
! *****************************************************************************
  SUBROUTINE pao_calc_CG_beta(pao, beta)
    TYPE(pao_env_type), POINTER              :: pao
    REAL(KIND=dp), INTENT(OUT)               :: beta

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_CG_beta', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: trace_G_mix, trace_G_new, &
                                                trace_G_prev

    beta = 0.0_dp !default

    IF(pao%istep > pao%cg_init_steps) THEN
       CALL cp_dbcsr_trace(pao%matrix_G, pao%matrix_G, trace_G_new)
       CALL cp_dbcsr_trace(pao%matrix_G_prev, pao%matrix_G_prev, trace_G_prev)
       CALL cp_dbcsr_trace(pao%matrix_G, pao%matrix_G_prev, trace_G_mix)
       IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| trace_G_new ", trace_G_new
       IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| trace_G_prev ", trace_G_prev
       IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| trace_G_mix ", trace_G_mix

       IF(trace_G_prev/=0.0_dp) THEN
          beta = (trace_G_new - trace_G_mix) / trace_G_prev  !Polak–Ribiere
       ENDIF

       IF(beta < 0.0_dp) THEN
          IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| resetting because beta < 0"
          beta = 0.0_dp
       ENDIF

       IF(pao%max_outer_pao>0 .AND. pao%istep-pao%last_reset > pao%max_outer_pao) THEN
          IF(pao%iw>0) WRITE(pao%iw,*) "PAO|CG| resetting because of MAX_OUTER_PAO"
          pao%last_reset = pao%istep
          beta = 0.0_dp
       ENDIF
    ENDIF

    IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| beta: ", beta

  END SUBROUTINE pao_calc_CG_beta


! *****************************************************************************
!> \brief Calculate the pre-optimization energy
!> \param ls_scf_env ...
!> \param energy ...
! *****************************************************************************
  SUBROUTINE pao_calc_trace_PH(ls_scf_env, energy)
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    REAL(KIND=dp), INTENT(OUT)               :: energy

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_trace_PH', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(KIND=dp)                            :: trace_PH

    CALL timeset(routineN,handle)
    energy = 0.0_dp
    DO ispin=1, ls_scf_env%nspins
       CALL cp_dbcsr_trace(ls_scf_env%matrix_p(ispin), ls_scf_env%matrix_ks(ispin), trace_PH)
       energy = energy + trace_PH
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE pao_calc_trace_PH


! *****************************************************************************
!> \brief Ensure that the number of electrons is correct.
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_check_trace_PS(ls_scf_env)
    TYPE(ls_scf_env_type)                    :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_check_trace_PS', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(KIND=dp)                            :: tmp, trace_PS
    TYPE(cp_dbcsr_type)                      :: matrix_S_desym

    CALL timeset(routineN,handle)
    CALL cp_dbcsr_init(matrix_S_desym)
    CALL cp_dbcsr_create(matrix_S_desym, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_desymmetrize(ls_scf_env%matrix_s, matrix_S_desym)

    trace_PS = 0.0_dp
    DO ispin=1, ls_scf_env%nspins
       CALL cp_dbcsr_trace(ls_scf_env%matrix_p(ispin), matrix_S_desym, tmp)
       trace_PS = trace_PS + tmp
    ENDDO

    CALL cp_dbcsr_release(matrix_S_desym)

    IF(ABS(ls_scf_env%nelectron_total - trace_PS) > 0.5)&
         CPABORT("Number of electrons wrong. Trace(PS) ="//cp_to_string(trace_PS))

    CALL timestop(handle)
  END SUBROUTINE pao_check_trace_PS


! *****************************************************************************
!> \brief Read primary density matrix from file.
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE pao_read_preopt_dm(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_read_preopt_dm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(KIND=dp)                            :: cs_pos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, rho_ao
    TYPE(cp_dbcsr_type)                      :: matrix_tmp
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    dft_control=dft_control,&
                    matrix_s=matrix_s,&
                    rho=rho,&
                    energy=energy)

    CALL qs_rho_get(rho, rho_ao=rho_ao)

    IF(dft_control%nspins/=1) CPABORT("open shell not yet implemented")

    DO ispin=1, dft_control%nspins
       CALL cp_dbcsr_init(matrix_tmp) ! read into tmp-matrix to preserve sparsity of rho_ao
       CALL cp_dbcsr_binary_read(pao%preopt_dm_file, matrix_new=matrix_tmp,&
                                 distribution=cp_dbcsr_distribution(matrix_s(1)%matrix))
       cs_pos = cp_dbcsr_checksum(matrix_tmp, pos=.TRUE.)
       IF (pao%iw>0) WRITE(pao%iw,'(T2,A,E20.8)') "Read restart DM "//TRIM(pao%preopt_dm_file)//&
                                                  " with checksum: ",cs_pos
       CALL cp_dbcsr_copy_into_existing(rho_ao(ispin)%matrix, matrix_tmp)
       CALL cp_dbcsr_release(matrix_tmp)
    ENDDO

    ! calculate corresponding ks matrix
    CALL qs_rho_update_rho(rho,qs_env=qs_env)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.)
    CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE.,&
                             just_energy=.FALSE., print_active=.TRUE.)
    IF (pao%iw>0) WRITE(pao%iw,*) "PAO| Quickstep energy from restart density:",energy%total

    CALL timestop(handle)

  END SUBROUTINE pao_read_preopt_dm


! *****************************************************************************
!> \brief Calculate density matrix using TRS4 purification
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_dm_trs4(qs_env, ls_scf_env)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_dm_trs4', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: project_name
    INTEGER                                  :: handle, ispin, &
                                                nelectron_spin_real, nspin
    LOGICAL                                  :: converged
    REAL(KIND=dp)                            :: homo_spin, lumo_spin, mu_spin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    project_name = logger%iter_info%project_name
    nspin = ls_scf_env%nspins

    CALL get_qs_env(qs_env,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s)

    ! rebuild S, S_inv, S_sqrt, and S_sqrt_inv in pao basis
    CALL cp_dbcsr_release(ls_scf_env%matrix_s_inv)
    CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt)
    CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt_inv)
    CALL ls_scf_init_matrix_s(matrix_s(1)%matrix, ls_scf_env)

    DO ispin=1, nspin
       CALL matrix_qs_to_ls(ls_scf_env%matrix_ks(ispin),matrix_ks(ispin)%matrix,&
                            ls_scf_env%ls_mstruct,covariant=.TRUE.)

       nelectron_spin_real = ls_scf_env%nelectron_spin(ispin)
       IF(ls_scf_env%nspins==1) nelectron_spin_real = nelectron_spin_real/2
       CALL density_matrix_trs4(ls_scf_env%matrix_p(ispin), ls_scf_env%matrix_ks(ispin),&
                                ls_scf_env%matrix_s_sqrt_inv,&
                                nelectron_spin_real, ls_scf_env%eps_filter, homo_spin, lumo_spin, mu_spin, &
                                dynamic_threshold=.FALSE., converged=converged,&
                                max_iter_lanczos=ls_scf_env%max_iter_lanczos,&
                                eps_lanczos=ls_scf_env%eps_lanczos)
       IF(.NOT.converged) CPABORT("TRS4 did not converge")
    ENDDO

    IF(nspin==1) CALL cp_dbcsr_scale(ls_scf_env%matrix_p(1), 2.0_dp)

    CALL timestop(handle)
  END SUBROUTINE pao_dm_trs4


! *****************************************************************************
!> \brief Helper routine, calculates partial derivative dE/dU
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_M ...
! *****************************************************************************
   SUBROUTINE pao_calc_outer_grad_lnv(qs_env,ls_scf_env,matrix_M)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    TYPE(cp_dbcsr_type)                      :: matrix_M

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_outer_grad_lnv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, rho_ao
    TYPE(cp_dbcsr_type) :: matrix_HB, matrix_HPS, matrix_M1, matrix_M2, &
      matrix_M3, matrix_NHB, matrix_NHBM2, matrix_NPA, matrix_NPAM1, &
      matrix_NSB, matrix_NSBM3, matrix_PA, matrix_PH, matrix_PHP, matrix_PSP, &
      matrix_SB, matrix_SP
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    ls_mstruct => ls_scf_env%ls_mstruct
    pao => ls_scf_env%pao_env

    CALL get_qs_env(qs_env,&
                    rho=rho,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    dft_control=dft_control)
    CALL qs_rho_get(rho, rho_ao=rho_ao)
    nspin = dft_control%nspins

    IF(nspin /= 1) CPABORT("open shell not yet implemented")
    !TODO: handle openshell case properly

    ! notation according to pao_math_lnv.pdf

    CALL cp_dbcsr_init(matrix_M)
    CALL cp_dbcsr_create(matrix_M, template=matrix_s(1)%matrix, matrix_type="N")

    !---------------------------------------------------------------------------
    ! calculate need products in pao basis
    CALL cp_dbcsr_init(matrix_PH)
    CALL cp_dbcsr_create(matrix_PH, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), ls_scf_env%matrix_ks(1),&
                          0.0_dp, matrix_PH)

    CALL cp_dbcsr_init(matrix_PHP)
    CALL cp_dbcsr_create(matrix_PHP, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_PH, ls_scf_env%matrix_p(1),&
                           0.0_dp, matrix_PHP)

    CALL cp_dbcsr_init(matrix_SP)
    CALL cp_dbcsr_create(matrix_SP, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s, ls_scf_env%matrix_p(1),&
                          0.0_dp, matrix_SP)

    IF(SIZE(ls_scf_env%matrix_p)==1) CALL cp_dbcsr_scale(matrix_SP, 0.5_dp)

    CALL cp_dbcsr_init(matrix_HPS)
    CALL cp_dbcsr_create(matrix_HPS, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, ls_scf_env%matrix_ks(1), matrix_SP,&
                          0.0_dp, matrix_HPS)

    CALL cp_dbcsr_init(matrix_PSP)
    CALL cp_dbcsr_create(matrix_PSP, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), matrix_SP,&
                           0.0_dp, matrix_PSP)

    !---------------------------------------------------------------------------
    ! M1 = dE_lnv / dP_pao
    CALL cp_dbcsr_init(matrix_M1)
    CALL cp_dbcsr_create(matrix_M1, template=ls_scf_env%matrix_s, matrix_type="N")

    CALL cp_dbcsr_multiply("N", "T", 3.0_dp, ls_scf_env%matrix_ks(1), matrix_SP,&
                           1.0_dp, matrix_M1)

    CALL cp_dbcsr_multiply("N", "N", 3.0_dp, matrix_SP, ls_scf_env%matrix_ks(1),&
                           1.0_dp, matrix_M1)

    CALL cp_dbcsr_multiply("N", "T", -2.0_dp, matrix_HPS, matrix_SP,&
                           1.0_dp, matrix_M1)

    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_SP, matrix_HPS, &
                           1.0_dp, matrix_M1)

    CALL cp_dbcsr_multiply("N", "T", -2.0_dp, matrix_SP, matrix_HPS,&
                           1.0_dp, matrix_M1)


    !---------------------------------------------------------------------------
    ! M2 = dE_lnv / dH
    CALL cp_dbcsr_init(matrix_M2)
    CALL cp_dbcsr_create(matrix_M2, template=ls_scf_env%matrix_s, matrix_type="N")

    CALL cp_dbcsr_add(matrix_M2, matrix_PSP, 1.0_dp, 3.0_dp)

    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PSP, matrix_SP,&
                           1.0_dp, matrix_M2)

    !---------------------------------------------------------------------------
    ! M3 = dE_lnv / dS
    CALL cp_dbcsr_init(matrix_M3)
    CALL cp_dbcsr_create(matrix_M3, template=ls_scf_env%matrix_s, matrix_type="N")

    CALL cp_dbcsr_add(matrix_M3, matrix_PHP, 1.0_dp, 3.0_dp)

    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PHP, matrix_SP,&
                           1.0_dp, matrix_M3)

    CALL cp_dbcsr_multiply("N", "T", -2.0_dp, matrix_PSP, matrix_PH,&
                           1.0_dp, matrix_M3)


    !---------------------------------------------------------------------------
    ! combine M1 with matrices from primary basis
    CALL cp_dbcsr_init(matrix_PA)
    CALL cp_dbcsr_create(matrix_PA, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, rho_ao(1)%matrix, ls_mstruct%matrix_A ,&
                           0.0_dp, matrix_PA)

    CALL cp_dbcsr_init(matrix_NPA)
    CALL cp_dbcsr_create(matrix_NPA, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N_inv, matrix_PA,&
                           0.0_dp, matrix_NPA)

    CALL cp_dbcsr_init(matrix_NPAM1)
    CALL cp_dbcsr_create(matrix_NPAM1, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NPA, matrix_M1,&
                           0.0_dp, matrix_NPAM1)

    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NPAM1, pao%matrix_Y,&
                           1.0_dp, matrix_M)


    !---------------------------------------------------------------------------
    ! combine M2 with matrices from primary basis
    CALL cp_dbcsr_init(matrix_HB)
    CALL cp_dbcsr_create(matrix_HB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ks(1)%matrix, ls_mstruct%matrix_B,&
                           0.0_dp, matrix_HB)

    CALL cp_dbcsr_init(matrix_NHB)
    CALL cp_dbcsr_create(matrix_NHB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, matrix_HB,&
                           0.0_dp, matrix_NHB)

    CALL cp_dbcsr_init(matrix_NHBM2)
    CALL cp_dbcsr_create(matrix_NHBM2, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NHB, matrix_M2,&
                           0.0_dp, matrix_NHBM2)

    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NHBM2, pao%matrix_Y,&
                           1.0_dp, matrix_M)


    !---------------------------------------------------------------------------
    ! combine M3 with matrices from primary basis
    CALL cp_dbcsr_init(matrix_SB)
    CALL cp_dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B,&
                           0.0_dp, matrix_SB)

    IF(SIZE(ls_scf_env%matrix_p)==1) CALL cp_dbcsr_scale(matrix_SB, 0.5_dp)

    CALL cp_dbcsr_init(matrix_NSB)
    CALL cp_dbcsr_create(matrix_NSB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, matrix_SB,&
                           0.0_dp, matrix_NSB)

    CALL cp_dbcsr_init(matrix_NSBM3)
    CALL cp_dbcsr_create(matrix_NSBM3, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NSB, matrix_M3,&
                           0.0_dp, matrix_NSBM3)

    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NSBM3, pao%matrix_Y,&
                           1.0_dp, matrix_M)

    IF(SIZE(ls_scf_env%matrix_p)==1) CALL cp_dbcsr_scale(matrix_M, 2.0_dp)

    !---------------------------------------------------------------------------
    ! cleanup
    CALL cp_dbcsr_release(matrix_PH)
    CALL cp_dbcsr_release(matrix_PHP)
    CALL cp_dbcsr_release(matrix_SP)
    CALL cp_dbcsr_release(matrix_HPS)
    CALL cp_dbcsr_release(matrix_PSP)
    CALL cp_dbcsr_release(matrix_M1)
    CALL cp_dbcsr_release(matrix_M2)
    CALL cp_dbcsr_release(matrix_M3)
    CALL cp_dbcsr_release(matrix_PA)
    CALL cp_dbcsr_release(matrix_NPA)
    CALL cp_dbcsr_release(matrix_NPAM1)
    CALL cp_dbcsr_release(matrix_HB)
    CALL cp_dbcsr_release(matrix_NHB)
    CALL cp_dbcsr_release(matrix_NHBM2)
    CALL cp_dbcsr_release(matrix_SB)
    CALL cp_dbcsr_release(matrix_NSB)
    CALL cp_dbcsr_release(matrix_NSBM3)

    CALL timestop(handle)
  END SUBROUTINE pao_calc_outer_grad_lnv


! *****************************************************************************
!> \brief Debugging routine for checking the analytic gradient.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_check_grad_full(pao, qs_env, ls_scf_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_check_grad_full', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, j, natoms
    INTEGER, DIMENSION(:), POINTER           :: blk_sizes_col, blk_sizes_row
    LOGICAL                                  :: found
    REAL(dp)                                 :: delta_max, energy1, energy2, &
                                                eps, Gij_num, symm
    REAL(dp), DIMENSION(:, :), POINTER       :: block_G, block_X, block_X0
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct

    IF(pao%check_grad_full_tol<0.0_dp) RETURN ! no checking

    CALL timeset(routineN,handle)

    ls_mstruct => ls_scf_env%ls_mstruct

    SELECT CASE(pao%parameterization)
      CASE(pao_exp_param)
        symm = -1.0_dp ! anti-symmetric
      CASE DEFAULT
        symm = 0.0_dp  ! no symmetry
    END SELECT

    CALL get_qs_env(qs_env, para_env=para_env, natom=natoms)

    eps = 1.0e-5_dp
    delta_max = 0.0_dp

    blk_sizes_col => cp_dbcsr_col_block_sizes(pao%matrix_X)
    blk_sizes_row => cp_dbcsr_row_block_sizes(pao%matrix_X)

    ! can not use an iterator here, because other DBCSR routines are called within loop.
    DO iatom=1, natoms
       IF(pao%iw>0) WRITE(pao%iw,*) 'PAO| checking full gradient of atom ',iatom
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_X, row=iatom, col=iatom, block=block_X, found=found)

       IF(ASSOCIATED(block_X)) THEN !only one node actually has the block
          CALL cp_dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
          CPASSERT(ASSOCIATED(block_G))
          ALLOCATE(block_X0(blk_sizes_row(iatom),blk_sizes_col(iatom)))
          block_X0 = block_X
       ENDIF

       DO i=1, blk_sizes_row(iatom)
          DO j=1, blk_sizes_col(iatom)
             IF(ASSOCIATED(block_X)) THEN
                block_X = block_X0
                block_X(i,j) = block_X(i,j) + eps
                IF(symm/=0.0) block_X(j,i) = block_X(j,i) + symm*eps
             ENDIF

             CALL pao_update_AB(pao, ls_mstruct)
             CALL pao_dm_trs4(qs_env, ls_scf_env)
             CALL pao_calc_trace_PH(ls_scf_env, energy1)

             IF(ASSOCIATED(block_X)) THEN
                block_X = block_X0
                block_X(i,j) = block_X(i,j) - eps
                IF(symm/=0.0) block_X(j,i) = block_X(j,i) - symm*eps
             ENDIF

             CALL pao_update_AB(pao, ls_mstruct)
             CALL pao_dm_trs4(qs_env, ls_scf_env)
             CALL pao_calc_trace_PH(ls_scf_env, energy2)

             IF(ASSOCIATED(block_X)) THEN
                Gij_num = (energy1 - energy2) / (2.0_dp*eps)
                delta_max = MAX(delta_max, ABS(Gij_num - block_G(i,j)))
             ENDIF
          ENDDO
       ENDDO

       IF(ASSOCIATED(block_X)) THEN
          block_X = block_X0
          DEALLOCATE(block_X0)
       ENDIF
    END DO

    CALL mp_max(delta_max, para_env%group)
    IF(pao%iw>0) WRITE(pao%iw,*) 'PAO| checked full gradient, max delta:', delta_max
    IF(delta_max > pao%check_grad_full_tol)  CALL cp_abort(__LOCATION__,&
       "Analytic and numeric gradients differ too much:"//cp_to_string(delta_max))

    CALL timestop(handle)
  END SUBROUTINE pao_check_grad_full

END MODULE pao_methods