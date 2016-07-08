MODULE tea_leaf_dpcg_module

  USE tea_leaf_dpcg_kernel_module
  USE tea_leaf_cheby_module
  USE tea_leaf_common_module
  USE tea_leaf_cg_module

  USE definitions_module
  use global_mpi_module
  USE update_halo_module

  IMPLICIT NONE

  LOGICAL :: inner_use_ppcg
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: inner_cg_alphas, inner_cg_betas
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: inner_ch_alphas, inner_ch_betas
  REAL(KIND=8) :: eigmin, eigmax, inner_ch_theta

CONTAINS

SUBROUTINE tea_leaf_dpcg_init_x0(solve_time, ppcg_inner_steps, ch_alphas, ch_betas, ch_theta)

  USE tea_leaf_ppcg_module

  IMPLICIT NONE

  INTEGER :: ppcg_inner_steps
  REAL(KIND=8) :: solve_time, ch_theta
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas

  REAL(KIND=8) :: tile_sum2

  INTEGER, PARAMETER :: level=1
  INTEGER :: t, it_count, info
  INTEGER :: fields(NUM_FIELDS)
  REAL(KIND=8) :: halo_time,timer,rro
!$ INTEGER :: OMP_GET_THREAD_NUM

  IF (.NOT. ALLOCATED(inner_cg_alphas)) THEN
    ALLOCATE(inner_cg_alphas(tl_ch_cg_presteps))
    ALLOCATE(inner_cg_betas (tl_ch_cg_presteps))
    ALLOCATE(inner_ch_alphas(tl_ppcg_inner_coarse_steps))
    ALLOCATE(inner_ch_betas (tl_ppcg_inner_coarse_steps))
  ENDIF

  CALL tea_leaf_dpcg_coarsen_matrix_level(level,solve_time)

  ! just use CG on the first one
  inner_use_ppcg = .FALSE.

!$OMP PARALLEL
!$OMP DO
    ! store the RHS in u which will then be copied into u0 
    DO t=1,tiles_per_task
      chunk(level+1)%tiles(t)%field%u = 0.0_8
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    CALL tea_leaf_dpcg_restrict_ZT_level(level,.TRUE.)

    CALL tea_leaf_dpcg_local_solve_level(level,                  &
                                         solve_time,             &
                                         tl_ch_cg_epslim,        &
                                         tl_ch_cg_presteps+1,    & ! need beta(tl_ch_cg_presteps)!=0
                                         it_count,               &
                                         inner_use_ppcg,         &
                                         inner_cg_alphas,        &
                                         inner_cg_betas,         &
                                         inner_ch_theta,         &
                                         inner_ch_alphas,        &
                                         inner_ch_betas)

    CALL tea_leaf_dpcg_subtract_z_level(level)

  ! for all subsequent steps, use ppcg if requested
  inner_use_ppcg = coarse_solve_ppcg

  IF (parallel%boss.AND.verbose_on) THEN
!$  IF (OMP_GET_THREAD_NUM().EQ.0) THEN
      do t=1,it_count-1
        write(g_out,*) "i=",t," alpha=",inner_cg_alphas(t)," beta=",inner_cg_betas(t)
      enddo
!$  ENDIF
  ENDIF
  CALL tea_calc_eigenvalues(inner_cg_alphas, inner_cg_betas, eigmin, eigmax, &
      coarse_solve_max_iters, it_count-1, info) ! we don't have the beta value for the last iteration
  IF (parallel%boss.AND.verbose_on) THEN
!$  IF (OMP_GET_THREAD_NUM().EQ.0) THEN
      write(g_out,*) "eigmin=",eigmin," eigmax=",eigmax
!$  ENDIF
  ENDIF
  !info = 0
  eigmin = eigmin * 0.95
  eigmax = eigmax * 1.05

  !! With jacobi preconditioner on
  !!eigmin = 0.1_8
  !eigmin = tl_ppcg_coarse_eigmin
  !!eigmax = 2.0_8
  !!New bound with the l1 Jacobi preconditioner
  !eigmax = 1.0_8

  IF (info .NE. 0) CALL report_error('tea_leaf_dpcg_init_x0', 'Error in calculating eigenvalues')

  CALL tea_calc_ch_coefs(inner_ch_alphas, inner_ch_betas, eigmin, eigmax, &
      inner_ch_theta, tl_ppcg_inner_coarse_steps)

  ! With jacobi preconditioner on
  !eigmin = 0.1_8
  eigmin = tl_ppcg_steps_eigmin
  !eigmax = 2.0_8
  !New bound with the l1 Jacobi preconditioner
  eigmax = 1.0_8

  CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
      ch_theta, ppcg_inner_steps)
  !write(6,*) maxval(abs(inner_ch_alphas)), maxval(abs(inner_ch_betas)), inner_ch_theta, &
  !  tl_ppcg_inner_coarse_steps, level

  fields = 0
  fields(FIELD_U) = 1

  ! update the halo for u prior to recalculating the residual
  IF (profiler_on) halo_time = timer()
  CALL update_halo(level, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  ! calc residual again, and do initial solve
  CALL tea_leaf_calc_residual(level)

  ! init z using the preconditioner - note that we don't need rro
  CALL tea_leaf_cg_init(level, ppcg_inner_steps, ch_alphas, ch_betas, ch_theta, solve_time, 1, rro)

  !write(6,*) "tea_leaf_run_ppcg_inner_steps:",ppcg_inner_steps,level, &
  !    maxval(abs(ch_alphas(1:ppcg_inner_steps))), &
  !    maxval(abs(ch_betas (1:ppcg_inner_steps))), ch_theta
  CALL tea_leaf_run_ppcg_inner_steps(level, ch_alphas, ch_betas, ch_theta, &
      ppcg_inner_steps, solve_time)

  CALL tea_leaf_dpcg_setup_and_solve_E_level(level,solve_time)

  CALL tea_leaf_dpcg_init_p(level)

END SUBROUTINE tea_leaf_dpcg_init_x0

SUBROUTINE tea_leaf_dpcg_setup_and_solve_E_level(level,solve_time)

  IMPLICIT NONE

  REAL(KIND=8) :: solve_time, tile_sum2

  INTEGER :: level, it_count, t

  CALL tea_leaf_dpcg_matmul_ZTA_level(level,solve_time)
  CALL tea_leaf_dpcg_restrict_ZT_level(level,.TRUE.)

  CALL tea_leaf_dpcg_local_solve_level(level,                  &
                                       solve_time,             &
                                       coarse_solve_eps,       &
                                       coarse_solve_max_iters, &
                                       it_count,               &
                                       inner_use_ppcg,         &
                                       inner_cg_alphas,        &
                                       inner_cg_betas,         &
                                       inner_ch_theta,         &
                                       inner_ch_alphas,        &
                                       inner_ch_betas)

  CALL tea_leaf_dpcg_prolong_Z_level(level)

END SUBROUTINE tea_leaf_dpcg_setup_and_solve_E_level

SUBROUTINE tea_leaf_dpcg_coarsen_matrix_level(level,solve_time)

  IMPLICIT NONE
  INTEGER :: t, err
  INTEGER :: level

  INTEGER :: sub_tile_dx, sub_tile_dy
  INTEGER :: sx, sy
  INTEGER :: xrem, yrem

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: jj,j_start,j_end,kk,k_start,k_end

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: solve_time
  REAL(KIND=8) :: timer,halo_time
  REAL(KIND=8) :: tile_size
  !REAL(KIND=8) :: kx_tot,ky_tot,Di_tot,Mi_tot
  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: kx_local, ky_local

  ! 3 different options for preconditioners
  INTEGER,PARAMETER::           TL_PREC_NONE       = 1 &
                               ,TL_PREC_JAC_DIAG   = 2 &
                               ,TL_PREC_JAC_BLOCK  = 3

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(Kx_local, Ky_local,sub_tile_dx,sub_tile_dy)
!$OMP DO
    DO t=1,tiles_per_task
      ! define the coarse Kx and Ky with unit scaling factor as the coefficient values are obtained by a Galerkin projection of the matrix
      chunk(level+1)%tiles(t)%field%rx = 1.0_8
      chunk(level+1)%tiles(t)%field%ry = 1.0_8
    ENDDO
!$OMP END DO

!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1)/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1)/chunk(level)%sub_tile_dims(2)

      xrem        =mod(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1,chunk(level)%sub_tile_dims(1))
      yrem        =mod(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1,chunk(level)%sub_tile_dims(2))

      chunk(level+1)%tiles(t)%field%vector_Kx = 0.0_8
      chunk(level+1)%tiles(t)%field%vector_Ky = 0.0_8

      kx_local = 0.0_8
      ky_local = 0.0_8

      CALL tea_leaf_dpcg_coarsen_matrix_kernel(chunk(level)%tiles(t)%field%x_min,     &
                                               chunk(level)%tiles(t)%field%x_max,     &
                                               chunk(level)%tiles(t)%field%y_min,     &
                                               chunk(level)%tiles(t)%field%y_max,     &
                                               chunk(level)%halo_exchange_depth,      &
                                               chunk(level)%tiles(t)%field%vector_Kx, &
                                               chunk(level)%tiles(t)%field%vector_Ky, &
                                               chunk(level)%sub_tile_dims(1),         &
                                               sub_tile_dx, xrem,                     &
                                               chunk(level)%sub_tile_dims(2),         &
                                               sub_tile_dy, yrem,                     &
                                               kx_local,                       &
                                               ky_local,                       &
                                               chunk(level)%tiles(t)%field%rx,        &
                                               chunk(level)%tiles(t)%field%ry)

      chunk(level+1)%tiles(t)%field%vector_Kx(1:chunk(level)%sub_tile_dims(1),1:chunk(level)%sub_tile_dims(2)) = kx_local
      chunk(level+1)%tiles(t)%field%vector_Ky(1:chunk(level)%sub_tile_dims(1),1:chunk(level)%sub_tile_dims(2)) = ky_local
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

!Need a depth chunk(level+1)%halo_exchange_depth halo exchange on Kx and Ky
!use custom comms for Kx, Ky (and Kz in 3D)
  fields=0
  fields(FIELD_KX)=1
  fields(FIELD_KY)=1
  IF (profiler_on) halo_time = timer()
  CALL update_halo(level+1, fields, chunk(level+1)%halo_exchange_depth)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_size,sub_tile_dx,sub_tile_dy,k,kk,k_start,k_end,j,jj,j_start,j_end)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1)/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1)/chunk(level)%sub_tile_dims(2)

      xrem        =mod(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1,chunk(level)%sub_tile_dims(1))
      yrem        =mod(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1,chunk(level)%sub_tile_dims(2))

      chunk(level+1)%tiles(t)%field%vector_Di = 0.0_8

      DO kk=1,chunk(level)%sub_tile_dims(2)
        k_start=min(kk-1,yrem)+(kk-1)*sub_tile_dy+chunk(level)%tiles(t)%field%y_min
        k_end  =min(kk  ,yrem)+kk    *sub_tile_dy+chunk(level)%tiles(t)%field%y_min-1
        DO jj=1,chunk(level)%sub_tile_dims(1)
          j_start=min(jj-1,xrem)+(jj-1)*sub_tile_dx+chunk(level)%tiles(t)%field%x_min
          j_end  =min(jj  ,xrem)+jj    *sub_tile_dx+chunk(level)%tiles(t)%field%x_min-1
          tile_size=(j_end-j_start+1)*(k_end-k_start+1)
          chunk(level+1)%tiles(t)%field%vector_Di(jj,kk) = &
            tile_size + &
            chunk(level+1)%tiles(t)%field%vector_Kx(jj    ,kk    ) + &
            chunk(level+1)%tiles(t)%field%vector_Ky(jj    ,kk    ) + &
            chunk(level+1)%tiles(t)%field%vector_Kx(jj + 1,kk    ) + &
            chunk(level+1)%tiles(t)%field%vector_Ky(jj    ,kk + 1)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

!Need a depth chunk(level+1)%halo_exchange_depth-1 halo exchange on Di
!use custom comms for Di
  IF (chunk(level+1)%halo_exchange_depth-1 > 0) THEN
    fields=0
    fields(FIELD_DI)=1
    IF (profiler_on) halo_time = timer()
    CALL update_halo(level+1, fields, chunk(level+1)%halo_exchange_depth-1)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)
  ENDIF

  IF (use_fortran_kernels) THEN
    IF (tl_preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
! The preconditioner construction must be called outside of a parallel region
! as the routine is threaded otherwise the result array is not shared correctly
      DO t=1,tiles_per_task
        CALL tea_block_init(chunk(level+1)%tiles(t)%field%x_min,     &
                            chunk(level+1)%tiles(t)%field%x_max,     &
                            chunk(level+1)%tiles(t)%field%y_min,     &
                            chunk(level+1)%tiles(t)%field%y_max,     &
                            chunk(level+1)%halo_exchange_depth,      &
                            chunk(level+1)%tiles(t)%field%tri_cp,    &
                            chunk(level+1)%tiles(t)%field%tri_bfp,   &
                            chunk(level+1)%tiles(t)%field%vector_Kx, &
                            chunk(level+1)%tiles(t)%field%vector_Ky, &
                            chunk(level+1)%tiles(t)%field%vector_Di, &
                            chunk(level+1)%tiles(t)%field%rx,        &
                            chunk(level+1)%tiles(t)%field%ry)
      ENDDO
    ELSE IF (tl_preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
! The preconditioner construction must be called outside of a parallel region
! as the routine is threaded otherwise the result array is not shared correctly
      DO t=1,tiles_per_task
        CALL tea_diag_init(chunk(level+1)%tiles(t)%field%x_min,     &
                           chunk(level+1)%tiles(t)%field%x_max,     &
                           chunk(level+1)%tiles(t)%field%y_min,     &
                           chunk(level+1)%tiles(t)%field%y_max,     &
                           chunk(level+1)%halo_exchange_depth,      &
                           chunk(level+1)%tiles(t)%field%vector_Mi, &
                           chunk(level+1)%tiles(t)%field%vector_Kx, &
                           chunk(level+1)%tiles(t)%field%vector_Ky, &
                           chunk(level+1)%tiles(t)%field%vector_Di, &
                           chunk(level+1)%tiles(t)%field%rx,        &
                           chunk(level+1)%tiles(t)%field%ry)
      ENDDO
    ENDIF
  ENDIF

END SUBROUTINE tea_leaf_dpcg_coarsen_matrix_level

SUBROUTINE tea_leaf_dpcg_matmul_ZTA_level(level,solve_time)

  IMPLICIT NONE

  INTEGER :: level
  REAL(KIND=8) :: solve_time

  INTEGER :: t, err

  INTEGER :: sub_tile_dx, sub_tile_dy
  INTEGER :: xrem, yrem

  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: ztaz

  REAL(KIND=8) :: halo_time,timer

  INTEGER :: fields(NUM_FIELDS)

  fields = 0
  fields(FIELD_Z) = 1

  IF (profiler_on) halo_time = timer()
  CALL update_halo(level, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  ! store the RHS in u which will then be copied into u0
!$OMP PARALLEL
!$OMP DO
  DO t=1,tiles_per_task
    chunk(level+1)%tiles(t)%field%u = 0.0_8
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztaz,sub_tile_dx,sub_tile_dy)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1)/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1)/chunk(level)%sub_tile_dims(2)

      xrem        =mod(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1,chunk(level)%sub_tile_dims(1))
      yrem        =mod(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1,chunk(level)%sub_tile_dims(2))

      ztaz = 0.0_8

      CALL tea_leaf_dpcg_matmul_ZTA_kernel(chunk(level)%tiles(t)%field%x_min, &
          chunk(level)%tiles(t)%field%x_max,                                  &
          chunk(level)%tiles(t)%field%y_min,                                  &
          chunk(level)%tiles(t)%field%y_max,                                  &
          chunk(level)%halo_exchange_depth,                                   &
          chunk(level)%sub_tile_dims(1),sub_tile_dx,xrem,                     &
          chunk(level)%sub_tile_dims(2),sub_tile_dy,yrem,                     &
          chunk(level)%tiles(t)%field%vector_z,                               &
          chunk(level)%tiles(t)%field%vector_Kx,                              &
          chunk(level)%tiles(t)%field%vector_Ky,                              &
          chunk(level)%tiles(t)%field%vector_Di,                              &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,  &
          ztaz)

      ! write back into the GLOBAL vector
      chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                      1:chunk(level)%sub_tile_dims(2)) = ztaz
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_matmul_ZTA_level

SUBROUTINE tea_leaf_dpcg_restrict_ZT_level(level,not_init)

  IMPLICIT NONE
  LOGICAL :: not_init
  INTEGER ::level

  INTEGER :: t, err
  INTEGER :: sub_tile_dx, sub_tile_dy
  INTEGER :: xrem, yrem
  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: ZTr

  IF (use_fortran_kernels) THEN
    IF (not_init) THEN
!$OMP PARALLEL PRIVATE(ZTr,sub_tile_dx,sub_tile_dy)
!$OMP DO
      DO t=1,tiles_per_task
        sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1)/chunk(level)%sub_tile_dims(1)
        sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1)/chunk(level)%sub_tile_dims(2)

        xrem        =mod(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1,chunk(level)%sub_tile_dims(1))
        yrem        =mod(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1,chunk(level)%sub_tile_dims(2))

        ztr = 0.0_8

        CALL tea_leaf_dpcg_restrict_ZT_kernel(chunk(level)%tiles(t)%field%x_min, &
                                              chunk(level)%tiles(t)%field%x_max, &
                                              chunk(level)%tiles(t)%field%y_min,           &
                                              chunk(level)%tiles(t)%field%y_max,           &
                                              chunk(level)%halo_exchange_depth,            &
                                              chunk(level)%sub_tile_dims(1),               &
                                              sub_tile_dx, xrem,                           &
                                              chunk(level)%sub_tile_dims(2),               &
                                              sub_tile_dy, yrem,                           &
                                              chunk(level)%tiles(t)%field%vector_r,        &
                                              ztr)

        chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                        1:chunk(level)%sub_tile_dims(2)) = &
        chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                        1:chunk(level)%sub_tile_dims(2)) - ztr
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ELSE
!$OMP PARALLEL PRIVATE(ZTr,sub_tile_dx,sub_tile_dy)
!$OMP DO
      DO t=1,tiles_per_task
        sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1)/chunk(level)%sub_tile_dims(1)
        sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1)/chunk(level)%sub_tile_dims(2)

        xrem        =mod(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1,chunk(level)%sub_tile_dims(1))
        yrem        =mod(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1,chunk(level)%sub_tile_dims(2))

        ztr = 0.0_8

        CALL tea_leaf_dpcg_restrict_ZT_kernel(chunk(level)%tiles(t)%field%x_min, &
                                              chunk(level)%tiles(t)%field%x_max, &
                                              chunk(level)%tiles(t)%field%y_min,           &
                                              chunk(level)%tiles(t)%field%y_max,           &
                                              chunk(level)%halo_exchange_depth,            &
                                              chunk(level)%sub_tile_dims(1),               &
                                              sub_tile_dx, xrem,                           &
                                              chunk(level)%sub_tile_dims(2),               &
                                              sub_tile_dy, yrem,                           &
                                              chunk(level)%tiles(t)%field%vector_r,        &
                                              ztr)

        chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                        1:chunk(level)%sub_tile_dims(2)) = - ztr
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDIF
  ENDIF

END SUBROUTINE tea_leaf_dpcg_restrict_ZT_level

SUBROUTINE tea_leaf_dpcg_prolong_Z_level(level)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t
  INTEGER :: sub_tile_dx, sub_tile_dy
  INTEGER :: xrem, yrem

  REAL(KIND=8), DIMENSION(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: t2_local

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(sub_tile_dx,sub_tile_dy,t2_local)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1)/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1)/chunk(level)%sub_tile_dims(2)

      xrem        =mod(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1,chunk(level)%sub_tile_dims(1))
      yrem        =mod(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1,chunk(level)%sub_tile_dims(2))

      t2_local = chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                                 1:chunk(level)%sub_tile_dims(2))
      CALL tea_leaf_dpcg_prolong_Z_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                          chunk(level)%tiles(t)%field%x_max,    &
                                          chunk(level)%tiles(t)%field%y_min,    &
                                          chunk(level)%tiles(t)%field%y_max,    &
                                          chunk(level)%halo_exchange_depth,     &
                                          chunk(level)%sub_tile_dims(1),        &
                                          sub_tile_dx, xrem,                    &
                                          chunk(level)%sub_tile_dims(2),        &
                                          sub_tile_dy, yrem,                    &
                                          chunk(level)%tiles(t)%field%vector_z, &
                                          t2_local)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_prolong_Z_level

SUBROUTINE tea_leaf_dpcg_subtract_z_level(level)

  IMPLICIT NONE
  INTEGER :: level
  INTEGER :: t

  INTEGER :: sub_tile_dx, sub_tile_dy
  INTEGER :: xrem, yrem

  REAL(KIND=8), DIMENSION(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: t2_local

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(sub_tile_dx,sub_tile_dy,t2_local)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1)/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1)/chunk(level)%sub_tile_dims(2)

      xrem        =mod(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+1,chunk(level)%sub_tile_dims(1))
      yrem        =mod(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+1,chunk(level)%sub_tile_dims(2))

      t2_local=chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                               1:chunk(level)%sub_tile_dims(2))
      CALL tea_leaf_dpcg_subtract_z_kernel(chunk(level)%tiles(t)%field%x_min, &
                                           chunk(level)%tiles(t)%field%x_max, &
                                           chunk(level)%tiles(t)%field%y_min, &
                                           chunk(level)%tiles(t)%field%y_max, &
                                           chunk(level)%halo_exchange_depth,  &
                                           chunk(level)%sub_tile_dims(1),     &
                                           sub_tile_dx, xrem,                 &
                                           chunk(level)%sub_tile_dims(2),     &
                                           sub_tile_dy, yrem,                 &
                                           chunk(level)%tiles(t)%field%u,     &
                                           t2_local)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_subtract_z_level

SUBROUTINE tea_leaf_dpcg_init_p(level)

  IMPLICIT NONE

  INTEGER :: t, level

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_init_p_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                       chunk(level)%tiles(t)%field%x_max,    &
                                       chunk(level)%tiles(t)%field%y_min,    &
                                       chunk(level)%tiles(t)%field%y_max,    &
                                       chunk(level)%halo_exchange_depth,     &
                                       chunk(level)%tiles(t)%field%vector_p, &
                                       chunk(level)%tiles(t)%field%vector_z)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_init_p

SUBROUTINE tea_leaf_dpcg_store_r(level)

  IMPLICIT NONE
  INTEGER :: t, level

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_store_r_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                        chunk(level)%tiles(t)%field%x_max,    &
                                        chunk(level)%tiles(t)%field%y_min,    &
                                        chunk(level)%tiles(t)%field%y_max,    &
                                        chunk(level)%halo_exchange_depth,     &
                                        chunk(level)%tiles(t)%field%vector_r, &
                                        chunk(level)%tiles(t)%field%vector_r_m1 )
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_store_r

SUBROUTINE tea_leaf_dpcg_calc_rrn(level, rrn)

  IMPLICIT NONE
  INTEGER :: t, level
  REAL(KIND=8) :: rrn, tile_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rrn)
!$OMP DO REDUCTION(+:rrn)
    DO t=1,tiles_per_task
      tile_rrn = 0.0_8

      CALL tea_leaf_dpcg_calc_rrn_kernel(chunk(level)%tiles(t)%field%x_min,       &
                                         chunk(level)%tiles(t)%field%x_max,       &
                                         chunk(level)%tiles(t)%field%y_min,       &
                                         chunk(level)%tiles(t)%field%y_max,       &
                                         chunk(level)%halo_exchange_depth,        &
                                         chunk(level)%tiles(t)%field%vector_r,    &
                                         chunk(level)%tiles(t)%field%vector_r_m1, &
                                         chunk(level)%tiles(t)%field%vector_z,    &
                                         tile_rrn)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_rrn

SUBROUTINE tea_leaf_dpcg_calc_p(level,beta)

  IMPLICIT NONE
  INTEGER :: t, level
  REAL(KIND=8) :: beta

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_calc_p_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                       chunk(level)%tiles(t)%field%x_max,    &
                                       chunk(level)%tiles(t)%field%y_min,    &
                                       chunk(level)%tiles(t)%field%y_max,    &
                                       chunk(level)%halo_exchange_depth,     &
                                       chunk(level)%tiles(t)%field%vector_p, &
                                       chunk(level)%tiles(t)%field%vector_z, &
                                       beta)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_p

SUBROUTINE tea_leaf_dpcg_local_solve_level(level,               &
                                           solve_time,          &
                                           eps_inner,           &
                                           inner_iters,         &
                                           n,                   &
                                           use_ppcg,            &
                                           inner_cg_alphas,     &
                                           inner_cg_betas,      &
                                           inner_ch_theta,      &
                                           inner_ch_alphas,     &
                                           inner_ch_betas)

  USE tea_leaf_cg_module
  USE tea_leaf_ppcg_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_THREAD_NUM
  LOGICAL :: use_ppcg
  INTEGER :: level, inner_iters
  INTEGER :: t,ppcg_inner_iters,n
  INTEGER :: fields(NUM_FIELDS)
  REAL(KIND=8) :: solve_time,eps_inner
  REAL(KIND=8), DIMENSION(inner_iters) :: inner_cg_alphas, inner_cg_betas
  REAL(KIND=8), DIMENSION(inner_iters) :: inner_ch_alphas, inner_ch_betas
  REAL(KIND=8) :: inner_ch_theta

  REAL(KIND=8) :: alpha, beta
  REAL(KIND=8) :: timer,halo_time,dot_product_time
  REAL(KIND=8) :: rro,initial_residual,pw,rrn,error,old_error
  REAL(KIND=8) :: normu,normz

!$OMP PARALLEL
!$OMP DO
  DO t=1,tiles_per_task
    chunk(level+1)%tiles(t)%field%u0 = chunk(level+1)%tiles(t)%field%u
    chunk(level+1)%tiles(t)%field%u  = 0.0_8
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  fields=0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(level+1, fields, 1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  CALL tea_leaf_calc_residual(level+1)
  CALL tea_leaf_calc_2norm(level+1, 1, initial_residual)

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(initial_residual)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  !write(6,*) "rnorm:",sqrt(initial_residual)

  ! All 3 of these solvers use the CG kernels
  CALL tea_leaf_cg_init(level+1, ppcg_inner_iters, inner_ch_alphas, inner_ch_betas, inner_ch_theta, solve_time, 1, rro)

  IF (use_ppcg) THEN
    CALL tea_leaf_run_ppcg_inner_steps(level+1, inner_ch_alphas, inner_ch_betas, inner_ch_theta, &
        tl_ppcg_inner_coarse_steps, solve_time)
    ppcg_inner_iters = ppcg_inner_iters + tl_ppcg_inner_coarse_steps

    IF (tl_ppcg_inner_coarse_steps >= 0) THEN
      CALL tea_leaf_ppcg_pupdate    (level+1)
      CALL tea_leaf_ppcg_calc_zrnorm(level+1, rro)
    ENDIF
  ENDIF

  ! and globally sum rro
  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(rro)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  initial_residual = rro
  initial_residual=SQRT(abs(initial_residual))
  old_error = initial_residual

  IF (parallel%boss.AND.verbose_on) THEN
!$  IF (OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(g_out,*)"Coarse solve - Initial residual ",initial_residual
!$  ENDIF
  ENDIF

  ! need to update p when using CG due to matrix/vector multiplication
  fields=0
  fields(FIELD_U) = 1
  fields(FIELD_P) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(level+1,fields,1)
  IF (profiler_on) solve_time=solve_time+(timer()-halo_time)

  fields=0
  fields(FIELD_P) = 1

  ppcg_inner_iters = 0
  DO n=1,inner_iters

    CALL tea_leaf_cg_calc_w(level+1, pw)
    !normu = 0.0_8; normz = 0.0_8
    !DO t=1,tiles_per_task
    !  normu = normu+sum(chunk(level+1)%tiles(t)%field%vector_w(1:chunk(level)%sub_tile_dims(1), &
    !                                                           1:chunk(level)%sub_tile_dims(2))**2)
    !  normz = normz+sum(chunk(level+1)%tiles(t)%field%vector_p(1:chunk(level)%sub_tile_dims(1), &
    !                                                           1:chunk(level)%sub_tile_dims(2))**2)
    !ENDDO
    !write(6,*) "normw,normp:",sqrt(normu),sqrt(normz)

    IF (profiler_on) dot_product_time=timer()
    CALL tea_allsum(pw)
    IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

    alpha = rro/pw
    if (n <= tl_ch_cg_presteps) inner_cg_alphas(n) = alpha
    !write(g_out,*) "alpha:",alpha,rro,pw
    if  (alpha /= alpha) then
      write(6,*) "alpha not finite:",n,alpha,rro,pw
      stop
    endif

    CALL tea_leaf_cg_calc_ur(level+1, alpha, rrn)
    IF (n == inner_iters) EXIT ! this is the last solution update, so no point doing any extra work

    ! not using rrn, so don't do a tea_allsum

    IF (use_ppcg) THEN
      !write(6,*) "tea_leaf_run_ppcg_inner_steps:",tl_ppcg_inner_coarse_steps,level+1, &
      !    maxval(abs(inner_ch_alphas(1:tl_ppcg_inner_coarse_steps))), &
      !    maxval(abs(inner_ch_betas (1:tl_ppcg_inner_coarse_steps))), inner_ch_theta
      CALL tea_leaf_run_ppcg_inner_steps(level+1, inner_ch_alphas, inner_ch_betas, inner_ch_theta, &
          tl_ppcg_inner_coarse_steps, solve_time)
      ppcg_inner_iters = ppcg_inner_iters + tl_ppcg_inner_coarse_steps
    ENDIF

    CALL tea_leaf_ppcg_calc_zrnorm(level+1, rrn)

    IF (profiler_on) dot_product_time=timer()
    CALL tea_allsum(rrn)
    IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

    beta = rrn/rro
    if (n <= tl_ch_cg_presteps) inner_cg_betas(n) = beta

    CALL tea_leaf_cg_calc_p(level+1, beta)

    ! updates u and possibly p
    IF (profiler_on) halo_time = timer()
    CALL update_halo(level+1,fields,1)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    error = rrn
    rro = rrn

    error=SQRT(abs(error))

    IF (parallel%boss.AND.verbose_on) THEN
!$    IF (OMP_GET_THREAD_NUM().EQ.0) THEN
        WRITE(g_out,*)"Coarse solve - Residual         ",error
!$    ENDIF
    ENDIF

    !write(6,*) "level  solve:",error,eps_inner,initial_residual
    IF (ABS(error) .LT. eps_inner*initial_residual) EXIT

    old_error = error

  ENDDO

END SUBROUTINE tea_leaf_dpcg_local_solve_level

END MODULE tea_leaf_dpcg_module

