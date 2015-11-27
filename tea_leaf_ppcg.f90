
MODULE tea_leaf_ppcg_module

  USE tea_leaf_ppcg_kernel_module
  USE tea_leaf_cheby_module
  USE definitions_module
  USE update_halo_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_ppcg_init_sd(level, theta)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t
  REAL(KIND=8) :: theta

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_ppcg_init_sd(chunk(level)%tiles(t)%field%x_min,&
          chunk(level)%tiles(t)%field%x_max,                              &
          chunk(level)%tiles(t)%field%y_min,                              &
          chunk(level)%tiles(t)%field%y_max,                              &
          chunk(level)%halo_exchange_depth,                               &
          chunk(level)%tiles(t)%field%vector_r,                           &
          chunk(level)%tiles(t)%field%vector_Kx,                        &
          chunk(level)%tiles(t)%field%vector_Ky,                        &
          chunk(level)%tiles(t)%field%vector_Di,                        &
          chunk(level)%tiles(t)%field%vector_sd,                          &
          chunk(level)%tiles(t)%field%vector_z,                          &
          chunk(level)%tiles(t)%field%tri_cp,                          &
          chunk(level)%tiles(t)%field%tri_bfp,                          &
          chunk(level)%tiles(t)%field%vector_Mi,                          &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,                          &
          theta, tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_init_sd

SUBROUTINE tea_leaf_ppcg_inner(level, ch_alphas, ch_betas, inner_step, bounds_extra)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t, inner_step, bounds_extra
  INTEGER :: x_min_bound, x_max_bound, y_min_bound, y_max_bound
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(x_min_bound, x_max_bound, y_min_bound, y_max_bound)
!$OMP DO
    DO t=1,tiles_per_task
      IF (chunk(level)%tiles(t)%tile_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE .OR. &
          chunk(level)%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
        x_min_bound = chunk(level)%tiles(t)%field%x_min - bounds_extra
      ELSE
        x_min_bound = chunk(level)%tiles(t)%field%x_min
      ENDIF

      IF (chunk(level)%tiles(t)%tile_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE .OR. &
          chunk(level)%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
        x_max_bound = chunk(level)%tiles(t)%field%x_max + bounds_extra
      ELSE
        x_max_bound = chunk(level)%tiles(t)%field%x_max
      ENDIF

      IF (chunk(level)%tiles(t)%tile_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE .OR. &
          chunk(level)%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
        y_min_bound = chunk(level)%tiles(t)%field%y_min - bounds_extra
      ELSE
        y_min_bound = chunk(level)%tiles(t)%field%y_min
      ENDIF

      IF (chunk(level)%tiles(t)%tile_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE .OR. &
          chunk(level)%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
        y_max_bound = chunk(level)%tiles(t)%field%y_max + bounds_extra
      ELSE
        y_max_bound = chunk(level)%tiles(t)%field%y_max
      ENDIF

      CALL tea_leaf_kernel_ppcg_inner(chunk(level)%tiles(t)%field%x_min,&
          chunk(level)%tiles(t)%field%x_max,                            &
          chunk(level)%tiles(t)%field%y_min,                            &
          chunk(level)%tiles(t)%field%y_max,                            &
          chunk(level)%halo_exchange_depth,                             &
          x_min_bound,                                    &
          x_max_bound,                                    &
          y_min_bound,                                    &
          y_max_bound,                                    &
          ch_alphas, ch_betas,                              &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,                                           &
          inner_step,                                     &
          chunk(level)%tiles(t)%field%u,                                &
          chunk(level)%tiles(t)%field%vector_r,                         &
          chunk(level)%tiles(t)%field%vector_Kx,                        &
          chunk(level)%tiles(t)%field%vector_Ky,                        &
          chunk(level)%tiles(t)%field%vector_Di,                        &
          chunk(level)%tiles(t)%field%vector_sd,                        &
          chunk(level)%tiles(t)%field%vector_z,                          &
          chunk(level)%tiles(t)%field%tri_cp,                          &
          chunk(level)%tiles(t)%field%tri_bfp,                          &
          chunk(level)%tiles(t)%field%vector_Mi,                          &
          tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_inner

SUBROUTINE tea_leaf_ppcg_inner_nouup(level, ch_alphas, ch_betas, inner_step, bounds_extra)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t, inner_step, bounds_extra
  INTEGER :: x_min_bound, x_max_bound, y_min_bound, y_max_bound
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(x_min_bound, x_max_bound, y_min_bound, y_max_bound)
!$OMP DO
    DO t=1,tiles_per_task
      IF (chunk(level)%tiles(t)%tile_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE .OR. &
          chunk(level)%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
        x_min_bound = chunk(level)%tiles(t)%field%x_min - bounds_extra
      ELSE
        x_min_bound = chunk(level)%tiles(t)%field%x_min
      ENDIF

      IF (chunk(level)%tiles(t)%tile_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE .OR. &
          chunk(level)%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
        x_max_bound = chunk(level)%tiles(t)%field%x_max + bounds_extra
      ELSE
        x_max_bound = chunk(level)%tiles(t)%field%x_max
      ENDIF

      IF (chunk(level)%tiles(t)%tile_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE .OR. &
          chunk(level)%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
        y_min_bound = chunk(level)%tiles(t)%field%y_min - bounds_extra
      ELSE
        y_min_bound = chunk(level)%tiles(t)%field%y_min
      ENDIF

      IF (chunk(level)%tiles(t)%tile_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE .OR. &
          chunk(level)%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
        y_max_bound = chunk(level)%tiles(t)%field%y_max + bounds_extra
      ELSE
        y_max_bound = chunk(level)%tiles(t)%field%y_max
      ENDIF

      CALL tea_leaf_kernel_ppcg_inner_nouup(chunk(level)%tiles(t)%field%x_min,&
          chunk(level)%tiles(t)%field%x_max,                            &
          chunk(level)%tiles(t)%field%y_min,                            &
          chunk(level)%tiles(t)%field%y_max,                            &
          chunk(level)%halo_exchange_depth,                             &
          x_min_bound,                                    &
          x_max_bound,                                    &
          y_min_bound,                                    &
          y_max_bound,                                    &
          ch_alphas, ch_betas,                              &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,                                           &
          inner_step,                                     &
          chunk(level)%tiles(t)%field%u,                                &
          chunk(level)%tiles(t)%field%vector_r,                         &
          chunk(level)%tiles(t)%field%vector_Kx,                        &
          chunk(level)%tiles(t)%field%vector_Ky,                        &
          chunk(level)%tiles(t)%field%vector_Di,                        &
          chunk(level)%tiles(t)%field%vector_sd,                        &
          chunk(level)%tiles(t)%field%vector_z,                          &
          chunk(level)%tiles(t)%field%tri_cp,                          &
          chunk(level)%tiles(t)%field%tri_bfp,                          &
          chunk(level)%tiles(t)%field%vector_Mi,                          &
          tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_inner_nouup

SUBROUTINE tea_leaf_ppcg_calc_zrnorm(level, rrn)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t
  REAL(KIND=8) :: rrn, tile_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rrn)
!$OMP DO REDUCTION(+:rrn)
    DO t=1,tiles_per_task
      tile_rrn = 0.0_8

      CALL tea_leaf_ppcg_calc_zrnorm_kernel(chunk(level)%tiles(t)%field%x_min, &
            chunk(level)%tiles(t)%field%x_max,                           &
            chunk(level)%tiles(t)%field%y_min,                           &
            chunk(level)%tiles(t)%field%y_max,                           &
            chunk(level)%halo_exchange_depth,                            &
            chunk(level)%tiles(t)%field%vector_z,                        &
            chunk(level)%tiles(t)%field%vector_r,                        &
            tl_preconditioner_type, tile_rrn)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_calc_zrnorm

SUBROUTINE tea_calc_ls_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                             theta, ppcg_inner_steps)

  INTEGER :: ppcg_inner_steps
  REAL(KIND=8), DIMENSION(ppcg_inner_steps) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax, theta

  ! TODO
  CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                         theta, ppcg_inner_steps)

END SUBROUTINE

SUBROUTINE tea_leaf_run_ppcg_inner_steps(level, ch_alphas, ch_betas, theta, &
    ppcg_inner_steps, solve_time)

  IMPLICIT NONE

  INTEGER :: level, fields(NUM_FIELDS)
  INTEGER :: ppcg_inner_steps, ppcg_cur_step
  REAL(KIND=8) :: theta
  REAL(KIND=8) :: halo_time, timer, solve_time
  REAL(KIND=8), DIMENSION(max(ppcg_inner_steps,0)) :: ch_alphas, ch_betas

  INTEGER(KIND=4) :: inner_step, bounds_extra

  IF (ppcg_inner_steps <= 0) RETURN

  !write(6,*) maxval(abs(ch_alphas)), maxval(abs(ch_betas)), theta, &
  !  ppcg_inner_steps, level
  !stop

  fields = 0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(level, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer() - halo_time)

  CALL tea_leaf_ppcg_init_sd(level, theta)

  ! inner steps
  DO ppcg_cur_step=1,ppcg_inner_steps,chunk(level)%halo_exchange_depth

    fields = 0
    fields(FIELD_SD) = 1
    fields(FIELD_R) = 1

    IF (profiler_on) halo_time = timer()
    CALL update_halo(level, fields,chunk(level)%halo_exchange_depth)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    inner_step = ppcg_cur_step

    fields = 0
    fields(FIELD_SD) = 1

    DO bounds_extra = chunk(level)%halo_exchange_depth-1, 0, -1
      CALL tea_leaf_ppcg_inner(level, ch_alphas, ch_betas, inner_step, bounds_extra)

      IF (profiler_on) halo_time = timer()
      CALL update_boundary(level, fields, 1)
      IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

      inner_step = inner_step + 1
      IF (inner_step .gt. ppcg_inner_steps) EXIT
    ENDDO
  ENDDO

  fields = 0
  fields(FIELD_P) = 1

END SUBROUTINE tea_leaf_run_ppcg_inner_steps

SUBROUTINE tea_leaf_run_ppcg_inner_steps_nouup(level, ch_alphas, ch_betas, theta, &
    ppcg_inner_steps, solve_time)

  IMPLICIT NONE

  INTEGER :: level, fields(NUM_FIELDS)
  INTEGER :: ppcg_inner_steps, ppcg_cur_step
  REAL(KIND=8) :: theta
  REAL(KIND=8) :: halo_time, timer, solve_time
  REAL(KIND=8), DIMENSION(max(ppcg_inner_steps,0)) :: ch_alphas, ch_betas

  INTEGER(KIND=4) :: inner_step, bounds_extra

  IF (ppcg_inner_steps <= 0) RETURN

  !write(6,*) maxval(abs(ch_alphas)), maxval(abs(ch_betas)), theta, &
  !  ppcg_inner_steps, level
  !stop

  fields = 0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(level, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer() - halo_time)

  CALL tea_leaf_ppcg_init_sd(level, theta)

  ! inner steps
  DO ppcg_cur_step=1,ppcg_inner_steps,chunk(level)%halo_exchange_depth

    fields = 0
    fields(FIELD_SD) = 1
    fields(FIELD_R) = 1

    IF (profiler_on) halo_time = timer()
    CALL update_halo(level, fields,chunk(level)%halo_exchange_depth)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    inner_step = ppcg_cur_step

    fields = 0
    fields(FIELD_SD) = 1

    DO bounds_extra = chunk(level)%halo_exchange_depth-1, 0, -1
      CALL tea_leaf_ppcg_inner_nouup(level, ch_alphas, ch_betas, inner_step, bounds_extra)

      IF (profiler_on) halo_time = timer()
      CALL update_boundary(level, fields, 1)
      IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

      inner_step = inner_step + 1
      IF (inner_step .gt. ppcg_inner_steps) EXIT
    ENDDO
  ENDDO

  fields = 0
  fields(FIELD_P) = 1

END SUBROUTINE tea_leaf_run_ppcg_inner_steps_nouup

END MODULE tea_leaf_ppcg_module

