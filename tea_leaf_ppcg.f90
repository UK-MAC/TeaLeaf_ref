
MODULE tea_leaf_ppcg_module

  USE tea_leaf_ppcg_kernel_module
  USE tea_leaf_cheby_module
  USE definitions_module
  USE update_halo_module
  
  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_ppcg_init_sd(theta)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: theta

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_ppcg_init_sd(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                              &
          chunk%tiles(t)%field%y_min,                              &
          chunk%tiles(t)%field%y_max,                              &
          halo_exchange_depth,                              &
          chunk%tiles(t)%field%vector_r,                           &
	  chunk%tiles(t)%field%vector_rtemp,                           &          
          chunk%tiles(t)%field%vector_Kx,                        &
          chunk%tiles(t)%field%vector_Ky,                        &
          chunk%tiles(t)%field%vector_sd,                          &
          chunk%tiles(t)%field%vector_z,                          &
          chunk%tiles(t)%field%vector_utemp,                          &          
          chunk%tiles(t)%field%tri_cp,                          &
          chunk%tiles(t)%field%tri_bfp,                          &
          chunk%tiles(t)%field%vector_Mi,                          &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,                          &
          theta, tl_preconditioner_type)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_ppcg_init_sd

! This routine is needed to initialise PPCG
! Move this to the ppcg kernel!

SUBROUTINE tea_leaf_ppcg_init(rro,ch_alphas,ch_betas,ppcg_inner_iters,theta,solve_time,step)

  IMPLICIT NONE

  INTEGER :: t, ppcg_inner_iters,step
  REAL(KIND=8) :: rro, tile_rro, theta, solve_time
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas 
  rro = 0.0_8

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      tile_rro = 0.0_8

      CALL tea_leaf_ppcg_init_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%tiles(t)%field%vector_p,                               &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_Mi,                              &
          chunk%tiles(t)%field%vector_z,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          tile_rro, tl_preconditioner_type,	&
          ppcg_inner_iters, ch_alphas, ch_betas, theta, solve_time,step)

      rro = rro + tile_rro
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_ppcg_init

SUBROUTINE tea_leaf_ppcg_calc_zrnorm(rrn)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: rrn, tile_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      tile_rrn = 0.0_8

      CALL tea_leaf_ppcg_calc_zrnorm_kernel(chunk%tiles(t)%field%x_min, &
            chunk%tiles(t)%field%x_max,                           &
            chunk%tiles(t)%field%y_min,                           &
            chunk%tiles(t)%field%y_max,                           &
            halo_exchange_depth,                           &
            chunk%tiles(t)%field%vector_z,                        &
            chunk%tiles(t)%field%vector_r,                        &
            tl_preconditioner_type, tile_rrn)

      rrn = rrn + tile_rrn
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_ppcg_calc_zrnorm

SUBROUTINE tea_calc_ls_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                             theta, ppcg_inner_steps)

  INTEGER :: ppcg_inner_steps
  REAL(KIND=8), DIMENSION(ppcg_inner_steps) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax, theta

  ! TODO. What is to do? JDS
  CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                         theta, ppcg_inner_steps)

END SUBROUTINE

! A new routine needed to update z in ppcg
SUBROUTINE tea_leaf_ppcg_update_z()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task


      CALL tea_leaf_ppcg_update_z_kernel(chunk%tiles(t)%field%x_min, &
					 chunk%tiles(t)%field%x_max, &
					 chunk%tiles(t)%field%y_min, &
					 chunk%tiles(t)%field%y_max, &
					 halo_exchange_depth,        &
					 chunk%tiles(t)%field%vector_z, &
					 chunk%tiles(t)%field%vector_utemp)

    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_ppcg_update_z

! This has been moved from tea_solve and been adapted slightly
SUBROUTINE tea_leaf_run_ppcg_inner_steps(ch_alphas, ch_betas, theta, &
    tl_ppcg_inner_steps, solve_time)

  IMPLICIT NONE

  INTEGER :: fields(NUM_FIELDS)
  INTEGER :: tl_ppcg_inner_steps, ppcg_cur_step
  REAL(KIND=8) :: theta
  REAL(KIND=8) :: halo_time, timer, solve_time
  INTEGER :: t, inner_step, bounds_extra
  INTEGER :: x_min_bound, x_max_bound, y_min_bound, y_max_bound
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas


  !INTEGER(KIND=4) :: inner_step, bounds_extra

  fields = 0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(fields,1)
  IF (profiler_on) solve_time = solve_time + (timer() - halo_time)

  CALL tea_leaf_ppcg_init_sd(theta)

  ! inner steps
  DO ppcg_cur_step=1,tl_ppcg_inner_steps,halo_exchange_depth

    fields = 0
    fields(FIELD_SD) = 1
    fields(FIELD_R) = 1

    IF (profiler_on) halo_time = timer()
    CALL update_halo(fields,halo_exchange_depth)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    inner_step = ppcg_cur_step

    fields = 0
    fields(FIELD_SD) = 1

    DO bounds_extra = halo_exchange_depth-1, 0, -1
      IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      IF (chunk%chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
        x_min_bound = chunk%tiles(t)%field%x_min
      ELSE
        x_min_bound = chunk%tiles(t)%field%x_min - bounds_extra
      ENDIF

      IF (chunk%chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
        x_max_bound = chunk%tiles(t)%field%x_max
      ELSE
        x_max_bound = chunk%tiles(t)%field%x_max + bounds_extra
      ENDIF

      IF (chunk%chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
        y_min_bound = chunk%tiles(t)%field%y_min
      ELSE
        y_min_bound = chunk%tiles(t)%field%y_min - bounds_extra
      ENDIF

      IF (chunk%chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
        y_max_bound = chunk%tiles(t)%field%y_max
      ELSE
        y_max_bound = chunk%tiles(t)%field%y_max + bounds_extra
      ENDIF

      CALL tea_leaf_kernel_ppcg_inner(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                            &
          chunk%tiles(t)%field%y_min,                            &
          chunk%tiles(t)%field%y_max,                            &
          halo_exchange_depth,                            &
          x_min_bound,                                    &
          x_max_bound,                                    &
          y_min_bound,                                    &
          y_max_bound,                                    &
          ch_alphas, ch_betas,                              &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,                                           &
          inner_step,                                     &
          chunk%tiles(t)%field%u,                                &
          chunk%tiles(t)%field%vector_r,                         &
          chunk%tiles(t)%field%vector_rtemp,                         &
          chunk%tiles(t)%field%vector_Kx,                        &
          chunk%tiles(t)%field%vector_Ky,                        &
          chunk%tiles(t)%field%vector_sd,                        &
          chunk%tiles(t)%field%vector_z,                          &
          chunk%tiles(t)%field%vector_utemp,                          &
          chunk%tiles(t)%field%tri_cp,                          &
          chunk%tiles(t)%field%tri_bfp,                          &
          chunk%tiles(t)%field%vector_Mi,                          &
          tl_preconditioner_type)
    ENDDO
  ENDIF
    ENDDO
  ENDDO

  fields = 0
  fields(FIELD_P) = 1
  
  ! Then we update z = utemp
  CALL tea_leaf_ppcg_update_z()
  
END SUBROUTINE tea_leaf_run_ppcg_inner_steps

! new routine to store the previous residual for use with FCG(1)

SUBROUTINE tea_leaf_ppcg_store_r()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task

      CALL tea_leaf_ppcg_store_r_kernel(chunk%tiles(t)%field%x_min,	&
					chunk%tiles(t)%field%x_max,	&
					chunk%tiles(t)%field%y_min,	&
					chunk%tiles(t)%field%y_max,	&
					halo_exchange_depth,		&
					chunk%tiles(t)%field%vector_r,	&
					chunk%tiles(t)%field%vector_rstore)

    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_ppcg_store_r

! New routine which does the residual compute for FCG(1)

SUBROUTINE tea_leaf_ppcg_calc_rrn(rrn)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: rrn, tile_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rrn)
!$OMP DO REDUCTION(+:rrn)
    DO t=1,tiles_per_task
      tile_rrn = 0.0_8

      CALL tea_leaf_ppcg_calc_rrn_kernel(chunk%tiles(t)%field%x_min,		&
					    chunk%tiles(t)%field%x_max,         &
					    chunk%tiles(t)%field%y_min,         &
					    chunk%tiles(t)%field%y_max,         &
					    halo_exchange_depth,                &
					    chunk%tiles(t)%field%vector_r,      &
					    chunk%tiles(t)%field%vector_rstore, &
					    chunk%tiles(t)%field%vector_z,      &
					    tile_rrn)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_calc_rrn

END MODULE tea_leaf_ppcg_module

