!Crown Copyright 2014 AWE.
!
! This file is part of TeaLeaf.
!
! TeaLeaf is free software: you can redistribute it and/or modIFy it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! TeaLeaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! TeaLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran heat conduction kernel
!>  @author Michael Boulton, Wayne Gaudin, Douglas Shanks
!>  @details Implicitly calculates the change in temperature using the PPCG method

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
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_ppcg_init_sd(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                              &
          chunk%tiles(t)%field%y_min,                              &
          chunk%tiles(t)%field%y_max,                              &
          chunk%halo_exchange_depth,                               &
          chunk%tiles(t)%field%vector_r,                           &
          chunk%tiles(t)%field%vector_Kx,                        &
          chunk%tiles(t)%field%vector_Ky,                        &
          chunk%tiles(t)%field%vector_Di,                        &
          chunk%tiles(t)%field%vector_sd,                          &
          chunk%tiles(t)%field%vector_z,                          &
          chunk%tiles(t)%field%vector_utemp,                          &
          chunk%tiles(t)%field%vector_rtemp,                          &
          chunk%tiles(t)%field%tri_cp,                          &
          chunk%tiles(t)%field%tri_bfp,                          &
          chunk%tiles(t)%field%vector_Mi,                          &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,                          &
          theta, tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_init_sd

! This routine is needed to initialise PPCG

SUBROUTINE tea_leaf_ppcg_init(ppcg_inner_iters, ch_alphas, ch_betas, theta, solve_time, step, rro)

  IMPLICIT NONE

  INTEGER :: ppcg_inner_iters,step
  REAL(KIND=8) :: rro,theta,solve_time
  REAL(KIND=8), DIMENSION(:) :: ch_alphas,ch_betas

  INTEGER :: t
  REAL(KIND=8) :: tile_rro

  rro = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rro) REDUCTION(+:rro)
!$OMP DO
    DO t=1,tiles_per_task
      tile_rro = 0.0_8

      CALL tea_leaf_ppcg_init_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          chunk%halo_exchange_depth,                                   &
          chunk%tiles(t)%field%vector_p,                               &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_Mi,                              &
          chunk%tiles(t)%field%vector_z,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%vector_Di,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          tile_rro, tl_preconditioner_type, &
          ppcg_inner_iters, ch_alphas, ch_betas, theta, solve_time, step)

      rro = rro + tile_rro
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_init

SUBROUTINE tea_leaf_ppcg_calc_zrnorm(rrn)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: rrn, tile_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rrn)
!$OMP DO REDUCTION(+:rrn)
    DO t=1,tiles_per_task
      tile_rrn = 0.0_8

      CALL tea_leaf_ppcg_calc_zrnorm_kernel(chunk%tiles(t)%field%x_min, &
            chunk%tiles(t)%field%x_max,                           &
            chunk%tiles(t)%field%y_min,                           &
            chunk%tiles(t)%field%y_max,                           &
            chunk%halo_exchange_depth,                            &
            chunk%tiles(t)%field%vector_z,                        &
            chunk%tiles(t)%field%vector_r,                        &
            tl_preconditioner_type, tile_rrn)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_calc_zrnorm

SUBROUTINE tea_leaf_ppcg_update_z()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_ppcg_update_z_kernel(chunk%tiles(t)%field%x_min, &
            chunk%tiles(t)%field%x_max,                           &
            chunk%tiles(t)%field%y_min,                           &
            chunk%tiles(t)%field%y_max,                           &
            chunk%halo_exchange_depth,                            &
            chunk%tiles(t)%field%vector_z,                        &
            chunk%tiles(t)%field%vector_utemp)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_update_z

SUBROUTINE tea_leaf_ppcg_pupdate

  IMPLICIT NONE

  INTEGER      :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_ppcg_pupdate_kernel(chunk%tiles(t)%field%x_min, &
            chunk%tiles(t)%field%x_max,                           &
            chunk%tiles(t)%field%y_min,                           &
            chunk%tiles(t)%field%y_max,                           &
            chunk%halo_exchange_depth,                            &
            chunk%tiles(t)%field%vector_z,                        &
            chunk%tiles(t)%field%vector_p)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_pupdate

SUBROUTINE tea_calc_ls_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                             theta, ppcg_inner_steps)

  INTEGER :: ppcg_inner_steps
  REAL(KIND=8), DIMENSION(ppcg_inner_steps) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax, theta

  ! TODO
  CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                         theta, ppcg_inner_steps)

END SUBROUTINE

SUBROUTINE tea_leaf_run_ppcg_inner_steps(ch_alphas, ch_betas, theta, &
    ppcg_inner_steps, solve_time)

  IMPLICIT NONE

  INTEGER :: fields(NUM_FIELDS)
  INTEGER :: ppcg_inner_steps, ppcg_cur_step
  REAL(KIND=8) :: theta
  REAL(KIND=8) :: halo_time, timer, solve_time
  REAL(KIND=8), DIMENSION(max(ppcg_inner_steps,0)) :: ch_alphas, ch_betas

  INTEGER :: t, inner_step, bounds_extra
  INTEGER :: x_min_bound, x_max_bound, y_min_bound, y_max_bound
!$  INTEGER :: OMP_GET_THREAD_NUM

  IF (ppcg_inner_steps < 0) RETURN

  CALL tea_leaf_ppcg_init_sd(theta)

  ! inner steps
  DO ppcg_cur_step=1,ppcg_inner_steps,chunk%halo_exchange_depth

    fields = 0
    fields(FIELD_SD) = 1
    fields(FIELD_R) = 1

    IF (profiler_on) halo_time = timer()
    CALL update_halo(fields,chunk%halo_exchange_depth)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    inner_step = ppcg_cur_step

    fields = 0
    fields(FIELD_SD) = 1

    IF (reflective_boundary) then

!$OMP PARALLEL PRIVATE(x_min_bound, x_max_bound, y_min_bound, y_max_bound)
!$OMP DO
    DO t=1,tiles_per_task

    DO bounds_extra = chunk%halo_exchange_depth-1, 0, -1

    IF (use_fortran_kernels) THEN
      IF (chunk%tiles(t)%tile_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE .OR. &
          chunk%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
        x_min_bound = chunk%tiles(t)%field%x_min - bounds_extra
      ELSE
        x_min_bound = chunk%tiles(t)%field%x_min
      ENDIF

      IF (chunk%tiles(t)%tile_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE .OR. &
          chunk%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
        x_max_bound = chunk%tiles(t)%field%x_max + bounds_extra
      ELSE
        x_max_bound = chunk%tiles(t)%field%x_max
      ENDIF

      IF (chunk%tiles(t)%tile_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE .OR. &
          chunk%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
        y_min_bound = chunk%tiles(t)%field%y_min - bounds_extra
      ELSE
        y_min_bound = chunk%tiles(t)%field%y_min
      ENDIF

      IF (chunk%tiles(t)%tile_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE .OR. &
          chunk%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
        y_max_bound = chunk%tiles(t)%field%y_max + bounds_extra
      ELSE
        y_max_bound = chunk%tiles(t)%field%y_max
      ENDIF

      CALL tea_leaf_kernel_ppcg_inner(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                            &
          chunk%tiles(t)%field%y_min,                            &
          chunk%tiles(t)%field%y_max,                            &
          chunk%halo_exchange_depth,                             &
          x_min_bound,                                    &
          x_max_bound,                                    &
          y_min_bound,                                    &
          y_max_bound,                                    &
          ch_alphas, ch_betas,                              &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,                                           &
          ppcg_cur_step + chunk%halo_exchange_depth-1 - bounds_extra, &
          chunk%tiles(t)%field%u,                                &
          chunk%tiles(t)%field%vector_r,                         &
          chunk%tiles(t)%field%vector_Kx,                        &
          chunk%tiles(t)%field%vector_Ky,                        &
          chunk%tiles(t)%field%vector_Di,                        &
          chunk%tiles(t)%field%vector_sd,                        &
          chunk%tiles(t)%field%vector_z,                          &
          chunk%tiles(t)%field%vector_utemp,                          &
          chunk%tiles(t)%field%vector_rtemp,                          &
          chunk%tiles(t)%field%tri_cp,                          &
          chunk%tiles(t)%field%tri_bfp,                          &
          chunk%tiles(t)%field%vector_Mi,                          &
          tl_preconditioner_type)

    ENDIF

!$    IF(OMP_GET_THREAD_NUM().EQ.0) THEN
      IF (profiler_on) halo_time = timer()
!$    ENDIF

    IF (reflective_boundary .EQV. .TRUE. .AND. ANY(chunk%chunk_neighbours .EQ. EXTERNAL_FACE)) THEN
      IF (use_fortran_kernels)THEN
        CALL update_halo_kernel(chunk%tiles(t)%field%x_min,          &
                                chunk%tiles(t)%field%x_max,          &
                                chunk%tiles(t)%field%y_min,          &
                                chunk%tiles(t)%field%y_max,          &
                                chunk%halo_exchange_depth,           &
                                chunk%chunk_neighbours,              &
                                chunk%tiles(t)%tile_neighbours,      &
                                chunk%tiles(t)%field%density,        &
                                chunk%tiles(t)%field%energy0,        &
                                chunk%tiles(t)%field%energy1,        &
                                chunk%tiles(t)%field%u,              &
                                chunk%tiles(t)%field%vector_p,       &
                                chunk%tiles(t)%field%vector_sd,      &
                                chunk%tiles(t)%field%vector_rtemp,      &
                                chunk%tiles(t)%field%vector_z,       &
                                chunk%tiles(t)%field%vector_kx,      &
                                chunk%tiles(t)%field%vector_ky,      &
                                chunk%tiles(t)%field%vector_di,      &
                                fields,                                     &
                                1)
      ENDIF
    ENDIF

!$    IF(OMP_GET_THREAD_NUM().EQ.0) THEN
      IF (profiler_on) profiler%halo_update = profiler%halo_update + (timer() - halo_time)
      IF (profiler_on) solve_time           = solve_time           + (timer() - halo_time)
!$    ENDIF

      IF (ppcg_cur_step + chunk%halo_exchange_depth-1 - bounds_extra .eq. ppcg_inner_steps) EXIT
      ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    else ! .not. reflective_boundary

!$OMP PARALLEL PRIVATE(x_min_bound, x_max_bound, y_min_bound, y_max_bound)
!$OMP DO
    DO t=1,tiles_per_task

    DO bounds_extra = chunk%halo_exchange_depth-1, 0, -1

  IF (use_fortran_kernels) THEN
      IF (chunk%tiles(t)%tile_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE .OR. &
          chunk%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
        x_min_bound = chunk%tiles(t)%field%x_min - bounds_extra
      ELSE
        x_min_bound = chunk%tiles(t)%field%x_min
      ENDIF

      IF (chunk%tiles(t)%tile_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE .OR. &
          chunk%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
        x_max_bound = chunk%tiles(t)%field%x_max + bounds_extra
      ELSE
        x_max_bound = chunk%tiles(t)%field%x_max
      ENDIF

      IF (chunk%tiles(t)%tile_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE .OR. &
          chunk%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
        y_min_bound = chunk%tiles(t)%field%y_min - bounds_extra
      ELSE
        y_min_bound = chunk%tiles(t)%field%y_min
      ENDIF

      IF (chunk%tiles(t)%tile_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE .OR. &
          chunk%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
        y_max_bound = chunk%tiles(t)%field%y_max + bounds_extra
      ELSE
        y_max_bound = chunk%tiles(t)%field%y_max
      ENDIF

      CALL tea_leaf_kernel_ppcg_inner(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                            &
          chunk%tiles(t)%field%y_min,                            &
          chunk%tiles(t)%field%y_max,                            &
          chunk%halo_exchange_depth,                             &
          x_min_bound,                                    &
          x_max_bound,                                    &
          y_min_bound,                                    &
          y_max_bound,                                    &
          ch_alphas, ch_betas,                              &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,                                           &
          ppcg_cur_step + chunk%halo_exchange_depth-1 - bounds_extra, &
          chunk%tiles(t)%field%u,                                &
          chunk%tiles(t)%field%vector_r,                         &
          chunk%tiles(t)%field%vector_Kx,                        &
          chunk%tiles(t)%field%vector_Ky,                        &
          chunk%tiles(t)%field%vector_Di,                        &
          chunk%tiles(t)%field%vector_sd,                        &
          chunk%tiles(t)%field%vector_z,                          &
          chunk%tiles(t)%field%vector_utemp,                          &
          chunk%tiles(t)%field%vector_rtemp,                          &
          chunk%tiles(t)%field%tri_cp,                          &
          chunk%tiles(t)%field%tri_bfp,                          &
          chunk%tiles(t)%field%vector_Mi,                          &
          tl_preconditioner_type)
  
  ENDIF

      IF (ppcg_cur_step + chunk%halo_exchange_depth-1 - bounds_extra .eq. ppcg_inner_steps) EXIT
    ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF !reflective_boundary
  ENDDO

  fields = 0
  fields(FIELD_P) = 1

  CALL tea_leaf_ppcg_update_z()

END SUBROUTINE tea_leaf_run_ppcg_inner_steps

SUBROUTINE tea_leaf_ppcg_store_r()

  IMPLICIT NONE
  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_ppcg_store_r_kernel(chunk%tiles(t)%field%x_min,    &
                                        chunk%tiles(t)%field%x_max,    &
                                        chunk%tiles(t)%field%y_min,    &
                                        chunk%tiles(t)%field%y_max,    &
                                        chunk%halo_exchange_depth,     &
                                        chunk%tiles(t)%field%vector_r, &
                                        chunk%tiles(t)%field%vector_r_store )
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_store_r

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

      CALL tea_leaf_ppcg_calc_rrn_kernel(chunk%tiles(t)%field%x_min,       &
                                         chunk%tiles(t)%field%x_max,       &
                                         chunk%tiles(t)%field%y_min,       &
                                         chunk%tiles(t)%field%y_max,       &
                                         chunk%halo_exchange_depth,        &
                                         chunk%tiles(t)%field%vector_r,    &
                                         chunk%tiles(t)%field%vector_r_store, &
                                         chunk%tiles(t)%field%vector_z,    &
                                         tile_rrn)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_calc_rrn

END MODULE tea_leaf_ppcg_module

