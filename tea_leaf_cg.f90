
MODULE tea_leaf_cg_module

  USE tea_leaf_cg_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_cg_init(level, ppcg_inner_iters, ch_alphas, ch_betas, theta, solve_time, rro)

  IMPLICIT NONE

  INTEGER :: level,ppcg_inner_iters
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

      CALL tea_leaf_cg_init_kernel(chunk(level)%tiles(t)%field%x_min, &
          chunk(level)%tiles(t)%field%x_max,                                  &
          chunk(level)%tiles(t)%field%y_min,                                  &
          chunk(level)%tiles(t)%field%y_max,                                  &
          chunk(level)%halo_exchange_depth,                                   &
          chunk(level)%tiles(t)%field%vector_p,                               &
          chunk(level)%tiles(t)%field%vector_r,                               &
          chunk(level)%tiles(t)%field%vector_Mi,                              &
          chunk(level)%tiles(t)%field%vector_z,                               &
          chunk(level)%tiles(t)%field%vector_Kx,                              &
          chunk(level)%tiles(t)%field%vector_Ky,                              &
          chunk(level)%tiles(t)%field%vector_Di,                              &
          chunk(level)%tiles(t)%field%tri_cp,   &
          chunk(level)%tiles(t)%field%tri_bfp,    &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,  &
          tile_rro, tl_preconditioner_type, &
          level, ppcg_inner_iters, ch_alphas, ch_betas, theta, solve_time)

      !write(6,'("tile_rro:",i3,4es25.18)') t,tile_rro,sum(chunk(level)%tiles(t)%field%vector_r**2), &
      !                                                sum(chunk(level)%tiles(t)%field%vector_z**2), &
      !                                                sum(chunk(level)%tiles(t)%field%vector_p**2)
      !write(6,'("tile_rro:",i3,5es25.18)') t,tile_rro,sum(chunk(level)%tiles(t)%field%vector_Mi**2), &
      !                                                sum(chunk(level)%tiles(t)%field%vector_Kx**2), &
      !                                                sum(chunk(level)%tiles(t)%field%vector_Ky**2), &
      !                                                sum(chunk(level)%tiles(t)%field%vector_Di**2)
      rro = rro + tile_rro
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_init

SUBROUTINE tea_leaf_cg_calc_w(level, pw)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t
  REAL(KIND=8) :: pw, tile_pw

  pw = 0.0_08

  IF (level > 1) THEN
  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_pw) REDUCTION(+:pw)
!$OMP DO
    DO t=1,tiles_per_task
      tile_pw = 0.0_8

      CALL tea_leaf_cg_calc_w_kernel_norxy(chunk(level)%tiles(t)%field%x_min,&
          chunk(level)%tiles(t)%field%x_max,                                         &
          chunk(level)%tiles(t)%field%y_min,                                         &
          chunk(level)%tiles(t)%field%y_max,                                         &
          chunk(level)%halo_exchange_depth,                                          &
          chunk(level)%tiles(t)%field%vector_p,                                      &
          chunk(level)%tiles(t)%field%vector_w,                                      &
          chunk(level)%tiles(t)%field%vector_Kx,                                     &
          chunk(level)%tiles(t)%field%vector_Ky,                                     &
          chunk(level)%tiles(t)%field%vector_Di,                                     &
          tile_pw)

      pw = pw + tile_pw
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF
  ELSE
  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_pw) REDUCTION(+:pw)
!$OMP DO
    DO t=1,tiles_per_task
      tile_pw = 0.0_8

      CALL tea_leaf_cg_calc_w_kernel(chunk(level)%tiles(t)%field%x_min,&
          chunk(level)%tiles(t)%field%x_max,                                         &
          chunk(level)%tiles(t)%field%y_min,                                         &
          chunk(level)%tiles(t)%field%y_max,                                         &
          chunk(level)%halo_exchange_depth,                                          &
          chunk(level)%tiles(t)%field%vector_p,                                      &
          chunk(level)%tiles(t)%field%vector_w,                                      &
          chunk(level)%tiles(t)%field%vector_Kx,                                     &
          chunk(level)%tiles(t)%field%vector_Ky,                                     &
          chunk(level)%tiles(t)%field%vector_Di,                                     &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,  &
          tile_pw)

      pw = pw + tile_pw
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_w

SUBROUTINE tea_leaf_cg_calc_ur(level, alpha, rrn)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t
  REAL(KIND=8) :: alpha, rrn, tile_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rrn) REDUCTION(+:rrn)
!$OMP DO
    DO t=1,tiles_per_task
      tile_rrn = 0.0_8

      CALL tea_leaf_cg_calc_ur_kernel(chunk(level)%tiles(t)%field%x_min,&
          chunk(level)%tiles(t)%field%x_max,                                          &
          chunk(level)%tiles(t)%field%y_min,                                          &
          chunk(level)%tiles(t)%field%y_max,                                          &
          chunk(level)%halo_exchange_depth,                                           &
          chunk(level)%tiles(t)%field%u,                                              &
          chunk(level)%tiles(t)%field%vector_p,                                       &
          chunk(level)%tiles(t)%field%vector_r,                                       &
          chunk(level)%tiles(t)%field%vector_Mi,                                      &
          chunk(level)%tiles(t)%field%vector_w,                                       &
          chunk(level)%tiles(t)%field%vector_z,                                       &
          chunk(level)%tiles(t)%field%tri_cp,   &
          chunk(level)%tiles(t)%field%tri_bfp,    &
          chunk(level)%tiles(t)%field%vector_Kx,                              &
          chunk(level)%tiles(t)%field%vector_Ky,                              &
          chunk(level)%tiles(t)%field%vector_Di,                              &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry, &
          alpha, tile_rrn, tl_preconditioner_type)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_ur

SUBROUTINE tea_leaf_cg_calc_p(level, beta)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t
  REAL(KIND=8) :: beta

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_calc_p_kernel(chunk(level)%tiles(t)%field%x_min,&
          chunk(level)%tiles(t)%field%x_max,                                         &
          chunk(level)%tiles(t)%field%y_min,                                         &
          chunk(level)%tiles(t)%field%y_max,                                         &
          chunk(level)%halo_exchange_depth,                                          &
          chunk(level)%tiles(t)%field%vector_p,                                      &
          chunk(level)%tiles(t)%field%vector_r,                                      &
          chunk(level)%tiles(t)%field%vector_z,                                      &
          beta, tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_p

END MODULE tea_leaf_cg_module

