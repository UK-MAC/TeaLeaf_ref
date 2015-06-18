
MODULE tea_leaf_ppcg_module

  USE tea_leaf_ppcg_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_ppcg_init_sd(rx, ry, theta)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry,rx,theta

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_ppcg_init_sd(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                              &
          chunk%tiles(t)%field%y_min,                              &
          chunk%tiles(t)%field%y_max,                              &
          halo_exchange_depth,                              &
          chunk%tiles(t)%field%vector_r,                           &
          chunk%tiles(t)%field%vector_Kx,                        &
          chunk%tiles(t)%field%vector_Ky,                        &
          chunk%tiles(t)%field%vector_sd,                          &
          chunk%tiles(t)%field%vector_z,                          &
          chunk%tiles(t)%field%tri_cp,                          &
          chunk%tiles(t)%field%tri_bfp,                          &
          chunk%tiles(t)%field%vector_Mi,                          &
          rx, ry,                          &
          theta, tl_preconditioner_type)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_init_sd

SUBROUTINE tea_leaf_ppcg_inner(rx, ry, ch_alphas, ch_betas, inner_step, bounds_extra)

  IMPLICIT NONE

  INTEGER :: t, inner_step, bounds_extra
  REAL(KIND=8) :: ry,rx
  INTEGER :: x_min_bound, x_max_bound, y_min_bound, y_max_bound
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(x_min_bound, x_max_bound, y_min_bound, y_max_bound)
!$OMP DO
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
          rx, ry,                                           &
          inner_step,                                     &
          chunk%tiles(t)%field%u,                                &
          chunk%tiles(t)%field%vector_r,                         &
          chunk%tiles(t)%field%vector_Kx,                        &
          chunk%tiles(t)%field%vector_Ky,                        &
          chunk%tiles(t)%field%vector_sd,                        &
          chunk%tiles(t)%field%vector_z,                          &
          chunk%tiles(t)%field%tri_cp,                          &
          chunk%tiles(t)%field%tri_bfp,                          &
          chunk%tiles(t)%field%vector_Mi,                          &
          tl_preconditioner_type)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_inner

SUBROUTINE tea_leaf_ppcg_calc_zrnorm(rrn)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: rrn, private_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(private_rrn)
!$OMP DO REDUCTION(+:rrn)
    DO t=1,tiles_per_task
      private_rrn = 0.0_8

      CALL tea_leaf_ppcg_calc_zrnorm_kernel(chunk%tiles(t)%field%x_min, &
            chunk%tiles(t)%field%x_max,                           &
            chunk%tiles(t)%field%y_min,                           &
            chunk%tiles(t)%field%y_max,                           &
            halo_exchange_depth,                           &
            chunk%tiles(t)%field%vector_z,                        &
            chunk%tiles(t)%field%vector_r,                        &
            tl_preconditioner_type, private_rrn)

      rrn = rrn + private_rrn
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_ppcg_calc_zrnorm

END MODULE tea_leaf_ppcg_module
