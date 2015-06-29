
MODULE tea_leaf_cheby_module

  USE tea_leaf_cheby_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_cheby_init(rx, ry, ch_alphas, ch_betas, max_cheby_iters, theta)

  IMPLICIT NONE

  INTEGER :: t, max_cheby_iters
  REAL(KIND=8) :: ry, rx, theta
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_cheby_init(chunk%tiles(t)%field%x_min,&
            chunk%tiles(t)%field%x_max,                          &
            chunk%tiles(t)%field%y_min,                          &
            chunk%tiles(t)%field%y_max,                          &
            halo_exchange_depth,                          &
            chunk%tiles(t)%field%u,                              &
            chunk%tiles(t)%field%u0,                             &
            chunk%tiles(t)%field%vector_p,                       &
            chunk%tiles(t)%field%vector_r,                       &
            chunk%tiles(t)%field%vector_Mi,                      &
            chunk%tiles(t)%field%vector_w,                       &
            chunk%tiles(t)%field%vector_z,                       &
            chunk%tiles(t)%field%vector_Kx,                      &
            chunk%tiles(t)%field%vector_Ky,                      &
            chunk%tiles(t)%field%tri_cp,   &
            chunk%tiles(t)%field%tri_bfp,    &
            ch_alphas, ch_betas, max_cheby_iters,           &
            rx, ry, theta, tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cheby_init

SUBROUTINE tea_leaf_cheby_iterate(rx, ry, ch_alphas, ch_betas, max_cheby_iters, cheby_calc_steps)

  IMPLICIT NONE

  INTEGER :: t, cheby_calc_steps, max_cheby_iters
  REAL(KIND=8) :: ry, rx
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_cheby_iterate(chunk%tiles(t)%field%x_min,&
                  chunk%tiles(t)%field%x_max,                       &
                  chunk%tiles(t)%field%y_min,                       &
                  chunk%tiles(t)%field%y_max,                       &
                  halo_exchange_depth,                       &
                  chunk%tiles(t)%field%u,                           &
                  chunk%tiles(t)%field%u0,                          &
                  chunk%tiles(t)%field%vector_p,                    &
                  chunk%tiles(t)%field%vector_r,                    &
                  chunk%tiles(t)%field%vector_Mi,                   &
                  chunk%tiles(t)%field%vector_w,                    &
                  chunk%tiles(t)%field%vector_z,                    &
                  chunk%tiles(t)%field%vector_Kx,                   &
                  chunk%tiles(t)%field%vector_Ky,                   &
                  chunk%tiles(t)%field%tri_cp,   &
                  chunk%tiles(t)%field%tri_bfp,    &
                  ch_alphas, ch_betas, max_cheby_iters,        &
                  rx, ry, cheby_calc_steps, tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cheby_iterate

END MODULE tea_leaf_cheby_module

