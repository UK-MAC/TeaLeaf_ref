
MODULE tea_leaf_cg_module

  USE tea_leaf_kernel_cg_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_init_cg(rx, ry, rro)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry, rx, rro

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_init_cg_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%tiles(t)%field%density,                                &
          chunk%tiles(t)%field%energy1,                                &
          chunk%tiles(t)%field%u,                                      &
          chunk%tiles(t)%field%vector_p,                               &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_Mi,                              &
          chunk%tiles(t)%field%vector_w,                               &
          chunk%tiles(t)%field%vector_z,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          rx, ry, rro, tl_preconditioner_type)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_init_cg

SUBROUTINE tea_leaf_cg_calc_w(rx, ry, pw)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry, rx, pw

  pw = 0.0_08

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_calc_w_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                         &
          chunk%tiles(t)%field%y_min,                                         &
          chunk%tiles(t)%field%y_max,                                         &
          halo_exchange_depth,                                         &
          chunk%tiles(t)%field%vector_p,                                      &
          chunk%tiles(t)%field%vector_w,                                      &
          chunk%tiles(t)%field%vector_Kx,                                     &
          chunk%tiles(t)%field%vector_Ky,                                     &
          rx, ry, pw)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_w

SUBROUTINE tea_leaf_cg_calc_ur(rx, ry, alpha, rrn)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry, rx, alpha, rrn

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_calc_ur_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                          &
          chunk%tiles(t)%field%y_min,                                          &
          chunk%tiles(t)%field%y_max,                                          &
          halo_exchange_depth,                                          &
          chunk%tiles(t)%field%u,                                              &
          chunk%tiles(t)%field%vector_p,                                       &
          chunk%tiles(t)%field%vector_r,                                       &
          chunk%tiles(t)%field%vector_Mi,                                      &
          chunk%tiles(t)%field%vector_w,                                       &
          chunk%tiles(t)%field%vector_z,                                       &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          rx, ry, &
          alpha, rrn, tl_preconditioner_type)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_ur

SUBROUTINE tea_leaf_cg_calc_p(rx, ry, beta)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry, rx, beta

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_calc_p_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                         &
          chunk%tiles(t)%field%y_min,                                         &
          chunk%tiles(t)%field%y_max,                                         &
          halo_exchange_depth,                                         &
          chunk%tiles(t)%field%vector_p,                                      &
          chunk%tiles(t)%field%vector_r,                                      &
          chunk%tiles(t)%field%vector_z,                                      &
          beta, tl_preconditioner_type)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_p

END MODULE tea_leaf_cg_module

