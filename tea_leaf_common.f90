
MODULE tea_leaf_common_module

  USE tea_leaf_common_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_init_common()

  IMPLICIT NONE

  INTEGER :: t

  INTEGER :: zero_boundary(4)

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      chunk%tiles(t)%field%rx = dt/(chunk%tiles(t)%field%celldx(chunk%tiles(t)%field%x_min)**2)
      chunk%tiles(t)%field%ry = dt/(chunk%tiles(t)%field%celldy(chunk%tiles(t)%field%y_min)**2)

      ! CG never needs matrix defined outside of boundaries, PPCG does
      IF (tl_use_cg) THEN
        zero_boundary = chunk%tiles(t)%tile_neighbours
      ELSE
        zero_boundary = chunk%chunk_neighbours
      ENDIF

      CALL tea_leaf_common_init_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%chunk_neighbours,                             &
          zero_boundary,                               &
          reflective_boundary,                                    &
          chunk%tiles(t)%field%density,                                &
          chunk%tiles(t)%field%energy1,                                &
          chunk%tiles(t)%field%u,                                      &
          chunk%tiles(t)%field%u0,                                      &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_w,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%vector_Mi,    &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          tl_preconditioner_type, coefficient)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_init_common

SUBROUTINE tea_leaf_calc_residual()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_calc_residual_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                        &
          chunk%tiles(t)%field%y_min,                        &
          chunk%tiles(t)%field%y_max,                        &
          halo_exchange_depth,                        &
          chunk%tiles(t)%field%u,                            &
          chunk%tiles(t)%field%u0,                           &
          chunk%tiles(t)%field%vector_r,                     &
          chunk%tiles(t)%field%vector_Kx,                    &
          chunk%tiles(t)%field%vector_Ky,                    &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry)
    ENDDO
  ENDIF

END SUBROUTINE

SUBROUTINE tea_leaf_calc_2norm(norm)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: norm, private_norm

  norm = 0.0_8

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      private_norm = 0.0_8

      CALL tea_leaf_calc_2norm_kernel(chunk%tiles(t)%field%x_min,        &
          chunk%tiles(t)%field%x_max,                                    &
          chunk%tiles(t)%field%y_min,                                    &
          chunk%tiles(t)%field%y_max,                                    &
          halo_exchange_depth,                                    &
          chunk%tiles(t)%field%vector_r,                                 &
          private_norm)

      norm = norm + private_norm
    ENDDO
  ENDIF

END SUBROUTINE

SUBROUTINE tea_leaf_finalise()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_finalise(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                           &
          chunk%tiles(t)%field%y_min,                           &
          chunk%tiles(t)%field%y_max,                           &
          halo_exchange_depth,                           &
          chunk%tiles(t)%field%energy1,                         &
          chunk%tiles(t)%field%density,                         &
          chunk%tiles(t)%field%u)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_finalise

END MODULE

