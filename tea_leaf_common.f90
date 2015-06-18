
MODULE tea_leaf_common_module

  USE tea_leaf_common_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_init_common(rx, ry)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry,rx

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO LASTPRIVATE(rx, ry)
    DO t=1,tiles_per_task
      rx = dt/(chunk%tiles(t)%field%celldx(chunk%tiles(t)%field%x_min)**2)
      ry = dt/(chunk%tiles(t)%field%celldy(chunk%tiles(t)%field%y_min)**2)

      CALL tea_leaf_common_init_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%chunk_neighbours,                             &
          chunk%tiles(t)%tile_neighbours,                               &
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
          rx, ry, tl_preconditioner_type, coefficient)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_init_common

SUBROUTINE tea_leaf_calc_residual(rx, ry)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry,rx

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
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
          rx, ry)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE

SUBROUTINE tea_leaf_calc_2norm(norm)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: norm, private_norm

  norm = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(private_norm)
!$OMP DO REDUCTION(+:norm)
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
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE

SUBROUTINE tea_leaf_finalise()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
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
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_finalise

END MODULE

