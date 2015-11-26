
MODULE tea_leaf_common_module

  USE tea_leaf_common_kernel_module
  USE definitions_module
  USE report_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_init_common(level)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t

  LOGICAL :: zero_boundary(4)

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(zero_boundary)
!$OMP DO
    DO t=1,tiles_per_task
      chunk(level)%tiles(t)%field%rx = dt/(chunk(level)%tiles(t)%field%celldx(chunk(level)%tiles(t)%field%x_min)**2)
      chunk(level)%tiles(t)%field%ry = dt/(chunk(level)%tiles(t)%field%celldy(chunk(level)%tiles(t)%field%y_min)**2)

      WHERE (chunk(level)%tiles(t)%tile_neighbours .EQ. EXTERNAL_FACE .AND. &
             chunk(level)%chunk_neighbours .EQ. EXTERNAL_FACE)
        zero_boundary = .TRUE.
      ELSE WHERE
        zero_boundary = .FALSE.
      END WHERE

      CALL tea_leaf_common_init_kernel(chunk(level)%tiles(t)%field%x_min, &
          chunk(level)%tiles(t)%field%x_max,                                  &
          chunk(level)%tiles(t)%field%y_min,                                  &
          chunk(level)%tiles(t)%field%y_max,                                  &
          chunk(level)%halo_exchange_depth,                                   &
          zero_boundary,                               &
          reflective_boundary,                                    &
          chunk(level)%tiles(t)%field%density,                                &
          chunk(level)%tiles(t)%field%energy1,                                &
          chunk(level)%tiles(t)%field%u,                                      &
          chunk(level)%tiles(t)%field%u0,                                      &
          chunk(level)%tiles(t)%field%vector_r,                               &
          chunk(level)%tiles(t)%field%vector_w,                               &
          chunk(level)%tiles(t)%field%vector_Kx,                              &
          chunk(level)%tiles(t)%field%vector_Ky,                              &
          chunk(level)%tiles(t)%field%vector_Di,                              &
          chunk(level)%tiles(t)%field%tri_cp,   &
          chunk(level)%tiles(t)%field%tri_bfp,    &
          chunk(level)%tiles(t)%field%vector_Mi,    &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,  &
          tl_preconditioner_type, coefficient)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_init_common

SUBROUTINE tea_leaf_calc_residual(level)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_calc_residual_kernel(                           &
          chunk(level)%tiles(t)%field%x_min,                        &
          chunk(level)%tiles(t)%field%x_max,                        &
          chunk(level)%tiles(t)%field%y_min,                        &
          chunk(level)%tiles(t)%field%y_max,                        &
          chunk(level)%halo_exchange_depth,                         &
          chunk(level)%tiles(t)%field%u,                            &
          chunk(level)%tiles(t)%field%u0,                           &
          chunk(level)%tiles(t)%field%vector_r,                     &
          chunk(level)%tiles(t)%field%vector_Kx,                    &
          chunk(level)%tiles(t)%field%vector_Ky,                    &
          chunk(level)%tiles(t)%field%vector_Di,                    &
          chunk(level)%tiles(t)%field%rx,                           &
          chunk(level)%tiles(t)%field%ry)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF
  !if (level > 1) then
  !  DO t=1,tiles_per_task
  !    write(6,*) t,sum(abs(chunk(level)%tiles(t)%field%vector_r**2))
  !  ENDDO
  !endif

END SUBROUTINE

SUBROUTINE tea_leaf_calc_2norm(level, norm_array, norm)

  IMPLICIT NONE

  INTEGER :: t, level, norm_array
  REAL(KIND=8) :: norm, tile_norm

  norm = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_norm)
!$OMP DO REDUCTION(+:norm)
    DO t=1,tiles_per_task
      tile_norm = 0.0_8

      ! 0 = u0.u0
      ! 1 = r.r
      ! XXX add some parameters in defintions.f90?
      IF (norm_array .EQ. 0) THEN
        CALL tea_leaf_calc_2norm_kernel(chunk(level)%tiles(t)%field%x_min,        &
            chunk(level)%tiles(t)%field%x_max,                                    &
            chunk(level)%tiles(t)%field%y_min,                                    &
            chunk(level)%tiles(t)%field%y_max,                                    &
            chunk(level)%halo_exchange_depth,                                     &
            chunk(level)%tiles(t)%field%u0,                                 &
            tile_norm)
      ELSE IF (norm_array .EQ. 1) THEN
        CALL tea_leaf_calc_2norm_kernel(chunk(level)%tiles(t)%field%x_min,        &
            chunk(level)%tiles(t)%field%x_max,                                    &
            chunk(level)%tiles(t)%field%y_min,                                    &
            chunk(level)%tiles(t)%field%y_max,                                    &
            chunk(level)%halo_exchange_depth,                                     &
            chunk(level)%tiles(t)%field%vector_r,                                 &
            tile_norm)
      ELSE
        CALL report_error("tea_leaf_common.f90", "Invalid value for norm_array")
      ENDIF

      norm = norm + tile_norm
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE

SUBROUTINE tea_leaf_finalise(level)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_finalise(chunk(level)%tiles(t)%field%x_min, &
          chunk(level)%tiles(t)%field%x_max,                           &
          chunk(level)%tiles(t)%field%y_min,                           &
          chunk(level)%tiles(t)%field%y_max,                           &
          chunk(level)%halo_exchange_depth,                            &
          chunk(level)%tiles(t)%field%energy1,                         &
          chunk(level)%tiles(t)%field%density,                         &
          chunk(level)%tiles(t)%field%u)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_finalise

END MODULE

