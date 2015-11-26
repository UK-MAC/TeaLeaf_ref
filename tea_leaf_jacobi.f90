
MODULE tea_leaf_jacobi_module

  USE tea_leaf_jacobi_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_jacobi_solve(level, error)

  IMPLICIT NONE

  INTEGER :: level, t
  REAL(KIND=8) :: error, tile_error

  error = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_error)
!$OMP DO REDUCTION(+:error)
    DO t=1,tiles_per_task
      tile_error = 0.0_8

      CALL tea_leaf_jacobi_solve_kernel(chunk(level)%tiles(t)%field%x_min,&
          chunk(level)%tiles(t)%field%x_max,                       &
          chunk(level)%tiles(t)%field%y_min,                       &
          chunk(level)%tiles(t)%field%y_max,                       &
          chunk(level)%halo_exchange_depth,                        &
          chunk(level)%tiles(t)%field%rx,                                          &
          chunk(level)%tiles(t)%field%ry,                                          &
          chunk(level)%tiles(t)%field%vector_Kx,                   &
          chunk(level)%tiles(t)%field%vector_Ky,                   &
          tile_error,                                       &
          chunk(level)%tiles(t)%field%u0,                          &
          chunk(level)%tiles(t)%field%u,                           &
          chunk(level)%tiles(t)%field%vector_r)

      error = error + tile_error
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_jacobi_solve

END MODULE tea_leaf_jacobi_module

