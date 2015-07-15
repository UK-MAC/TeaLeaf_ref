
MODULE tea_leaf_jacobi_module

  USE tea_leaf_jacobi_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_jacobi_solve(error)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: error, tile_error

  error = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL NUM_THREADS(tiles_per_task) PRIVATE(tile_error)
!$OMP DO REDUCTION(+:error)
    DO t=1,tiles_per_task
      tile_error = 0.0_8

      CALL tea_leaf_jacobi_solve_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                       &
          chunk%tiles(t)%field%y_min,                       &
          chunk%tiles(t)%field%y_max,                       &
          halo_exchange_depth,                       &
          chunk%tiles(t)%field%rx,                                          &
          chunk%tiles(t)%field%ry,                                          &
          chunk%tiles(t)%field%vector_Kx,                   &
          chunk%tiles(t)%field%vector_Ky,                   &
          tile_error,                                       &
          chunk%tiles(t)%field%u0,                          &
          chunk%tiles(t)%field%u,                           &
          chunk%tiles(t)%field%vector_r)

      error = error + tile_error
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_jacobi_solve

END MODULE tea_leaf_jacobi_module

