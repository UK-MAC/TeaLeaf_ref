
MODULE tea_leaf_jacobi_module

  USE tea_leaf_jacobi_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_jacobi_solve(rx, ry, error)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry,rx, error, private_error

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(private_error)
!$OMP DO REDUCTION(+:error)
    DO t=1,tiles_per_task
      CALL tea_leaf_jacobi_solve_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                       &
          chunk%tiles(t)%field%y_min,                       &
          chunk%tiles(t)%field%y_max,                       &
          halo_exchange_depth,                       &
          rx,                                          &
          ry,                                          &
          chunk%tiles(t)%field%vector_Kx,                   &
          chunk%tiles(t)%field%vector_Ky,                   &
          private_error,                                       &
          chunk%tiles(t)%field%u0,                          &
          chunk%tiles(t)%field%u,                           &
          chunk%tiles(t)%field%vector_r)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_jacobi_solve

END MODULE tea_leaf_jacobi_module

