
MODULE tea_leaf_dpcg_module

  USE tea_leaf_dpcg_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_dpcg_init()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_zero_t2_kernel(chunk%def%x_min,    &
          chunk%def%x_max,                                  &
          chunk%def%y_min,                                  &
          chunk%def%y_max,                                  &
          halo_exchange_depth,                              &
          chunk%def%t2)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  CALL tea_leaf_dpcg_solve_E()

END SUBROUTINE tea_leaf_dpcg_init

SUBROUTINE tea_leaf_dpcg_solve_E

  IMPLICIT NONE

  INTEGER :: t

  REAL(KIND=8) :: ztr

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztr)
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_sum_r_kernel(chunk%tiles(t)%field%x_min,       &
          chunk%tiles(t)%field%x_max,                                   &
          chunk%tiles(t)%field%y_min,                                   &
          chunk%tiles(t)%field%y_max,                                   &
          halo_exchange_depth,                                          &
          ztr,                                                          &
          e,                                                          &
          chunk%tiles(t)%field%vector_kx, &
          chunk%tiles(t)%field%vector_ky, &
          chunk%tiles(t)%field%vector_r)
      
      ! write back into the GLOBAL vector
      chunk%def%t1(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = ztr
      chunk%def%def_e(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = e
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_solve_E

END MODULE tea_leaf_dpcg_module


