
MODULE tea_leaf_dpcg_module

  USE tea_leaf_dpcg_kernel_module
  USE definitions_module

  use global_mpi_module
  USE tea_leaf_common_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_dpcg_init_x0()

  IMPLICIT NONE

  INTEGER :: t

  ! get E
  CALL tea_leaf_dpcg_sum_matrix_rows()

  ! TODO do initial bit
  ! may not need to do ?

  ! initial solve
  CALL tea_leaf_dpcg_setup_and_solve_E()

  CALL tea_leaf_dpcg_init_p()

END SUBROUTINE tea_leaf_dpcg_init_x0

SUBROUTINE tea_leaf_dpcg_sum_matrix_rows()

  IMPLICIT NONE
  INTEGER :: t, err
  REAL(KIND=8) :: E_local

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(E_local)
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_sum_matrix_rows_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%vector_kx, &
          chunk%tiles(t)%field%vector_ky, &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          E_local)

      chunk%def%def_e(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = E_local
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk%def%def_e, size(chunk%def%def_e), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_sum_matrix_rows

SUBROUTINE tea_leaf_dpcg_restrict_ZT()

  IMPLICIT NONE
  INTEGER :: t, err
  REAL(KIND=8) :: ZTr

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ZTr)
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_restrict_ZT_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%vector_r,    &
          ztr)

      chunk%def%t1(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = ztr
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk%def%t1, size(chunk%def%t1), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_restrict_ZT

SUBROUTINE tea_leaf_dpcg_setup_and_solve_E

  IMPLICIT NONE

  INTEGER :: t, err

  CALL tea_leaf_dpcg_matmul_ZTA()

  CALL tea_leaf_dpcg_local_solve(   &
      chunk%def%x_min, &
      chunk%def%x_max,                                  &
      chunk%def%y_min,                                  &
      chunk%def%y_max,                                  &
      halo_exchange_depth,                                  &
      chunk%def%t2,                               &
      chunk%def%t1,                               &
      chunk%def%def_e,                               &
      chunk%def%def_p,                               &
      chunk%def%def_r,                               &
      chunk%def%def_Mi,                               &
      chunk%def%def_w,                               &
      chunk%def%def_z)

  CALL tea_leaf_dpcg_prolong_Z()

END SUBROUTINE tea_leaf_dpcg_setup_and_solve_E

SUBROUTINE tea_leaf_dpcg_matmul_ZTA()

  IMPLICIT NONE

  INTEGER :: t, err
  REAL(KIND=8) :: ztaz

  ! always done first
  CALL tea_leaf_dpcg_restrict_ZT()

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztaz)
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_matmul_ZTA_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_Mi,                              &
          chunk%tiles(t)%field%vector_z,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          ztaz,                     &
          tl_preconditioner_type)

      ! write back into the GLOBAL vector
      chunk%def%t1(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = &
        ztaz + chunk%def%t2(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2))
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk%def%t2, size(chunk%def%t1), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_matmul_ZTA

SUBROUTINE tea_leaf_dpcg_prolong_Z()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_prolong_Z_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%vector_z, &
          chunk%def%t2(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)))
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_prolong_Z

SUBROUTINE tea_leaf_dpcg_init_p()

  IMPLICIT NONE

  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_init_p_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                         &
          chunk%tiles(t)%field%y_min,                                         &
          chunk%tiles(t)%field%y_max,                                         &
          halo_exchange_depth,                                         &
          chunk%tiles(t)%field%vector_p,                                      &
          chunk%tiles(t)%field%vector_z)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_init_p

SUBROUTINE tea_leaf_dpcg_store_r()

  IMPLICIT NONE
  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_store_r_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%vector_r, &
          chunk%tiles(t)%field%vector_r_m1 )
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_store_r

SUBROUTINE tea_leaf_dpcg_calc_rrn(rrn)

  IMPLICIT NONE
  INTEGER :: t
  REAL(KIND=8) :: rrn, tile_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rrn)
!$OMP DO REDUCTION(+:rrn)
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_calc_rrn_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%vector_r, &
          chunk%tiles(t)%field%vector_r_m1, &
          chunk%tiles(t)%field%vector_z, &
          tile_rrn)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_rrn

SUBROUTINE tea_leaf_dpcg_calc_p(beta)

  IMPLICIT NONE
  INTEGER :: t
  REAL(KIND=8) :: beta

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_calc_p_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%vector_p, &
          chunk%tiles(t)%field%vector_r, &
          beta )
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_p

END MODULE tea_leaf_dpcg_module

