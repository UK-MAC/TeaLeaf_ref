
MODULE tea_leaf_dpcg_module

  USE tea_leaf_dpcg_kernel_module
  USE tea_leaf_cheby_module
  USE tea_leaf_common_module

  USE definitions_module
  use global_mpi_module
  USE update_halo_module

  IMPLICIT NONE

  LOGICAL :: inner_use_ppcg
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: inner_cg_alphas, inner_cg_betas
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: inner_ch_alphas, inner_ch_betas
  REAL(KIND=8) :: eigmin, eigmax, theta

CONTAINS

SUBROUTINE tea_leaf_dpcg_init_x0(solve_time)

  IMPLICIT NONE

  REAL(KIND=8) :: solve_time

  INTEGER :: it_count, info
  INTEGER :: fields(NUM_FIELDS)
  REAL(KIND=8) :: halo_time,timer

  IF (.NOT. ALLOCATED(inner_cg_alphas)) THEN
    ALLOCATE(inner_cg_alphas(coarse_solve_max_iters))
    ALLOCATE(inner_cg_betas (coarse_solve_max_iters))
    ALLOCATE(inner_ch_alphas(coarse_solve_max_iters))
    ALLOCATE(inner_ch_betas (coarse_solve_max_iters))
  ENDIF

  CALL tea_leaf_dpcg_coarsen_matrix()

  CALL tea_leaf_dpcg_restrict_ZT()

  ! just use CG on the first one
  inner_use_ppcg = .FALSE.

  ! FIXME if initial residual is very small, not enough steps to provide an
  ! accurate guess for the eigenvalues (if diagonal scaling on the coarse
  ! grid correction is disablee). Need to run CG for at least ~30 steps to
  ! get a good guess

  CALL tea_leaf_dpcg_local_solve(   &
      chunk%def%x_min, &
      chunk%def%x_max,                                  &
      chunk%def%y_min,                                  &
      chunk%def%y_max,                                  &
      halo_exchange_depth,                                  &
      chunk%def%t2,                               &
      chunk%def%t1,                               &
      chunk%def%def_Kx, &
      chunk%def%def_Ky, &
      chunk%def%def_di, &
      chunk%def%def_p,                               &
      chunk%def%def_r,                               &
      chunk%def%def_Mi,                               &
      chunk%def%def_w,                               &
      chunk%def%def_z, &
      chunk%def%def_sd, &
      coarse_solve_eps, &
      coarse_solve_max_iters,                          &
      it_count,         &
      0.0_8,            &
      inner_use_ppcg,       &
      inner_cg_alphas, inner_cg_betas,      &
      inner_ch_alphas, inner_ch_betas       &
      )

  ! add back onto the fine grid
  CALL tea_leaf_dpcg_add_z()

  ! for all subsequent steps, use ppcg
  inner_use_ppcg = .TRUE.

  !CALL tea_calc_eigenvalues(inner_cg_alphas, inner_cg_betas, eigmin, eigmax, &
  !    max_iters, it_count, info)
  info = 0

  ! With jacobi preconditioner on
  eigmin = 0.01_8
  eigmax = 2.0_8

  IF (info .NE. 0) CALL report_error('tea_leaf_dpcg_init_x0', 'Error in calculating eigenvalues')

  CALL tea_calc_ch_coefs(inner_ch_alphas, inner_ch_betas, eigmin, eigmax, &
      theta, it_count)

  fields = 0
  fields(FIELD_U) = 1

  ! update the halo for u prior to recalculating the residual
  IF (profiler_on) halo_time = timer()
  CALL update_halo(fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  ! calc residual again, and do initial solve
  CALL tea_leaf_calc_residual()

  CALL tea_leaf_dpcg_setup_and_solve_E(solve_time)

  CALL tea_leaf_dpcg_init_p()

END SUBROUTINE tea_leaf_dpcg_init_x0

SUBROUTINE tea_leaf_dpcg_coarsen_matrix()

  IMPLICIT NONE
  INTEGER :: t, err

  REAL(KIND=8) :: kx_local, ky_local, tile_size

  chunk%def%def_Kx = 0.0_8
  chunk%def%def_Ky = 0.0_8
  chunk%def%def_di = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(Kx_local, Ky_local)
!$OMP DO
    DO t=1,tiles_per_task
      kx_local = 0.0_8
      ky_local = 0.0_8

      CALL tea_leaf_dpcg_coarsen_matrix_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          kx_local,                             &
          ky_local,                             &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry)

      chunk%def%def_kx(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = kx_local
      chunk%def%def_ky(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = ky_local
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk%def%def_kx, size(chunk%def%def_kx), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)
  CALL MPI_Allreduce(MPI_IN_PLACE, chunk%def%def_ky, size(chunk%def%def_ky), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_size)
!$OMP DO
    DO t=1,tiles_per_task
      tile_size = chunk%tiles(t)%x_cells*chunk%tiles(t)%y_cells

      chunk%def%def_di(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = &
          tile_size + &
          chunk%def%def_kx(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) + &
          chunk%def%def_ky(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) + &
          chunk%def%def_kx(chunk%tiles(t)%def_tile_coords(1) + 1, chunk%tiles(t)%def_tile_coords(2)) + &
          chunk%def%def_ky(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2) + 1)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk%def%def_di, size(chunk%def%def_di), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_coarsen_matrix

SUBROUTINE tea_leaf_dpcg_setup_and_solve_E(solve_time)

  IMPLICIT NONE

  REAL(KIND=8) :: solve_time

  INTEGER :: it_count

  CALL tea_leaf_dpcg_matmul_ZTA(solve_time)
  CALL tea_leaf_dpcg_restrict_ZT()

  CALL tea_leaf_dpcg_local_solve(   &
      chunk%def%x_min, &
      chunk%def%x_max,                                  &
      chunk%def%y_min,                                  &
      chunk%def%y_max,                                  &
      halo_exchange_depth,                                  &
      chunk%def%t2,                               &
      chunk%def%t1,                               &
      chunk%def%def_Kx, &
      chunk%def%def_Ky, &
      chunk%def%def_di, &
      chunk%def%def_p,                               &
      chunk%def%def_r,                               &
      chunk%def%def_Mi,                               &
      chunk%def%def_w,                               &
      chunk%def%def_z, &
      chunk%def%def_sd, &
      coarse_solve_eps, &
      coarse_solve_max_iters,                          &
      it_count,         &
      theta,            &
      inner_use_ppcg,       &
      inner_cg_alphas, inner_cg_betas,      &
      inner_ch_alphas, inner_ch_betas       &
      )

  CALL tea_leaf_dpcg_prolong_Z()

END SUBROUTINE tea_leaf_dpcg_setup_and_solve_E

SUBROUTINE tea_leaf_dpcg_matmul_ZTA(solve_time)

  IMPLICIT NONE

  REAL(KIND=8) :: solve_time

  INTEGER :: t, err
  REAL(KIND=8) :: ztaz,halo_time,timer

  INTEGER :: fields(NUM_FIELDS)

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztaz)
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_solve_z_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_z,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%vector_Mi,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  chunk%def%t1 = 0.0_8

  fields = 0
  fields(FIELD_Z) = 1

  IF (profiler_on) halo_time = timer()
  CALL update_halo(fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  fields = 0
  fields(FIELD_P) = 1

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztaz)
!$OMP DO
    DO t=1,tiles_per_task
      ztaz = 0.0_8

      CALL tea_leaf_dpcg_matmul_ZTA_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%tiles(t)%field%vector_z,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          ztaz)

      ! write back into the GLOBAL vector
      chunk%def%t1(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = ztaz
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk%def%t1, size(chunk%def%t1), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_matmul_ZTA

SUBROUTINE tea_leaf_dpcg_restrict_ZT()

  IMPLICIT NONE
  INTEGER :: t, err
  REAL(KIND=8) :: ZTr

  chunk%def%t2 = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ZTr)
!$OMP DO
    DO t=1,tiles_per_task
      ztr = 0.0_8

      CALL tea_leaf_dpcg_restrict_ZT_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%vector_r,    &
          ztr)

      chunk%def%t2(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)) = ztr
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk%def%t2, size(chunk%def%t2), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_restrict_ZT

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
      tile_rrn = 0.0_8

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
          chunk%tiles(t)%field%vector_z, &
          beta)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_p

SUBROUTINE tea_leaf_dpcg_calc_zrnorm(rro)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: rro, tile_rro

  rro = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rro)
!$OMP DO REDUCTION(+:rro)
    DO t=1,tiles_per_task
      tile_rro = 0.0_8

      CALL tea_leaf_dpcg_calc_zrnorm_kernel(chunk%tiles(t)%field%x_min, &
            chunk%tiles(t)%field%x_max,                           &
            chunk%tiles(t)%field%y_min,                           &
            chunk%tiles(t)%field%y_max,                           &
            halo_exchange_depth,                           &
            chunk%tiles(t)%field%vector_z,                        &
            chunk%tiles(t)%field%vector_r,                        &
            tile_rro)

      rro = rro + tile_rro
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_zrnorm

SUBROUTINE tea_leaf_dpcg_add_z()

  IMPLICIT NONE
  INTEGER :: t

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_add_z_kernel(chunk%tiles(t)%field%x_min,    &
          chunk%tiles(t)%field%x_max,           &
          chunk%tiles(t)%field%y_min,           &
          chunk%tiles(t)%field%y_max,           &
          halo_exchange_depth,                  &
          chunk%tiles(t)%field%u, &
          chunk%def%t2(chunk%tiles(t)%def_tile_coords(1), chunk%tiles(t)%def_tile_coords(2)))
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_add_z

SUBROUTINE tea_leaf_dpcg_local_solve(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           u,                      &
                           u0,                      &
                           def_kx,                      &
                           def_ky,                      &
                           def_di,                      &
                           p,                      &
                           r,                      &
                           Mi,                     &
                           w,                     &
                           z,       &
                           sd,       &
                           eps, &
                           inner_iters,         &
                           it_count,    &
                           theta,       &
                           use_ppcg,    &
                           inner_cg_alphas, inner_cg_betas,     &
                           inner_ch_alphas, inner_ch_betas)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: u, u0, def_kx, def_ky, def_di, p, r, Mi, w, z, sd

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: it_count

  REAL(KIND=8) :: rro, smvp, initial_residual, eps
  REAL(KIND=8) ::  alpha, beta, pw, rrn

  INTEGER :: inner_iters, inner_step
  LOGICAL :: use_ppcg
  REAL(KIND=8), DIMENSION(inner_iters) :: inner_cg_alphas, inner_cg_betas
  REAL(KIND=8), DIMENSION(inner_iters) :: inner_ch_alphas, inner_ch_betas
  REAL(KIND=8) :: theta

  rro = 0.0_8
  initial_residual = 0.0_8
  pw = 0.0_8

  rrn = 1e10

  it_count = 0

!$OMP PARALLEL private(alpha, beta, smvp, inner_step)
  IF (use_ppcg) THEN
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        ! t1 = t1 - t2
        u0(j, k) = u0(j, k) - u(j, k)
      ENDDO
    ENDDO
!$OMP END DO
  ELSE
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        ! first step - t1 = t2
        u0(j, k) = u(j, k)
      ENDDO
    ENDDO
!$OMP END DO
  ENDIF

!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      ! zero for approximate solve
      u(j, k) = 0.0_8
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO REDUCTION(+:initial_residual)
  DO k=y_min, y_max
    DO j=x_min, x_max
      smvp = def_di(j, k)*u(j, k)             &
        - (def_ky(j, k+1)*u(j, k+1) + def_ky(j, k)*u(j, k-1))  &
        - (def_kx(j+1, k)*u(j+1, k) + def_kx(j, k)*u(j-1, k))

      Mi(j, k) = 1.0_8/def_di(j, k)

      r(j, k) = u0(j, k) - smvp
      z(j, k) = r(j, k)*Mi(j, k)
      p(j, k) = z(j, k)

      initial_residual = initial_residual + r(j, k)*p(j, k)
    ENDDO
  ENDDO
!$OMP END DO

!$OMP SINGLE
    rro = initial_residual
    initial_residual = sqrt(abs(initial_residual))
!$OMP END SINGLE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO WHILE ((sqrt(abs(rrn)) .gt. eps*initial_residual) .and. (it_count < inner_iters))

!$OMP BARRIER

!$OMP SINGLE
    pw = 0.0_8
    rrn = 0.0_8
!$OMP END SINGLE

!$OMP DO REDUCTION(+:pw)
    DO k=y_min,y_max
      DO j=x_min,x_max
        smvp = def_di(j, k)*p(j, k)             &
          - (def_ky(j, k+1)*p(j, k+1) + def_ky(j, k)*p(j, k-1))  &
          - (def_kx(j+1, k)*p(j+1, k) + def_kx(j, k)*p(j-1, k))

        w(j, k) = smvp
        pw = pw + smvp*p(j, k)
      ENDDO
    ENDDO
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    alpha = rro/pw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        u(j, k) = u(j, k) + alpha*p(j, k)
        r(j, k) = r(j, k) - alpha*w(j, k)
        z(j, k) = r(j, k)*Mi(j, k)
      ENDDO
    ENDDO
!$OMP END DO

    IF (use_ppcg) THEN
!$OMP DO
      DO k=y_min,y_max
        DO j=x_min,x_max
          sd(j, k) = z(j, k)/theta
        ENDDO
      ENDDO
!$OMP END DO

      DO inner_step=1,10
!$OMP DO
        DO k=y_min,y_max
          DO j=x_min,x_max
            smvp = def_di(j, k)*sd(j, k)             &
              - (def_ky(j, k+1)*sd(j, k+1) + def_ky(j, k)*sd(j, k-1))  &
              - (def_kx(j+1, k)*sd(j+1, k) + def_kx(j, k)*sd(j-1, k))

            r(j, k) = r(j, k) - smvp
            z(j, k) = r(j, k)*Mi(j, k)

            u(j, k) = u(j, k) + sd(j, k)
          ENDDO
        ENDDO
!$OMP END DO
!$OMP DO
        DO k=y_min,y_max
          DO j=x_min,x_max
            sd(j, k) = inner_ch_alphas(inner_step)*sd(j, k) + inner_ch_betas(inner_step)*z(j, k)
          ENDDO
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF

!$OMP DO REDUCTION(+:rrn)
    DO k=y_min,y_max
      DO j=x_min,x_max
        rrn = rrn + r(j, k)*z(j, k)
      ENDDO
    ENDDO
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    beta = rrn/rro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        p(j, k) = z(j, k) + beta*p(j, k)
      ENDDO
    ENDDO
!$OMP END DO

!$OMP SINGLE
    rro = rrn
    it_count = it_count + 1
    inner_cg_alphas(it_count) = alpha
    inner_cg_betas(it_count) = beta
!$OMP END SINGLE

  ENDDO

!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_local_solve

END MODULE tea_leaf_dpcg_module

