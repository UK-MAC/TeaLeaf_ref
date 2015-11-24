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

  REAL(KIND=8) :: solve_time,tile_sum2

  INTEGER, PARAMETER :: level=1
  INTEGER :: t, it_count, info
  INTEGER :: fields(NUM_FIELDS)
  REAL(KIND=8) :: halo_time,timer

  IF (.NOT. ALLOCATED(inner_cg_alphas)) THEN
    ALLOCATE(inner_cg_alphas(coarse_solve_max_iters))
    ALLOCATE(inner_cg_betas (coarse_solve_max_iters))
    ALLOCATE(inner_ch_alphas(coarse_solve_max_iters))
    ALLOCATE(inner_ch_betas (coarse_solve_max_iters))
  ENDIF

  !CALL tea_leaf_dpcg_coarsen_matrix()
  CALL tea_leaf_dpcg_coarsen_matrix_level(level,solve_time)

  ! just use CG on the first one
  inner_use_ppcg = .FALSE.

  !chunk(level)%def%t1 = 0.0_8
  !CALL tea_leaf_dpcg_restrict_ZT(.TRUE.)

  !tile_sum2 = sum(chunk(level)%def%t2**2)
  !IF (parallel%boss) write(6,"(a17,es25.18)") "in -serial solve:",sqrt(tile_sum2)

!$OMP PARALLEL
!$OMP DO
  ! store the RHS in u which will then be copied into u0 
  DO t=1,tiles_per_task
    chunk(level+1)%tiles(t)%field%u = 0.0_8
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
  CALL tea_leaf_dpcg_restrict_ZT_level(level,.TRUE.)

  !tile_sum2 = 0.0_8
  !DO t=1,tiles_per_task
  !  !write(6,*) "Tile:",t
  !  !write(6,"(7es12.5)") chunk(level+1)%tiles(t)%field%u
  !  tile_sum2 = tile_sum2+sum(chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
  !                                                            1:chunk(level)%sub_tile_dims(2))**2)
  !ENDDO
  !CALL tea_allsum(tile_sum2)
  !IF (parallel%boss) write(6,"(a17,es25.18)") "in -tiled  solve:",sqrt(tile_sum2)

  !CALL tea_leaf_dpcg_local_solve(   &
  !    chunk(level)%def%x_min, &
  !    chunk(level)%def%x_max,                                  &
  !    chunk(level)%def%y_min,                                  &
  !    chunk(level)%def%y_max,                                  &
  !    halo_exchange_depth,                                  &
  !    chunk(level)%def%t2,                               &
  !    chunk(level)%def%t1,                               &
  !    chunk(level)%def%def_Kx, &
  !    chunk(level)%def%def_Ky, &
  !    chunk(level)%def%def_di, &
  !    chunk(level)%def%def_p,                               &
  !    chunk(level)%def%def_r,                               &
  !    chunk(level)%def%def_Mi,                               &
  !    chunk(level)%def%def_w,                               &
  !    chunk(level)%def%def_z, &
  !    chunk(level)%def%def_sd, &
  !    coarse_solve_eps, &
  !    coarse_solve_max_iters,                          &
  !    it_count,         &
  !    0.0_8,            &
  !    inner_use_ppcg,       &
  !    inner_cg_alphas, inner_cg_betas,      &
  !    inner_ch_alphas, inner_ch_betas       &
  !    )

  !!write(6,"(12es12.5)") chunk(level)%def%t2
  !tile_sum2 = sum(chunk(level)%def%t2**2)
  !IF (parallel%boss) write(6,"(a17,es25.18)") "out-serial solve:",sqrt(tile_sum2)

  CALL tea_leaf_dpcg_local_solve_level(level,                  &
                                       solve_time,             &
                                       coarse_solve_eps,       &
                                       coarse_solve_max_iters, &
                                       it_count,               &
                                       inner_use_ppcg,         &
                                       theta,                  &
                                       inner_ch_alphas,        &
                                       inner_ch_betas)

  !tile_sum2 = 0.0_8
  !DO t=1,tiles_per_task
  !  !write(6,*) "Tile:",t
  !  !write(6,"(7es12.5)") chunk(level+1)%tiles(t)%field%u
  !  tile_sum2 = tile_sum2+sum(chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
  !                                                            1:chunk(level)%sub_tile_dims(2))**2)
  !ENDDO
  !CALL tea_allsum(tile_sum2)
  !IF (parallel%boss) write(6,"(a17,es25.18)") "out-tiled  solve:",sqrt(tile_sum2)

  ! add back onto the fine grid
  !CALL tea_leaf_dpcg_subtract_z()

!!$OMP PARALLEL
!!$OMP DO
  ! zero the coarse grid solution u so that we don't change the answer 
  !DO t=1,tiles_per_task
  !  chunk(level+1)%tiles(t)%field%u = 0.0_8
  !ENDDO
!!$OMP END DO
!!$OMP END PARALLEL
  CALL tea_leaf_dpcg_subtract_z_level(level)

  ! for all subsequent steps, use ppcg
  !inner_use_ppcg = .TRUE.

  !CALL tea_calc_eigenvalues(inner_cg_alphas, inner_cg_betas, eigmin, eigmax, &
  !    coarse_solve_max_iters, it_count, info)
  info = 0

  ! With jacobi preconditioner on
  eigmin = 0.1_8
  eigmax = 2.0_8

  IF (info .NE. 0) CALL report_error('tea_leaf_dpcg_init_x0', 'Error in calculating eigenvalues')

  CALL tea_calc_ch_coefs(inner_ch_alphas, inner_ch_betas, eigmin, eigmax, &
      theta, coarse_solve_max_iters)

  fields = 0
  fields(FIELD_U) = 1

  ! update the halo for u prior to recalculating the residual
  IF (profiler_on) halo_time = timer()
  CALL update_halo(level, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  ! calc residual again, and do initial solve
  CALL tea_leaf_calc_residual(level)

  !CALL tea_leaf_dpcg_setup_and_solve_E(solve_time)
  CALL tea_leaf_dpcg_setup_and_solve_E_level(level,solve_time)

  CALL tea_leaf_dpcg_init_p()

END SUBROUTINE tea_leaf_dpcg_init_x0

SUBROUTINE tea_leaf_dpcg_setup_and_solve_E(solve_time)

  IMPLICIT NONE

  REAL(KIND=8) :: solve_time, tile_sum1, tile_sum2

  INTEGER :: level, it_count, t

  level=1

  CALL tea_leaf_dpcg_matmul_ZTA(solve_time)
  CALL tea_leaf_dpcg_restrict_ZT(.TRUE.)

  !write(6,"(12es12.5)") chunk(level)%def%t2
  tile_sum1 = sum(chunk(level)%def%t2**2)
  IF (parallel%boss) write(6,"(a17,es25.18)") "in -serial solve:",sqrt(tile_sum1)

  CALL tea_leaf_dpcg_local_solve(   &
      chunk(level)%def%x_min, &
      chunk(level)%def%x_max,                                  &
      chunk(level)%def%y_min,                                  &
      chunk(level)%def%y_max,                                  &
      halo_exchange_depth,                                  &
      chunk(level)%def%t2,                               &
      chunk(level)%def%t1,                               &
      chunk(level)%def%def_Kx, &
      chunk(level)%def%def_Ky, &
      chunk(level)%def%def_di, &
      chunk(level)%def%def_p,                               &
      chunk(level)%def%def_r,                               &
      chunk(level)%def%def_Mi,                               &
      chunk(level)%def%def_w,                               &
      chunk(level)%def%def_z, &
      chunk(level)%def%def_sd, &
      coarse_solve_eps, &
      coarse_solve_max_iters,                          &
      it_count,         &
      theta,            &
      inner_use_ppcg,       &
      inner_cg_alphas, inner_cg_betas,      &
      inner_ch_alphas, inner_ch_betas       &
      )

  !write(6,"(12es12.5)") chunk(level)%def%t2
  tile_sum2 = sum(chunk(level)%def%t2**2)
  IF (parallel%boss) write(6,"(a17,es25.18)") "out-serial solve:",sqrt(tile_sum2)

  CALL tea_leaf_dpcg_setup_and_solve_E_level(level,solve_time)

  tile_sum2 = 0.0_8
  DO t=1,tiles_per_task
    !write(6,*) "Tile:",t
    !write(6,"(7es12.5)") chunk(level+1)%tiles(t)%field%u
    tile_sum2 = tile_sum2+sum(chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                                              1:chunk(level)%sub_tile_dims(2))**2)
  ENDDO
  CALL tea_allsum(tile_sum2)
  IF (parallel%boss) write(6,"(a17,es25.18)") "out-tiled  solve:",sqrt(tile_sum2)

  !CALL tea_leaf_dpcg_prolong_Z()

END SUBROUTINE tea_leaf_dpcg_setup_and_solve_E

SUBROUTINE tea_leaf_dpcg_setup_and_solve_E_level(level,solve_time)

  IMPLICIT NONE

  REAL(KIND=8) :: solve_time, tile_sum2

  INTEGER :: level, it_count, t

  CALL tea_leaf_dpcg_matmul_ZTA_level(level,solve_time)
  CALL tea_leaf_dpcg_restrict_ZT_level(level,.TRUE.)

  !tile_sum2 = 0.0_8
  !DO t=1,tiles_per_task
  !  !write(6,*) "Tile:",t
  !  !write(6,"(7es12.5)") chunk(level+1)%tiles(t)%field%u
  !  tile_sum2 = tile_sum2+sum(chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
  !                                                            1:chunk(level)%sub_tile_dims(2))**2)
  !ENDDO
  !CALL tea_allsum(tile_sum2)
  !IF (parallel%boss) write(6,"(a17,es25.18)") "in -tiled  solve:",sqrt(tile_sum2)

  CALL tea_leaf_dpcg_local_solve_level(level,                  &
                                       solve_time,             &
                                       coarse_solve_eps,       &
                                       coarse_solve_max_iters, &
                                       it_count,               &
                                       inner_use_ppcg,         &
                                       theta,                  &
                                       inner_ch_alphas,        &
                                       inner_ch_betas)

  CALL tea_leaf_dpcg_prolong_Z_level(level)

END SUBROUTINE tea_leaf_dpcg_setup_and_solve_E_level

SUBROUTINE tea_leaf_dpcg_coarsen_matrix()

  IMPLICIT NONE
  INTEGER :: t, err
  INTEGER, PARAMETER :: level=1

  INTEGER :: sub_tile_dx, sub_tile_dy

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: jj,j_start,j_end,kk,k_start,k_end
  REAL(KIND=8) :: tile_size
  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: kx_local, ky_local

  chunk(level)%def%def_Kx = 0.0_8
  chunk(level)%def%def_Ky = 0.0_8
  chunk(level)%def%def_di = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(Kx_local, Ky_local,sub_tile_dx,sub_tile_dy)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
                   chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
                   chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      kx_local = 0.0_8
      ky_local = 0.0_8

      CALL tea_leaf_dpcg_coarsen_matrix_kernel(chunk(level)%tiles(t)%field%x_min,     &
                                               chunk(level)%tiles(t)%field%x_max,     &
                                               chunk(level)%tiles(t)%field%y_min,     &
                                               chunk(level)%tiles(t)%field%y_max,     &
                                               halo_exchange_depth,            &
                                               chunk(level)%tiles(t)%field%vector_Kx, &
                                               chunk(level)%tiles(t)%field%vector_Ky, &
                                               chunk(level)%sub_tile_dims(1),         &
                                               sub_tile_dx,                    &
                                               chunk(level)%sub_tile_dims(2),         &
                                               sub_tile_dy,                    &
                                               kx_local,                       &
                                               ky_local,                       &
                                               chunk(level)%tiles(t)%field%rx,        &
                                               chunk(level)%tiles(t)%field%ry)

      chunk(level)%def%def_kx(chunk(level)%tiles(t)%def_tile_coords(1): &
                       chunk(level)%tiles(t)%def_tile_coords(1)+chunk(level)%sub_tile_dims(1)-1, &
                       chunk(level)%tiles(t)%def_tile_coords(2): &
                       chunk(level)%tiles(t)%def_tile_coords(2)+chunk(level)%sub_tile_dims(2)-1) = kx_local
      chunk(level)%def%def_ky(chunk(level)%tiles(t)%def_tile_coords(1): &
                       chunk(level)%tiles(t)%def_tile_coords(1)+chunk(level)%sub_tile_dims(1)-1, &
                       chunk(level)%tiles(t)%def_tile_coords(2): &
                       chunk(level)%tiles(t)%def_tile_coords(2)+chunk(level)%sub_tile_dims(2)-1) = ky_local
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk(level)%def%def_kx, size(chunk(level)%def%def_kx), &
    MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)
  CALL MPI_Allreduce(MPI_IN_PLACE, chunk(level)%def%def_ky, size(chunk(level)%def%def_ky), &
    MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_size,sub_tile_dx,sub_tile_dy,k,kk,k_start,k_end,j,jj,j_start,j_end)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
                   chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
                   chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      DO kk=1,chunk(level)%sub_tile_dims(2)
        k_start=chunk(level)%tiles(t)%field%y_min+(kk-1)*sub_tile_dy
        k_end  =min(k_start+sub_tile_dy-1,chunk(level)%tiles(t)%field%y_max)
        DO jj=1,chunk(level)%sub_tile_dims(1)
          j_start=chunk(level)%tiles(t)%field%x_min+(jj-1)*sub_tile_dx
          j_end  =min(j_start+sub_tile_dx-1,chunk(level)%tiles(t)%field%x_max)
          tile_size=(j_end-j_start+1)*(k_end-k_start+1)
          chunk(level)%def%def_di(chunk(level)%tiles(t)%def_tile_coords(1)+jj-1, &
                                  chunk(level)%tiles(t)%def_tile_coords(2)+kk-1) = &
            tile_size + &
            chunk(level)%def%def_kx(chunk(level)%tiles(t)%def_tile_coords(1)+jj-1    , &
                                    chunk(level)%tiles(t)%def_tile_coords(2)+kk-1    ) + &
            chunk(level)%def%def_ky(chunk(level)%tiles(t)%def_tile_coords(1)+jj-1    , &
                                    chunk(level)%tiles(t)%def_tile_coords(2)+kk-1    ) + &
            chunk(level)%def%def_kx(chunk(level)%tiles(t)%def_tile_coords(1)+jj-1 + 1, &
                                    chunk(level)%tiles(t)%def_tile_coords(2)+kk-1    ) + &
            chunk(level)%def%def_ky(chunk(level)%tiles(t)%def_tile_coords(1)+jj-1    , &
                                    chunk(level)%tiles(t)%def_tile_coords(2)+kk-1 + 1)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

  CALL MPI_Allreduce(MPI_IN_PLACE, chunk(level)%def%def_di, size(chunk(level)%def%def_di), &
     MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

  !write(6,*) "Deflation matrix:",shape(chunk(level)%def%def_kx),shape(chunk(level)%def%def_ky)
  !write(6,*) "Kx:"
  !write(6,"(12es12.5)") chunk(level)%def%def_kx
  !write(6,*) "Ky:"
  !write(6,"(12es12.5)") chunk(level)%def%def_ky

END SUBROUTINE tea_leaf_dpcg_coarsen_matrix

SUBROUTINE tea_leaf_dpcg_coarsen_matrix_level(level,solve_time)

  IMPLICIT NONE
  INTEGER :: t, err
  INTEGER :: level

  INTEGER :: sub_tile_dx, sub_tile_dy

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: jj,j_start,j_end,kk,k_start,k_end

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: solve_time
  REAL(KIND=8) :: timer,halo_time
  REAL(KIND=8) :: tile_size
  !REAL(KIND=8) :: kx_tot,ky_tot,Di_tot,Mi_tot
  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: kx_local, ky_local

  ! 3 different options for preconditioners
  INTEGER,PARAMETER::           TL_PREC_NONE       = 1 &
                               ,TL_PREC_JAC_DIAG   = 2 &
                               ,TL_PREC_JAC_BLOCK  = 3

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(Kx_local, Ky_local,sub_tile_dx,sub_tile_dy)
!$OMP DO
    DO t=1,tiles_per_task
      ! define the coarse Kx and Ky with unit scaling factor
      chunk(level+1)%tiles(t)%field%rx = 1.0_8
      chunk(level+1)%tiles(t)%field%ry = 1.0_8
    ENDDO
!$OMP END DO

!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
                   chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
                   chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      chunk(level+1)%tiles(t)%field%vector_Kx = 0.0_8
      chunk(level+1)%tiles(t)%field%vector_Ky = 0.0_8

      kx_local = 0.0_8
      ky_local = 0.0_8

      CALL tea_leaf_dpcg_coarsen_matrix_kernel(chunk(level)%tiles(t)%field%x_min,     &
                                               chunk(level)%tiles(t)%field%x_max,     &
                                               chunk(level)%tiles(t)%field%y_min,     &
                                               chunk(level)%tiles(t)%field%y_max,     &
                                               halo_exchange_depth,            &
                                               chunk(level)%tiles(t)%field%vector_Kx, &
                                               chunk(level)%tiles(t)%field%vector_Ky, &
                                               chunk(level)%sub_tile_dims(1),         &
                                               sub_tile_dx,                    &
                                               chunk(level)%sub_tile_dims(2),         &
                                               sub_tile_dy,                    &
                                               kx_local,                       &
                                               ky_local,                       &
                                               chunk(level)%tiles(t)%field%rx,        &
                                               chunk(level)%tiles(t)%field%ry)

      chunk(level+1)%tiles(t)%field%vector_Kx(1:chunk(level)%sub_tile_dims(1),1:chunk(level)%sub_tile_dims(2)) = kx_local
      chunk(level+1)%tiles(t)%field%vector_Ky(1:chunk(level)%sub_tile_dims(1),1:chunk(level)%sub_tile_dims(2)) = ky_local
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF
  !kx_tot=0.0_8; ky_tot=0.0_8
  !DO t=1,tiles_per_task
  !  kx_tot=kx_tot+sum(chunk(level+1)%tiles(t)%field%vector_Kx**2)
  !  ky_tot=ky_tot+sum(chunk(level+1)%tiles(t)%field%vector_Ky**2)
  !ENDDO
  !call tea_allsum(kx_tot); call tea_allsum(ky_tot)
  !write(6,*) "Kx_tot,Ky_tot:",kx_tot,ky_tot

!Need a depth one halo exchange on Kx and Ky
!use custom comms for Kx, Ky (and Kz in 3D)
  fields=0
  fields(FIELD_KX)=1
  fields(FIELD_KY)=1
  IF (profiler_on) halo_time = timer()
  CALL update_halo(level+1, fields, 1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_size,sub_tile_dx,sub_tile_dy,k,kk,k_start,k_end,j,jj,j_start,j_end)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
                   chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
                   chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      chunk(level+1)%tiles(t)%field%vector_Di = 0.0_8

      DO kk=1,chunk(level)%sub_tile_dims(2)
        k_start=chunk(level)%tiles(t)%field%y_min+(kk-1)*sub_tile_dy
        k_end  =min(k_start+sub_tile_dy-1,chunk(level)%tiles(t)%field%y_max)
        DO jj=1,chunk(level)%sub_tile_dims(1)
          j_start=chunk(level)%tiles(t)%field%x_min+(jj-1)*sub_tile_dx
          j_end  =min(j_start+sub_tile_dx-1,chunk(level)%tiles(t)%field%x_max)
          tile_size=(j_end-j_start+1)*(k_end-k_start+1)
          chunk(level+1)%tiles(t)%field%vector_Di(jj,kk) = &
            tile_size + &
            chunk(level+1)%tiles(t)%field%vector_Kx(jj    ,kk    ) + &
            chunk(level+1)%tiles(t)%field%vector_Ky(jj    ,kk    ) + &
            chunk(level+1)%tiles(t)%field%vector_Kx(jj + 1,kk    ) + &
            chunk(level+1)%tiles(t)%field%vector_Ky(jj    ,kk + 1)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF
  !Di_tot=0.0_8
  !DO t=1,tiles_per_task
  !  Di_tot=Di_tot+sum(chunk(level+1)%tiles(t)%field%vector_Di**2)
  !ENDDO
  !call tea_allsum(Di_tot)
  !write(6,*) "Di_tot:",Di_tot

  IF (use_fortran_kernels) THEN
    IF (tl_preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
! The preconditioner construction must be called outside of a parallel region
! as the routine is threaded otherwise the result array is not shared correctly
      DO t=1,tiles_per_task
        CALL tea_block_init(chunk(level+1)%tiles(t)%field%x_min,     &
                            chunk(level+1)%tiles(t)%field%x_max,     &
                            chunk(level+1)%tiles(t)%field%y_min,     &
                            chunk(level+1)%tiles(t)%field%y_max,     &
                            halo_exchange_depth,                     &
                            chunk(level+1)%tiles(t)%field%tri_cp,    &
                            chunk(level+1)%tiles(t)%field%tri_bfp,   &
                            chunk(level+1)%tiles(t)%field%vector_Kx, &
                            chunk(level+1)%tiles(t)%field%vector_Ky, &
                            chunk(level+1)%tiles(t)%field%vector_Di, &
                            chunk(level+1)%tiles(t)%field%rx,        &
                            chunk(level+1)%tiles(t)%field%ry)
      ENDDO
    ELSE IF (tl_preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
! The preconditioner construction must be called outside of a parallel region
! as the routine is threaded otherwise the result array is not shared correctly
      DO t=1,tiles_per_task
        CALL tea_diag_init(chunk(level+1)%tiles(t)%field%x_min,     &
                           chunk(level+1)%tiles(t)%field%x_max,     &
                           chunk(level+1)%tiles(t)%field%y_min,     &
                           chunk(level+1)%tiles(t)%field%y_max,     &
                           halo_exchange_depth,                     &
                           chunk(level+1)%tiles(t)%field%vector_Mi, &
                           chunk(level+1)%tiles(t)%field%vector_Kx, &
                           chunk(level+1)%tiles(t)%field%vector_Ky, &
                           chunk(level+1)%tiles(t)%field%vector_Di, &
                           chunk(level+1)%tiles(t)%field%rx,        &
                           chunk(level+1)%tiles(t)%field%ry)
      ENDDO
    ENDIF
  ENDIF
  !Mi_tot=0.0_8
  !DO t=1,tiles_per_task
  !  Mi_tot=Mi_tot+sum(chunk(level+1)%tiles(t)%field%vector_Mi**2)
  !ENDDO
  !call tea_allsum(Mi_tot)
  !write(6,*) "Mi_tot:",Mi_tot

  !DO t=1,tiles_per_task
  !  write(6,*) "Tile:",t,size(chunk(level+1)%tiles(t)%field%vector_Kx),size(chunk(level+1)%tiles(t)%field%vector_Ky)
  !  write(6,*) "Kx:"
  !  write(6,"(7es12.5)") chunk(level+1)%tiles(t)%field%vector_Kx
  !  write(6,*) "Ky:"
  !  write(6,"(7es12.5)") chunk(level+1)%tiles(t)%field%vector_Ky
  !ENDDO

END SUBROUTINE tea_leaf_dpcg_coarsen_matrix_level

SUBROUTINE tea_leaf_dpcg_matmul_ZTA(solve_time)

  IMPLICIT NONE

  REAL(KIND=8) :: solve_time
  INTEGER, PARAMETER :: level=1

  INTEGER :: t, err

  INTEGER :: sub_tile_dx, sub_tile_dy

  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: ztaz

  REAL(KIND=8) :: halo_time,timer

  INTEGER :: fields(NUM_FIELDS)

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztaz)
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_solve_z_kernel(chunk(level)%tiles(t)%field%x_min,     &
                                        chunk(level)%tiles(t)%field%x_max,     &
                                        chunk(level)%tiles(t)%field%y_min,     &
                                        chunk(level)%tiles(t)%field%y_max,     &
                                        halo_exchange_depth,            &
                                        chunk(level)%tiles(t)%field%vector_r,  &
                                        chunk(level)%tiles(t)%field%vector_z,  &
                                        chunk(level)%tiles(t)%field%vector_Kx, &
                                        chunk(level)%tiles(t)%field%vector_Ky, &
                                        chunk(level)%tiles(t)%field%vector_Di, &
                                        chunk(level)%tiles(t)%field%vector_Mi, &
                                        chunk(level)%tiles(t)%field%tri_cp,    &
                                        chunk(level)%tiles(t)%field%tri_bfp,   &
                                        chunk(level)%tiles(t)%field%rx,        &
                                        chunk(level)%tiles(t)%field%ry,        &
                                        tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  fields = 0
  fields(FIELD_Z) = 1

  IF (profiler_on) halo_time = timer()
  CALL update_halo(level, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  fields = 0
  fields(FIELD_P) = 1

  chunk(level)%def%t1 = 0.0_8
  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztaz,sub_tile_dx,sub_tile_dy)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
        chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
        chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      ztaz = 0.0_8

      CALL tea_leaf_dpcg_matmul_ZTA_kernel(chunk(level)%tiles(t)%field%x_min, &
          chunk(level)%tiles(t)%field%x_max,                                  &
          chunk(level)%tiles(t)%field%y_min,                                  &
          chunk(level)%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                         &
          chunk(level)%sub_tile_dims(1),sub_tile_dx,                          &
          chunk(level)%sub_tile_dims(2),sub_tile_dy,                          &
          chunk(level)%tiles(t)%field%vector_z,                               &
          chunk(level)%tiles(t)%field%vector_Kx,                              &
          chunk(level)%tiles(t)%field%vector_Ky,                              &
          chunk(level)%tiles(t)%field%vector_Di,                              &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,  &
          ztaz)

      ! write back into the GLOBAL vector
      chunk(level)%def%t1(chunk(level)%tiles(t)%def_tile_coords(1): &
                          chunk(level)%tiles(t)%def_tile_coords(1)+chunk(level)%sub_tile_dims(1)-1, &
                          chunk(level)%tiles(t)%def_tile_coords(2): &
                          chunk(level)%tiles(t)%def_tile_coords(2)+chunk(level)%sub_tile_dims(2)-1) = ztaz
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  !CALL MPI_Allreduce(MPI_IN_PLACE, chunk(level)%def%t1, size(chunk(level)%def%t1), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_matmul_ZTA

SUBROUTINE tea_leaf_dpcg_matmul_ZTA_level(level,solve_time)

  IMPLICIT NONE

  INTEGER :: level
  REAL(KIND=8) :: solve_time

  INTEGER :: t, err

  INTEGER :: sub_tile_dx, sub_tile_dy

  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: ztaz

  REAL(KIND=8) :: halo_time,timer

  INTEGER :: fields(NUM_FIELDS)

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztaz)
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_solve_z_kernel(chunk(level)%tiles(t)%field%x_min,     &
                                        chunk(level)%tiles(t)%field%x_max,     &
                                        chunk(level)%tiles(t)%field%y_min,     &
                                        chunk(level)%tiles(t)%field%y_max,     &
                                        halo_exchange_depth,            &
                                        chunk(level)%tiles(t)%field%vector_r,  &
                                        chunk(level)%tiles(t)%field%vector_z,  &
                                        chunk(level)%tiles(t)%field%vector_Kx, &
                                        chunk(level)%tiles(t)%field%vector_Ky, &
                                        chunk(level)%tiles(t)%field%vector_Di, &
                                        chunk(level)%tiles(t)%field%vector_Mi, &
                                        chunk(level)%tiles(t)%field%tri_cp,    &
                                        chunk(level)%tiles(t)%field%tri_bfp,   &
                                        chunk(level)%tiles(t)%field%rx,        &
                                        chunk(level)%tiles(t)%field%ry,        &
                                        tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  fields = 0
  fields(FIELD_Z) = 1

  IF (profiler_on) halo_time = timer()
  CALL update_halo(level, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  fields = 0
  fields(FIELD_P) = 1

  ! store the RHS in u which will then be copied into u0
!$OMP PARALLEL
!$OMP DO
  DO t=1,tiles_per_task
    chunk(level+1)%tiles(t)%field%u = 0.0_8
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(ztaz,sub_tile_dx,sub_tile_dy)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
        chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
        chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      ztaz = 0.0_8

      CALL tea_leaf_dpcg_matmul_ZTA_kernel(chunk(level)%tiles(t)%field%x_min, &
          chunk(level)%tiles(t)%field%x_max,                                  &
          chunk(level)%tiles(t)%field%y_min,                                  &
          chunk(level)%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                         &
          chunk(level)%sub_tile_dims(1),sub_tile_dx,                          &
          chunk(level)%sub_tile_dims(2),sub_tile_dy,                          &
          chunk(level)%tiles(t)%field%vector_z,                               &
          chunk(level)%tiles(t)%field%vector_Kx,                              &
          chunk(level)%tiles(t)%field%vector_Ky,                              &
          chunk(level)%tiles(t)%field%vector_Di,                              &
          chunk(level)%tiles(t)%field%rx,  &
          chunk(level)%tiles(t)%field%ry,  &
          ztaz)

      ! write back into the GLOBAL vector
      chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                      1:chunk(level)%sub_tile_dims(2)) = ztaz
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  !CALL MPI_Allreduce(MPI_IN_PLACE, chunk(level)%def%t1, size(chunk(level)%def%t1), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_matmul_ZTA_level

SUBROUTINE tea_leaf_dpcg_restrict_ZT(not_init)

  IMPLICIT NONE
  LOGICAL :: not_init
  INTEGER :: t, err
  INTEGER, PARAMETER :: level=1

  INTEGER :: sub_tile_dx, sub_tile_dy
  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: ZTr

  IF (use_fortran_kernels) THEN
    IF (not_init) THEN
!$OMP PARALLEL PRIVATE(ZTr,sub_tile_dx,sub_tile_dy)
!$OMP DO
      DO t=1,tiles_per_task
        sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
          chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
        sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
          chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

        ztr = 0.0_8

        CALL tea_leaf_dpcg_restrict_ZT_kernel(chunk(level)%tiles(t)%field%x_min, &
                                              chunk(level)%tiles(t)%field%x_max, &
                                              chunk(level)%tiles(t)%field%y_min,           &
                                              chunk(level)%tiles(t)%field%y_max,           &
                                              halo_exchange_depth,                  &
                                              chunk(level)%sub_tile_dims(1),               &
                                              sub_tile_dx,                          &
                                              chunk(level)%sub_tile_dims(2),               &
                                              sub_tile_dy,                          &
                                              chunk(level)%tiles(t)%field%vector_r,        &
                                              ztr)

        chunk(level)%def%t1(chunk(level)%tiles(t)%def_tile_coords(1): &
                            chunk(level)%tiles(t)%def_tile_coords(1)+chunk(level)%sub_tile_dims(1)-1, &
                            chunk(level)%tiles(t)%def_tile_coords(2): &
                            chunk(level)%tiles(t)%def_tile_coords(2)+chunk(level)%sub_tile_dims(2)-1) = &
        chunk(level)%def%t1(chunk(level)%tiles(t)%def_tile_coords(1): &
                            chunk(level)%tiles(t)%def_tile_coords(1)+chunk(level)%sub_tile_dims(1)-1, &
                            chunk(level)%tiles(t)%def_tile_coords(2): &
                            chunk(level)%tiles(t)%def_tile_coords(2)+chunk(level)%sub_tile_dims(2)-1) - ztr
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ELSE
!$OMP PARALLEL PRIVATE(ZTr,sub_tile_dx,sub_tile_dy)
!$OMP DO
      DO t=1,tiles_per_task
        sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
          chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
        sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
          chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

        ztr = 0.0_8

        CALL tea_leaf_dpcg_restrict_ZT_kernel(chunk(level)%tiles(t)%field%x_min, &
                                              chunk(level)%tiles(t)%field%x_max, &
                                              chunk(level)%tiles(t)%field%y_min,           &
                                              chunk(level)%tiles(t)%field%y_max,           &
                                              halo_exchange_depth,                  &
                                              chunk(level)%sub_tile_dims(1),               &
                                              sub_tile_dx,                          &
                                              chunk(level)%sub_tile_dims(2),               &
                                              sub_tile_dy,                          &
                                              chunk(level)%tiles(t)%field%vector_r,        &
                                              ztr)

        chunk(level)%def%t1(chunk(level)%tiles(t)%def_tile_coords(1): &
                            chunk(level)%tiles(t)%def_tile_coords(1)+chunk(level)%sub_tile_dims(1)-1, &
                            chunk(level)%tiles(t)%def_tile_coords(2): &
                            chunk(level)%tiles(t)%def_tile_coords(2)+chunk(level)%sub_tile_dims(2)-1) = - ztr
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDIF
  ENDIF

  CALL MPI_Allreduce(chunk(level)%def%t1, chunk(level)%def%t2, size(chunk(level)%def%t2), &
    MPI_DOUBLE_PRECISION, MPI_SUM, mpi_cart_comm, err)

END SUBROUTINE tea_leaf_dpcg_restrict_ZT

SUBROUTINE tea_leaf_dpcg_restrict_ZT_level(level,not_init)

  IMPLICIT NONE
  LOGICAL :: not_init
  INTEGER ::level

  INTEGER :: t, err
  INTEGER :: sub_tile_dx, sub_tile_dy
  REAL(KIND=8),dimension(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: ZTr

  IF (use_fortran_kernels) THEN
    IF (not_init) THEN
!$OMP PARALLEL PRIVATE(ZTr,sub_tile_dx,sub_tile_dy)
!$OMP DO
      DO t=1,tiles_per_task
        sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
          chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
        sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
          chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

        ztr = 0.0_8

        CALL tea_leaf_dpcg_restrict_ZT_kernel(chunk(level)%tiles(t)%field%x_min, &
                                              chunk(level)%tiles(t)%field%x_max, &
                                              chunk(level)%tiles(t)%field%y_min,           &
                                              chunk(level)%tiles(t)%field%y_max,           &
                                              halo_exchange_depth,                  &
                                              chunk(level)%sub_tile_dims(1),               &
                                              sub_tile_dx,                          &
                                              chunk(level)%sub_tile_dims(2),               &
                                              sub_tile_dy,                          &
                                              chunk(level)%tiles(t)%field%vector_r,        &
                                              ztr)

        chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                        1:chunk(level)%sub_tile_dims(2)) = &
        chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                        1:chunk(level)%sub_tile_dims(2)) - ztr
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ELSE
!$OMP PARALLEL PRIVATE(ZTr,sub_tile_dx,sub_tile_dy)
!$OMP DO
      DO t=1,tiles_per_task
        sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
          chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
        sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
          chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

        ztr = 0.0_8

        CALL tea_leaf_dpcg_restrict_ZT_kernel(chunk(level)%tiles(t)%field%x_min, &
                                              chunk(level)%tiles(t)%field%x_max, &
                                              chunk(level)%tiles(t)%field%y_min,           &
                                              chunk(level)%tiles(t)%field%y_max,           &
                                              halo_exchange_depth,                  &
                                              chunk(level)%sub_tile_dims(1),               &
                                              sub_tile_dx,                          &
                                              chunk(level)%sub_tile_dims(2),               &
                                              sub_tile_dy,                          &
                                              chunk(level)%tiles(t)%field%vector_r,        &
                                              ztr)

        chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                        1:chunk(level)%sub_tile_dims(2)) = - ztr
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDIF
  ENDIF

END SUBROUTINE tea_leaf_dpcg_restrict_ZT_level

SUBROUTINE tea_leaf_dpcg_prolong_Z()

  IMPLICIT NONE

  INTEGER :: t
  INTEGER, PARAMETER :: level=1

  INTEGER :: sub_tile_dx, sub_tile_dy

  REAL(KIND=8), DIMENSION(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: t2_local

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(sub_tile_dx,sub_tile_dy,t2_local)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
        chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
        chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      t2_local = chunk(level)%def%t2(chunk(level)%tiles(t)%def_tile_coords(1): &
                                     chunk(level)%tiles(t)%def_tile_coords(1)+chunk(level)%sub_tile_dims(1)-1, &
                                     chunk(level)%tiles(t)%def_tile_coords(2): &
                                     chunk(level)%tiles(t)%def_tile_coords(2)+chunk(level)%sub_tile_dims(2)-1)
      CALL tea_leaf_dpcg_prolong_Z_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                          chunk(level)%tiles(t)%field%x_max,    &
                                          chunk(level)%tiles(t)%field%y_min,    &
                                          chunk(level)%tiles(t)%field%y_max,    &
                                          halo_exchange_depth,           &
                                          chunk(level)%sub_tile_dims(1),        &
                                          sub_tile_dx,                   &
                                          chunk(level)%sub_tile_dims(2),        &
                                          sub_tile_dy,                   &
                                          chunk(level)%tiles(t)%field%vector_z, &
                                          t2_local)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_prolong_Z

SUBROUTINE tea_leaf_dpcg_prolong_Z_level(level)

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: t
  INTEGER :: sub_tile_dx, sub_tile_dy

  REAL(KIND=8), DIMENSION(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: t2_local

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(sub_tile_dx,sub_tile_dy,t2_local)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
        chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
        chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      t2_local = chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                                 1:chunk(level)%sub_tile_dims(2))
      CALL tea_leaf_dpcg_prolong_Z_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                          chunk(level)%tiles(t)%field%x_max,    &
                                          chunk(level)%tiles(t)%field%y_min,    &
                                          chunk(level)%tiles(t)%field%y_max,    &
                                          halo_exchange_depth,                  &
                                          chunk(level)%sub_tile_dims(1),        &
                                          sub_tile_dx,                          &
                                          chunk(level)%sub_tile_dims(2),        &
                                          sub_tile_dy,                          &
                                          chunk(level)%tiles(t)%field%vector_z, &
                                          t2_local)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_prolong_Z_level

SUBROUTINE tea_leaf_dpcg_subtract_z()

  IMPLICIT NONE
  INTEGER :: t
  INTEGER, PARAMETER :: level=1

  INTEGER :: sub_tile_dx, sub_tile_dy

  REAL(KIND=8), DIMENSION(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: t2_local

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(sub_tile_dx,sub_tile_dy,t2_local)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
        chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
        chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      t2_local=chunk(level)%def%t2(chunk(level)%tiles(t)%def_tile_coords(1): &
                                   chunk(level)%tiles(t)%def_tile_coords(1)+chunk(level)%sub_tile_dims(1)-1, &
                                   chunk(level)%tiles(t)%def_tile_coords(2): &
                                   chunk(level)%tiles(t)%def_tile_coords(2)+chunk(level)%sub_tile_dims(2)-1)
      CALL tea_leaf_dpcg_subtract_z_kernel(chunk(level)%tiles(t)%field%x_min, &
                                           chunk(level)%tiles(t)%field%x_max, &
                                           chunk(level)%tiles(t)%field%y_min, &
                                           chunk(level)%tiles(t)%field%y_max, &
                                           halo_exchange_depth,        &
                                           chunk(level)%sub_tile_dims(1),     &
                                           sub_tile_dx,                &
                                           chunk(level)%sub_tile_dims(2),     &
                                           sub_tile_dy,                &
                                           chunk(level)%tiles(t)%field%u,     &
                                           t2_local)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_subtract_z

SUBROUTINE tea_leaf_dpcg_subtract_z_level(level)

  IMPLICIT NONE
  INTEGER :: level
  INTEGER :: t

  INTEGER :: sub_tile_dx, sub_tile_dy

  REAL(KIND=8), DIMENSION(chunk(level)%sub_tile_dims(1), chunk(level)%sub_tile_dims(2)) :: t2_local

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(sub_tile_dx,sub_tile_dy,t2_local)
!$OMP DO
    DO t=1,tiles_per_task
      sub_tile_dx=(chunk(level)%tiles(t)%field%x_max-chunk(level)%tiles(t)%field%x_min+ &
        chunk(level)%sub_tile_dims(1))/chunk(level)%sub_tile_dims(1)
      sub_tile_dy=(chunk(level)%tiles(t)%field%y_max-chunk(level)%tiles(t)%field%y_min+ &
        chunk(level)%sub_tile_dims(2))/chunk(level)%sub_tile_dims(2)

      t2_local=chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
                                               1:chunk(level)%sub_tile_dims(2))
      CALL tea_leaf_dpcg_subtract_z_kernel(chunk(level)%tiles(t)%field%x_min, &
                                           chunk(level)%tiles(t)%field%x_max, &
                                           chunk(level)%tiles(t)%field%y_min, &
                                           chunk(level)%tiles(t)%field%y_max, &
                                           halo_exchange_depth,        &
                                           chunk(level)%sub_tile_dims(1),     &
                                           sub_tile_dx,                &
                                           chunk(level)%sub_tile_dims(2),     &
                                           sub_tile_dy,                &
                                           chunk(level)%tiles(t)%field%u,     &
                                           t2_local)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_subtract_z_level

SUBROUTINE tea_leaf_dpcg_init_p()

  IMPLICIT NONE

  INTEGER :: t, level

  level=1
  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_init_p_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                       chunk(level)%tiles(t)%field%x_max,    &
                                       chunk(level)%tiles(t)%field%y_min,    &
                                       chunk(level)%tiles(t)%field%y_max,    &
                                       halo_exchange_depth,           &
                                       chunk(level)%tiles(t)%field%vector_p, &
                                       chunk(level)%tiles(t)%field%vector_z)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_init_p

SUBROUTINE tea_leaf_dpcg_store_r()

  IMPLICIT NONE
  INTEGER :: t, level

  level=1
  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_store_r_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                        chunk(level)%tiles(t)%field%x_max,    &
                                        chunk(level)%tiles(t)%field%y_min,    &
                                        chunk(level)%tiles(t)%field%y_max,    &
                                        halo_exchange_depth,           &
                                        chunk(level)%tiles(t)%field%vector_r, &
                                        chunk(level)%tiles(t)%field%vector_r_m1 )
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_store_r

SUBROUTINE tea_leaf_dpcg_calc_rrn(rrn)

  IMPLICIT NONE
  INTEGER :: t, level
  REAL(KIND=8) :: rrn, tile_rrn

  rrn = 0.0_8

  level=1
  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rrn)
!$OMP DO REDUCTION(+:rrn)
    DO t=1,tiles_per_task
      tile_rrn = 0.0_8

      CALL tea_leaf_dpcg_calc_rrn_kernel(chunk(level)%tiles(t)%field%x_min,       &
                                         chunk(level)%tiles(t)%field%x_max,       &
                                         chunk(level)%tiles(t)%field%y_min,       &
                                         chunk(level)%tiles(t)%field%y_max,       &
                                         halo_exchange_depth,              &
                                         chunk(level)%tiles(t)%field%vector_r,    &
                                         chunk(level)%tiles(t)%field%vector_r_m1, &
                                         chunk(level)%tiles(t)%field%vector_z,    &
                                         tile_rrn)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_rrn

SUBROUTINE tea_leaf_dpcg_calc_p(beta)

  IMPLICIT NONE
  INTEGER :: t, level
  REAL(KIND=8) :: beta

  level=1
  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_dpcg_calc_p_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                       chunk(level)%tiles(t)%field%x_max,    &
                                       chunk(level)%tiles(t)%field%y_min,    &
                                       chunk(level)%tiles(t)%field%y_max,    &
                                       halo_exchange_depth,           &
                                       chunk(level)%tiles(t)%field%vector_p, &
                                       chunk(level)%tiles(t)%field%vector_z, &
                                       beta)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_p

SUBROUTINE tea_leaf_dpcg_calc_zrnorm(rro)

  IMPLICIT NONE

  INTEGER :: t, level
  REAL(KIND=8) :: rro, tile_rro

  rro = 0.0_8

  level=1
  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rro)
!$OMP DO REDUCTION(+:rro)
    DO t=1,tiles_per_task
      tile_rro = 0.0_8

      CALL tea_leaf_dpcg_calc_zrnorm_kernel(chunk(level)%tiles(t)%field%x_min,    &
                                            chunk(level)%tiles(t)%field%x_max,    &
                                            chunk(level)%tiles(t)%field%y_min,    &
                                            chunk(level)%tiles(t)%field%y_max,    &
                                            halo_exchange_depth,           &
                                            chunk(level)%tiles(t)%field%vector_z, &
                                            chunk(level)%tiles(t)%field%vector_r, &
                                            tile_rro)

      rro = rro + tile_rro
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_dpcg_calc_zrnorm

SUBROUTINE tea_leaf_dpcg_local_solve(x_min,               &
                                     x_max,               &
                                     y_min,               &
                                     y_max,               &
                                     halo_exchange_depth, &
                                     u,                   &
                                     u0,                  &
                                     def_kx,              &
                                     def_ky,              &
                                     def_di,              &
                                     p,                   &
                                     r,                   &
                                     Mi,                  &
                                     w,                   &
                                     z,                   &
                                     sd,                  &
                                     eps,                 &
                                     inner_iters,         &
                                     it_count,            &
                                     theta,               &
                                     use_ppcg,            &
                                     inner_cg_alphas,     &
                                     inner_cg_betas,      &
                                     inner_ch_alphas,     &
                                     inner_ch_betas)

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
  REAL(KIND=8), PARAMETER :: omega=1.0_8

  rro = 0.0_8
  initial_residual = 0.0_8
  pw = 0.0_8

  rrn = 1e10

  it_count = 0

!$OMP PARALLEL private(alpha, beta, smvp, inner_step)
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      ! copy the RHS vector
      u0(j, k) = u(j, k)
      ! zero for approximate solve
      u(j, k) = 0.0_8
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO REDUCTION(+:initial_residual)
  DO k=y_min, y_max
    DO j=x_min, x_max
      !smvp = (    def_di(j  , k  )*u(j  , k  ) &
      !         - (def_ky(j  , k+1)*u(j  , k+1) + def_ky(j  , k  )*u(j  , k-1)) &
      !         - (def_kx(j+1, k  )*u(j+1, k  ) + def_kx(j  , k  )*u(j-1, k  )) )

      Mi(j, k) = omega/def_di(j, k)

      r(j, k) = u0(j, k) !- smvp
      z(j, k) = r(j, k)*Mi(j, k)
      p(j, k) = z(j, k)

      initial_residual = initial_residual + r(j, k)*p(j, k)
    ENDDO
  ENDDO
!$OMP END DO
  write(6,*) "rnorm:",sqrt(sum(r**2))

  write(6,*) "rro:",rro
!$OMP SINGLE
    rro = initial_residual
    initial_residual = sqrt(abs(initial_residual))
!$OMP END SINGLE
    write(6,*) "initial_residual:",initial_residual

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
      smvp = (    def_di(j  , k  )*p(j  , k  ) &
               - (def_ky(j  , k+1)*p(j  , k+1) + def_ky(j  , k  )*p(j  , k-1)) &
               - (def_kx(j+1, k  )*p(j+1, k  ) + def_kx(j  , k  )*p(j-1, k  )) )

        w(j, k) = smvp
        pw = pw + smvp*p(j, k)
      ENDDO
    ENDDO
!$OMP END DO
    write(6,*) "normw,normp:",sqrt(sum(w**2)),sqrt(sum(p**2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    alpha = rro/pw
    write(6,*) "serial:",rro,pw

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

      DO inner_step=1,min(10,coarse_solve_max_iters)
!$OMP DO
        DO k=y_min,y_max
          DO j=x_min,x_max
            smvp = (    def_di(j  , k  )*sd(j  , k  ) &
                     - (def_ky(j  , k+1)*sd(j  , k+1) + def_ky(j  , k  )*sd(j  , k-1)) &
                     - (def_kx(j+1, k  )*sd(j+1, k  ) + def_kx(j  , k  )*sd(j-1, k  )) )


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
    if (parallel%boss) then
      !if (it_count == 1) write(6,*) use_ppcg
      write(6,'("serial iteration:",i3," alpha=",es20.13," beta=",es20.13," rrn=",es20.13," norm u:",es20.13," norm z:",es20.13)') &
        it_count+1,alpha,beta,rrn,sqrt(sum(u**2)),sqrt(sum(z**2))
    endif

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

    !write(6,*) "serial solve:",sqrt(abs(rrn)),eps,initial_residual
  ENDDO

!$OMP END PARALLEL

  if (parallel%boss) then
    !if (it_count == 1) write(6,*) use_ppcg
    write(6,'("serial iteration:",i3," alpha=",es20.13," beta=",es20.13," rrn=",es20.13," norm u:",es20.13," norm z:",es20.13)') &
      it_count,alpha,beta,rrn,sqrt(sum(u**2)),sqrt(sum(z**2))
  endif
  if (use_ppcg) write(6,*) theta,inner_ch_alphas(1:10),inner_ch_betas(1:10)
  if (it_count == 0) stop

END SUBROUTINE tea_leaf_dpcg_local_solve

SUBROUTINE tea_leaf_dpcg_local_solve_level(level,               &
                                           solve_time,          &
                                           eps_inner,           &
                                           inner_iters,         &
                                           n,                   &
                                           use_ppcg,            &
                                           theta,               &
                                           inner_ch_alphas,     &
                                           inner_ch_betas)

  USE tea_leaf_cg_module
  USE tea_leaf_ppcg_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_THREAD_NUM
  LOGICAL :: use_ppcg
  INTEGER :: level, inner_iters
  INTEGER :: t,ppcg_inner_iters,n
  INTEGER :: fields(NUM_FIELDS)
  REAL(KIND=8) :: solve_time,eps_inner
  REAL(KIND=8), DIMENSION(inner_iters) :: inner_ch_alphas, inner_ch_betas
  REAL(KIND=8) :: theta

  REAL(KIND=8) :: alpha, beta
  REAL(KIND=8) :: timer,halo_time,dot_product_time
  REAL(KIND=8) :: rro,initial_residual,pw,rrn,error,old_error
  REAL(KIND=8) :: normu,normz

!$OMP PARALLEL
!$OMP DO
  DO t=1,tiles_per_task
    chunk(level+1)%tiles(t)%field%u0 = chunk(level+1)%tiles(t)%field%u
    chunk(level+1)%tiles(t)%field%u  = 0.0_8
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  fields=0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(level+1, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  CALL tea_leaf_calc_residual(level+1)
  CALL tea_leaf_calc_2norm(level+1, 1, initial_residual)

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(initial_residual)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  IF (parallel%boss.AND.verbose_on) THEN
!$  IF (OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(g_out,*)"Coarse solve - Initial residual ",initial_residual
!$  ENDIF
  ENDIF
  !write(6,*) "rnorm:",sqrt(initial_residual)

  ! All 3 of these solvers use the CG kernels
  CALL tea_leaf_cg_init(level+1,rro)
  !write(6,*) "rro:",rro

  ! and globally sum rro
  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(rro)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  initial_residual = rro
  initial_residual=SQRT(initial_residual)
  old_error = initial_residual
  !write(6,*) "initial_residual:",initial_residual

  ! need to update p when using CG due to matrix/vector multiplication
  fields=0
  fields(FIELD_U) = 1
  fields(FIELD_P) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(level+1,fields,1)
  IF (profiler_on) solve_time=solve_time+(timer()-halo_time)

  fields=0
  fields(FIELD_P) = 1

  ppcg_inner_iters = 0
  DO n=1,max_iters

    CALL tea_leaf_cg_calc_w(level+1, pw)
    !normu = 0.0_8; normz = 0.0_8
    !DO t=1,tiles_per_task
    !  normu = normu+sum(chunk(level+1)%tiles(t)%field%vector_w(1:chunk(level)%sub_tile_dims(1), &
    !                                                           1:chunk(level)%sub_tile_dims(2))**2)
    !  normz = normz+sum(chunk(level+1)%tiles(t)%field%vector_p(1:chunk(level)%sub_tile_dims(1), &
    !                                                           1:chunk(level)%sub_tile_dims(2))**2)
    !ENDDO
    !write(6,*) "normw,normp:",sqrt(normu),sqrt(normz)

    IF (profiler_on) dot_product_time=timer()
    CALL tea_allsum(pw)
    IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

    alpha = rro/pw
    !write(6,*) "level:",rro,pw
    if  (alpha /= alpha) then
      write(6,*) "alpha not finite:",n,alpha
      stop
    endif

    CALL tea_leaf_cg_calc_ur(level+1, alpha, rrn)

    ! not using rrn, so don't do a tea_allsum

    IF (use_ppcg) THEN
      CALL tea_leaf_run_dpcg_inner_steps(level+1, inner_ch_alphas, inner_ch_betas, theta, &
          tl_ppcg_inner_steps, solve_time)
      ppcg_inner_iters = ppcg_inner_iters + tl_ppcg_inner_steps
    ENDIF

    CALL tea_leaf_ppcg_calc_zrnorm(level+1, rrn)

    IF (profiler_on) dot_product_time=timer()
    CALL tea_allsum(rrn)
    IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

    beta = rrn/rro
    !normu = 0.0_8; normz = 0.0_8
    !DO t=1,tiles_per_task
    !  normu = normu+sum(chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
    !                                                    1:chunk(level)%sub_tile_dims(2))**2)
    !  normz = normz+sum(chunk(level+1)%tiles(t)%field%vector_z(1:chunk(level)%sub_tile_dims(1), &
    !                                                           1:chunk(level)%sub_tile_dims(2))**2)
    !ENDDO
    !CALL tea_allsum(normu)
    !CALL tea_allsum(normz)
    !if (parallel%boss) &
    !  write(6,'("level  iteration:",i3," alpha=",es20.13," beta=",es20.13," rrn=",es20.13," norm u:",es20.13," norm z:",es20.13)') &
    !    n,alpha,beta,rrn,sqrt(normu),sqrt(normz)

    CALL tea_leaf_cg_calc_p(level+1, beta)

    ! updates u and possibly p
    IF (profiler_on) halo_time = timer()
    CALL update_halo(level+1,fields,1)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    error = rrn
    rro = rrn

    error=SQRT(error)

    IF (parallel%boss.AND.verbose_on) THEN
!$    IF (OMP_GET_THREAD_NUM().EQ.0) THEN
        WRITE(g_out,*)"Coarse solve - Residual ",error
!$    ENDIF
    ENDIF

    !write(6,*) "level  solve:",error,eps_inner,initial_residual
    IF (ABS(error) .LT. eps_inner*initial_residual) EXIT

    old_error = error

  ENDDO

  !normu = 0.0_8; normz = 0.0_8
  !DO t=1,tiles_per_task
  !  normu = normu+sum(chunk(level+1)%tiles(t)%field%u(1:chunk(level)%sub_tile_dims(1), &
  !                                                    1:chunk(level)%sub_tile_dims(2))**2)
  !  normz = normz+sum(chunk(level+1)%tiles(t)%field%vector_z(1:chunk(level)%sub_tile_dims(1), &
  !                                                           1:chunk(level)%sub_tile_dims(2))**2)
  !ENDDO
  !CALL tea_allsum(normu)
  !CALL tea_allsum(normz)
  !if (parallel%boss) &
  !  write(6,'("level  iteration:",i3," alpha=",es20.13," beta=",es20.13," rrn=",es20.13," norm u:",es20.13," norm z:",es20.13)') &
  !    n,alpha,beta,rrn,sqrt(normu),sqrt(normz)

END SUBROUTINE tea_leaf_dpcg_local_solve_level

SUBROUTINE tea_leaf_run_dpcg_inner_steps(level, ch_alphas, ch_betas, theta, &
    tl_ppcg_inner_steps, solve_time)

  USE tea_leaf_ppcg_module, only : tea_leaf_ppcg_init_sd,tea_leaf_ppcg_inner

  IMPLICIT NONE

  INTEGER :: level, fields(NUM_FIELDS)
  INTEGER :: tl_ppcg_inner_steps, ppcg_cur_step
  REAL(KIND=8) :: theta
  REAL(KIND=8) :: halo_time, timer, solve_time
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas

  INTEGER(KIND=4) :: inner_step, bounds_extra

  fields = 0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(level, fields,1)
  IF (profiler_on) solve_time = solve_time + (timer() - halo_time)

  CALL tea_leaf_ppcg_init_sd(level, theta)

  ! inner steps
  DO ppcg_cur_step=1,tl_ppcg_inner_steps,halo_exchange_depth

    fields = 0
    fields(FIELD_SD) = 1
    fields(FIELD_R) = 1

    IF (profiler_on) halo_time = timer()
    CALL update_halo(level, fields,halo_exchange_depth)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    inner_step = ppcg_cur_step

    fields = 0
    fields(FIELD_SD) = 1

    DO bounds_extra = halo_exchange_depth-1, 0, -1
      CALL tea_leaf_ppcg_inner(level, ch_alphas, ch_betas, inner_step, bounds_extra)

      IF (profiler_on) halo_time = timer()
      CALL update_boundary(level, fields, 1)
      IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

      inner_step = inner_step + 1
      IF (inner_step .gt. tl_ppcg_inner_steps) EXIT
    ENDDO
  ENDDO

  fields = 0
  fields(FIELD_P) = 1

END SUBROUTINE tea_leaf_run_dpcg_inner_steps

END MODULE tea_leaf_dpcg_module

