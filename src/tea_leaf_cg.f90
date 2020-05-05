
MODULE tea_leaf_cg_module

  USE tea_leaf_cg_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_cg_init(rro)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: rro,tile_rro

  rro = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rro) REDUCTION(+:rro)
!$OMP DO
    DO t=1,tiles_per_task
      tile_rro = 0.0_8

      CALL tea_leaf_cg_init_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          chunk%halo_exchange_depth,                                   &
          chunk%tiles(t)%field%vector_p,                               &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_Mi,                              &
          chunk%tiles(t)%field%vector_z,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%vector_Di,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          tile_rro, tl_preconditioner_type)

      rro = rro + tile_rro
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_init

SUBROUTINE tea_leaf_cg_calc_w( pw)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: pw, tile_pw

  pw = 0.0_08

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_pw) REDUCTION(+:pw)
!$OMP DO
    DO t=1,tiles_per_task
      tile_pw = 0.0_8

      CALL tea_leaf_cg_calc_w_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                         &
          chunk%tiles(t)%field%y_min,                                         &
          chunk%tiles(t)%field%y_max,                                         &
          chunk%halo_exchange_depth,                                          &
          chunk%tiles(t)%field%vector_p,                                      &
          chunk%tiles(t)%field%vector_w,                                      &
          chunk%tiles(t)%field%vector_Kx,                                     &
          chunk%tiles(t)%field%vector_Ky,                                     &
          chunk%tiles(t)%field%vector_Di,                                     &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry,  &
          tile_pw)

      pw = pw + tile_pw
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF


END SUBROUTINE tea_leaf_cg_calc_w

SUBROUTINE tea_leaf_cg_calc_ur( alpha, rrn)

  IMPLICIT NONE
  
  INTEGER :: t
  REAL(KIND=8) :: alpha, rrn, tile_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(tile_rrn) REDUCTION(+:rrn)
!$OMP DO
    DO t=1,tiles_per_task
      tile_rrn = 0.0_8

      CALL tea_leaf_cg_calc_ur_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                          &
          chunk%tiles(t)%field%y_min,                                          &
          chunk%tiles(t)%field%y_max,                                          &
          chunk%halo_exchange_depth,                                           &
          chunk%tiles(t)%field%u,                                              &
          chunk%tiles(t)%field%vector_p,                                       &
          chunk%tiles(t)%field%vector_r,                                       &
          chunk%tiles(t)%field%vector_Mi,                                      &
          chunk%tiles(t)%field%vector_w,                                       &
          chunk%tiles(t)%field%vector_z,                                       &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%vector_Di,                              &
          chunk%tiles(t)%field%rx,  &
          chunk%tiles(t)%field%ry, &
          alpha, tile_rrn, tl_preconditioner_type)

      rrn = rrn + tile_rrn
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_ur

SUBROUTINE tea_leaf_cg_calc_p( beta)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: beta

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_calc_p_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                         &
          chunk%tiles(t)%field%y_min,                                         &
          chunk%tiles(t)%field%y_max,                                         &
          chunk%halo_exchange_depth,                                          &
          chunk%tiles(t)%field%vector_p,                                      &
          chunk%tiles(t)%field%vector_r,                                      &
          chunk%tiles(t)%field%vector_z,                                      &
          beta, tl_preconditioner_type)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_p

END MODULE tea_leaf_cg_module

