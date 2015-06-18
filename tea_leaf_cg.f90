
MODULE tea_leaf_cg_module

  USE tea_leaf_cg_kernel_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_cg_init(rx, ry, rro)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry, rx, rro, private_rro

  rro = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(private_rro)
!$OMP DO REDUCTION(+:rro)
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_init_kernel(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%tiles(t)%field%density,                                &
          chunk%tiles(t)%field%energy1,                                &
          chunk%tiles(t)%field%u,                                      &
          chunk%tiles(t)%field%vector_p,                               &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_Mi,                              &
          chunk%tiles(t)%field%vector_w,                               &
          chunk%tiles(t)%field%vector_z,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          rx, ry, private_rro, tl_preconditioner_type)

      rro = rro + private_rro
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_init

SUBROUTINE tea_leaf_cg_calc_w(rx, ry, pw)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry, rx, pw, private_pw

  pw = 0.0_08

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(private_pw)
!$OMP DO REDUCTION(+:pw)
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_calc_w_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                         &
          chunk%tiles(t)%field%y_min,                                         &
          chunk%tiles(t)%field%y_max,                                         &
          halo_exchange_depth,                                         &
          chunk%tiles(t)%field%vector_p,                                      &
          chunk%tiles(t)%field%vector_w,                                      &
          chunk%tiles(t)%field%vector_Kx,                                     &
          chunk%tiles(t)%field%vector_Ky,                                     &
          rx, ry, private_pw)

      pw = pw + private_pw
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_w

SUBROUTINE tea_leaf_cg_calc_ur(rx, ry, alpha, rrn)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry, rx, alpha, rrn, private_rrn

  rrn = 0.0_8

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(private_rrn)
!$OMP DO REDUCTION(+:rrn)
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_calc_ur_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                          &
          chunk%tiles(t)%field%y_min,                                          &
          chunk%tiles(t)%field%y_max,                                          &
          halo_exchange_depth,                                          &
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
          rx, ry, &
          alpha, private_rrn, tl_preconditioner_type)

      rrn = rrn + private_rrn
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_ur

SUBROUTINE tea_leaf_cg_calc_p(rx, ry, beta)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: ry, rx, beta

  IF (use_fortran_kernels) THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL tea_leaf_cg_calc_p_kernel(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                                         &
          chunk%tiles(t)%field%y_min,                                         &
          chunk%tiles(t)%field%y_max,                                         &
          halo_exchange_depth,                                         &
          chunk%tiles(t)%field%vector_p,                                      &
          chunk%tiles(t)%field%vector_r,                                      &
          chunk%tiles(t)%field%vector_z,                                      &
          beta, tl_preconditioner_type)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE tea_leaf_cg_calc_p

END MODULE tea_leaf_cg_module

