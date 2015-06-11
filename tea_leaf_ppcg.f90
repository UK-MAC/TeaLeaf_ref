MODULE tea_leaf_kernel_ppcg_module

USE tea_leaf_kernel_common_module
USE tea_leaf_kernel_cheby_module

IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_kernel_ppcg_init_sd(x_min,             &
                                        x_max,             &
                                        y_min,             &
                                        y_max,             &
                                        halo_exchange_depth,             &
                                        r,                 &
                                        kx,                 &
                                        ky,                 &
                                        sd,                &
                                        z,                &
                                        cp,                &
                                        bfp,                &
                                        Mi,                &
                                        rx, ry,             &
                                        theta,             &
                                        preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, sd, kx, ky , z, Mi
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: theta, theta_r, rx, ry

  INTEGER :: j,k

  theta_r = 1.0_8/theta

!$OMP PARALLEL

  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

    IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
      CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, cp, bfp, Kx, Ky, rx, ry)
    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
      CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, Mi, Kx, Ky, rx, ry)
    ENDIF

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k) = z(j, k)*theta_r
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k) = r(j, k)*theta_r
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_init_sd

SUBROUTINE tea_leaf_kernel_ppcg_inner(x_min,             &
                                      x_max,             &
                                      y_min,             &
                                      y_max,             &
                                      halo_exchange_depth,             &
                                      x_min_bound,      &
                                      x_max_bound,      &
                                      y_min_bound,      &
                                      y_max_bound,      &
                                      alpha,             &
                                      beta,              &
                                      rx, ry,            &
                                      inner_step,       &
                                      u,                 &
                                      r,                 &
                                      Kx,                &
                                      Ky,                &
                                      sd,                &
                                      z,                &
                                      cp,                &
                                      bfp,                &
                                      Mi,                &
                                      preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u, r, Kx, Ky, sd , z, Mi
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  INTEGER(KIND=4) :: j,k
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, inner_step

!$OMP PARALLEL PRIVATE(smvp)

!$OMP DO
    DO k=y_min_bound,y_max_bound
        DO j=x_min_bound,x_max_bound
            smvp = (1.0_8                                           &
                + ry*(Ky(j, k+1) + Ky(j, k))                        &
                + rx*(Kx(j+1, k) + Kx(j, k)))*sd(j, k)              &
                - ry*(Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
                - rx*(Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))

            r(j, k) = r(j, k) - smvp
            u(j, k) = u(j, k) + sd(j, k)
        ENDDO
    ENDDO
!$OMP END DO

    IF (preconditioner_type .NE. TL_PREC_NONE) THEN
      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                               r, z, cp, bfp, Kx, Ky, rx, ry)
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
        CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                               r, z, Mi, Kx, Ky, rx, ry)
      ENDIF
  
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*z(j, k)
          ENDDO
      ENDDO
!$OMP END DO
    ELSE
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*r(j, k)
          ENDDO
      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth,             &
                          z, r,               &
                          preconditioner_type,    &
                          norm)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, z
  REAL(KIND=8) :: norm
  integer :: j, k

  norm = 0.0_8

!$OMP PARALLEL
  IF (preconditioner_type .NE. TL_PREC_NONE) THEN
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + z(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  ELSE
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + r(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL

end SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel

SUBROUTINE tea_calc_ls_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                             theta, ppcg_inner_steps)

  INTEGER :: ppcg_inner_steps
  REAL(KIND=8), DIMENSION(ppcg_inner_steps) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax, theta

  ! TODO
  CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                         theta, ppcg_inner_steps)

END SUBROUTINE

END MODULE


