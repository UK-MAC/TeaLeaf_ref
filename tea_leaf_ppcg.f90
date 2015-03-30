MODULE tea_leaf_kernel_ppcg_module

USE tea_leaf_kernel_common_module
USE tea_leaf_kernel_cheby_module

IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_kernel_ppcg_init_sd(x_min,             &
                                        x_max,             &
                                        y_min,             &
                                        y_max,             &
                                        r,                 &
                                        kx,                 &
                                        ky,                 &
                                        sd,                &
                                        z,                &
                                        cp,                &
                                        bfp,                &
                                        rx, ry,             &
                                        theta,             &
                                        preconditioner_on)

  IMPLICIT NONE

  LOGICAL :: preconditioner_on
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: r, sd, kx, ky
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: z, cp, bfp
  REAL(KIND=8) :: theta, theta_r, rx, ry

  INTEGER :: j,k

  theta_r = 1.0_8/theta

!$OMP PARALLEL
  if (preconditioner_on) then

    call tea_block_solve(x_min, x_max, y_min, y_max,             &
                        r, z,                 &
                        cp,                     &
                        bfp,                     &
                        Kx, Ky, rx, ry)

!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            sd(j, k) = z(j, k)*theta_r
        ENDDO
    ENDDO
!$OMP END DO
  else
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            sd(j, k) = r(j, k)*theta_r
        ENDDO
    ENDDO
!$OMP END DO
  endif
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_init_sd

SUBROUTINE tea_leaf_kernel_ppcg_inner(x_min,             &
                                      x_max,             &
                                      y_min,             &
                                      y_max,             &
                                      step,              &
                                      alpha,             &
                                      beta,              &
                                      rx, ry,            &
                                      u,                 &
                                      r,                 &
                                      Kx,                &
                                      Ky,                &
                                      sd,                &
                                      z,                &
                                      cp,                &
                                      bfp,                &
                                      preconditioner_on)
  IMPLICIT NONE

  LOGICAL :: preconditioner_on
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u, r, Kx, Ky, sd
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: z, cp, bfp
  INTEGER(KIND=4) :: j,k, step, upper_k, ko
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry

!$OMP PARALLEL
!$OMP DO private(upper_k)
    !DO k=y_min,y_max
    DO ko=y_min, y_max, kstep
        upper_k = min(ko+kstep - 1, y_max)
        do k=ko,upper_k
        DO j=x_min,x_max
            smvp = (1.0_8                                           &
                + ry*(Ky(j, k+1) + Ky(j, k))                        &
                + rx*(Kx(j+1, k) + Kx(j, k)))*sd(j, k)              &
                - ry*(Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
                - rx*(Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))
            r(j, k) = r(j, k) - smvp

            u(j, k) = u(j, k) + sd(j, k)
        ENDDO
        ENDDO
    ENDDO
!$OMP END DO

  if (preconditioner_on) then

    call tea_block_solve(x_min, x_max, y_min, y_max,             &
                        r, z,                 &
                        cp,                     &
                        bfp,                     &
                        Kx, Ky, rx, ry)

!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            sd(j, k) = alpha(step)*sd(j, k) + beta(step)*z(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  else
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            sd(j, k) = alpha(step)*sd(j, k) + beta(step)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  endif
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          z, r,               &
                          preconditioner_on,    &
                          norm)

  IMPLICIT NONE

  LOGICAL :: preconditioner_on
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: z
  REAL(KIND=8) :: norm
  integer :: j, k

  norm = 0.0_8

!$OMP PARALLEL
  if (preconditioner_on) then
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + z(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  else
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + r(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  endif
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


