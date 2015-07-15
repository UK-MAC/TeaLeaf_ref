!Crown Copyright 2014 AWE.
!
! This file is part of TeaLeaf.
!
! TeaLeaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! TeaLeaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! TeaLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran heat conduction kernel
!>  @author Michael Boulton, Wayne Gaudin
!>  @details Implicitly calculates the change in temperature using accelerated Chebyshev method

MODULE tea_leaf_cheby_kernel_module

  USE tea_leaf_common_kernel_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_kernel_cheby_init(x_min,  &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           halo_exchange_depth,             &
                           u,                 &
                           u0,                &
                           p,                 &
                           r,                 &
                           Mi,                &
                           w,                 &
                           z,                 &
                           Kx,                &
                           Ky,                &
                           cp,                     &
                           bfp,                     &
                           rx,                &
                           ry,                &
                           theta,             &
                           preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u, u0, p , w , r, Mi, z , Kx, Ky
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER :: j,k
  REAL(KIND=8) ::  rx, ry, theta

!$OMP PARALLEL NUM_THREADS(INNER_NUM_THREADS)
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
                - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
                - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))
            r(j, k) = u0(j, k) - w(j, k)
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
    DO k=y_min,y_max
      DO j=x_min,x_max
        p(j, k) = z(j, k)/theta
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = r(j, k)/theta
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP DO
  DO k=y_min,y_max
      DO j=x_min,x_max
          u(j, k) = u(j, k) + p(j, k)
      ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE tea_leaf_kernel_cheby_iterate(x_min, &
                           x_max,               &
                           y_min,               &
                           y_max,               &
                           halo_exchange_depth,               &
                           u,                   &
                           u0,                  &
                           p,                   &
                           r,                   &
                           Mi,                  &
                           w                ,   &
                           z,                   &
                           Kx,                  &
                           Ky,                  &
                           cp,   &
                           bfp,    &
                           ch_alphas,           &
                           ch_betas,            &
                           max_cheby_iters,     &
                           rx,                  &
                           ry,                  &
                           step,                &
                           preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u, u0, p , w , r, Mi, z , Kx, Ky
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER :: j,k

    REAL(KIND=8) ::  rx, ry

    INTEGER :: step, max_cheby_iters
    REAL(KIND=8), DIMENSION(max_cheby_iters) :: ch_alphas, ch_betas

!$OMP PARALLEL NUM_THREADS(INNER_NUM_THREADS)
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
                - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
                - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))
            r(j, k) = u0(j, k) - w(j, k)
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
    DO k=y_min,y_max
      DO j=x_min,x_max
        p(j, k) = ch_alphas(step)*p(j, k) + ch_betas(step)*z(j, k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        p(j, k) = ch_alphas(step)*p(j, k) + ch_betas(step)*r(j, k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
        u(j, k) = u(j, k) + p(j, k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_cheby_iterate

END MODULE

