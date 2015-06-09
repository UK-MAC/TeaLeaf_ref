!cROWn Copyright 2014 AWE.
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
!>  @details Implicitly calculates the change in temperature using CG method

MODULE tea_leaf_kernel_cg_module

USE tea_leaf_kernel_common_module

CONTAINS

SUBROUTINE tea_leaf_kernel_init_cg_fortran(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           density,                &
                           energy,                 &
                           u,                      &
                           p,                      &
                           r,                      &
                           Mi,                     &
                           w,                      &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           cp,                     &
                           bfp,                    &
                           rx,                     &
                           ry,                     &
                           rro,                    &
                           preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u, r, w, Kx, Ky, z, Mi, density, energy, p
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: j,k

  REAL(kind=8) :: rro
  REAL(KIND=8) ::  rx, ry

  rro = 0.0_8

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k) = 0.0_8
      z(j, k) = 0.0_8
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
            p(j, k) = z(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = r(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP DO REDUCTION(+:rro)
  DO k=y_min,y_max
    DO j=x_min,x_max
      rro = rro + r(j, k)*p(j, k);
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_init_cg_fortran

SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_w(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   halo_exchange_depth,             &
                                                   p,                 &
                                                   w,                 &
                                                   Kx,                &
                                                   Ky,                &
                                                   rx,                &
                                                   ry,                &
                                                   pw                 )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: w, Kx, Ky, p

    REAL(KIND=8) ::  rx, ry

    INTEGER(KIND=4) :: j,k
    REAL(kind=8) :: pw

!$OMP PARALLEL
!$OMP DO REDUCTION(+:pw)
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))*p(j, k)             &
                - ry*(Ky(j, k+1)*p(j, k+1) + Ky(j, k)*p(j, k-1))  &
                - rx*(Kx(j+1, k)*p(j+1, k) + Kx(j, k)*p(j-1, k))

            pw = pw + w(j, k)*p(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_w

SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_ur(x_min,             &
                                                    x_max,             &
                                                    y_min,             &
                                                    y_max,             &
                                                    halo_exchange_depth,             &
                                                    u,                 &
                                                    p,                 &
                                                    r,                 &
                                                    Mi,                &
                                                    w,                 &
                                                    z,                 &
                                                    cp,                     &
                                                    bfp,                     &
                                                    Kx, &
                                                    Ky, &
                                                    rx, &
                                                    ry, &
                                                    alpha,             &
                                                    rrn,               &
                                                    preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u, r, Mi, w, z, Kx, Ky, p
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: rx, ry

    INTEGER(KIND=4) :: j,k
    REAL(kind=8) :: alpha, rrn

!$OMP PARALLEL
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            u(j, k) = u(j, k) + alpha*p(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT

  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

!$OMP DO
    DO k=y_min, y_max
      DO j=x_min,x_max
        r(j, k) = r(j, k) - alpha*w(j, k)
      ENDDO
    ENDDO
!$OMP END DO

    IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
      CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, cp, bfp, Kx, Ky, rx, ry)
    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
      CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, Mi, Kx, Ky, rx, ry)
    ENDIF

!$OMP DO REDUCTION(+:rrn)
    DO k=y_min,y_max
        DO j=x_min,x_max
            rrn = rrn + r(j, k)*z(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  ELSE
!$OMP DO REDUCTION(+:rrn)
    DO k=y_min,y_max
        DO j=x_min,x_max
            r(j, k) = r(j, k) - alpha*w(j, k)
            rrn = rrn + r(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_ur

SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_p(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   halo_exchange_depth,             &
                                                   p,                 &
                                                   r,                 &
                                                   z,                 &
                                                   beta,              &
                                                   preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: z, r, p

  INTEGER(KIND=4) :: j,k
  REAL(kind=8) :: beta

!$OMP PARALLEL
  IF (preconditioner_type .NE. TL_PREC_NONE) THEN
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = z(j, k) + beta*p(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = r(j, k) + beta*p(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_p

END MODULE

