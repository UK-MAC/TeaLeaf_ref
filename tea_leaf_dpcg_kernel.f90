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

MODULE tea_leaf_dpcg_kernel_module

  USE tea_leaf_common_kernel_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_dpcg_zero_t2_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           t2)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: t2

  INTEGER(KIND=4) :: j,k

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      t2(j, k) = 0.0_8
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_zero_t2_kernel

SUBROUTINE tea_leaf_dpcg_sum_r_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           rx, ry, &
                           ztr, &
                           e, &
                           kx, &
                           ky, &
                           r)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, kx, ky

  INTEGER(KIND=4) :: j,k

  REAL(kind=8) :: rx, ry
  REAL(kind=8) :: ztr, e

  ztr = 0.0_8
  e = 0.0_8

!$OMP PARALLEL REDUCTION(+:ztr, e)
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      ztr = ztr + r(j, k)
      e = e + (1.0_8                            &
            + (ry*Ky(j, k+1) + ry*Ky(j, k))     &
            + (rx*Kx(j+1, k) + rx*Kx(j, k)))    &
            - (ry*Ky(j, k+1) + ry*Ky(j, k))     &
            - (rx*Kx(j+1, k) + rx*Kx(j, k))
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_sum_r_kernel

SUBROUTINE tea_leaf_dpcg_local_solve(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           u,                      &
                           u0,                      &
                           def_e,                      &
                           p,                      &
                           r,                      &
                           Mi,                     &
                           w,                     &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           cp,                     &
                           bfp,                    &
                           rx,                     &
                           ry,                     &
                           preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: r, Kx, Ky, z, Mi, p, w, u, u0, def_e
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: j,k, lcount

  REAL(kind=8) :: rro, smvp
  REAL(KIND=8) ::  rx, ry
  REAL(KIND=8) ::  alpha, beta, pw, rrn

  rro = 0.0_8
  rrn = 0.0_8
  pw = 0.0_8

!$OMP PARALLEL private(alpha, beta)

!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      u(j, k) = u0(j, k)
      p(j, k) = 0.0_8
      z(j, k) = 0.0_8
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      IF (j .GT. x_min) Kx(j,k)=(def_e(j-1,k  ) + def_e(j,k))/(2.0_8*def_e(j-1,k  )*def_e(j,k))
      IF (k .GT. y_min) Ky(j,k)=(def_e(j  ,k-1) + def_e(j,k))/(2.0_8*def_e(j  ,k-1)*def_e(j,k))
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
    DO k=y_min, y_max
      DO j=x_min, x_max
        smvp = (1.0_8                                         &
            + (Ky(j, k+1) + Ky(j, k))                      &
            + (Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
            - (Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
            - (Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))
        r(j, k) = u0(j, k) - smvp
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
!$OMP END DO
  ELSE
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        p(j, k) = r(j, k)
      ENDDO
    ENDDO
!$OMP END DO
  ENDIF
!$OMP DO REDUCTION(+:rro)
  DO k=y_min,y_max
    DO j=x_min,x_max
      rro = rro + r(j, k)*p(j, k);
    ENDDO
  ENDDO
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO lcount=1,20

!$OMP DO REDUCTION(+:pw)
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + (Ky(j, k+1) + Ky(j, k))                      &
                + (Kx(j+1, k) + Kx(j, k)))*p(j, k)             &
                - (Ky(j, k+1)*p(j, k+1) + Ky(j, k)*p(j, k-1))  &
                - (Kx(j+1, k)*p(j+1, k) + Kx(j, k)*p(j-1, k))

            pw = pw + w(j, k)*p(j, k)
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
      ENDDO
    ENDDO
!$OMP END DO

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  beta = rrn/rro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!$OMP SINGLE
  rro = rrn
!$OMP END SINGLE

ENDDO

!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_local_solve

SUBROUTINE tea_leaf_dpcg_add_t2_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           tile_t2, &
                           u)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u

  INTEGER(KIND=4) :: j,k

  REAL(kind=8) :: tile_t2

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      u(j, k) = u(j, k) + tile_t2
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_add_t2_kernel

SUBROUTINE tea_leaf_dpcg_solve_z_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           r,                      &
                           Mi,                     &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           cp,                     &
                           bfp,                    &
                           rx,                     &
                           ry,                     &
                           ztaz,                    &
                           ztr,                     &
                           preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, z, Mi
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) ::  rx, ry
  REAL(kind=8) :: ztr, ztaz

  ztr = 0.0_8
  ztaz = 0.0_8

!$OMP PARALLEL REDUCTION(+:ztr, ztaz)
  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

    IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
      CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, cp, bfp, Kx, Ky, rx, ry)
    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
      CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, Mi, Kx, Ky, rx, ry)
    ENDIF

  ELSE

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        z(j, k) = r(j, k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      ztaz = ztaz + z(j, k)*((1.0_8                            &
            + (ry*Ky(j, k+1) + ry*Ky(j, k))     &
            + (rx*Kx(j+1, k) + rx*Kx(j, k)))    &
            - (ry*Ky(j, k+1) + ry*Ky(j, k))     &
            - (rx*Kx(j+1, k) + rx*Kx(j, k)))
      ztr = ztr + r(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_solve_z_kernel

SUBROUTINE tea_leaf_cg_init_zp_kernel(x_min,             &
                                       x_max,             &
                                       y_min,             &
                                       y_max,             &
                                       halo_exchange_depth,             &
                                       p,                 &
                                       z,                 &
                                       tile_t2)


  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: z, p

  INTEGER(KIND=4) :: j,k
  REAL(kind=8) :: tile_t2

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      z(j, k) = z(j, k) - tile_t2
      p(j, k) = z(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_cg_init_zp_kernel

END MODULE

