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

SUBROUTINE tea_leaf_dpcg_restrict_ZT_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    r , &
    ZTr )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: ZTr
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r

  ZTr = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION (+:ZTr)
  DO k=y_min,y_max
    DO j=x_min,x_max
      ZTr = ZTr + r(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_restrict_ZT_kernel

SUBROUTINE tea_leaf_dpcg_sum_matrix_rows_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    kx, ky , &
    rx, ry, E_local )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: rx, ry, E_local
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: kx, ky

  E_local = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION (+:E_local)
  DO k=y_min,y_max
    DO j=x_min,x_max
      E_local = E_local + (1.0_8                &
            + (ry*Ky(j, k+1) + ry*Ky(j, k))     &
            + (rx*Kx(j+1, k) + rx*Kx(j, k)))    &
            - (ry*Ky(j, k+1) + ry*Ky(j, k))     &
            - (rx*Kx(j+1, k) + rx*Kx(j, k))
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_sum_matrix_rows_kernel

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
                           z)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: r, z, Mi, p, w, u, u0, def_e

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: it_count

  REAL(kind=8) :: rro, smvp
  REAL(KIND=8) ::  alpha, beta, pw, rrn

  rro = 0.0_8
  pw = 0.0_8

  rrn = 1e10

  it_count = 0

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

!$OMP DO REDUCTION(+:rro)
  DO k=y_min, y_max
    DO j=x_min, x_max
      smvp = (1.0_8                                         &
          + (def_e(j, k+1) + def_e(j, k))                      &
          + (def_e(j+1, k) + def_e(j, k)))*u(j, k)             &
          + (def_e(j, k+1)*u(j, k+1) + def_e(j, k)*u(j, k-1))  &
          + (def_e(j+1, k)*u(j+1, k) + def_e(j, k)*u(j-1, k))

      r(j, k) = u0(j, k) - smvp
      p(j, k) = r(j, k)

      rro = rro + r(j, k)*p(j, k);
    ENDDO
  ENDDO
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO WHILE ((abs(sqrt(rrn)) .gt. 1e-10*rro) .and. (it_count < 20))

!$OMP SINGLE
    pw = 0.0_8
!$OMP END SINGLE

!$OMP DO REDUCTION(+:pw)
    DO k=y_min,y_max
      DO j=x_min,x_max
        smvp = (1.0_8                                         &
            + (def_e(j, k+1) + def_e(j, k))                      &
            + (def_e(j+1, k) + def_e(j, k)))*p(j, k)             &
            + (def_e(j, k+1)*p(j, k+1) + def_e(j, k)*p(j, k-1))  &
            + (def_e(j+1, k)*p(j+1, k) + def_e(j, k)*p(j-1, k))

        w(j, k) = smvp
        pw = pw + smvp*p(j, k)
      ENDDO
    ENDDO
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    alpha = rro/pw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP SINGLE
    rrn = 0.0_8
!$OMP END SINGLE

!$OMP DO REDUCTION(+:rrn)
    DO k=y_min,y_max
      DO j=x_min,x_max
        u(j, k) = u(j, k) + alpha*p(j, k)
        r(j, k) = r(j, k) - alpha*w(j, k)
        rrn = rrn + r(j, k)*r(j, k)
      ENDDO
    ENDDO
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    beta = rrn/rro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        p(j, k) = r(j, k) + beta*p(j, k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

!$OMP SINGLE
    rro = rrn
    it_count = it_count + 1
!$OMP END SINGLE

  ENDDO

!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_local_solve

SUBROUTINE tea_leaf_dpcg_init_p_kernel(x_min,             &
                                       x_max,             &
                                       y_min,             &
                                       y_max,             &
                                       halo_exchange_depth,             &
                                       p,                 &
                                       z)


  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: z, p

  INTEGER(KIND=4) :: j,k

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k) = z(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_init_p_kernel

SUBROUTINE tea_leaf_dpcg_matmul_ZTA_kernel(x_min,  &
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
                           preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, z, Mi
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) ::  rx, ry
  REAL(kind=8) :: ztaz

  ztaz = 0.0_8

!$OMP PARALLEL REDUCTION(+:ztaz)
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
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_matmul_ZTA_kernel

SUBROUTINE tea_leaf_dpcg_store_r_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    r, r_m1 )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, r_m1

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      r_m1(j, k) = r(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_store_r_kernel

SUBROUTINE tea_leaf_dpcg_calc_rrn_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    r, r_m1 , z, &
    rrn )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: rrn
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, r_m1, z

  rrn = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION (+:rrn)
  DO k=y_min,y_max
    DO j=x_min,x_max
      rrn = rrn + (r(j, k) - r_m1(j, k))*z(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_calc_rrn_kernel

SUBROUTINE tea_leaf_dpcg_calc_p_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    p, r , &
    beta )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: beta
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: p, r

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k) = r(j, k) + beta*p(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_calc_p_kernel

SUBROUTINE tea_leaf_dpcg_prolong_Z_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    z , &
    tile_t2 )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: tile_t2
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: z

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      z(j, k) = z(j, k) - tile_t2
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_prolong_Z_kernel

SUBROUTINE tea_leaf_dpcg_calc_zrnorm_kernel(x_min, &
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
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, z
  REAL(KIND=8) :: norm
  integer :: j, k

  norm = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + z(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

end SUBROUTINE tea_leaf_dpcg_calc_zrnorm_kernel

END MODULE

