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
!>  @details Implicitly calculates the change in temperature using deflated CG method

MODULE tea_leaf_dpcg_kernel_module

  USE tea_leaf_common_kernel_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_dpcg_coarsen_matrix_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    kx, ky, &
    kx_local, ky_local, &
    rx, ry)

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: rx, ry, kx_local, ky_local
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: kx, ky

  kx_local = 0.0_8
  ky_local = 0.0_8

!$OMP PARALLEL REDUCTION(+:kx_local, ky_local)
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      kx_local = kx_local + rx*Kx(j, k)
      ky_local = ky_local + ry*Ky(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_coarsen_matrix_kernel

SUBROUTINE tea_leaf_dpcg_add_z_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    u, &
    tile_t2 )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: tile_t2
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      u(j, k) = u(j, k) + tile_t2
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_add_z_kernel

SUBROUTINE tea_leaf_dpcg_restrict_ZT_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    r , &
    ZTr )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: ZTr
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r

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

SUBROUTINE tea_leaf_dpcg_solve_z_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           r,                      &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           Mi,                     &
                           cp,                     &
                           bfp,                    &
                           rx,                     &
                           ry,                     &
                           preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, z, Mi
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) :: rx, ry

!$OMP PARALLEL
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
!$OMP END DO
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_solve_z_kernel

SUBROUTINE tea_leaf_dpcg_matmul_ZTA_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           rx,                     &
                           ry,                     &
                           ztaz,                    &
                           preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, z, Mi
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) :: rx, ry
  REAL(kind=8) :: ztaz

  ztaz = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION(+:ztaz)
  DO k=y_min,y_max
    DO j=x_min,x_max
      ztaz = ztaz + (1.0_8                                  &
          + ry*(Ky(j, k+1) + Ky(j, k))                      &
          + rx*(Kx(j+1, k) + Kx(j, k)))*z(j, k)             &
          - ry*(Ky(j, k+1)*z(j, k+1) + Ky(j, k)*z(j, k-1))  &
          - rx*(Kx(j+1, k)*z(j+1, k) + Kx(j, k)*z(j-1, k))
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_matmul_ZTA_kernel

SUBROUTINE tea_leaf_dpcg_matmul_ZTA_kernel1(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           rx,                     &
                           ry,                     &
                           ztaz,                    &
                           preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, z, Mi
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) :: rx, ry
  REAL(kind=8) :: ztaz
  REAL(kind=8), DIMENSION(x_min:x_max) :: ztazj ! store the values on the row to enable vectorization of the loop

  ztaz = 0.0_8

!$OMP PARALLEL PRIVATE(ztazj)
!$OMP DO REDUCTION(+:ztaz)
  DO k=y_min,y_max
    DO j=x_min,x_max
      ztazj(j) = (1.0_8                                     &
          + ry*(Ky(j, k+1) + Ky(j, k))                      &
          + rx*(Kx(j+1, k) + Kx(j, k)))*z(j, k)             &
          - ry*(Ky(j, k+1)*z(j, k+1) + Ky(j, k)*z(j, k-1))  &
          - rx*(Kx(j+1, k)*z(j+1, k) + Kx(j, k)*z(j-1, k))
    ENDDO
    DO j=x_min,x_max
      ztaz = ztaz + ztazj(j)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_matmul_ZTA_kernel1

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
                          :: r, z, Mi, p, w, u, u0, def_kx, def_ky, def_di, sd

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

!$OMP BARRIER
!$OMP MASTER
    rro = initial_residual
    initial_residual = sqrt(abs(initial_residual))
!$OMP END MASTER
!$OMP BARRIER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO WHILE ((sqrt(abs(rrn)) .gt. eps*initial_residual) .and. (it_count < inner_iters))

!$OMP BARRIER
!$OMP MASTER
    pw = 0.0_8
!$OMP END MASTER
!$OMP BARRIER

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

!$OMP BARRIER
!$OMP MASTER
    rrn = 0.0_8
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        u(j, k) = u(j, k) + alpha*p(j, k)
        r(j, k) = r(j, k) - alpha*w(j, k)
        z(j, k) = r(j, k)*Mi(j, k)
      ENDDO
    ENDDO
!$OMP END DO

!    IF (use_ppcg) THEN
!      DO inner_step=1,10
!!$OMP DO
!        DO k=y_min,y_max
!          DO j=x_min,x_max
!            sd(j, k) = z(j, k)/theta
!          ENDDO
!        ENDDO
!!$OMP END DO
!!$OMP DO
!        DO k=y_min,y_max
!          DO j=x_min,x_max
!            smvp = def_di(j, k)*sd(j, k)             &
!              - (def_ky(j, k+1)*sd(j, k+1) + def_ky(j, k)*sd(j, k-1))  &
!              - (def_kx(j+1, k)*sd(j+1, k) + def_kx(j, k)*sd(j-1, k))
!
!            r(j, k) = r(j, k) - smvp
!            z(j, k) = r(j, k)*Mi(j, k)
!
!            u(j, k) = u(j, k) + sd(j, k)
!          ENDDO
!        ENDDO
!!$OMP END DO
!!$OMP DO
!        DO k=y_min,y_max
!          DO j=x_min,x_max
!            sd(j, k) = inner_ch_alphas(inner_step)*sd(j, k) + inner_ch_betas(inner_step)*z(j, k)
!          ENDDO
!        ENDDO
!!$OMP END DO
!      ENDDO
!    ENDIF

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

!$OMP BARRIER
!$OMP MASTER
    rro = rrn
    it_count = it_count + 1
    inner_cg_alphas(it_count) = alpha
    inner_cg_betas(it_count) = beta
!$OMP END MASTER
!$OMP BARRIER

  ENDDO

!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_local_solve

SUBROUTINE tea_leaf_dpcg_init_p_kernel(x_min,               &
                                       x_max,               &
                                       y_min,               &
                                       y_max,               &
                                       halo_exchange_depth, &
                                       p,                   &
                                       z)


  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
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

SUBROUTINE tea_leaf_dpcg_store_r_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    r, r_m1 )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, r_m1

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
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, r_m1, z

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
    p, r, z, Kx, Ky, cp, bfp, &
    rx, ry, beta, preconditioner_type )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  INTEGER :: preconditioner_type
  REAL(KIND=8) :: rx, ry, beta
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: p, r, Kx, Ky, Mi, z
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

!$OMP PARALLEL
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        p(j, k) = z(j, k) + beta*p(j, k)
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
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: z

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
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
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

END SUBROUTINE tea_leaf_dpcg_calc_zrnorm_kernel

END MODULE

