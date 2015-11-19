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
    nx, dx, ny, dy, kx_local, ky_local, &
    rx, ry)

  IMPLICIT NONE
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  INTEGER(KIND=4) :: nx, dx, ny, dy
  REAL(KIND=8) :: rx, ry
  REAL(KIND=8), DIMENSION(nx,ny) :: kx_local, ky_local
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: kx, ky

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: jj,j_start,j_end,kk,k_start,k_end

  kx_local = 0.0_8
  ky_local = 0.0_8

  DO kk=1,ny
    k_start=y_min+(kk-1)*dy
    k_end  =min(k_start+dy-1,y_max)
!$OMP PARALLEL REDUCTION(+:kx_local, ky_local) PRIVATE(j,jj,j_start,j_end)
!$OMP DO
    DO k=k_start,k_end
      DO jj=1,nx
        j_start=x_min+(jj-1)*dx
        j_end  =min(j_start+dx-1,x_max)
        DO j=j_start,j_end
          kx_local(jj,kk) = kx_local(jj,kk) + rx*Kx(j, k)
          ky_local(jj,kk) = ky_local(jj,kk) + ry*Ky(j, k)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDDO

END SUBROUTINE tea_leaf_dpcg_coarsen_matrix_kernel

SUBROUTINE tea_leaf_dpcg_prolong_Z_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    nx, dx, ny, dy, &
    z , &
    tile_t2 )

  IMPLICIT NONE
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  INTEGER(KIND=4) :: nx, dx, ny, dy
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: z
  REAL(KIND=8), DIMENSION(nx, ny) :: tile_t2

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: jj,j_start,j_end,kk,k_start,k_end

  DO kk=1,ny
    k_start=y_min+(kk-1)*dy
    k_end  =min(k_start+dy-1,y_max)
!$OMP PARALLEL PRIVATE(j,jj,j_start,j_end)
!$OMP DO
    DO k=k_start,k_end
      DO jj=1,nx
        j_start=x_min+(jj-1)*dx
        j_end  =min(j_start+dx-1,x_max)
        DO j=j_start,j_end
          z(j, k) = z(j, k) - tile_t2(jj,kk)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDDO

END SUBROUTINE tea_leaf_dpcg_prolong_Z_kernel

SUBROUTINE tea_leaf_dpcg_subtract_z_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    nx, dx, ny, dy, &
    u, &
    tile_t2 )

  IMPLICIT NONE
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  INTEGER(KIND=4) :: nx, dx, ny, dy
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u
  REAL(KIND=8), DIMENSION(nx, ny) :: tile_t2

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: jj,j_start,j_end,kk,k_start,k_end

  DO kk=1,ny
    k_start=y_min+(kk-1)*dy
    k_end  =min(k_start+dy-1,y_max)
!$OMP PARALLEL PRIVATE(j,jj,j_start,j_end)
!$OMP DO
    DO k=k_start,k_end
      DO jj=1,nx
        j_start=x_min+(jj-1)*dx
        j_end  =min(j_start+dx-1,x_max)
        DO j=j_start,j_end
          u(j, k) = u(j, k) - tile_t2(jj,kk)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDDO

END SUBROUTINE tea_leaf_dpcg_subtract_z_kernel

SUBROUTINE tea_leaf_dpcg_restrict_ZT_kernel(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    nx, dx, ny, dy, &
    r , &
    ZTr )

  IMPLICIT NONE
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  INTEGER(KIND=4) :: nx, dx, ny, dy
  REAL(KIND=8), DIMENSION(nx, ny) :: ZTr
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: jj,j_start,j_end,kk,k_start,k_end

  ZTr = 0.0_8

  DO kk=1,ny
    k_start=y_min+(kk-1)*dy
    k_end  =min(k_start+dy-1,y_max)
!$OMP PARALLEL
!$OMP DO REDUCTION (+:ZTr) PRIVATE(j,jj,j_start,j_end)
    DO k=k_start,k_end
      DO jj=1,nx
        j_start=x_min+(jj-1)*dx
        j_end  =min(j_start+dx-1,x_max)
        DO j=j_start,j_end
          ZTr(jj,kk) = ZTr(jj,kk) + r(j, k)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDDO

END SUBROUTINE tea_leaf_dpcg_restrict_ZT_kernel

SUBROUTINE tea_leaf_dpcg_matmul_ZTA_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,    &
                           nx, dx, ny, dy,         &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           rx,                     &
                           ry,                     &
                           ztaz)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth

  INTEGER(KIND=4) :: nx, dx, ny, dy
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: Kx, Ky, z
  REAL(KIND=8) :: rx, ry
  REAL(KIND=8), DIMENSION(nx, ny) :: ztaz

  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: jj,j_start,j_end,kk,k_start,k_end

  ztaz = 0.0_8

  DO kk=1,ny
    k_start=y_min+(kk-1)*dy
    k_end  =min(k_start+dy-1,y_max)
!$OMP PARALLEL
!$OMP DO REDUCTION (+:ztaz) PRIVATE(j,jj,j_start,j_end)
    DO k=k_start,k_end
      DO jj=1,nx
        j_start=x_min+(jj-1)*dx
        j_end  =min(j_start+dx-1,x_max)
        DO j=j_start,j_end
          ztaz(jj,kk) = ztaz(jj,kk) + (1.0_8                &
          + ry*(Ky(j, k+1) + Ky(j, k))                      &
          + rx*(Kx(j+1, k) + Kx(j, k)))*z(j, k)             &
          - ry*(Ky(j, k+1)*z(j, k+1) + Ky(j, k)*z(j, k-1))  &
          - rx*(Kx(j+1, k)*z(j+1, k) + Kx(j, k)*z(j-1, k))
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDDO

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
                           ztaz)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: Kx, Ky, z

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
                             r, z, Mi)
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
    p, z, &
    beta )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: beta
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: p, z

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

SUBROUTINE tea_leaf_dpcg_calc_zrnorm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth,             &
                          z, r,               &
                          norm)

  IMPLICIT NONE

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

