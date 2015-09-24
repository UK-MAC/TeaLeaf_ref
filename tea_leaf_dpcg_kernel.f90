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
                           ztr, &
                           e, &
                           kx, &
                           ky, &
                           r)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, kx, ky

  INTEGER(KIND=4) :: j,k

  REAL(kind=8) :: ztr, e

  ztr = 0.0_8
  e = 0.0_8

!$OMP PARALLEL REDUCTION(+:ztr, e)
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      ztr = ztr + r(j, k)
      e = e + (1.0_8                    &
            + (Ky(j, k+1) + Ky(j, k))   &
            + (Kx(j+1, k) + Kx(j, k)))  &
            + (Ky(j, k+1) + Ky(j, k))   &
            + (Kx(j+1, k) + Kx(j, k))
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_dpcg_sum_r_kernel

END MODULE

