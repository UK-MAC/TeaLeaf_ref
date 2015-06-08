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
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Implicitly calculates the change in temperature using a Jacobi iteration

MODULE tea_leaf_kernel_jacobi_module

CONTAINS

SUBROUTINE tea_leaf_kernel_jacobi_solve(x_min,       &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           halo_exchange_depth,             &
                           rx,                &
                           ry,                &
                           Kx,                &
                           Ky,                &
                           error,             &
                           u0,                &
                           u1,                &
                           un)


  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: un , u1, u0 , Kx , Ky

  REAL(KIND=8) :: ry,rx, error

  INTEGER(KIND=4) :: j,k

  error = 0.0_8

!$OMP PARALLEL
!$OMP DO
    DO k=y_min, y_max
      DO j=x_min, x_max
        un(j,k) = u1(j,k)
      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO REDUCTION(+:error)
    DO k=y_min, y_max
      DO j=x_min, x_max
        u1(j,k) = (u0(j,k) + rx*(Kx(j+1,k  )*un(j+1,k  ) + Kx(j  ,k  )*un(j-1,k  ))  &
                           + ry*(Ky(j  ,k+1)*un(j  ,k+1) + Ky(j  ,k  )*un(j  ,k-1))) &
                             /(1.0_8 &
                                + rx*(Kx(j,k)+Kx(j+1,k)) &
                                + ry*(Ky(j,k)+Ky(j,k+1)))

        error = error +  ABS(u1(j,k)-un(j,k))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_jacobi_solve

END MODULE tea_leaf_kernel_jacobi_module

