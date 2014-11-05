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

MODULE tea_leaf_kernel_module

CONTAINS

SUBROUTINE tea_leaf_kernel_init(x_min,        &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           celldx,            &
                           celldy,            &
                           volume,            &
                           density,           &
                           energy,            &
                           u0,                &
                           u1,                &
                           un,                &
                           Kx_tmp,            &
                           Ky_tmp,            &
                           Kx,                &
                           Ky,                &
                           coef)

  IMPLICIT NONE

   INTEGER         ::            CONDUCTIVITY        = 1 &
                                ,RECIP_CONDUCTIVITY  = 2

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: un
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Kx_tmp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Ky_tmp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Ky

  INTEGER(KIND=4) :: coef

  INTEGER(KIND=4) :: j,k,n


! CALC DIFFUSION COEFFICIENT
!$OMP PARALLEL
  IF(coef .EQ. RECIP_CONDUCTIVITY) THEN
!$OMP DO 
    DO k=y_min-1,y_max+2
      DO j=x_min-1,x_max+2
         Kx_tmp(j  ,k  )=1.0_8/density(j  ,k  )
         Ky_tmp(j  ,k  )=1.0_8/density(j  ,k  )
      ENDDO
    ENDDO
!$OMP END DO
  ELSE IF(coef .EQ. CONDUCTIVITY) THEN
!$OMP DO
    DO k=y_min-1,y_max+2
      DO j=x_min-1,x_max+2
         Kx_tmp(j  ,k  )=density(j  ,k  )
         Ky_tmp(j  ,k  )=density(j  ,k  )
      ENDDO
    ENDDO
!$OMP END DO
  ENDIF

!$OMP DO
  DO k=y_min,y_max+1
    DO j=x_min,x_max+1
         Kx(j,k)=(Kx_tmp(j-1,k  )+Kx_tmp(j,k))/(2.0_8*Kx_tmp(j-1,k  )*Kx_tmp(j,k))
         Ky(j,k)=(Ky_tmp(j  ,k-1)+Ky_tmp(j,k))/(2.0_8*Ky_tmp(j,  k-1)*Ky_tmp(j,k))
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO 
  DO k=y_min-1, y_max+1
    DO j=x_min-1, x_max+1
      u0(j,k) =  energy(j,k) * density(j,k)
    ENDDO
  ENDDO
!$OMP END DO

  ! INITIAL GUESS
!$OMP DO
  DO k=y_min-1, y_max+1
    DO j=x_min-1, x_max+1
      u1(j,k) = u0(j,k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE tea_leaf_kernel_init

SUBROUTINE tea_leaf_kernel_solve(x_min,       &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           rx,                &
                           ry,                &
                           Kx,                &
                           Ky,                &
                           error,             &
                           u0,                &
                           u1,                &
                           un)


  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: un
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u1, u0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Ky

  REAL(KIND=8) :: ry,rx, error

  INTEGER(KIND=4) :: j,k

  error = 0.0_8

!$OMP PARALLEL
!$OMP DO
    DO k=y_min-1, y_max+1
      DO j=x_min-1, x_max+1
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

END SUBROUTINE tea_leaf_kernel_solve

! Finalise routine is used by both implementations
SUBROUTINE tea_leaf_kernel_finalise(x_min,    &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           energy,            &
                           density,           &
                           u)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density

  INTEGER(KIND=4) :: j,k

!$OMP PARALLEL DO 
  DO k=y_min, y_max
    DO j=x_min, x_max
      energy(j,k) = u(j,k) / density(j,k)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE tea_leaf_kernel_finalise

SUBROUTINE tea_leaf_calc_residual(x_min,       &
                                  x_max,       &
                                  y_min,       &
                                  y_max,       &
                                  u ,          &
                                  u0,          &
                                  r,           &
                                  Kx,          &
                                  Ky,          &
                                  rx, ry       )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u0, u, r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Ky

  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: j,k

!$OMP PARALLEL
!$OMP DO private(smvp)
    DO k=y_min, y_max
      DO j=x_min, x_max
        smvp = (1.0_8                                         &
            + ry*(Ky(j, k+1) + Ky(j, k))                      &
            + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
            - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
            - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))
        r(j, k) = u0(j, k) - smvp
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_calc_residual

END MODULE tea_leaf_kernel_module

