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

IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_kernel_init_cg_fortran(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           density,                &
                           energy,                 &
                           u,                      &
                           p,                      & ! 1
                           r,                      & ! 2
                           Mi,                     & ! 3
                           w,                      & ! 4
                           z,                      & ! 5
                           Kx,                     & ! 6
                           Ky,                     & ! 7
                           rx,                     &
                           ry,                     &
                           rro,                    &
                           coef)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: p
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Mi
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: w
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: z
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky

  INTEGER(KIND=4) :: coef
  INTEGER(KIND=4) :: j,k,n

  REAL(kind=8) :: rro
  REAL(KIND=8) ::  rx, ry

   INTEGER         ::            CONDUCTIVITY        = 1 &
                                ,RECIP_CONDUCTIVITY  = 2

  rro = 0.0_8
  p = 0.0_8
  r = 0.0_8

!$OMP PARALLEL
!$OMP DO 
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+2
      u(j,k) = energy(j,k)*density(j,k)
    ENDDO
  ENDDO
!$OMP END DO

  IF(coef .EQ. RECIP_CONDUCTIVITY) THEN
!$OMP DO 
    ! use w as temp val
    DO k=y_min-1,y_max+1
      DO j=x_min-1,x_max+1
         w(j  ,k  )=1.0_8/density(j  ,k  )
      ENDDO
    ENDDO
!$OMP END DO
  ELSE IF(coef .EQ. CONDUCTIVITY) THEN
!$OMP DO
    DO k=y_min-1,y_max+1
      DO j=x_min-1,x_max+1
         w(j  ,k  )=density(j  ,k  )
      ENDDO
    ENDDO
!$OMP END DO
  ENDIF

!$OMP DO
   DO k=y_min,y_max+1
     DO j=x_min,x_max+1
          Kx(j,k)=(w(j-1,k  ) + w(j,k))/(2.0_8*w(j-1,k  )*w(j,k))
          Ky(j,k)=(w(j  ,k-1) + w(j,k))/(2.0_8*w(j  ,k-1)*w(j,k))
     ENDDO
   ENDDO
!$OMP END DO

!$OMP DO REDUCTION(+:rro)
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
                - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
                - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))

            r(j, k) = u(j, k) - w(j, k)

#if defined(USE_PRECONDITIONER)
            ! inverse diagonal used as preconditioner
            Mi(j, k) = (1.0_8                                     &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))
            Mi(j, k) = 1.0_8/Mi(j, k)

            z(j, k) = Mi(j, k)*r(j, k)
            p(j, k) = z(j, k)
#else
            p(j, k) = r(j, k)
#endif

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
                                                   p,                 &
                                                   w,                 &
                                                   Kx,                &
                                                   Ky,                &
                                                   rx,                &
                                                   ry,                &
                                                   pw                 )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: p
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: w
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky

    REAL(KIND=8) ::  rx, ry

    INTEGER(KIND=4) :: j,k,n
    REAL(kind=8) :: pw

    pw = 0.0_08

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
                                                    u,                 &
                                                    p,                 &
                                                    r,                 &
                                                    Mi,                &
                                                    w,                 &
                                                    z,                 &
                                                    alpha,             &
                                                    rrn                )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: p
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Mi
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: w
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: z

    INTEGER(KIND=4) :: j,k,n
    REAL(kind=8) :: alpha, rrn

    rrn = 0.0_08

!$OMP PARALLEL
!$OMP DO REDUCTION(+:rrn)
    DO k=y_min,y_max
        DO j=x_min,x_max
            u(j, k) = u(j, k) + alpha*p(j, k)
            r(j, k) = r(j, k) - alpha*w(j, k)
#if defined(USE_PRECONDITIONER)
            z(j, k) = Mi(j, k)*r(j, k)
            rrn = rrn + r(j, k)*z(j, k)
#else
            rrn = rrn + r(j, k)*r(j, k)
#endif
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_ur

SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_p(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   p,                 &
                                                   r,                 &
                                                   z,                 &
                                                   beta               )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: p
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: z

    REAL(kind=8) :: error

    INTEGER(KIND=4) :: j,k,n
    REAL(kind=8) :: beta

!$OMP PARALLEL
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
#if defined(USE_PRECONDITIONER)
            p(j, k) = z(j, k) + beta*p(j, k)
#else
            p(j, k) = r(j, k) + beta*p(j, k)
#endif
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_p

END MODULE
