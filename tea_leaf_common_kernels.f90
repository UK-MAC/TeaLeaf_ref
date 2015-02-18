MODULE tea_leaf_kernel_common_module

CONTAINS

SUBROUTINE tea_leaf_kernel_init_common(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           density,                &
                           energy,                 &
                           u,                      &
                           u0,                     &
                           r,                      &
                           w,                      &
                           Kx,                     &
                           Ky,                     &
                           rx,                     &
                           ry,                     &
                           coef)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density, energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u, u0, r, w, Kx, Ky

  INTEGER(KIND=4) :: coef
  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) ::  rx, ry

   INTEGER         ::            CONDUCTIVITY        = 1 &
                                ,RECIP_CONDUCTIVITY  = 2

!$OMP PARALLEL
!$OMP DO 
  DO k=y_min, y_max
    DO j=x_min, x_max
      u(j,k) = energy(j,k)*density(j,k)
      u0(j,k) = energy(j,k)*density(j,k)
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

!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
                - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
                - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))

            r(j, k) = u(j, k) - w(j, k)
            !r(j, k) = u(j, k)! This is required to make a zero initial guess to match petsc errant behaviour
                              ! Only works one timestep is run
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_init_common

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

!$OMP PARALLEL
!$OMP DO
  DO k=y_min, y_max
    DO j=x_min, x_max
      energy(j,k) = u(j,k) / density(j,k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

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

SUBROUTINE tea_leaf_calc_2norm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          arr,               &
                          norm)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: arr
  REAL(KIND=8) :: norm
  integer :: j, k

  norm = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + arr(j, k)*arr(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

end SUBROUTINE tea_leaf_calc_2norm_kernel

END MODULE tea_leaf_kernel_common_module
