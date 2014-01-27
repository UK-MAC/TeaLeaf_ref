MODULE tea_leaf_kernel_module

CONTAINS

SUBROUTINE tea_leaf_kernel_init(x_min,             &
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
                           heat_capacity,     &
                           Kx_tmp,            &
                           Ky_tmp,            &
                           Kx,                &
                           Ky,                &
                           coef)

! clover_module used for coefficient constants
  USE clover_module
  USE report_module

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: u0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: heat_capacity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx_tmp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky_tmp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: un

  INTEGER(KIND=4) :: coef

  INTEGER(KIND=4) :: j,k,n


! CALC DIFFUSION COEFFICIENT
!$OMP PARALLEL
  IF(coef .EQ. RECIP_CONDUCTIVITY) THEN
!$OMP DO 
    DO k=y_min-1,y_max+1
      DO j=x_min-1,x_max+1
         Kx_tmp(j  ,k  )=1.0_8/density(j  ,k  )
         Kx_tmp(j+1,k  )=1.0_8/density(j+1,k  )
         Ky_tmp(j  ,k  )=1.0_8/density(j  ,k  )
         Ky_tmp(j  ,k+1)=1.0_8/density(j  ,k+1)
      ENDDO
    ENDDO
!$OMP END DO
  ELSE IF(coef .EQ. CONDUCTIVITY) THEN
!$OMP DO
    DO k=y_min-1,y_max+1
      DO j=x_min-1,x_max+1
         Kx_tmp(j  ,k  )=density(j  ,k  )
         Kx_tmp(j+1,k  )=density(j+1,k  )
         Ky_tmp(j  ,k  )=density(j  ,k  )
         Ky_tmp(j  ,k+1)=density(j  ,k+1)
      ENDDO
    ENDDO
!$OMP END DO
  ELSE
    CALL report_error('tea_leaf', 'unknown coefficient option')
  ENDIF

!$OMP DO
  DO k=y_min,y_max+1
    DO j=x_min,x_max+1
         Kx(j,k)=(Kx_tmp(j-1,k  )+Kx_tmp(j,k  ))/(2.0_8*Kx_tmp(j-1,k)*Kx_tmp(j,k))
         Ky(j,k)=(Ky_tmp(j,k-1)+Ky_tmp(j,k))/(2.0_8*Ky_tmp(j,k-1)*Ky_tmp(j,k))
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
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: u0, un
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky

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

!$OMP DO REDUCTION(MAX:error)
    DO k=y_min, y_max
      DO j=x_min, x_max
        u1(j,k) = (u0(j,k) + Kx(j+1,k)*rx*un(j+1,k) + Kx(j,k)*rx*un(j-1,k) &
                           + Ky(j,k+1)*ry*un(j,k+1) + Ky(j,k)*ry*un(j,k-1)) &
                             /(1.0_8+2.0_8*(0.5_8*(Kx(j,k)+Kx(j+1,k)))*rx &
                                +2.0_8*(0.5_8*(Ky(j,k)+Ky(j+1,k)))*ry)

        error = MAX(error, ABS(u1(j,k)-u0(j,k)))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve

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

END MODULE tea_leaf_kernel_module
