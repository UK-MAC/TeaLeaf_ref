MODULE flux_calc_kernel_module

CONTAINS

SUBROUTINE flux_calc_kernel(x_min,x_max,y_min,y_max,dt,              &
                            xarea,                           &
                            yarea,                           &
                            xvel0,                           &
                            yvel0,                           &
                            xvel1,                           &
                            yvel1,                           &
                            vol_flux_x,                      &
                            vol_flux_y                       )

  IMPLICIT NONE

  INTEGER       :: x_min, x_max, y_min, y_max
  REAL(KIND=8) :: dt
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: xarea
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: yarea
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1,yvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y

  INTEGER :: j,k

!$OMP PARALLEL

!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max+1 
      vol_flux_x(j,k)=0.25*dt*xarea(j,k)                  &
                     *(xvel0(j,k)+xvel0(j,k+1)+xvel1(j,k)+xvel1(j,k+1))
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max+1
    DO j=x_min,x_max
      vol_flux_y(j,k)=0.25*dt*yarea(j,k)                  &
                     *(yvel0(j,k)+yvel0(j+1,k)+yvel1(j,k)+yvel1(j+1,k))
    ENDDO
  ENDDO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE flux_calc_kernel

END MODULE flux_calc_kernel_module
