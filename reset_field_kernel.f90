MODULE reset_field_kernel_module

CONTAINS

SUBROUTINE reset_field_kernel(x_min,x_max,y_min,y_max,    &
                              density0,           &
                              density1,           &
                              energy0,            &
                              energy1,            &
                              xvel0,              &
                              xvel1,              &
                              yvel0,              &
                              yvel1)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1,energy1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1,yvel1

  INTEGER :: j,k

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
     DO j=x_min,x_max
        density0(j,k)=density1(j,k)
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max
     DO j=x_min,x_max
        energy0(j,k)=energy1(j,k)
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max+1
     DO j=x_min,x_max+1
        xvel0(j,k)=xvel1(j,k)
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max+1
     DO j=x_min,x_max+1
        yvel0(j,k)=yvel1(j,k)
     ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE reset_field_kernel

END MODULE reset_field_kernel_module
