MODULE set_field_kernel_module

CONTAINS

SUBROUTINE set_field_kernel(x_min,x_max,y_min,y_max,    &
                              density0,           &
                              density1,           &
                              energy0,            &
                              energy1)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1,energy1

  INTEGER :: j,k

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
     DO j=x_min,x_max
        density1(j,k)=density0(j,k)
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max
     DO j=x_min,x_max
        energy1(j,k)=energy0(j,k)
     ENDDO
  ENDDO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE set_field_kernel

END MODULE set_field_kernel_module
