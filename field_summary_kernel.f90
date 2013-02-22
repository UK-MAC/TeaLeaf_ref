MODULE field_summary_kernel_module

CONTAINS

SUBROUTINE field_summary_kernel(x_min,x_max,y_min,y_max, &
                                volume,                  &
                                density0,                &
                                energy0,                 &
                                pressure,                &
                                xvel0,                   &
                                yvel0,                   &
                                vol,mass,ie,ke,press     )

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
  REAL(KIND=8) :: vol,mass,ie,ke,press

  INTEGER      :: j,k,jv,kv
  REAL(KIND=8) :: vsqrd,cell_vol,cell_mass

  vol=0.0
  mass=0.0
  ie=0.0
  ke=0.0
  press=0.0

!$OMP PARALLEL
!$OMP DO PRIVATE(vsqrd,cell_vol,cell_mass) REDUCTION(+ : vol,mass,press,ie,ke)
  DO k=y_min,y_max
    DO j=x_min,x_max
      vsqrd=0.0
      DO kv=k,k+1
        DO jv=j,j+1
          vsqrd=vsqrd+0.25*(xvel0(jv,kv)**2+yvel0(jv,kv)**2)
        ENDDO
      ENDDO
      cell_vol=volume(j,k)
      cell_mass=cell_vol*density0(j,k)
      vol=vol+cell_vol
      mass=mass+cell_mass
      ie=ie+cell_mass*energy0(j,k)
      ke=ke+cell_mass*0.5*vsqrd
      press=press+cell_vol*pressure(j,k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE field_summary_kernel

END MODULE field_summary_kernel_module
