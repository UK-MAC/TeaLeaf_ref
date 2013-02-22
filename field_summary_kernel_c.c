#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void field_summary_kernel(x_min,x_max,y_min,y_max,
                                volume,
                                density0,
                                energy0,
                                pressure,
                                xvel0,
                                yvel0,
                                vol,mass,ie,ke,press)
{
  int x_min,x_max,y_min,y_max;
  double volume;
  double density0,energy0;
  double pressure;
  double xvel0,yvel0;
  double vol,mass,ie,ke,press;

  int j,k,jv,kv;
  double vsqrd,cell_vol,cell_mass;

  vol=0.0;
  mass=0.0;
  ie=0.0;
  ke=0.0;;
  press=0.0;

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

}
