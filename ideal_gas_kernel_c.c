#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void ideal_gas_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                      double *density,
                      double *energy,
                      double *pressure,
                      double *soundspeed)
{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;

  int j,k;

  double sound_speed_squared,v,pressurebyenergy,pressurebyvolume;
  
#pragma omp parallel
 {
#pragma omp for private(v,pressurebyenergy,pressurebyvolume,sound_speed_squared)
  for (k=y_min;k<=y_max;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max;j++) {								 
      v=1.0/density[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
      pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=(1.4-1.0)*density[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                   *energy[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
      pressurebyenergy=(1.4-1.0)*density[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
      pressurebyvolume=-density[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]*pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
      sound_speed_squared=v*v*(pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]*pressurebyenergy-pressurebyvolume);
      soundspeed[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=sqrt(sound_speed_squared);
    }
  }

 }

}
