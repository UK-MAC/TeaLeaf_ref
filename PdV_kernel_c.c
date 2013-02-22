#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void pdv_kernel_c_(int *prdct,
                   int *xmin,int *xmax,int *ymin,int *ymax,
                double *dtbyt,
                double *xarea,
                double *yarea,
                double *volume,
                double *density0,
                double *density1,
                double *energy0,
                double *energy1,
                double *pressure,
                double *viscosity,
                double *xvel0,
                double *xvel1,
                double *yvel0,
                double *yvel1,
                double *volume_change)
{
  int predict=*prdct;
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  double dt=*dtbyt;

  int j,k;
  double recip_volume,energy_change,min_cell_volume,right_flux,left_flux,top_flux,bottom_flux,total_flux;
  
#pragma omp parallel
 {

  if(predict==0) {
    
#pragma omp for private(right_flux,left_flux,top_flux,bottom_flux,total_flux,min_cell_volume,energy_change,recip_volume)
    for (k=y_min;k<=y_max;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max;j++) {
        left_flux=  (xarea[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                   *(xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt*0.5;
        right_flux= (xarea[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])
                                   *(xvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt*0.5;
        bottom_flux=(yarea[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)])
                                   *(yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt*0.5;
        top_flux=   (yarea[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)])
                                   *(yvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt*0.5;

        total_flux=right_flux-left_flux+top_flux-bottom_flux;

        volume_change[FTNREF2D(j  ,k  ,x_max,x_min,y_min)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                         /(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+total_flux);

        min_cell_volume=MIN(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+right_flux-left_flux+top_flux-bottom_flux
                           ,MIN(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+right_flux-left_flux
                           ,volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+top_flux-bottom_flux));

        recip_volume=1.0/volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];

        energy_change=(pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]/density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                     +viscosity[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]/density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)])
                      *total_flux*recip_volume;

        energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]-energy_change;

        density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           *volume_change[FTNREF2D(j  ,k  ,x_max,x_min,y_min)];

      }
    }
  }
  else{
#pragma omp for private(right_flux,left_flux,top_flux,bottom_flux,total_flux,min_cell_volume,energy_change,recip_volume)
    for (k=y_min;k<=y_max;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max;j++) {
        left_flux=  (xarea[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                   *(xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +xvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel1[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt;
        right_flux= (xarea[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])
                                   *(xvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                                    +xvel1[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel1[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt;
        bottom_flux=(yarea[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)])
                                   *(yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel1[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt;
        top_flux=   (yarea[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)])
                                   *(yvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel1[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel1[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt;

        total_flux=right_flux-left_flux+top_flux-bottom_flux;

        volume_change[FTNREF2D(j  ,k  ,x_max,x_min,y_min)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                         /(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+total_flux);

        min_cell_volume=MIN(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+right_flux-left_flux+top_flux-bottom_flux
                           ,MIN(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+right_flux-left_flux
                           ,volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+top_flux-bottom_flux));

        recip_volume=1.0/volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];

        energy_change=(pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]/density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                     +viscosity[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]/density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)])
                      *total_flux*recip_volume;

        energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]-energy_change;

        density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           *volume_change[FTNREF2D(j  ,k  ,x_max,x_min,y_min)];

      }
    }

  }

 }

}

