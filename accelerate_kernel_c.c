#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void accelerate_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                          double *dbyt,
                          double *xarea,
                          double *yarea,
                          double *volume,
                          double *density0,
                          double *pressure,
                          double *viscosity,
                          double *xvel0,
                          double *yvel0,
                          double *xvel1,
                          double *yvel1,
                          double *stepbymass)
{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  double dt=*dbyt;

  int j,k,err;
  double nodal_mass;

#pragma omp parallel
 {

#pragma omp for private(nodal_mass)
  for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max+1;j++) {
      nodal_mass=(density0[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]*volume[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]
                 +density0[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]*volume[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]
                 +density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]*volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                 +density0[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]*volume[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)])
                 *0.25;
      stepbymass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.5*dt/nodal_mass;
    }
  }

#pragma omp for
  for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max+1;j++) {
      xvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           -stepbymass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           *(xarea[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           *(pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]-pressure[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)])
                                           +xarea[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                           *(pressure[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]-pressure[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]));
    }
  }

#pragma omp for
  for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max+1;j++) {
      yvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           -stepbymass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           *(yarea[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                           *(pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]-pressure[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)])
                                           +yarea[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                           *(pressure[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]-pressure[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]));
    }
  }

#pragma omp for
  for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max+1;j++) {
      xvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=xvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           -stepbymass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           *(xarea[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           *(viscosity[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]-viscosity[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)])
                                           +xarea[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                           *(viscosity[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]-viscosity[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]));
    }
  }

#pragma omp for
  for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max+1;j++) {
      yvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=yvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           -stepbymass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                           *(yarea[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                           *(viscosity[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]-viscosity[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)])
                                           +yarea[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                           *(viscosity[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]-viscosity[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]));
    }
  }

 }

}
