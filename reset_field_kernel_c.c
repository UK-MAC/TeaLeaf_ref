#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void reset_field_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                        double *density0,
                        double *density1,
                        double *energy0,
                        double *energy1,
                        double *xvel0,
                        double *xvel1,
                        double *yvel0,
                        double *yvel1)
{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;

  int j,k;
  
#pragma omp parallel
 {
#pragma omp for
  for (k=y_min;k<=y_max;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max;j++) {
      density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
    }
  }
#pragma omp for
  for (k=y_min;k<=y_max;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max;j++) {
      energy0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
    }
  }
#pragma omp for
  for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max+1;j++) {
      xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=xvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
    }
  }
#pragma omp for
  for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max+1;j++) {
      yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=yvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
    }
  }

 }

}
