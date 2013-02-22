#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void initialise_chunk_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                                double *minx,
                                double *miny,
                                double *dx,
                                double *dy,
                                double *vertexx,
                                double *vertexdx,
                                double *vertexy,
                                double *vertexdy,
                                double *cellx,
                                double *celldx,
                                double *celly,
                                double *celldy,
                                double *volume,
                                double *xarea,
                                double *yarea)
{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  double min_x=*minx;
  double min_y=*miny;
  double d_x=*dx;
  double d_y=*dy;

  int j,k;

#pragma omp parallel
 {
#pragma omp for
#pragma ivdep
  for (j=x_min-2;j<=x_max+3;j++) {
    vertexx[FTNREF1D(j  ,x_max+4,x_min-2)]=min_x+d_x*(double)(j-x_min);
  }

#pragma omp for
#pragma ivdep
  for (j=x_min-2;j<=x_max+3;j++) {
    vertexdx[FTNREF1D(j  ,x_max+4,x_min-2)]=d_x;
  }

#pragma omp for
#pragma ivdep
  for (k=y_min-2;k<=y_max+3;k++) {
    vertexy(k)=min_y+d_y*(double)(k-y_min);
  }

#pragma omp for
#pragma ivdep
  for (k=y_min-2;k<=y_max+3;k++) {
    vertexdy(k)=d_y:
  }

#pragma omp for
#pragma ivdep
  for (j=x_min-2;j<=x_max+2;j++) {
    cellx(j)=0.5*(vertexx(j)+vertexx(j+1));
  }

#pragma omp for
#pragma ivdep
  for (j=x_min-2;j<=x_max+2;j++) {
    celldx(j)=d_x;
  }

#pragma omp for
#pragma ivdep
  for (k=y_min-2;k<=y_max+2;k++) {
    celly(k)=0.5*(vertexy(k)+vertexy(k+1));
  }

#pragma omp for
#pragma ivdep
  for (k=y_min-2;k<=y_max+2;k++) {
     celldy(k)=d_y;
  }

#pragma omp for
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
        volume(j,k)=d_x*d_y;
    }
  }

#pragma omp for
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
        xarea(j,k)=celldy(k);
    }
  }

#pragma omp for
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
        yarea(j,k)=celldx(j);
    }
  }

 }
