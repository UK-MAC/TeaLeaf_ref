#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void viscosity_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                      double *celldx,
                      double *celldy,
                      double *density0,
                      double *pressure,
                      double *viscosity,
                      double *xvel0,
                      double *yvel0)
{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;

  int j,k;
  double ugrad,vgrad,grad2,pgradx,pgrady,pgradx2,pgrady2,grad
        ,ygrad,pgrad,xgrad,div,strain2,limiter;
	
#pragma omp parallel	
 {
#pragma omp for private(ugrad,vgrad,div,strain2,pgradx,pgrady,pgradx2,pgrady2,limiter,pgrad,xgrad,ygrad,grad,grad2)
  for (k=y_min;k<=y_max;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max;j++) {

      ugrad=(xvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
            +xvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
           -(xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
            +xvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]);

      vgrad=(yvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
            +yvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
           -(yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
            +yvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]);

      div=(celldx[FTNREF1D(j,x_min-2)]*(ugrad)
          +celldy[FTNREF1D(k,y_min-2)]*(vgrad));

      strain2=0.5*(xvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                  +xvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                  -xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                  -xvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])/celldy[FTNREF1D(k,y_min-2)]
             +0.5*(yvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                  +yvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                  -yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                  -yvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)])/celldx[FTNREF1D(j,x_min-2)];

      pgradx=(pressure[FTNREF2D(j+1,k  ,x_max+4,x_min-2,y_min-2)]
             -pressure[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)])
            /(celldx[FTNREF1D(j,x_min-2)]+celldx[FTNREF1D(j+1,x_min-2)]);
      pgrady=(pressure[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
             -pressure[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)])
            /(celldy[FTNREF1D(k,y_min-2)]+celldy[FTNREF1D(k+1,y_min-2)]);

      pgradx2 = pgradx*pgradx;
      pgrady2 = pgrady*pgrady;

      limiter = ((0.5*(ugrad)/celldx[FTNREF1D(j,x_min-2)])*pgradx2+(0.5*(vgrad)/celldy[FTNREF1D(k,y_min-2)])*pgrady2+strain2*pgradx*pgrady)
              /MAX(pgradx2+pgrady2,1.0e-16);

      pgradx = SIGN(MAX(1.0e-16,fabs(pgradx)),pgradx);
      pgrady = SIGN(MAX(1.0e-16,fabs(pgrady)),pgrady);
      pgrad = sqrt(pgradx*pgradx+pgrady*pgrady);
      xgrad = fabs(celldx[FTNREF1D(j,x_min-2)]*pgrad/pgradx);
      ygrad = fabs(celldy[FTNREF1D(k,y_min-2)]*pgrad/pgrady);
      grad  = MIN(xgrad,ygrad);
      grad2 = grad*grad;

      if(limiter>0.0 || div>=0.0){
        viscosity[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=0.0;
      }
      else{
        viscosity[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=2.0*density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]*grad2*limiter*limiter;
      }

    }
  }

 }

}
