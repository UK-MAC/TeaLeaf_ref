/*Crown Copyright 2014 AWE.
*
* This file is part of TeaLeaf.
*
* TeaLeaf is free software: you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the
* Free Software Foundation, either version 3 of the License, or (at your option)
* any later version.
*
* TeaLeaf is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public License along with
* TeaLeaf. If not, see http://www.gnu.org/licenses/. */

/**
 *  @brief C kernel to update the external halo cells in a chunk.
 *  @author David Beckingsale, Wayne Gaudin
 *  @details Updates halo cells for the required fields at the required depth
 *  for any halo cells that lie on an external boundary. The location and type
 *  of data governs how this is carried out. External boundaries are always
 *  reflective.
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void update_halo_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                        int *lft,int *bttm,int *rght,int *tp,
                        int *lft_bndry,int *bttm_bndry,int *rght_bndry,int *tp_bndry,
                        int *chunk_neighbours,
                        double *density,
                        double *energy0,
                        double *energy1,
                        double *u,
                        double *p,
                        double *sd,
                        int *fields,
                        int *dpth)
{

  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  int left=*lft;
  int bottom=*bttm;
  int right=*rght;
  int top=*tp;
  int left_boundary=*lft_bndry;
  int bottom_boundary=*bttm_bndry;
  int right_boundary=*rght_bndry;
  int top_boundary=*tp_bndry;
  int depth=*dpth;

  /* These need to be kept consistent with the data module to avoid use statement */
  int CHUNK_LEFT=1,CHUNK_RIGHT=2,CHUNK_BOTTOM=3,CHUNK_TOP=4,EXTERNAL_FACE=-1;

  int FIELD_DENSITY    = 1;
  int FIELD_ENERGY0    = 3;
  int FIELD_ENERGY1    = 4;
  int FIELD_U          =16;
  int FIELD_P          =17;
  int FIELD_SD         =17;
  int NUM_FIELDS       =18;

  int j,k;

  /* Update values in external halo cells based on depth and fields requested */

#pragma omp parallel
 {
  if(fields[FTNREF1D(FIELD_DENSITY,1)]==1) {
    if(chunk_neighbours[FTNREF1D(CHUNK_BOTTOM,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          density[FTNREF2D(j  ,1-k,x_max+4,x_min-2,y_min-2)]=density[FTNREF2D(j  ,0+k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_TOP,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          density[FTNREF2D(j  ,y_max+k,x_max+4,x_min-2,y_min-2)]=density[FTNREF2D(j  ,y_max+1-k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_LEFT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          density[FTNREF2D(1-j,k,x_max+4,x_min-2,y_min-2)]=density[FTNREF2D(0+j,k,x_max+4,x_min-2,y_min-2)];
        }
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_RIGHT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          density[FTNREF2D(x_max+j,k,x_max+4,x_min-2,y_min-2)]=density[FTNREF2D(x_max+1-j,k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }
  }

  if(fields[FTNREF1D(FIELD_ENERGY0,1)]==1) {
    if(chunk_neighbours[FTNREF1D(CHUNK_BOTTOM,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          energy0[FTNREF2D(j  ,1-k,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(j  ,0+k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_TOP,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          energy0[FTNREF2D(j  ,y_max+k,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(j  ,y_max+1-k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_LEFT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          energy0[FTNREF2D(1-j,k,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(0+j,k,x_max+4,x_min-2,y_min-2)];
        }
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_RIGHT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          energy0[FTNREF2D(x_max+j,k,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(x_max+1-j,k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }
  }

  if(fields[FTNREF1D(FIELD_ENERGY1,1)]==1) {
    if(chunk_neighbours[FTNREF1D(CHUNK_BOTTOM,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          energy1[FTNREF2D(j  ,1-k,x_max+4,x_min-2,y_min-2)]=energy1[FTNREF2D(j  ,0+k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_TOP,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          energy1[FTNREF2D(j  ,y_max+k,x_max+4,x_min-2,y_min-2)]=energy1[FTNREF2D(j  ,y_max+1-k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_LEFT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          energy1[FTNREF2D(1-j,k,x_max+4,x_min-2,y_min-2)]=energy1[FTNREF2D(0+j,k,x_max+4,x_min-2,y_min-2)];
        }
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_RIGHT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          energy1[FTNREF2D(x_max+j,k,x_max+4,x_min-2,y_min-2)]=energy1[FTNREF2D(x_max+1-j,k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }
  }

  if(fields[FTNREF1D(FIELD_U,1)]==1) {
    if(chunk_neighbours[FTNREF1D(CHUNK_BOTTOM,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          u[FTNREF2D(j  ,1-k,x_max+4,x_min-2,y_min-2)]=u[FTNREF2D(j  ,0+k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_TOP,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          u[FTNREF2D(j  ,y_max+k,x_max+4,x_min-2,y_min-2)]=u[FTNREF2D(j  ,y_max+1-k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_LEFT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          u[FTNREF2D(1-j,k,x_max+4,x_min-2,y_min-2)]=u[FTNREF2D(0+j,k,x_max+4,x_min-2,y_min-2)];
        }
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_RIGHT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          u[FTNREF2D(x_max+j,k,x_max+4,x_min-2,y_min-2)]=u[FTNREF2D(x_max+1-j,k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }
  }

  if(fields[FTNREF1D(FIELD_P,1)]==1) {
    if(chunk_neighbours[FTNREF1D(CHUNK_BOTTOM,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          p[FTNREF2D(j  ,1-k,x_max+4,x_min-2,y_min-2)]=p[FTNREF2D(j  ,0+k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_TOP,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          p[FTNREF2D(j  ,y_max+k,x_max+4,x_min-2,y_min-2)]=p[FTNREF2D(j  ,y_max+1-k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_LEFT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          p[FTNREF2D(1-j,k,x_max+4,x_min-2,y_min-2)]=p[FTNREF2D(0+j,k,x_max+4,x_min-2,y_min-2)];
        }
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_RIGHT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          p[FTNREF2D(x_max+j,k,x_max+4,x_min-2,y_min-2)]=p[FTNREF2D(x_max+1-j,k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }
  }

  if(fields[FTNREF1D(FIELD_SD,1)]==1) {
    if(chunk_neighbours[FTNREF1D(CHUNK_BOTTOM,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          sd[FTNREF2D(j  ,1-k,x_max+4,x_min-2,y_min-2)]=sd[FTNREF2D(j  ,0+k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_TOP,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (j=x_min-depth;j<=x_max+depth;j++) {
#pragma ivdep
        for (k=1;k<=depth;k++) {
          sd[FTNREF2D(j  ,y_max+k,x_max+4,x_min-2,y_min-2)]=sd[FTNREF2D(j  ,y_max+1-k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_LEFT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          sd[FTNREF2D(1-j,k,x_max+4,x_min-2,y_min-2)]=sd[FTNREF2D(0+j,k,x_max+4,x_min-2,y_min-2)];
        }
      }
    }

    if(chunk_neighbours[FTNREF1D(CHUNK_RIGHT,1)]==EXTERNAL_FACE) {
#pragma omp for private(j,k)
      for (k=y_min-depth;k<=y_max+depth;k++) {
#pragma ivdep
        for (j=1;j<=depth;j++) {
          sd[FTNREF2D(x_max+j,k,x_max+4,x_min-2,y_min-2)]=sd[FTNREF2D(x_max+1-j,k,x_max+4,x_min-2,y_min-2)];
	}
      }
    }
  }


 }

}
