/************************************************************************
* bilinear interpolator modified from GMTSAR                            *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 01/05/14                                                     *
************************************************************************/
#include <math.h>
#include <stdio.h>

double bilinear(double *xys, float *s_in, int xdims, int ydims)
{
	double dx, dy, sout;
	int k00, k01, k10, k11;
	int i0, j0;

	/* compute the residual offsets */
	j0 = (int)floor(xys[0]);
	i0 = (int)floor(xys[1]);
        dx = xys[0] - (double)j0;
        dy = xys[1] - (double)i0;
	if(dx < 0. || dx > 1. || dy < 0. || dy > 1) fprintf(stderr," dx or dy out of bounds %f %f \n",dx,dy);

	/* compute the indices of the 4 corners */
       
	k00 = xdims*i0     + j0;
	k01 = xdims*i0     + (j0+1);
	k10 = xdims*(i0+1) + j0;
	k11 = xdims*(i0+1) + (j0+1);

	/* do the interpolation if all 4 corners are within the bounds of the slave array */

	if(i0 < 0 || i0 >= (ydims-1) || j0 < 0 || j0 >= (xdims-1)) {
	  sout = 0;
	}
	else {
        sout = s_in[k00] * (1.0 - dy) * (1.0 - dx)
             + s_in[k10] * (dy)       * (1.0 - dx)
             + s_in[k01] * (1.0 - dy) * (dx) 
             + s_in[k11] * (dy)       * (dx);
	}
	return(sout);
}
