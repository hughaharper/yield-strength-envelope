/***************************************************************************/
/* grd_interp reads a files of age and curvature as well as a              */
/* grid of values such as flexural rigidity where the horizontal axis      */
/* is age and the vertical axis is curvature.  An output file is created   */
/* to match the age and curvature of the input files.                      */
/***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  01/02/14                                                      *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 ***************************************************************************/

# include "gmt.h"
# include "netcdf.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

double bilinear(double *, float *, int , int );

int main (int argc, char **argv)
{
	int ibufsizer,ibufsizea;
	int ii, jj, kk, xdimm, ydimm;
	int nage, ncurv;
	double age0, agef, dage, curv0, curvf, dcurv;
        double xys[2];
	char *file1, *file2, *file3, *file4;
	int argc2 = 1;
	char *argv2[2] = {"dummy",0};
	float *rig_in, *age, *curv, *rig_out;
	struct GRD_HEADER grd_rig, grd_age, grd_curv;

/* execute GMT_begin */
   	argc2 = GMT_begin (argc2, argv2);

/* get the information from the command line */
	if(argc < 5){
		printf("\nUsage: grd_interp rigid_in.grd age.grd curv.grd rigid_out.grd \n \n");
		printf("   rigid_in.grd   - master data to be interpolated created by curv2rigid \n");
		printf("   age.grd        - grid of ages \n");
		printf("   curv.grd       - grid of curvature \n");
		printf("   rigid_out.grd  - output rigidity grid \n \n");
		exit(-1);
	}
	file1 = argv[1];
	file2 = argv[2];
	file3 = argv[3];
	file4 = argv[4];

/* read the header of the rigid_in, age, and curvature files*/
	if (GMT_read_grd_info (file1, &grd_rig)) {
		fprintf (stderr, "Error opening file %s\n", file1);
		exit (EXIT_FAILURE);
	}
	if (GMT_read_grd_info (file2, &grd_age)) {
		fprintf (stderr, "Error opening file %s\n", file1);
		exit (EXIT_FAILURE);
	}
	if (GMT_read_grd_info (file3, &grd_curv)) {
		fprintf (stderr, "Error opening file %s\n", file2);
		exit (EXIT_FAILURE);
	}

/* make sure the dimensions of the age and curvature files match */
        if(grd_age.nx != grd_curv.nx || grd_age.ny != grd_curv.ny) {
		fprintf (stderr, "file dimensions do not match \n");
		exit (EXIT_FAILURE);
	}

/* allocate the memory for the four files */
        ibufsizer = grd_rig.nx * grd_rig.ny;
        ibufsizea = grd_age.nx * grd_age.ny;
	if((rig_in = (float *) malloc(ibufsizer*sizeof(float))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for rigidity file.\n");
		exit(-1);
        }
	if((age = (float *) malloc(ibufsizea*sizeof(float))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for age file.\n");
		exit(-1);
        }
	if((curv = (float *) malloc(ibufsizea*sizeof(float))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for curvature file.\n");
		exit(-1);
        }
	if((rig_out = (float *) malloc(ibufsizea*sizeof(float))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for rigidity output file.\n");
		exit(-1);
        }

/* read the three input grids */
	if(GMT_read_grd(file1, &grd_rig, rig_in, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) {
		fprintf (stderr, "Error reading file %s\n", file1);
		exit (EXIT_FAILURE);
	}
	if(GMT_read_grd(file2, &grd_age, age, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) {
		fprintf (stderr, "Error reading file %s\n", file1);
		exit (EXIT_FAILURE);
	}
	if(GMT_read_grd(file3, &grd_curv, curv, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) {
		fprintf (stderr, "Error reading file %s\n", file2);
		exit (EXIT_FAILURE);
	}

/* get the parameters to map age an curv onto the rig_in grid */
	age0=grd_rig.x_min;
	agef=grd_rig.x_max;
	curv0=grd_rig.y_min;
	curvf=grd_rig.y_max;
        nage=grd_rig.nx;
        ncurv=grd_rig.ny;
	dage=(agef-age0)/(nage-1);
        dcurv=(curvf-curv0)/(ncurv-1);

/* make sure this grid is pizel registered */
	if(grd_rig.node_offset != 1) {
		fprintf(stderr,"/n **** error -  grid is node registered /n");
		exit(-1);
	}

/* create the interpolated grid */
	ydimm = grd_age.ny;
	xdimm = grd_age.nx;
	for(ii=0; ii<ydimm; ii++) {
		for(jj=0; jj<xdimm; jj++) {
			kk = jj + xdimm * ii;
/* convert the age and curvature coordinates to x and y pixels */
			xys[0] = (age[kk]-age0)/dage;
			xys[1] = (curvf-curv[kk])/dcurv;
                        rig_out[kk]=(float) bilinear(xys,rig_in,nage,ncurv);
		}
	}

/* write the interpolated grid */
	GMT_write_grd(file4, &grd_age, rig_out, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE);		
	free(rig_out);
	free(rig_in);
	free(age);
	free(curv);
   return(0);
}
