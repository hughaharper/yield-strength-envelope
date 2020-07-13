/*******************************************************************************
 * Compute temperature and yield strength envelopes for the oceanic lithosphere*
 *******************************************************************************/

#include "litho.h"

void set_litho_defaults_(Litho *);
void print_litho_defaults_(Litho *);
double temp_plt_(Litho *, double *, double *);
double depth_sflr_(Litho *, double *);
double pressure_(Litho *, double *, double *, unsigned int *);
double byerlee_(Litho *, double *, double *, unsigned int *);
double ductile_(Litho *, double *, unsigned int *);
double elbendstress_(Litho *, double *, double *);

int main (int argc, char **argv)
{
	unsigned int verbose = 0;
        int i,j,nz=1000;
	double age,curv,z,dz,dsf,temp,obp,zn,zp,zmt,bystr,dustr,ystr;
	double ystrp,ystrm,moment;

        /* variables for root finding algorithm */
        double ival1, ival2, wi, itol=1.0e-4, stol=1.e4;
        double zmloc, z1loc, z2loc;
        double ebstrm, ebstr1, ebstr2;
        double efstrm, efstr1, efstr2;
        double npstr, npstrm, npstr1, npstr2;
        
        int maxiter=40;
        
        /* switches */
        unsigned int wcsw = 1; /* do not include water column overburden = 0, include it = 1 */
        unsigned int bysw; /* compression = 0, tension = 1 */
        unsigned int dusw; /* dorn law = 0, power law = 1 */

        Litho l;		
        Litho *lptr = &l; /* pointer to litho structure */

	/* get the information from the command line */
	    
        if(argc < 3){
                printf("\n Usage: test_litho age curvature \n \n");
                exit(-1);
        }

        age = atof(argv[1]);
        curv = atof(argv[2]);
        
        if (fabs(curv) == 0.0){
                printf("\n Please enter a non-zero value for curvature \n");
                exit(-1);
        } 
        
        /* insert some error handling for zero values of either curvature or age */ 

        set_litho_defaults_(lptr);
        if(verbose == 1) print_litho_defaults_(lptr);
        
        /* set value for in-plane stress (may be taken from user input later) */
        npstr=0.0;

	/* find depth of seafloor from thermal subsidence */
        dsf=depth_sflr_(lptr,&age);
        
        /* compute profile of ductile yielding to estimate "mechanical thickness" */
        zp = lptr->dp;
        dz = zp/nz;
        dusw = 0;
        zmt = zp;
        
        
        for(j=0;j<nz;j++){
        	 z=((double)j+0.5)*dz;
        	 temp=temp_plt_(lptr,&z,&age);
        	 dustr=ductile_(lptr,&temp,&dusw);
        	 if (dustr <= 2.e8 && dusw != 1) {
        			dusw = 1;
        			dustr=ductile_(lptr,&temp,&dusw);
        	 } 
        	 if ( dustr <= 5.e7) {
        			zmt=z;
        			break;
        	 } 
        
        }
        
        zn=0.5*zmt;
        if(verbose == 1) printf("age, Hm %f %f \n",age,zmt);
        
        /* set counter for root finding algorithm to zero */
        i=0;
        /* set initial boundaries of interval */
        ival1=1.5*dz;
        ival2=zmt-(1.5*dz);
        wi=ival2-ival1;
        
        do {
            zn=0.5*(ival1+ival2);
            npstrm=0.0;
            npstr1=0.0;
            npstr2=0.0;
            /* switch for function computing ductile deformation induced stresses */
            dusw = 0;
            /* inner loop computes stress profile */
        	for(j=0;j<nz;j++){
        		z=((double)j+0.5)*dz;
        		z1loc=(z-ival1);
        		z2loc=(z-ival2);
        		zmloc=(z-zn);
        		temp=temp_plt_(lptr,&z,&age);
			if( verbose == 1) printf(" age, depth, temp, obp, %f %f %f %f \n",age,z,temp,obp);
        		obp=pressure_(lptr,&z,&dsf,&wcsw);
			if( verbose == 1) printf(" age, depth, temp, obp, %f %f %f %f \n",age,z,temp,obp);
        		dustr=ductile_(lptr,&temp,&dusw);
        		ebstrm=elbendstress_(lptr,&zmloc,&curv);
        		ebstr1=elbendstress_(lptr,&z1loc,&curv);
        		ebstr2=elbendstress_(lptr,&z2loc,&curv);
				         		 
        		if (dustr <= 2.e8 && dusw != 1) {
        				dusw = 1;
        				dustr=ductile_(lptr,&temp,&dusw);
        		}  
        		 
        		if(curv*zmloc >= 0.0){
        			  bysw = 1;
    				  bystr = byerlee_(lptr,&z,&obp,&bysw);
    				  ystr = fmin(bystr,dustr);
        			  efstrm = fmin(ystr,ebstrm);
        		}
        		else {
        			  bysw = 0;
        			  bystr = byerlee_(lptr,&z,&obp,&bysw);
        			  ystr = fmax(bystr,-dustr);
        			  efstrm = fmax(ystr,ebstrm);
        		}
        		if(curv*z1loc >= 0.0){
        			  bysw = 1;
    				  bystr = byerlee_(lptr,&z,&obp,&bysw);
    				  ystr = fmin(bystr,dustr);
        			  efstr1 = fmin(ystr,ebstr1);
        		}
        		else {
        			  bysw = 0;
        			  bystr = byerlee_(lptr,&z,&obp,&bysw);
        			  ystr = fmax(bystr,-dustr);
        			  efstr1 = fmax(ystr,ebstr1);
        		}
        		if(curv*z2loc >= 0.0){
        			  bysw = 1;
    				  bystr = byerlee_(lptr,&z,&obp,&bysw);
    				  ystr = fmin(bystr,dustr);
        			  efstr2 = fmin(ystr,ebstr2);
        		}
        		else {
        			  bysw = 0;
        			  bystr = byerlee_(lptr,&z,&obp,&bysw);
        			  ystr = fmax(bystr,-dustr);
        			  efstr2 = fmax(ystr,ebstr2);
        		}
        		 	 
        	npstrm += efstrm*dz;
                npstr1 += efstr1*dz;
                npstr2 += efstr2*dz;               
        	}
        	
        	npstrm -= npstr;
		npstr1 -= npstr;
		npstr2 -= npstr;
            
            if(npstr1*npstr2 >= 0.0) {
                  fprintf(stderr,"Nodal plane not in specified interval. \n");
                  break;
            }
            
            if(fabs(npstrm) < stol) {
                  fprintf(stderr,"Tolerance value for stress difference reached. \n");
                  break;     
            }
        						  
        	i++;
        	
        	if(npstrm*npstr2 < 0) {
                 ival1=zn;        
        	}
        	else {
        	     ival2=zn;
        	}
        	
        	wi=ival2-ival1;
        	
        	if(wi < itol) {
        	     fprintf(stderr,"Tolerance value for width of interval reached. \n");
        	     break;
        	}			  
        	     
        	} while(i < maxiter);
        	
        	/* just print stress profile */
        	
        	dusw=0;
        	npstrm=0.0;
                moment=0.0;
        	for(j=0;j<nz;j++) {
                 z=((double)j+0.5)*dz;
                 zmloc=(z-zn);
        		 temp=temp_plt_(lptr,&z,&age);
        		 obp=pressure_(lptr,&z,&dsf,&wcsw);
        		 dustr=ductile_(lptr,&temp,&dusw);
        		 
        		 ebstrm=elbendstress_(lptr,&zmloc,&curv);
        		 
        		 if (dustr <= 2.e8 && dusw != 1) {
        				dusw = 1;
        				dustr=ductile_(lptr,&temp,&dusw);
        		 } 
			bysw = 1;
                        bystr = byerlee_(lptr,&z,&obp,&bysw);
    			ystrp = fmin(bystr,dustr);
			bysw = 0;
                        bystr = byerlee_(lptr,&z,&obp,&bysw);
    			ystrm = fmax(bystr,-dustr);
        		 
        		if(curv*zmloc >= 0.0){
        			bysw = 1;
    				bystr = byerlee_(lptr,&z,&obp,&bysw);
    				ystr = fmin(bystr,dustr);
        			efstrm = fmin(ystr,ebstrm);
        		}
        		else {
        			bysw = 0;
        			bystr = byerlee_(lptr,&z,&obp,&bysw);
        			ystr = fmax(bystr,-dustr);
        			efstrm = fmax(ystr,ebstrm);
        		}
			moment=moment+efstrm*zmloc*dz;
        		npstrm += efstrm*dz;
        		 
                fprintf(stdout,"%lf %lf %lf %lf %lf %lf \n",z*1.e-3,temp,obp*1.e-6,ystrp*1.e-6,ystrm*1.e-6,efstrm*1.e-6);
        	      
        	} 
        	
            fprintf(stderr,"end load, bending moment, zmloc %lf  %lf %lf \n",npstrm*1.e-6,moment*1.e-16,zn);
        
	return(EXIT_SUCCESS);
}
