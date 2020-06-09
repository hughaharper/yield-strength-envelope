#include "litho.h"

/*------------------------------------------------------------*/
/* computes yield strength according to form of Byerlee's law */
/* (values from Mueller and Philllips, 1995)                  */
/*------------------------------------------------------------*/

double temp_plt_(Litho *l, double *z, double *age);
double depth_sflr_(Litho *l, double *age);
double pressure_(Litho *l, double *z, double *dsf);
double byerlee_(Litho *l, double *z, double *obp, unsigned int *bysw);
double ductile_(Litho *l, double *temp, unsigned int *dusw);
double elbendstress_(Litho *l, double *zloc, double *curv);

double yse_moment_old_(Litho *l, double *age, double *curv, double *ival1, \
                   double *ival2, double *npstr, unsigned int *wcsw, \
                   unsigned int *znsw, int *rfiter)
{

     unsigned int verbose=0;

     /* variables for root finding algorithm */
     double wi, itol=1.0e2, stol=1.e6;
     double zmloc, z1loc, z2loc,zrt;
     double ebstrm, ebstr1, ebstr2;
     double efstrm, efstr1, efstr2;
     double npstrm, npstr1, npstr2;
     double z,zn,dz,dsf,temp,obp,zp,zmt,bystr,dustr,ystr,bendmo;
     double curvl, agel;


     int i, j;
     int maxiter = 40 , nz = 100;
     
     unsigned int bysw; /* compression = 0, tension = 1 */
     unsigned int dusw; /* dorn law = 0, power law = 1 */
     
     dz=l->dp/nz;
     curvl = *curv;
     agel = *age;
     
     /* find depth of seafloor from thermal subsidence,
        set to zero if water column overburden is not to be included */

     if(*wcsw == 1) dsf=depth_sflr_(l,age);
     else dsf =0.0;
     
     /* set counter for root finding algorithm to zero */
     i=0;
		
	 npstrm=0.0;
	 npstr1=0.0;
	 npstr2=0.0;
		

     wi=*ival2-*ival1;
	 zn=0.5*(*ival1+*ival2);
		
	 /* switch for function computing ductile deformation induced stresses */
	 dusw = 0;
	 

		
	 /* this loop computes stress profile */
	 for(j=0;j<nz;j++){
		 z=((double)j+0.5)*dz;
		 
		 z1loc=(z-*ival1);
		 z2loc=(z-*ival2);
		 zmloc=(z - zn);

		 temp=temp_plt_(l,&z,&agel);
		 obp=pressure_(l,&z,&dsf);
		 dustr=ductile_(l,&temp,&dusw);
		 

		 
		 ebstr1=elbendstress_(l,&z1loc,&curvl);
		 ebstr2=elbendstress_(l,&z2loc,&curvl);
						 
		 if (dustr <= 2.e8 && dusw != 1) {
				dusw = 1;
				dustr=ductile_(l,&temp,&dusw);
		 }  
		 
		 if(*curv*z1loc >= 0.0){
			  bysw = 1;
			  bystr = byerlee_(l,&z,&obp,&bysw);
			  ystr = fmin(bystr,dustr);
			  efstr1 = fmin(ystr,ebstr1);
		 }
		 else {
			  bysw = 0;
			  bystr = byerlee_(l,&z,&obp,&bysw);
			  ystr = fmax(bystr,-dustr);
			  efstr1 = fmax(ystr,ebstr1);
		}
		 
		if(*curv*z2loc >= 0.0){
			  bysw = 1;
			  bystr = byerlee_(l,&z,&obp,&bysw);
			  ystr = fmin(bystr,dustr);
			  efstr2 = fmin(ystr,ebstr2);
		}
		else {
			  bysw = 0;
			  bystr = byerlee_(l,&z,&obp,&bysw);
			  ystr = fmax(bystr,-dustr);
			  efstr2 = fmax(ystr,ebstr2);
		}
			 
		npstr1 += efstr1*dz;
		npstr2 += efstr2*dz;               

	}

	npstr1 -= *npstr;
	npstr2 -= *npstr;

                  
	if(npstr1*npstr2 >= 0.0) {
		 if(verbose == 1) fprintf(stderr,"Nodal plane not in specified interval. \n");
	     *znsw = 0;
		 return(EXIT_FAILURE);
	}
		
	if(npstr1 < 0) {
		 zrt=*ival1;
		 wi=*ival2-*ival1;
	}
	else {
		 zrt=*ival2;
		 wi=*ival1-*ival2;         
	}
            
            
	do {         
	
		 wi = 0.5*wi;
		 zn = zrt + wi;
		 
		 /* switch for function computing ductile deformation induced stresses */
		 dusw = 0;
		 npstrm=0.0;
		 
		 /* inner loop computes stress profile */
		 for(j=0;j<nz;j++){
			 z=((double)j+0.5)*dz;
			 
			 zmloc=z-zn;

			 temp=temp_plt_(l,&z,&agel);
			 obp=pressure_(l,&z,&dsf);
			 dustr=ductile_(l,&temp,&dusw);
			 
			 ebstrm=elbendstress_(l,&zmloc,&curvl);
							 
			 if (dustr <= 2.e8 && dusw != 1) {
					dusw = 1;
					dustr=ductile_(l,&temp,&dusw);
			 }  
			 
			 if(*curv*zmloc >= 0.0){
				  bysw = 1;
				  bystr = byerlee_(l,&z,&obp,&bysw);
				  ystr = fmin(bystr,dustr);
				  efstrm = fmin(ystr,ebstrm);
			 }
			 else {
				  bysw = 0;
				  bystr = byerlee_(l,&z,&obp,&bysw);
				  ystr = fmax(bystr,-dustr);
				  efstrm = fmax(ystr,ebstrm);
			 }
				 
			npstrm += efstrm*dz;
		}
		
		npstrm -= *npstr;
		
		if(npstrm <= 0.0) zrt = zn;
		
		
		if(fabs(npstrm) < stol) {
			  if(verbose == 1) fprintf(stderr,"Tolerance value for stress difference reached. \n");
			  break;     
		}
		
		if(wi < itol) {
			 if(verbose == 1) fprintf(stderr,"Tolerance value for width of interval reached. \n");
			 
			 break;
		}		
							  
		i++;
 
	} while(i < maxiter);
	
    *znsw = 1;
	*rfiter = i;

	dusw=0;
	npstrm=0.0;
	bendmo=0.0;
	
	// fprintf(stderr,"Nodal depth: %lf km \n",zn*1e-3);
	
	l->zn = zn;
	
	
	for(j=0;j<nz;j++) {
		 
		 z=((double)j+0.5)*dz;
		   
		 zmloc=(z-zn);

		 temp=temp_plt_(l,&z,&agel);
		 obp=pressure_(l,&z,&dsf);
		 dustr=ductile_(l,&temp,&dusw);
		 
		 ebstrm=elbendstress_(l,&zmloc,&curvl);
		 
		 if (dustr <= 2.e8 && dusw != 1) {
				dusw = 1;
				dustr=ductile_(l,&temp,&dusw);
		 } 
		 
		 if(*curv*zmloc >= 0.0){
			  bysw = 1;
			  bystr = byerlee_(l,&z,&obp,&bysw);
			  ystr = fmin(bystr,dustr);
			  efstrm = fmin(ystr,ebstrm);
		 }
		 else {
			  bysw = 0;
			  bystr = byerlee_(l,&z,&obp,&bysw);
			  ystr = fmax(bystr,-dustr);
			  efstrm = fmax(ystr,ebstrm);
		 }
		 
		 npstrm += efstrm*dz;
		 bendmo += efstrm*zmloc*dz;
		 // fprintf(stdout,"%lf %lf %lf %lf %lf %lf \n",z*1.e-3,temp,obp*1.e-6,ystr*1.e-6,ebstrm*1e-6,efstrm*1.e-6); 
    } 
    return bendmo;
}
