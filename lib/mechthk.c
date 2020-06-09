#include "litho.h"

/*------------------------------------------------------------*/
/* computes mechanical thickness of a plate given its age     */
/* assumes plate-cooling model                                */
/* (following McNutt and Menard, 1982)                        */
/*------------------------------------------------------------*/

double temp_plt_(Litho *l, double *z, double *age);
double ductile_(Litho *l, double *temp, unsigned int *dusw);

double mechthk_(Litho *l, double *age)
{
     double z,temp,dustr,mthk,agel;

     int i, nz;
     double dz = 1.e2;
     
     unsigned int dusw = 0; /* dorn law = 0, power law = 1 */
     
     agel = *age;
	 nz=(int)(floor(l->dp/dz));
     z = 0.5*dz;
     
     for(i=1;i<nz;i++){
          z+=dz;
 
          temp=temp_plt_(l,&z,&agel);
          dustr=ductile_(l,&temp,&dusw);
 
          if (dustr <= 2.e8 && dusw != 1) {
		       dusw = 1;
		       dustr=ductile_(l,&temp,&dusw);
          } 
 
          if ( dustr <= 5.e7) {
		       mthk=z;
		       break;
          } 
     }
     
     l->zmt = mthk;                

     return mthk;
}
