%
%  script to prepare three files used in curv2rigid.
%  age.grd spans 2 to 160 Ma
%  curv.grd spans -1.e-5 to 1.e-5
%
nj=1024;
ni=2048;
age0=2.;
agef=160.;
dage=(agef-age0)/(nj-1);
age=(0:nj-1)*dage+age0;
curv0=-.9e-5;
curvf= .9e-5;
dcurv=(curvf-curv0)/(ni-1);
curv=(0:ni-1)*dcurv+curv0;
[A,C]=meshgrid(age,curv);
%M=8000*sqrt(A);
%mt0=min(min(M));
%mtf=max(max(M));
%
%   now write the three files
%
cdfwrite(A,[0,nj,0,ni,age0,agef,1,1,1],'age_sub.grd');
cdfwrite(C,[0,nj,0,ni,curv0,curvf,1,1,1],'curv_sub.grd');
%cdfwrite(M,[1,nj,1,ni,mt0,mtf,1,1,1],'mt.grd');
