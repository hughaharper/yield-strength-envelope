function[Z,D] = cdfread(name)
%
% CDFREAD Read matrix from a GMT nertcdf grd-file
%
%	Z = CDFREAD('filename') will return Z, the matrix stored in the
%	GMT grdfile format.
%
%	[Z,D] = CDFREAD('filename') will also return an array D that
%	contains (xmin, xmax, ymin, ymax, zmin, zmax, format, xinc, yinc)
%	for this data set.  This code assumes pixel registration.
%	
%	See also CDFWRITE
  if(exist(name,'file')) 
%
%   first read the header information
%
	finfo = ncinfo(name);
	xnm=finfo.Variables(1).Name;
	nx=finfo.Variables(1).Size;
	ynm=finfo.Variables(2).Name;
	ny=finfo.Variables(2).Size;
	znm=finfo.Variables(3).Name;
	X = ncread(name,xnm);
	Y = ncread(name,ynm);
	Z = ncread(name,znm)';
	dx=abs(X(2)-X(1));
	dy=abs(Y(2)-Y(1));
	d1=min(X)-dx/2.;
	d2=max(X)+dx/2.;
	d3=min(Y)-dy/2.;
	d4=max(Y)+dy/2.;
	d5=min(min(Z));
	d6=max(max(Z));
	d7=1;
%
	D=[d1,d2,d3,d4,d5,d6,d7,dx,dy];
  else
	D=0.;
	Z=0.;
	display(' ')
	display(' error: file not found ')
	display(' ')
  end
  end
