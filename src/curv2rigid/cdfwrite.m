function cdfwrite(Z,D,filename,title)
%
% CDFWRITE Write matrix to a GMT netcdf grd-file
%
%	CDFWRITE(Z, D, 'filename') will write the matrix Z using the
%	GMT grdfile format.  The array D contains (xmin, xmax, ymin,
%	ymax, zmin, zmax, format, xinc, yinc) for this data set.  This code
%       assumes pixel registration (format - 1).
%
%	See also CDFREAD
%
  if(exist(filename,'file'))
        display(' ')
        display(' error: file exists ')
        display(' ')
  else
	dim = size(Z);
	nx = dim(2);
	ny = dim(1);
	nccreate(filename,'x','Dimensions',{'x',nx},'Datatype','single','Format','classic')
	nccreate(filename,'y','Dimensions',{'y',ny},'Datatype','single','Format','classic')
 	nccreate(filename,'z','Dimensions',{'x',nx,'y',ny},'Datatype','single','Format','classic')
	xmin=D(1);
    xmax=D(2);
	ncwriteatt(filename,'x','actual_range',[xmin xmax]);
	xmin=xmin+D(8)/2.;
    xmax=xmax-D(8)/2.;
	xvec=[xmin:D(8):xmax];
	length(xvec);
	% ncwrite(filename,'x',xvec);

	ymin=D(3);
    ymax=D(4);
	ncwriteatt(filename,'y','actual_range',[ymin ymax]);
	ymin=ymin+D(9)/2.;
    ymax=ymax-D(9)/2.;
	yvec=[ymin:D(9):ymax];
	length(yvec);
	% ncwrite(filename,'y',yvec);

%
% write out the grid
%
	zmin=min(min(Z));
	zmax=max(max(Z));
	ncwrite(filename,'z',Z');
	ncwriteatt(filename,'z','actual_range',[zmin zmax]);
  	ncwriteatt(filename,'z','node_offset',1);
  end
end
