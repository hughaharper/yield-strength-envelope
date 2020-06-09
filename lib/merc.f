c
      subroutine merc(rln,rlt,xj,yi,icall)
c
c  routine to compute the pixel locations xj, yi associated with the
c  mercator projection of rln, rlt.
c
c  input
c   rln   -   longitude (deg)
c   rlt   -   latitude (deg)
c   icall -   0-set up grid parameters
c             1-calculate index
c
c output
c   xj     -   column of matrix for rln 
c   yi     -   row of matrix for rlt
c
      implicit real*8 (a-h,o-z)
      common /grid/nlt,nln,dlt,dln,iproj
      common /bounds/rlt0,rltf,rln0,rlnf,rlnm
      data rad /.0174533/
      save arg,rad
c
c  if icall equals 0 setup the mapping and find the boundaries
c
      if(icall.eq.0) then
      rlnf=rln0+dln*nln
c
c  check to see if the left side of the box is negative
c  and add 360. if it's true.
c
      if(rln0.lt.0.) then
      rln0=rln0+360.
      rlnf=rlnf+360.
      endif
c
c  compute the maximum latitude
c
      arg=log(tan(rad*(45.+rlt0/2.)))
      arg2=rad*dln*nlt+arg
      term=exp(arg2)
      rltf=2.*atan(term)/rad-90.
c
c  print corners of area
c
c     write(*,903)
c 903 format('  corners of area  ')
c     write(*,904)rltf,rln0,rltf,rlnf
c     write(*,904)rlt0,rln0,rlt0,rlnf
c 904 format(2f7.2,6x,2f7.2)
      return
      endif
c
c compute the indices of the point
c
      if(icall.eq.1) then
        rln1=rln
        arg1=log(tan(rad*(45.+rlt0/2.)))
        arg2=log(tan(rad*(45.+rlt/2.)))
        yi=nlt-(arg2-arg1)/(dln*rad)
c       if(yi.lt.0.or.yi.gt.nlt) yi=-1
  20    continue
        xj=(rln1-rln0)/dln
        xj2=xj
c
c  check to see if the point lies to the left of the box
c
        if(xj2.lt.0) then
        rln1=rln1+360.
        if(rln1.le.rlnf)go to 20
        endif
c       if(xj.lt.0.or.xj.gt.nln) xj=-1
      else
c
c  compute latitude and longitude
c
        if(yi.lt.0.or.yi.gt.nlt) then
        rlt=-999.
        return
        endif
        if(xj.lt.0.or.xj.gt.nln) then
        rln=-999.
        return
        endif
        arg1=rad*dln*(nlt-yi)
        arg2=log(tan(rad*(45.+rlt0/2.)))
        term=exp(arg1+arg2)
        rlt=2.*atan(term)/rad-90.
        rln=rln0+dln*xj
      endif
      return
      end
