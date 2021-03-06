c
      program load
c
c  USAGE: load rln0 rlt0 segs.dat fz.grd 
c 
c  EXAMPLE:
c******************************************************************
c  load -70 5 segs.dat fz.grd
c******************************************************************
c  Trench elements are stored in the file segs.dat
c******************************************************************
c   rln1           rln2            rlt1          rlt2            Vo      Mo
c -68.2788067317 -67.2119819319   19.3195117256 19.3279862117   -1.e15 -1.e17
c -67.2119819319 -65.806047249    19.3279862117 19.3014021804   -1.e15 -1.e17
c -65.806047249  -64.619540884    19.3014021804 19.2798227214   -1.e15 -1.e17
c
c*****************   main program   ************************
c
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
c
c  User - change the parameters on the lines below to meet your needs
c*******************************************************************
c     parameter(ni=2048,nj=2048,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
c     parameter(ni=1800,nj=1500,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
c      parameter(ni=1587,nj=1800,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
c      parameter(ni=2,nj=2048,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
c      parameter(ni=581,nj=1200,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
      parameter(ni=906,nj=2040,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
c      parameter(ni=2029,nj=1920,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
c
c  parameter   description
c  ni          # of rows in load grid
c  nj          # of columns in load grid
c
c*******************************************************************
      real*8 kx,ky,kr2
      real*8 rln0,rlt0,dlt,dln,rland,rdum
      real*8 x0,y0,dx,dy
      character*80 fzdisp
      character*80 cseg,crln0,crlt0
c
      common /grid/nlt,nln,dlt,dln,iproj
      common /bounds/rlt0,rltf,rln0,rlnf,rlnm
c
      real*4 fz(nj,ni)
c
c set the parameters for the mercator projection
c
      nlt=ni
      nln=nj
      dln=1./60.
c
c  zero the array fz
c
      do 30 i=1,ni
      do 30 j=1,nj
      fz(j,i)=0.
  30  continue
c
c   get values from command line
c
      narg = iargc()
      if(narg.ne.4) then
       write(*,'(a)')'  '
       write(*,'(a)')
     & 'Usage: load rln0 rlt0 segs.dat fz.grd '
       write(*,'(a)')
       write(*,'(a)')
     & '       rln0 - longitude of origin '
       write(*,'(a)')
     & '       rlt0 - latitude of origin '
       write(*,'(a)')
     & '       segs.dat- trench segment file'
       write(*,'(a)')
     & '       fz.grd - output files of vertical load '
       write(*,'(a)')'  '
       stop
      else 
        call getarg(1,crln0)
        call getarg(2,crlt0)
        call getarg(3,cseg)
        read(crln0,*)rln0
        read(crlt0,*)rlt0
        call getarg(4,fzdisp)
        nc=index(fzdisp,' ')
        fzdisp(nc:nc)=char(0)
      endif
c
c   setup the mapping using a Mercator projection assuming a 1 minute grid
c
      nlt=ni
      nln=nj
      dln=1/60.
      call merc(rln,rlt,xj,yi,0)
c              
c  compute the height and width of the area in km
c
      rlat2=abs(rlt0+rltf)/2
      xscl=cos(rlat2/57.29578)
      write(*,*)rlt0,rltf,rlat2,xscl
      dx=xscl*111.*dln
      dy=dx
      width=nj*dx
      height=abs(ni*dy)
c
c  open the input file and load the force arrays
c
      open(unit=5,file=cseg,status='old')
c
c   generate the load from the segments
c   set sigma equal to one pixel
c
      sig=2.
      rland=0.
c
   50 read(5,*,end=999)rln1,rln2,rlt1,rlt2,F1,F2
      call merc(rln1,rlt1,x1,y1,1)
      call merc(rln2,rlt2,x2,y2,1)
      F1=F1/(dx*dx)
      F2=F2/(dx*dx*dx)
c
c  figure out the part of the array to add the new forces
c
      ipd=3*sig+2
      jx0=min(x1,x2)-ipd
      jxf=max(x1,x2)+.5+ipd
      iy0=min(y1,y2)-ipd
      iyf=max(y1,y2)+.5+ipd
      do 100 iy=iy0,iyf
      do 200 jx=jx0,jxf
        if(jx.le.0.or.jx.gt.nj) go to 200
        if(iy.le.0.or.iy.gt.ni) go to 200
        x=jx
        y=iy
        call moment(x,y,x1,y1,x2,y2,dx,sig,F1,F2,fzz)
c
c sum the forces and convert from km**2 to m**2
c
        fz(jx,iy)=fz(jx,iy)+fzz*1.e-6
c
c keep track of the max value
c
        rland=max(rland,fz(jx,iy))
 200  continue
 100  continue
      go to 50
 999  close(unit=5)
c
c  write fz.grd 
c
      rland=rland*1.01
      rdum=rland
c
c   write the grd files
c
      dlt=(rltf-rlt0)/ni
      x0=0.
      y0=0.
      call writegrd(fz,nj,ni,y0,x0,dy,dx,
     +              rland,rdum,fzdisp,fzdisp)
      stop
      end
