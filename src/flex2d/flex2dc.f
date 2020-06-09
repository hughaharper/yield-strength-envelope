c
      program flex2d
c
c*****************   main program   ************************
c
c Program to compute flexure of a thin plate due to a surface
c load on a plate of variable rigidity.  Follows the methods of 
c Sandwell 1984 by recasting the 2D flexure equation as a Fredholm
c integral equation of the second kind.  This integral equation is
c then solved iteratively via the method of successive approximation.
c (possibly making use of the successive over relaxation technique
c to more quickly achieve convergence).  Once the flexure surface
c has been satisfactorily determined, gravity from a deflected
c moho is computed.
c
c***********************************************************
c
      real*8 kx,ky,kh2,kh,beta
      real*8 x0,y0,dx,dy
      real*8 rland,rdum,trans
c
c  change ni and nj as needed
c
      parameter(ni=2048,nj=2048,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
      character*80 cvload,crig,cisurf,cwout
      character*80 cgrav,ccurv,title
      
      real*4 vload(nj*ni),rigid(nj*ni),curv(nj*ni)
      real*4 fz(nj,ni)
      real*4 D(nj,ni),b(nj,ni),bxy(nj,ni)
      real*4 wold(nj,ni),wnew(nj,ni)
      real*4 gz(nj,ni),dp(nj,ni)
      real*4 grav(nj,ni)
      real*4 dtx(nj,ni), dtxy(nj,ni), dty(nj,ni) 
      real*4 wtx(nj,ni), wtxy(nj,ni), wty(nj,ni)
      real*4 xwind(nj),ywind(ni)
      real*4 sxwind(nj),sywind(ni)
      real*4 swx, xwy, rigsub, tefac, tesub, temean, tesum
      
      real*8 vlmean, vlsum
      
      complex*8 fkz(nj2,ni),gkz(nj2,ni)
      complex*8 wconstk(nj2,ni),wpertk(nj2,ni)
      complex*8 woldk(nj2,ni),wnewk(nj2,ni)
      complex*8 Dk(nj2,ni),bk(nj2,ni),bxyk(nj2,ni)
      complex*8 gravk(nj2,ni)
      complex*8 dpk(nj2,ni) 
      complex*8 dtxk(nj2,ni), dtxyk(nj2,ni), dtyk(nj2,ni)
      complex*8 wtk(nj2,ni) 
      complex*8 wtxk(nj2,ni), wtxyk(nj2,ni), wtyk(nj2,ni)
      complex*8 fzval,gzval,tkval,grval,wctemp,wdtemp
      complex*8 gsurf,gmoho
      dimension n(2)
      complex*8 work(nwork)
      equivalence (fz(1,1),fkz(1,1))
      equivalence (gz(1,1),gkz(1,1))
      equivalence (grav(1,1),gravk(1,1))
      equivalence (D(1,1),Dk(1,1))
      equivalence (dp(1,1),dpk(1,1))
      equivalence (b(1,1),bk(1,1))
      equivalence (bxy(1,1),bxyk(1,1))
      equivalence (wold(1,1),woldk(1,1))
      equivalence (wnew(1,1),wnewk(1,1))
      equivalence (dtx(1,1),dtxk(1,1))
      equivalence (dtxy(1,1),dtxyk(1,1))
      equivalence (dty(1,1),dtyk(1,1))
      equivalence (wtx(1,1),wtxk(1,1))
      equivalence (wtxy(1,1),wtxyk(1,1))
      equivalence (wty(1,1),wtyk(1,1))
c
      pi=acos(-1.)
c
c  zero the arrays fz and gz
c
      do 30 i=1,ni
      do 30 j=1,nj
      fz(j,i)=0.
      gz(j,i)=0.
  30  continue
c
c  set the dimensions for fourt
c
      n(1)=nj
      n(2)=ni
c
c   get values from command line
c
      narg = iargc()
      if(narg.lt.6) then
        write(*,'(a)')'  '
        write(*,'(a)')
     & 'Usage: flex2dc vload.grd rig.grd isurf w.grd curv.grd grav.grd'
        write(*,'(a)')
        write(*,'(a)')
     &  '       vload.grd - surface topography load [m]'
        write(*,'(a)')
     &  '       rig.grd   - variable rigidity parameter (D, Te, or age)'
        write(*,'(a)')
     &  '       isurf     - which surface type to model?'
        write(*,'(a)')
     &  '                 - (0) land topography load'
        write(*,'(a)')
     &  '                 - (1) ocean topography load'
        write(*,'(a)')
     &  '                 - (2) tidal loading of floating ice'
        write(*,'(a)')
     &  '       w.grd     - output flexural surface [m]'
        write(*,'(a)')
     &  '       curv.grd  - output plate curvature [1/m]'
        write(*,'(a)')
     &  '       grav.grd  - output gravity value [mGal]'
        write(*,'(a)')
        stop
      else 
        call getarg(1,cvload)
        call getarg(2,crig)
        call getarg(3,cisurf)
        call getarg(4,cwout)
        call getarg(5,ccurv)
        call getarg(6,cgrav)
        read(cisurf,*)isurf
      endif
c
c  set the constants for the surface type
c
      if(isurf.eq.0) then
c
c  land topography load
c
        rhow=0.0
        rhoc=2800.
        rhom=3300.
        young=7.e10
c
c   seafloor topography load
c
      else if(isurf.eq.1) then
        rhow=1025.
        rhoc=2800.
        rhom=3300.
        young=6.5e10
c
c  ice parameters for tidal loading
c
      else if(isurf.eq.2) then
        rhow=0.
        rhoc=1025
        rhom=1025.
        young=3.e9
      else
        write(*,*)'Error: Incorrect surface type!'
        stop
      endif
c
c  set some constants
c  the sor = 0.80 based on trial runs
c
      sor=0.80
      Gcns=2.*pi*6.673e-11
      g=9.81
      rnu=0.25
c
c   set mean ocean depth and moho depth
c
      zs=4000.
      zc=zs+7000.
c
c   read the grd files
c
      call readgrd(vload,nj1,ni1,y0,x0,
     +            dy,dx,rdum,title,trim(cvload)//char(0))
      if(ni1.gt.ni.or.nj1.gt.nj) then
        print *,ni1,ni,nj1,nj,y0,x0,dy,dx,rdum
        write(*,'(a)')' recompile program to increase model dimensions'
        stop
      endif

      call readgrd(rigid,nj1,ni1,y0,x0,
     +            dy,dx,rdum,title,trim(crig)//char(0))
      if(ni1.gt.ni.or.nj1.gt.nj) then
        print *,ni1,ni,nj1,nj,y0,x0,dy,dx,rdum
        write(*,'(a)')' recompile program to increase model dimensions'
        stop
      endif
      
c
c  compute the windows for the full grid
c
      nsigy=ni/8.0
      do 70 i=1,ni
      if(i.lt.nsigy) then
       ywind(i)=0.5*(1.-cos(pi*(i-1)/nsigy))
      else if(i.gt.(ni-nsigy)) then
       ywind(i)=0.5*(1.-cos(pi*(ni-i)/nsigy))
      else
       ywind(i)=1.
      endif
   70 continue
      nsigx=nj/8.0
      do 80 j=1,nj
      if(j.lt.nsigx) then
       xwind(j)=0.5*(1.-cos(pi*(j-1)/nsigx))
      else if(j.gt.(nj-nsigx)) then
       xwind(j)=0.5*(1.-cos(pi*(nj-j)/nsigx))
      else
       xwind(j)=1.
      endif
   80 continue

c              
c  compute the height and width of the area in m
c  and also the sigma smoothing factor for the perturbation
c
      width=nj*dx*1.e3
      height=abs(ni*dy)*1.e3
      sigpert=4.*dx*1.e3

c
c  compute the windows for subgrid
c
      swy=ni1/8.
      do 104 i=1,ni1
      if(i.lt.swy) then
       sywind(i)=0.5*(1.-cos(pi*(i-1)/swy))
      else if(i.gt.(ni1-swy)) then
       sywind(i)=0.5*(1.-cos(pi*(ni1-i)/swy))
      else
       sywind(i)=1.
      endif
  104 continue
      swx=nj1/8.
      do 105 j=1,nj1
      if(j.lt.swx) then
       sxwind(j)=0.5*(1.-cos(pi*(j-1)/swx))
      else if(j.gt.(nj1-swx)) then
       sxwind(j)=0.5*(1.-cos(pi*(nj1-j)/swx))
      else
       sxwind(j)=1.
      endif
  105 continue

c
c   compute mean of the subgrid
c

      Dsum=0.
      rigsub = 0.
      tesub = 0.
      tesum = 0.
      do 106 i=1,ni1
      do 106 j=1,nj1
         k=nj1*(i-1)+j
         Dsum=Dsum+rigid(k)
  106 continue
      Dmean=Dsum/(ni1*nj1)

c
c   remove the mean from the subgrid and reapply the window
c
      tefac = young/(12*(1-(rnu*rnu)))
      rigsub = 0.
      do 107 i=1,ni1
      do 107 j=1,nj1   
         k=nj1*(i-1)+j

         if((rigid(k)-Dmean).lt.0) then
          rigsub=-sqrt(sqrt(abs(rigid(k)-Dmean)))*sywind(i)*sxwind(j)
         else
          rigsub=sqrt(sqrt(abs(rigid(k)-Dmean)))*sywind(i)*sxwind(j)
         endif        
         rigid(k) = rigsub**4.0 + Dmean
         
         
         
         if (rigid(k) < 0) then
          write(*,108) Dmean, rigid(k)
  108     format('Dmean =',e11.4,', rigidity =',e11.4)                
         endif
  107 continue
c
c  combine the subgrid and mean and recompute the mean
c
      Dsum=0.
      val=0.
      do 115 i=1,ni
      do 115 j=1,nj
         k=nj1*(i-1)+j
         if(i.le.ni1.and.j.le.nj1) then
           D(j,i)=rigid(k)
         else
           D(j,i)=Dmean
         endif
  115 continue

c
c   low pass filter the complete rigidity grid to eliminate sharp steps
c
      call fourt(D,n,2,-1,0,work,nwork)
      do 116 i=1,ni
      ky=-(i-1)/height
      if(i.ge.ni2) ky= (ni-i+1)/height
      do 116 j=1,nj2
      kx=(j-1)/width
      kh2=kx*kx+ky*ky
      kh=sqrt(kh2)
      beta=2*pi*kh
c
c apply a Gaussian filter and normalize by number of grid points
c also zero all the nyquist values
c
      sig=40.*dx*1.e3
      Dk(j,i)=exp(-beta*beta*sig*sig)*Dk(j,i)/(ni*nj)
      if(i.eq.ni2.or.j.eq.nj2)Dk(j,i)=cmplx(0.,0.)
      dpk(j,i)=Dk(j,i)
 116  continue

      D0=real(Dk(1,1))     
      Dk(1,1)=cmplx(0.,0.)
      call fourt(D,n,2,1,-1,work,nwork)

c     call writegrd(D,nj,ni,y0,x0,dx,dy,rland,rdum,
c    +              trim(cgrav)//char(0),trim(cgrav)//char(0))
c
c  load the vertical force subgrid 
c
      vlsum=0.0
      vlmean=0.0
      do 117 i=1,ni1
      do 117 j=1,nj1
        k=nj1*(i-1)+j
        vlsum=vlsum+vload(k)
  117 continue
  
      vlmean=vlsum/(nj1*ni1)


	  vlsum=0.0
      do 118 i=1,ni1
      do 118 j=1,nj1
        k=nj1*(i-1)+j
        vload(k)=(vload(k)-vlmean)*sywind(i)*sxwind(j)
c	vload(k)=(vload(k))*sywind(i)*sxwind(j)
	    vlsum=vlsum+vload(k)
  118 continue
  
  	  vlmean=vlsum/(nj1*ni1)
c
c generate load (demeaned) and apply window
c generate varying part of D: D(x)=D0+D'(x)
c generate starting guess for wold
c need to reverse the sign of fz so positive load produces pos. topo
c
      call fourt(wold,n,2,-1,0,work,nwork)
      
      do 120 i=1,ni
      do 120 j=1,nj
       k=nj1*(i-1)+j
       fz(j,i)=0.0
       if(i.le.ni1.and.j.le.nj1) then
       fz(j,i)=-(vload(k))*xwind(j)*ywind(i)/(ni*nj)
c       fz(j,i)=-(vload(k)-vlmean)*xwind(j)*ywind(i)/(ni*nj)
       endif
       D(j,i)=(D(j,i))*xwind(j)*ywind(i)
       woldk(j,i)=cmplx(0.,0.)
  120 continue
c
c  prepare the fourier transform of dp
c
      dpk(1,1)=cmplx(0.,0.)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prep iteration by defining convergence tolerance
c counting variables, and FFTing wold and fz
c
      tolerance=1.e-3
      icount=0
      icountmax=60
      wcmag=0.

      call fourt(fz,n,2,-1,0,work,nwork)
      fkz(1,1)=cmplx(0.,0.)
c     
c to test constant D solution, skip over perturbation stuff      
c      goto 899
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  begin Loop 1: b(k)=k^2*w(k)
c
  409 continue
      do 255 i=1,ni
      ky=-(i-1)/height
      if(i.ge.ni2) ky= (ni-i+1)/height
      do 255 j=1,nj2
      kx=(j-1)/width
      
      kh2=kx*kx+ky*ky
      bk(j,i)=kh2*woldk(j,i)/(ni*nj)
 255  continue

      do 257 i=1,ni
      ky=-(i-1)/height
      if(i.ge.ni2) ky= (ni-i+1)/height
      do 257 j=1,nj2
      kx=(j-1)/width

      if(icount.eq.0)then
      dtxk(j,i)=(kx**2)*dpk(j,i)/(ni*nj)
      dtxyk(j,i)=(kx*ky)*dpk(j,i)/(ni*nj)
      dtyk(j,i)=(ky**2)*dpk(j,i)/(ni*nj)
      endif
    
      wtxk(j,i)=(ky**2)*woldk(j,i)/(ni*nj)
      wtxyk(j,i)=(kx*ky)*woldk(j,i)/(ni*nj)
      wtyk(j,i)=(kx**2)*woldk(j,i)/(ni*nj)
        
 257  continue

c
c IFFT to get b(x) and all the cross terms
c
      call fourt(b,n,2,1,-1,work,nwork)
      if(icount.eq.0)then
      call fourt(dtx,n,2,1,-1,work,nwork)
      call fourt(dtxy,n,2,1,-1,work,nwork)
      call fourt(dty,n,2,1,-1,work,nwork)
      endif
      call fourt(wtx,n,2,1,-1,work,nwork)
      call fourt(wtxy,n,2,1,-1,work,nwork)      
      call fourt(wty,n,2,1,-1,work,nwork)      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  begin Loop 2: b(x)*D'(x)
c
      do 260 i=1,ni
      do 260 j=1,nj
c construct new b matrix and all the cross terms here      
       b(j,i)=b(j,i)*D(j,i)
       bxy(j,i)=dtx(j,i)*wtx(j,i)
       bxy(j,i)=bxy(j,i)-2*dtxy(j,i)*wtxy(j,i)
       bxy(j,i)=bxy(j,i)+dty(j,i)*wty(j,i)
  260 continue
c
c FFT to get back to b(k): this is the value of the Fredholm integral
c
      call fourt(b,n,2,-1,0,work,nwork)
      call fourt(bxy,n,2,-1,0,work,nwork)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  begin Loop 3: wnew(k)=T[wold(k)], and calculate convergence criteria
c
      wdiff=0.
      wnmag=0.
      womag=0.
      
      do 275 i=1,ni
      ky=-(i-1)/height
      if(i.ge.ni2) ky= (ni-i+1)/height
      do 275 j=1,nj2
      kx=(j-1)/width
      
      kh2=kx*kx+ky*ky
      kh=sqrt(kh2)
      beta=2*pi*kh
      den=D0*beta**4 + (rhom-rhow)*g
      
      if(icount.eq.0)then
        wconstk(j,i)=fkz(j,i)/den
    	wcmag=wcmag+(real(wconstk(j,i))**2+aimag(wconstk(j,i))**2)**.5
        wconstk(1,1)=cmplx(0.0)
      endif

c
c compute the perturbation to the depth 
c kill off the shortest wavelengths to aid convergence
c
      filt=exp(-kh2*sigpert*sigpert)
      wpertk(j,i)=filt*(2*pi)**4*(kh2*bk(j,i)-bxyk(j,i)*(1-rnu))/den
      wpertk(1,1)=cmplx(0.,0.)
      if(i.eq.ni2.or.j.eq.nj2)wpertk(j,i)=cmplx(0.,0.)
      if(icount.eq.0) then
        wnewk(j,i)=wconstk(j,i)-wpertk(j,i)
      else
        wnewk(j,i)=sor*(wconstk(j,i)-wpertk(j,i)) + (1.-sor)*woldk(j,i)
      endif
      

      
      wdtemp=wnewk(j,i)-woldk(j,i)
      wdiff=wdiff+(real(wdtemp)**2+aimag(wdtemp)**2)**.5
      wnmag=wnmag+(real(wnewk(j,i))**2+aimag(wnewk(j,i))**2)**.5
      womag=womag+(real(woldk(j,i))**2+aimag(woldk(j,i))**2)**.5
      
 275  continue

      wnewk(1,1)=cmplx(0.,0.)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c test for convergrnce and either loop back or go forward
c
      icount=icount+1
      err=wdiff/wcmag
c      write(*,697)icount,wnmag,err
      write(*,698)icount,icount,wnmag,icount,icount-1,wdiff,err
c      write(*,699)icount,wdiff,wcmag,err
  697 format(i3,': |wn| =',e11.4,', err =',e11.4)
  698 format(i3,': |w',i2,'| =',e11.4,', |w',i2,'-w',i2,'| =',e11.4,
     & ', err =',e11.4)
  699 format(i3,': wdiff =',e11.4,', wconst =',e11.4,', error =',e11.4)
      
      if(err.gt.tolerance.and.icount.lt.icountmax) then
        do 475 i=1,ni
        do 475 j=1,nj2
          woldk(j,i)=wnewk(j,i)
 475    continue

c
c   go back and do another iteration
c
      goto 409
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c DEBUGGING TEST: CAN I GET THE RIGHT ANSWER WITH CONSTANT D?
c
c  899 continue
c      do 888 i=1,ni
c      ky=-(i-1)/height
c      if(i.ge.ni2) ky= (ni-i+1)/height
c      do 888 j=1,nj2
c      kx=(j-1)/width
c      
c      kh2=kx*kx+ky*ky
c      kh=sqrt(kh2)
c      beta=2*pi*kh
c      wnewk(j,i)=-fkz(j,i)/(D0*beta**4 + (rhom-rhow)*g)
c 888  continue
c      write(*,*)D0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c
c compute gravity and curvature
c
      do 777 i=1,ni
        ky=-(i-1)/height
      if(i.ge.ni2) ky= (ni-i+1)/height
      do 777 j=1,nj2
        kx=(j-1)/width
        gzval=wnewk(j,i)
        kh2=kx*kx+ky*ky
        kh=sqrt(kh2)
        beta=2*pi*kh
c
c   upward continue each density interface
c
      gsurf=Gcns*gzval*(rhoc-rhow)*exp(-beta*zs)
	  gmoho=Gcns*gzval*(rhom-rhoc)*exp(-beta*zc)
      gravk(j,i)=(gsurf+gmoho)*1.e5
c
c   compute the plate curvature
c
        woldk(j,i)=-beta**2*wnewk(j,i)
	
 777  continue
c
c IFFT to get wnew(x), and grav(x), and curvature
c
      wnewk(1,1)=cmplx(0.0,0.0) 

      call fourt(grav,n,2,1,-1,work,nwork)
      call fourt(wnew,n,2,1,-1,work,nwork)
      call fourt(wold,n,2,1,-1,work,nwork)
c  899 continue
c
c  put the results into the subgrids for output
c
      do 901 i=1,ni1
      do 901 j=1,nj1
         k=nj1*(i-1)+j
         vload(k)=wnew(j,i)
         rigid(k)=grav(j,i)
         curv(k)=wold(j,i)
 901  continue
c
c  write out flexure and gravity grd files
c  some parameters must be real*8
c
      rland=9998.
      rdum=9999.
      
      call writegrd(vload,nj1,ni1,y0,x0,dy,dx,rland,rdum,
     +              trim(cwout)//char(0),trim(cwout)//char(0))
      call writegrd(rigid,nj1,ni1,y0,x0,dy,dx,rland,rdum,
     +              trim(cgrav)//char(0),trim(cgrav)//char(0))
      call writegrd(curv,nj1,ni1,y0,x0,dy,dx,rland,rdum,
     +              trim(ccurv)//char(0),trim(ccurv)//char(0))
      
      stop
      end
