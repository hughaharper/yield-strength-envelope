c
      subroutine moment(x,y,x1,y1,x2,y2,dxx,sig,F1,F2,fzz)
c
c Routine to make line loads and line moments into a regular grid by approximating
c the lines as Gaussians and the moments derivatives of Gaussians.
c
      implicit real*8(a-h,o-z)
      real*8 L
c
c       create the force vector for this moment
c
c       input:
c
c       x, y    - computation point      (km)
c       x1,y1   - start position         (km)
c       x2,y2   - end position           (km)
c       dxx     - grid spacing           (km)
c       sig     - source width in pixels (try 1)
c       F1      - vertical load
c       F2      - vertical moment
c
c       output:
c
c       fzz     - z-component of force
c
c  set the x and y taper lengths
c
      pi = 3.14159265358979
      sigy=sig
c
c  sigx  should be = 2 or perhaps 4
c
      sigx=2.0
c
c  compute the length and orientation of the fault
c
      Dx=(x2-x1) 
      Dy=(y2-y1) 
      L=sqrt(Dx*Dx+Dy*Dy) 
      theta=atan2(Dy,Dx) 
      st=sin(theta) 
      ct=cos(theta) 

c
c  rotate the vector into the model space
c
      x0=x-x1 
      y0=y-y1 
      xp=x0*ct+y0*st 
      yp=-x0*st+y0*ct 
 
c  compute the h and dh/dx functions
c
      h=exp(-0.5*((yp/sigy)**2))/(sqrt(2*pi)*sigy)
      dh=-yp*exp(-0.5*((yp/sigy)**2))/(sqrt(2*pi)*sigy**3)
c
c   compute the g and dgdx functions
c
      if(xp .gt. -sigx .and. xp .lt. sigx) then
         g=0.5*(1-cos(pi*(xp+sigx)/(2*sigx))) 
         dg=pi*sin(pi*(xp+sigx)/(2*sigx))/(4*sigx)
      else if(xp .ge. sigx .and. xp .le. L-sigx) then
         g=1. 
         dg=0.
      else if(xp .gt. L-sigx .and. xp .lt. L+sigx) then
         g=0.5*(1-cos(pi*(L+sigx-xp)/(2*sigx))) 
         dg=-pi*sin(pi*(L+sigx-xp)/(2*sigx))/(4*sigx)
      else
         g=0. 
         dg=0.
      endif
c
c compute the sum of the two 
c
      fzz=g*(h*F1+dh*F2)
c
      return
      end
