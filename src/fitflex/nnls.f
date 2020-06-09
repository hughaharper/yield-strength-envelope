c=======================================================================
c          UNIT N:  Solve Nonnegative Least Squares Problem
c=======================================================================
      subroutine nnls(a, mdim, m, n, b, x, rnorm, w, zz, index, mode)
      implicit double precision (a-h,o-z)                               *doub
c$$$$ calls h12,g1,g2
c  Non-negative least-squares fitting routine.  Given by Lawson & Hanson
c  in  `Solving Least Squares Problems ,' Prentice-Hall 1974.
c  Given  matrix  a  (m by n), and an m-vector  b  , finds the n-vector
c  x  that solves the problem 2-norm(a*x - b) minimum with x.ge.0.
c
c  a(),mdim,m,n  a  is the array, with mrows,  n  cols.  mdim  is the
c     actual dimension given  a  in the calling program (mdim.ge.m)
c     a  is destroyed by program.
c  b()  on entry must contain m-vector  b.  Destroyed by program
c  x()  solution vector  x.  Need not be initialized.
c  rnorm  minimum two-norm.
c  w()  an array of working space.  Must be dimensioned at least n
c  zz()  another working array.  Must be dimensioned at least m.
c  index  another working array. Must be dimensioned at least  n.
c  mode  a flag indicateng outcome of subroutine.
c     mode=1   means routine worked properly.
c     mode=2   means dimensions of arrays were bad (m.le0 .or. n.le.0)
c     mode=3   means iteration count exceeded (more than 3*n iterations)
c
      dimension a(mdim,n),b(m),x(n),w(n),zz(m),index(n),dummy(1)
c
      zero=0.
      one=1.
      two=2.
      factor=.01
c
      mode=1
      if (m.gt.0 .and. n.gt.0) go to 10
      mode=2
      return
 10   iter=0
      itmax=3*n
c
c
c  Initialize arrays index  and  x.
      do 20 i=1,n
      x(i)=zero
 20   index(i)=i
c
      iz2=n
      iz1=1
      nsetp=0
      npp1=1
c     Main loop begins here.
 30   continue
c
c  Quit if all coeffs are already in solution. Or if m cols of a have
c  been triangularized.
      if (iz1.gt.iz2 .or. nsetp.ge.m) go to 350
c
c
c  Compute components of dual vector w.
      do 50 iz=iz1,iz2
      j=index(iz)
      sm=zero
      do 40 l=npp1,m
 40   sm=a(l,j)*b(l) + sm
 50   w(j)=sm
c  Find largest +ve w(j).
 60   wmax=zero
      do 70 iz=iz1,iz2
      j=index(iz)
      if (w(j).le.wmax) go to 70
      wmax=w(j)
      izmax=iz
 70   continue
 
c
c
c
c  If wmax.le.0 go to termination (kuhn-tucker conditions ok).
      if (wmax) 350,350,80
 80   iz=izmax
      j=index(iz)
c
c
  
c
c  Sign of w(j) is ok for j to be moved to set  p.
c  Begin Householder trans. Check new diagonal element to avoid near ld.
      asave=a(npp1,j)

      call h12(1, npp1, npp1+1, m, a(1,j), 1, up, dummy, 1, 1, 0)

       
      unorm=zero
      if (nsetp.eq.0) go to 100
      do 90 l=1,nsetp
 90   unorm=a(l,j)**2 + unorm
 100  unorm=sqrt(unorm)
      temp=unorm + abs(a(npp1,j))*factor
      if (temp-unorm) 130,130,110
c  Col j is sufficiently indep. Copy  b  to zz. update  z  and solve
c  for  ztest (=proposed new val for x(j)).
c
 110  do 120 l=1,m
 120  zz(l)=b(l)
      call h12(2, npp1, npp1+1, m, a(1,j), 1, up, zz, 1, 1, 1)
      ztest=zz(npp1)/a(npp1,j)
c
c
c              See if  ztest is +ve.
      if (ztest) 130,130,140
c
c  Reject  j  as a candidate to be moved from set  z  to set  p.
c  Restore a(npp1,j), set w(j)=0. And loop back to test dual
c  coeffs again.
c
 130  a(npp1,j)=asave
      w(j)=zero
      go to 60

c
c  The index j=index(iz) has been selected to be moved ftom set  z
c  to set  p .  Update b,  update indices, apply householder trans
c  to cols in new set z.  Zero subdiagonal  elts in col  j, set
c  w(j)=0.
c
 140  do 150 l=1,m
 150  b(l)=zz(l)
c
      index(iz)=index(iz1)
      index(iz1)=j
      iz1=1 + iz1
      nsetp=npp1
      npp1=1 + npp1
c
      if (iz1.gt.iz2) go to 170
      do 160 jz=iz1,iz2
      jj=index(jz)
 160  call h12(2, nsetp, npp1, m, a(1,j), 1, up, a(1,jj), 1, mdim, 1)
 170  continue
c
      if  (nsetp.eq.m) go to 190
      do 180 l=npp1,m
 180  a(l,j)=zero
 190  continue
c
      w(j)=zero

c  Solve triangulate system. Store sol in zz,temporarily.
c
      next=1
      go to 400
 200  continue
c
c  Secondary loop begins here.
 210  iter=1 + iter
      if (iter.le.itmax) go to 220
      mode=3
c  Write statement here deleted*********************
      go to 350
 220  continue
c
c
c  See if all new constrained coeffs are feasible. If not, find alpha.
      alpha=two
      do 240 ip=1,nsetp
      l=index(ip)
      if (zz(ip)) 230,230,240
c
 230  t=-x(l)/(zz(ip) - x(l))
      if (alpha.le.t) go to 240
      alpha=t
      jj=ip
 240  continue
c  If all new constrained coeffs are feasible, alpha is still 2. If so,
c  exit into main loop.
c
      if (alpha.eq.two) go to 330
c  Otherwise alpha will be in (0,1) to interpolate between old x
c  and new zz.
c
      do 250 ip=1,nsetp
      l=index(ip)
 250  x(l)=alpha*(zz(ip)-x(l)) + x(l)
c
c  Modify  a  bb  and the index arrays to move coeff  i  from  set  p
c   to set  z.
c
      i=index(jj)
 260  x(i)=zero
c
      if (jj.eq.nsetp) go to 290
      jj=1 + jj
      do 280 j=jj,nsetp
      ii=index(j)
      index(j-1)=ii
      call g1(a(j-1,ii), a(j,ii), cc, ss, a(j-1,ii))
      a(j,ii)=zero
      do 270 l=1,n
      if (l.ne.ii) call g2(cc, ss, a(j-1,l), a(j,l))
 270  continue
 280  call g2(cc, ss, b(j-1), b(j))
 290  npp1=nsetp
      nsetp=nsetp - 1
      iz1=iz1 -     1
      index(iz1)=i
c  All coeffs in set p should be feasible.  If they are not, this is due
c  to round-off.  Non-positive ones set to 0 and moved to set z.
c
c
      do 300 jj=1,nsetp
      i=index(jj)
      if (x(i)) 260,260,300
 300  continue
c
c
c  Copy  b  to  zz  and solve again and loop back.
      do 310 i=1,m
 310  zz(i)=b(i)
      next=2
      go to 400
 320  continue
      go to 210
c
c  End of secondary loop.
 330  do 340 ip=1,nsetp
      i=index(ip)
 340  x(i)=zz(ip)
c  All new coeffs are +ve.  Loop back to beginning.
      go to 30
c  End of main loop.  Terminating section.
c
c
 350  sm=zero
      if (npp1.gt.m) go to 370
      do 360 i=npp1,m
 360  sm=b(i)**2 + sm
      go to 390
 370  do 380 j=1,n
 380  w(j)=zero
 390  rnorm=sqrt(sm)
      return
c
c  Solution of linear triangular system.
c
 400  do 430 l=1,nsetp
      ip=nsetp + 1 - l
      if (l.eq.1) go to 420
      do 410 ii=1,ip
 410  zz(ii)=zz(ii) - a(ii,jj)*zz(ip+1)
 420  jj=index(ip)
 430  zz(ip)=zz(ip)/a(ip,jj)
      go to (200, 320), next
c ******** format deleted for print statement here ****************
      end
c_______________________________________________________________________
      subroutine g2(cosa, sina, x, y)
      implicit double precision (a-h,o-z)                               *doub
c$$$$ calls no other routines
c  Rotates a 2-vector (x,y) by angle a. Over-writes original vector.
c
      z=cosa*x + sina*y
      y=cosa*y - sina*x
      x=z
      return
      end
c_______________________________________________________________________
      subroutine g1(a, b, cosa, sina, sig)
      implicit double precision (a-h,o-z)                               *doub
c$$$$ calls no other routines
c
c  From Lawson & Hanson `Solving Least Squares Problems' 1974 p 309.
c  Computes an orthogonal matrix to rotate (a,b) into (sqrt(a*a+b*b),0)
c  matrix in form (cosa, sina/-sina, cosa) and sig holds magnitude of
c  (a,b).  Sig is allowed to overwrite a or b in calling routine.
c
      zero=0.0
      one=1.0
      if (abs(a) .le. abs(b)) go to 10
      xr=b/a
      yr=sqrt(one + xr**2)
      cosa=sign(one/yr, a)
      sina=cosa*xr
      sig=abs(a)*yr
      return
 10   if (b) 20,30,20
 20   xr=a/b
      yr=sqrt(one + xr**2)
      sina=sign(one/yr, b)
      cosa=sina*xr
      sig=abs(b)*yr
      return
 30   sig=zero
      cosa=zero
      sina=one
      return
      end
c_______________________________________________________________________
      subroutine h12(mode, lpivot, l1, m, u, iue, up, c, ice, icv, ncv)
      implicit double precision (a-h,o-z)                               *doub
c$$$$ calls no other routines
c  Construction and application of householder transformation for nnls
c  from  h12  in Lawson & Hanson `Solving Least Squares Problems'
c  - Prentice-Hall 1974 (pp308,309).  Double precision throughout.
      dimension u(iue,m),c(*)
      if (0.ge.lpivot .or. lpivot.ge.l1 .or. l1.gt.m) return
      cl=abs(u(1,lpivot))
      if (mode.eq.2) go to 60
c  Construct transformation.
      do 10 j=l1,m
 10   cl=max(abs(u(1,j)),cl)
      if (cl) 130,130,20
 20   clinv=1.0/cl
      sm=(u(1,lpivot)*clinv)**2
      do 30 j=l1,m
 30   sm=(u(1,j)*clinv)**2 + sm
      sm1=sm
      cl=sqrt(sm1)*cl
      if (u(1,lpivot)) 50,50,40
 40   cl=-cl
 50   up=u(1,lpivot) - cl
      u(1,lpivot)=cl
      go to 70
c  Apply transformation.
 60   if (cl) 130,130,70
 70   if (ncv.le.0) return
      b=up*u(1,lpivot)
      if (b) 80,130,130
 80   b=1.0/b
      i2=1-icv+ice*(lpivot-1)
      incr=ice*(l1-lpivot)
      do 120 j=1,ncv
      i2=icv + i2
      i3=i2 + incr
      i4=i3
      sm=c(i2)*up
      do 90 i=l1,m
      sm=c(i3)*u(1,i) + sm
 90   i3=ice + i3
      if (sm) 100,120,100
 100  sm=b*sm
      c(i2)=sm*up + c(i2)
      do 110 i=l1,m
      c(i4)=sm*u(1,i) + c(i4)
 110  i4=ice + i4
 120  continue
 130  return
      end
c_______________________________________________________________________
