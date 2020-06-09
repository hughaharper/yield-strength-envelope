      program fitflexnp

c=======================================================================
c
c     Program fitflex.f to invert for 
c	  at the trench axis by calling nnls
c
c     This fortran program reads the input and G-matrix and 
c     (possibly) constraints into arrays, then calls nnls to solve for 
c	  applied moment and vertical load.
c
c     This program solve for the uncertainties by simulating noise 
c     based on the data variance. $numsim times of realization gives 
c     the uncertainties of the model parameters.
c
c     Input: 
c		- marine gravity anomaly values from satellite altimetry (mGal)       
c		- multibeam bathymetry (m)
c    	both data sets can possibly we weighted according to some scheme
c		perhaps distance from closest point on trench axis
c
c     G-matrix:
c       predicted deflection sampled at the same data points as 
c		the gravity and bathymetry
c       
c     Constraints:
c       may possibly be added later
c
c     Output:
c       slip rate on individual fault segments with uncertainties
c     
c     Date: September 2014
c	  Code originally written by X. Tong to find fault slip rates 
c	  modified by E. S. Garcia to solve for flexural loading
c
c=======================================================================


c     Define variables
      implicit none
      
      integer numbt, numgr, numts, nc, nc1, nc2, m, n
      integer i, j, k, p, q, numsimu, narg, mode
      integer lwork, info
      integer rank
c      
c     geol + constraints 
c      parameter (numbt=0,numgr=0,numts=56,nc1=30,nc2=10,
c     + m=96, n=56, nc=40) 
c     m = nc1 + nc2 + numts
c
c  CHANGE THIS PART AS NEEDED
c
c	  since there are no constraints or nuisance parameters,
c	  m=numbt+numgr	 
c	  n=numts*2
c
      parameter (numbt=183955,numgr=183955,numts=11,nc1=30,
     + nc2=9, m=367910, n=22, nc=39)

c     matrix A is m row n column
c     n is unknowns, m is knowns, 
c     Warning: I assume m>n, it's overdetermined problem
     
      parameter (lwork=n+m*64) 

c	  *** turn off random realizations for now
      parameter (numsimu=10)
c
      double precision sl, wtopo, wgrav, wc
      double precision xg, yg, corr, xs, ys, ul, vl, zl
      
      double precision xbb(numbt), ybb(numbt)
      double precision btdat(numbt), sdbtdat(numbt), btcorr(numbt)
      double precision bgm(numbt,numts), bgv(numbt,numts)
      
      double precision xgr(numgr), ygr(numgr)
      double precision grdat(numgr), sdgrdat(numgr), grcorr(numgr)
      double precision ggm(numgr,numts), ggv(numgr,numts)
      
      double precision cmatr(nc,numts), losgs(numgr,numts) 
      double precision amatr(m,n),bvec(m),xvec(n),wvec(n),rvec(m)
      double precision zvec(m), invec(n), atmp(m,n)
      double precision rnorm
      double precision DNRM2
      double precision fs(numts), sfs(numts)
      double precision A(m,n),B(m,1),work(lwork)
      double precision S(n)

      double precision u1, u2, nlos, nug, nvg, tmps
      double precision sr(n, numsimu), tmpsum(n), fsr(n), fssr(n) 
      double precision tl(n)
      double precision ptopo(numbt),pgrav(numgr)
      double precision sumgr,sumtp,chigrav,chibath,chiboth
      double precision wrmsbath,wrmsgrav,sumsdtp,sumsdgr
      double precision shiftconst      

      integer JPVT(n)
      character*20 inbath, inbathfnct, grav, gravfnct, csl, cmat
      character*20 topofs, sarfs, cwbath, cwgrav, cwc
      
c     parse the command line
      narg = iargc()
      if (narg .ne. 6) then 
        write(*,*) ' '
        write(*,*)
     +  'Usage: fitflex bath.dat bath.fnct weight.bath 
     +                  grav.dat grav.fnct weight.grav' 
        write(*,*) ' '
        write(*,*) 
c		corr might be correction due to seafloor subsidence
     +  ' bath.dat - x_pos (km) y_pos (km) bath (m) bath_stdev (m) corr'
        write(*,*) ' '
        write(*,*)
     +  ' bath.fnct - a list of the green function file for bath data'
        write(*,*) ' '
        write(*,*)
     +  ' green function file for bath - x_pos(km)  y_pos(km) bath (m)'
        write(*,*) ' '
        write(*,*)
     +  ' weight.bath - weight for bath data set'
        write(*,*) ' '
        write(*,*)
     +  ' grav.dat - x_pos (km)  y_pos (km)  grav (mGal) grav_stdev (mGal)
     +               corr'
        write(*,*) ' '
        write(*,*)
     +  ' grav.fnct - a list of the green function file for grav data'
        write(*,*) ' '
        write(*,*)
     +  ' green function file for grav - x_pos(km) y_pos(km) 
     +                                   grav (mGal)'
        write(*,*) ' '
        write(*,*) 
     +  ' weight.grav - weight for grav data set'
        write(*,*) ' '
        stop
      else
        call getarg(1, inbath)
        call getarg(2, inbathfnct)
        call getarg(3, cwbath)
        call getarg(4, grav)
        call getarg(5, gravfnct)
        call getarg(6, cwgrav)
        
        read(cwbath, *) wtopo
        read(cwgrav, *) wgrav
      endif
      
      print *, ' '
      print *, 'input: '
      print *, inbath
      print *, inbathfnct
      print *, grav
      print *, gravfnct
      print *, 'topo weight: ', wtopo
      print *, 'grav weight: ', wgrav
      print *, ' '
 
c     Initialize
      do 1 j=1, m
        do 2 i=1, n
          amatr(j,i)=0
          tmpsum(i)=0
2       enddo
        bvec(j)=0
1     enddo
   
      do 3 j=1, numbt
        ptopo(j)=0

3     enddo
      do 4 j=1, numgr
        pgrav(j)=0
4     enddo

      sumgr=0.
      sumsdgr=0.
      sumtp=0.
      sumsdtp=0. 
      shiftconst=100
     
c     Read in the topo data points
      open(unit=5, file=inbath, status='old')
      do 10 j = 1, numbt
        read(5,*) xbb(j), ybb(j), btdat(j), sdbtdat(j), btcorr(j)
10    enddo
      close(5)

c     Read in the G-matrix for the bath data
      open(unit=15, file=inbathfnct, status='old')
      do 20 j = 1, numts
        read(15,*) topofs
        open(unit=25, file=topofs, status='old')
        
        do 30 i = 1, numbt
          read(25,*) xbb(i), ybb(i), bgm(i,j), bgv(i,j)
30      enddo
        
        close(25)
20    enddo
      close(15)

c     Read in the grav data points
      open(unit=35, file=grav, status='old')
      do 40 i = 1, numgr
        read(35,*) xgr(i), ygr(i), grdat(i), sdgrdat(i), grcorr(i)
40    enddo
      close(35)

c     Read in the G-matrix for the grav data
      open(unit=45, file=gravfnct, status='old')
      do 50 j = 1, numts
        read(45,*) sarfs
        open(unit=55, file=sarfs, status='old')

        do 60 i = 1, numgr
          read(55,*) xgr(i), ygr(i), ggm(i,j), ggv(i,j)
60      enddo

        close(55)
50    enddo
      close(45)

      
c     Form the inversion matrix
c     Form the A-matrix 

c     Use the topo data
      do 80 i=1, numbt
        do 90 j=1, numts
c     Weighting by sigma
          amatr(i,j)=bgm(i,j)/(sdbtdat(i)*wtopo)
          amatr(i,j+numts)=bgv(i,j)/(sdbtdat(i)*wtopo)
90      enddo
80    enddo

c     Use the grav data
      do 100 i=1, numgr
        do 110 j=1, numts
          amatr(i+numbt,j)=ggm(i,j)/(sdgrdat(i)*wgrav)
          amatr(i+numbt,j+numts)=ggv(i,j)/(sdgrdat(i)*wgrav)
110     enddo
100   enddo

c     Form the b-vector and add random noise to it
c     Start loop here
c	  *** Turn off 400 loop for now
c     *** do 400 p=1, numsimu 

c     Use the topo data
        do 300 i=1, numbt

c     generate random numbers with uniform distribution
c     then transform it into normal distribution
c     Box-Muller transform    
c         u1 = rand()
c         u2 = rand()
           
c         print *, u1, u2
c         nug=sug(j)*sqrt(-2.0*log(u1))*cos(2*3.1416*u2)
c         nvg=svg(j)*sqrt(-2.0*log(u1))*sin(2*3.1416*u2)
c     Weighting by sigma

          bvec(i)=btdat(i)/(sdbtdat(i)*wtopo)
          
300     enddo

c     Use the grav data
        do 310 i=1, numgr

c          u1=rand()
c          u2=rand()
c          print *, u1, u2
c          nlos=slos(j)*sqrt(-2.0*log(u1))*cos(2*3.1416*u2)

          bvec(i+numbt)=grdat(i)/(sdgrdat(i)*wgrav)
          
310     enddo

c     print input parameters
        print *, 'paramters: '
        print *, 'numgr=', numgr
        print *, 'numbt=', numbt
        print *, 'm=', m
        print *, 'numts=', numts
        print *, 'n=', n
        print *, ' '
        print *, 'start - least square problem' 
      
c     copy the amatr and bvec to A and B
c     as they will be destroyed

        do 250 i=1, m
          do 240 j=1, n
            A(i,j)=amatr(i,j)
c            atmp(i,j)=amatr(i,j)
c            print *, A(i,j)
240       enddo
          B(i,1)=bvec(i)
c          rvec(i)=bvec(i)
c        print *, B(i,1)
250     enddo

c     Call nnls
        print *, 'calling nnls'
        call nnls(A,m,m,n,B,xvec,rnorm,wvec,zvec,
     + invec, mode)
       
        print *, 'end - nnls'
        print *, ' '

c     Print output
        print *, 'print output...'
        print *, 'mode=', mode
        print *, 'rnorm=', rnorm
        print *, ' '
      
c     store the model parameters in an array
        do 410 q=1, n
          tl(q)=xvec(q)
410     enddo

     
c     compute the predicted topo data 
c	  *** something clumsy about the following loops:
      do 500 i=1, numbt
        do 510 j=1, numts
          ptopo(i)=bgm(i,j)*tl(j)+ptopo(i)
          ptopo(i)=bgv(i,j)*tl(j+numts)+ptopo(i)
510     enddo  
500   enddo

c     compute the predicted grav data
c	  *** something clumsy about the following loops:
      do 520 i=1, numgr
        do 530 j=1, numts
          pgrav(i)=ggm(i,j)*tl(j)+pgrav(i)
          pgrav(i)=ggv(i,j)*tl(j+numts)+pgrav(i)
530     enddo
520   enddo

c     compute the chi square misfit for the topo
      do 540 i=1,numbt
        sumtp=sumtp+((ptopo(i)-btdat(i))/sdbtdat(i))**2
	    sumsdtp=sumsdtp+(1.0/sdbtdat(i))**2
540   enddo
      chibath=sumtp/(numbt-numts)
      wrmsbath=sqrt(sumtp/sumsdtp)

c     compute the chi square misfit for the grav
      do 550 i=1,numgr
        sumgr=sumgr+((pgrav(i)-grdat(i))/sdgrdat(i))**2
		sumsdgr=sumsdgr+(1.0/sdgrdat(i))**2
550   enddo
      chigrav=sumgr/(numgr-numts)
      wrmsgrav=sqrt(sumgr/sumsdgr)
     
      chiboth=(sumgr+sumtp)/(numgr+numbt-numts)     

      print *,'Chi Square misfit for topo=', chibath 
      print *,'Weighted RMS misfit for topo=', wrmsbath
      print *,'Chi Square misfit for grav=', chigrav
      print *,'Weighted RMS misfit for grav=', wrmsgrav
      print *,'Chi Square misfit for both=', chiboth 

c     write the topo into a file
      open(unit=95, file='pred.topo.dat')
      do 560 i=1, numbt
        write(95,19) btdat(i),ptopo(i),btdat(i)-ptopo(i),sdbtdat(i)
560   enddo
      close(95)
19    format(8(F7.2))

c     write the gravity into a file
      open(unit=105, file='pred.grav.dat')
      do 570 i=1, numgr
        write(105,29) grdat(i),pgrav(i),grdat(i)-pgrav(i),sdgrdat(i)
570   enddo
      close(105)
29    format(4(F7.2))

c     write the model vector into a file
	  write(*,*) 'Writing out model vector, should have lines: ',n
      open(unit=85, file='lq.dat')
      do 140 i=1, n
        write(85,*) i, tl(i)
140   enddo
      close(85)
      
88    stop
      end
