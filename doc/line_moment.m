%
%  demonstrate equivalence of inverse fourier method and the analytic solution
%
% set some constants
%
  clear
  rm=3300;
  rw=1000;
  g=9.8;
  E=6.e10;
  nu=0.25;
  h=30000;
  D=E*h^3/(12.*(1-nu^2));
  alph=(4*D/((rm-rw)*g))^.25;
  V0=1.;
  cnst=V0*alph^3/(8.*D);  % constant provided in T&S
  cnstm=cnst/alph;
%
%  set the x and the wavenumber
%
  N = 2048;             % number of points in the series should be a power of 2
  L = 3.e6;             % set the overall length of the profile
  dx= L/N;              % this is the data spacing
  x = dx*(1:N)-L/2;     % generate x-vector and put the origin arbitrarily at L/4.
  k = -(N/2):(N/2-1);   % generate the integer part pf the wavenumber vector 
  k = 2.*pi*k./L;       % convert from integer to radians per meter
  ks= ifftshift(k);     % move the zero wavenumber to the center of the array
%
%  place a 1 at the center of the array
%
   w=zeros(1,length(x));
   w(N/2)=1;
   cw=fft(w);
%
%  compute W(k)
%
  demon=D*(ks.^4)+(rm-rw)*g;
  cw=cw./demon/dx;            % divide by dx to get the units correct
  w1=real(ifft(cw))/cnst;      % shift the origin to the center of the array
                              % divide by the constant from T&S
%
% now compute analytic solution
%
  xa=abs(x)/alph;
  wm=exp(-xa).*(cos(xa)+sin(xa)); % don't multiply by the constant from T&S
%
%  plot the line load
%
  figure(1)
  clf
  subplot(2,1,1),plot(x/alph,w1,x/alph,wm)
  axis([-20,20,-.1,1.1])
  legend('analytic','FFT')
  xlabel('x/a')
  ylabel('deflection')
  subplot(2,1,2),plot(x/alph,wm-w1)
  axis([-20,20,-3.e-7,3.e-7])
  xlabel('x/a')
  ylabel('analytic - FFT')
%
%  now compute the response to a bending moment using two approaches
%
  w=zeros(1,length(x));
  w(N/2)=1;
  w(N/2+1)=-1;
  cw=fft(w);
  cw=cw./demon/dx;                 % divide by dx to get the units correct
  w2=real(ifft(cw))/(cnstm*2*dx) ;  % shift the origin to the center of the array
                                   % divide by the constant from T&S
%
% now compute analytic solution
%
  xa=(x-dx/2)/alph;
  wm=-exp(-abs(xa)).*sin(xa);   % don't multiply by the constant from T&S
  figure(2)
  clf
  subplot(2,1,1),plot(x/alph,w2,x/alph,wm)
  axis([-20,20,-.4,.4])
  legend('analytic','FFT')
  xlabel('x/a')
  ylabel('deflection')
  subplot(2,1,2),plot(x/alph,wm-w2)
  axis([-20,20,-5.e-5,5.e-5])
  xlabel('x/a')
  ylabel('analytic - FFT')
%
%  now make a trench simulation
%
  w=zeros(1,length(x));
  w(N/2)=1;
  w(N/2+1)=-1;
  w=w/(2*dx/alph);                 % need to scale the moment by 2*dx/alpha
  w(N/2)=w(N/2)-1.;
  cw=fft(w)/cnst;
  cw=cw./demon/dx;                 % divide by dx to get the units correct
  w3=real(ifft(cw));
  figure(3)
  clf
  subplot(2,1,1),plot(x/alph,w3)
  axis([1,10,-1.,.2])
  xlabel('x/a')
  ylabel('deflection')
