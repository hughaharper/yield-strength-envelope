nz100=load('mocurv_100Myr_nz100.txt'); 
nz1000=load('mocurv_100Myr_nz1000.txt');
tol1km=load('mocurv_100Myr_tol1km.txt');

nz100c=nz100(:,1); nz100m=nz100(:,2);
nz1000c=nz1000(:,1); nz1000m=nz1000(:,2);
tol1kmc=tol1km(:,1); tol1kmm=tol1km(:,2);

figure(1)
plot(nz100c,nz100m,nz1000c,nz1000m,tol1kmc,tol1kmm);
grid()
xlim([-10 10]);
ylim([-4 6]);