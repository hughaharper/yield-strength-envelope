% MATLAB script to test code 

age = input('Enter plate age in MYr: ','s');
file = input('Enter output file name: ','s');

runcmd=strcat('mocurv', 32, age, 32, '>', 32,file);

system(runcmd);

dat=load(file);
curv=-dat(:,1);
moment=dat(:,2);



figure(1)
plot(curv,moment);
grid;
xlabel('Curvature (10^{-6} m^{-1})');
ylabel('Moment (10^{17} N)');
xlim([-10 10]);
ylim([-4.0 6.0]);
title('Moment-Curvature Relationship (No Water Column Overburden)');
