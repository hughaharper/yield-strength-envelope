% MATLAB script to test code 

age = input('Enter plate age in MYr: ','s');
curv = input('Enter curvature in units of 10^-7 m^-1: ');
file = input('Enter output file name: ','s');

curvin = num2str(curv*1.e-7);


runcmd=strcat('test_litho', 32, age, 32, curvin, 32,'>', 32,file);

system(runcmd);

dat=load(file);
z=-dat(:,1);
temp=dat(:,2);
pres=dat(:,3);
estr=dat(:,4);
ystr=dat(:,5);
dstr=dat(:,6);


figure(1)

subplot(4,2,1)
plot(temp,z);
xlabel('Temperature (deg C)');
ylabel('Depth (km)');
title('Temperature Profile');

subplot(4,2,2)
plot(pres,z);
xlabel('Overburden Pressure (MPa)');
ylabel('Depth (km)');
title('Pressure Profile');

subplot(4,2,3:4)
plot(estr,z);
xlabel('Differential Stress (MPa)');
ylabel('Depth (km)');
ylim([-80 0]);
title('Yield Stress Profile');


subplot(4,2,5:6)
plot(ystr,z);
xlabel('Differential Stress (MPa)');
ylabel('Depth (km)');
ylim([-80 0]);
title('Elastic Stress Profile');

subplot(4,2,7:8)
plot(dstr,z);
xlabel('Differential Stress (MPa)');
ylabel('Depth (km)');
ylim([-80 0]);
title('Differential Stress Profile');