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
ystrp=dat(:,4);
ystrm=dat(:,5);
estr=dat(:,6);


figure()

% subplot(2,2,1)
% plot(temp,z);
% xlabel('Temperature (deg C)');
% ylabel('Depth (km)');
% title('Temperature Profile');
% 
% subplot(2,2,2)
% plot(pres,z);
% xlabel('Overburden Pressure (MPa)');
% ylabel('Depth (km)');
% title('Pressure Profile');

% subplot(2,2,3:4)
plot(ystrp,z,'linewidth',2); hold on;
plot(ystrm,z,'linewidth',1);
plot(estr,z,'linewidth',1);
xlabel('Differential Stress (MPa)');
ylabel('Depth (km)');
ylim([-80 0]);
title_string = sprintf('Yield Stress Profile, Age = %s, T_{lith} = 125500 km',age);
title(title_string);
legend('Tension','Compression');

