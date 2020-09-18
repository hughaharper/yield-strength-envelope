clear; close all;
% script to drive plate_cooling and plot the isotherms
age = 20; % Myr, max age of model

runcmd=strcat('plate_cooling', 32, num2str(age), 32,'>', 32, 'out.temp');
system(runcmd);

dat=load('out.temp');
ages = dat(2:end,1);
depths = (-1*dat(1,2:end))/1000;
temps = dat(2:end,2:end)';

% plot profiles...
figure(1);
plot(temps(:,1),depths); hold on;
plot(temps(:,20),depths);
plot(temps(:,50),depths);
xlabel('Temp (C)');
ylabel('Depth (km)');

figure(2);
[X,Z] = meshgrid(ages,depths);
contour(X,Z,temps,[250,500,750,1000,1250],'k-','Showtext','on'); hold on;
yline(-7,'k--','linewidth',0.75);
xlabel('Age (Ma)');
ylabel('Depth (km)');

plotfixer;