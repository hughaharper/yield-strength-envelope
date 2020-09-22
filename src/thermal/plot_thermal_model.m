clear;
% script to drive plate_cooling and plot the isotherms
age = 10; % Myr, max age of model
rate = 0.01; % rate in m / yr
% add in a flag for thermal model used
% 1 is plate model, 2 is sleep model
model=2;

switch model
    case 1
        runcmd=strcat('plate_cooling', 32, num2str(age), 32,'>', 32, 'out.temp');
    case 2
        runcmd=strcat('sleep_cooling', 32, num2str(age), 32, num2str(rate) , 32, '>', 32, 'out.temp');
    otherwise
        runcmd=strcat('plate_cooling', 32, num2str(age), 32,'>', 32, 'out.temp');
end

system(runcmd);

dat=load('out.temp');
ages = dat(2:end,1);
depths = (-1*dat(1,2:end))/1000;
temps = dat(2:end,2:end)';

% plot profiles...
figure();
plot(temps(:,1),depths); hold on;
plot(temps(:,20),depths);
plot(temps(:,50),depths);
xlabel('Temp (C)');
ylabel('Depth (km)');

figure();
[X,Z] = meshgrid(ages,depths);
contour(X(1:400,:),Z(1:400,:),temps(1:400,:),[100 300 600 800 1000 1185],'k-','Showtext','on'); hold on;
yline(-7,'k--','linewidth',0.75);
xlabel('Age (Ma)');
ylabel('Depth (km)');

plotfixer;