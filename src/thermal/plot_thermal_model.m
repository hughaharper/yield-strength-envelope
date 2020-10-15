clear;
% script to drive plate_cooling and plot the isotherms
age = 10; % Myr, max age of model
rate = 0.015; % rate in m / yr
% add in a flag for thermal model used
% 1 is plate model, 2 is sleep model
model=2;

switch model
    case 1
        runcmd=strcat('plate_cooling', 32, num2str(age), 32,'>', 32, 'out.temp');
    case 2
        runcmd=strcat('sleep_cooling', 32, num2str(rate) , 32, '>', 32, 'out.temp');
    case 3
        runcmd=strcat('sleep_modified', 32, num2str(rate) , 32, '>', 32, 'out.temp');
    otherwise
        runcmd=strcat('plate_cooling', 32, num2str(age), 32,'>', 32, 'out.temp');
end

system(runcmd);

dat=load('out.temp');
dist = dat(2:end,1)/1000;
depths = (-1*dat(1,2:end))/1000;
temps = dat(2:end,2:end)';

% plot profiles...
figure();
plot(temps(:,1),depths); hold on;
plot(temps(:,20),depths);
plot(temps(:,end),depths);
ylim([-30 0]);
xlabel('Temp (C)');
ylabel('Depth (km)');

figure();
[X,Z] = meshgrid(dist,depths);
contour(X(:,:),Z(:,:),temps(:,:),[0 100 300 600 800 1000 1185],'k-','Showtext','on'); hold on;
yline(-7,'k--','linewidth',0.75);
xlabel('Distance (km)');
ylabel('Depth (km)');

plotfixer;