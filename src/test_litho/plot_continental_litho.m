% script to drive plot_ocean_litho
clear;

file = 'out.temp';

runcmd=strcat('cont_litho_yse',32,'>', 32,file);

system(runcmd);

dat=load(file);
z=-dat(:,1);
temp=dat(:,2);
pres=dat(:,3);
ystrp=dat(:,4);
ystrm=dat(:,5);

figure()
plot(ystrp,z,'linewidth',1.5); hold on;
plot(ystrm,z,'linewidth',1.5);
xline(0,'k-','linewidth',1);
yline(-40,'k--','linewidth',0.5);

xlabel('Differential Stress (MPa)');
ylabel('Depth (km)');
ylim([-100 0]);
title_string = sprintf('Yield Stress Profile, Age = 100 Ma, T_{lith} = 125000 km');
title(title_string);

plotfixer;
