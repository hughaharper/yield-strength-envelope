% script to drive plot_ocean_litho
clear; close all;

file = 'out.temp';

ages = [0.5 1 1.5 2 2.5 3 4];
n = length(ages);

for i=1:n

    runcmd=strcat('ocean_litho_yse', 32, num2str(ages(i)), 32, '>', 32, file);

    system(runcmd);

    dat=load(file);
    z=-dat(:,1);
    temp=dat(:,2);
    pres=dat(:,3);
    ystrp=dat(:,4);
    ystrm=dat(:,5);

    figure(1)
    subplot(2,n,i);
    plot(ystrp,z,'linewidth',1.5); hold on;
    plot(ystrm,z,'linewidth',1.5);
    xline(0,'k-','linewidth',1);
    yline(-6,'k--','linewidth',0.5);
    
    if i==1
        xlabel('Differential Stress (MPa)');
        ylabel('Depth (km)');
    end
    ylim([-20 0]);
    xlim([-1000 1000]);
    title_string = sprintf('%.1f Myr',ages(i));
    title(title_string);
    
    subplot(2,n,[n+1:n+2]);
    plot(temp,z); hold on;
    
    if i==1
        xline(1350,'k-','linewidth',1);
        yline(-6,'k--','linewidth',0.5);
        xlabel('Temperature (C)');
        ylabel('Depth (km)');
    end
    
    strength(i) = trapz(-z*1e3,(ystrp - ystrm)*1e6);
end

subplot(2,n,[n+4:2*n]);
semilogy(ages,strength,'ko-','linewidth',1.5,'markersize',0.5);
ylim([1e11 1e13]);
xlabel('Age (Myr)');
ylabel('Strength N m^{-1}');

plotfixer;
