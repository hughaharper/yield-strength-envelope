% script to drive plot_ocean_litho
clear; close all;
%age = input('Enter plate age in MYr: ','s');
%file = input('Enter output file name: ','s');

% loop through some ages
%age = [144 100 64 36 16 4];
age = [8 4 2 1 0.5 0.1];

for ii = 1:length(age)
    
    runcmd=strcat('ocean_litho_yse', 32, num2str(age(ii)), 32,'>', 32, 'out.temp');

    system(runcmd);

    dat=load('out.temp');
    z=-dat(:,1);
    temp=dat(:,2);
    pres=dat(:,3);
    ystrp=dat(:,4);
    ystrm=dat(:,5);

    figure(1)
    subplot(1,length(age),ii);
    plot(ystrp,z,'linewidth',1.5); hold on;
    plot(ystrm,z,'linewidth',1.5);
    
    xline(0,'k-','linewidth',1);
    
    if ii == 1
        xlabel('Differential Stress (MPa)');
        ylabel('Depth (km)');
    end
    
    ylim([-30 0]);
    xlim([-1000 1000]);
    title_string = sprintf('%.1f Ma ',age(ii));
    title(title_string);
    
    % integrate the YSE through depth to get a total strength in N/m
    strength(ii) = trapz(-z*1e3,(ystrp - ystrm)*1e6);
    
end %for ages
figure(2)
semilogy(age,strength,'ko-','linewidth',1.5);
ylim([1e11 1e13]);
xlabel('Age (Ma)'); ylabel('Strength N m^{-1}');

plotfixer;