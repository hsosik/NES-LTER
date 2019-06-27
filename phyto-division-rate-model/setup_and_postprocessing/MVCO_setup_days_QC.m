% A simple script that creates a matlab movie form (results in a structure)
% to screen output from setup_days_all.m for model input days

filelist=dir([modelinputpath 'day*data.mat']);

figure, set(gcf,'Position',  [ 82         428        1524         521])
count=0;

for j=1:length(filelist)
    
    count=count+1;
    
    filename=filelist(j).name;
    eval(['load ' modelinputpath filename])
    
    clf
    subplot(1,3,2,'replace')
    imagesc([1 25],[1 57],Vhists)
    set(gca,'YDir','normal');
    %pcolor(Vhists), hading flat
    caxis([0 0.10])
    colorbar
    xlabel('Hours After Dawn','fontsize',14)
    ylabel('Cell Volume','fontsize',14)
    title({['Day: ' filename(4:9) ' : ' datestr(str2num(filename(4:9)))]; 'Vhists'},'fontsize',14)
    
    subplot(133)
    imagesc([1 25],[1 57],N_dist)
    set(gca,'YDir','normal');
    %         pcolor(N_dist)
    %         shading flat
    colorbar %caxis([0 10000])
    xlabel('Hours After Dawn','fontsize',14)
    ylabel('Cell Volume','fontsize',14)
    t1=title('N_dist','fontsize',14);
    set(t1,'interpreter','none')
    
    subplot(131)
    plot(Edata(:,1), Edata(:,2),'.-','color',[0 0.3 0.8],'markersize',10,'linewidth',1.5)
    xlabel('Hours After Dawn','fontsize',14)
    ylabel('Radiation','fontsize',14)
    ylim([0 1000])
    
    F1(count) = getframe(gcf);
end

notes='to play, type implay(syn_dist20XX) in MATLAB window';
eval(['syn_dist' num2str(year2do) '=F1;'])
eval(['save ' modelinputpath 'setup_days' num2str(year2do) '.mat syn_dist' num2str(year2do) ' notes'])

close all
clear F1
