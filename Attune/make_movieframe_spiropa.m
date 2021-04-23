% This function is to be called by the moviemaker for each fcs file
% it generates 11 scatter plots of useful channel data. using labels
% taken directly from

% This version based off the plotting options in
% assign_class_spiropaRB1904_v2.m

%For now, quality control thresholds are still calculated from fcsdata, at
%some point we may want to make sure that that is only calculated one
%place in the pipeline


function [] = make_movieframe_spiropa(fcsdat, fcshdr, class)

%First, Quality Control Flags
%check flowrate data to flag plots for low quality
t = find(fcsdat(:,12)>200 & fcsdat(:,3)>200);
QC_flowrate(1) = (median(fcsdat(t,3)./fcsdat(t,12)));
QC_flowrate(2) = (std(fcsdat(t,3)./fcsdat(t,12)));

QC_flag = 0; %default bad
if (QC_flowrate(2)<2 & QC_flowrate(1)<1.5)
    QC_flag = 1; %set to good
end

%parameter numbers for main polygons
npar_eukX = 15; %BL3-H
npar_eukY = 19; %GL1-H
npar_synX = 11; %FSC-H
npar_synY = 19; %GL1-H

%preset subplot locations 
figure(1), clf
a1 = 3; a2 = 5;
ax1 = subplot(a1,a2,1);
ax2 = subplot(a1,a2,2);
ax3 = subplot(a1,a2,3);
ax4 = subplot(a1,a2,4);
ax5 = subplot(a1,a2,5);
ax6 = subplot(a1,a2,6);
ax7 = subplot(a1,a2,7);
ax8 = subplot(a1,a2,8);
ax9 = subplot(a1,a2,9);
ax10 = subplot(a1,a2,10);
ax11 = subplot(a1,a2,11);

%axis limits
mylim1 = [1e1 1.1e6 1e1 1.1e6];
mylim2 = [1e3 1.1e6 1e3 1.1e6];

%plot plot plot
make_plot(ax1,npar_eukX,npar_eukY,mylim1)
title(ax1, fcshdr.filename, 'interpreter', 'none')
make_plot(ax2,npar_synX,npar_synY,mylim1)
make_plot(ax3,19,10,mylim1)
title(ax3, [fcshdr.date, ' ', fcshdr.starttime], 'interpreter', 'none')
make_plot(ax4,12,19,mylim1)
make_plot(ax5,17,15,mylim1) %BL3H v GL3H
if QC_flag == 0
    title(ax5, 'POOR QUALITY flag tripped', 'color', 'r', 'fontweight', 'bold')
end
make_plot(ax6,12,15,mylim1)
make_plot(ax7,11,15, mylim1)
make_plot(ax8,17,8, mylim1)
make_plot(ax9,18,19, mylim2)
make_plot(ax10,11,12, mylim1)
make_plot(ax11,12,14, mylim1)

%make legend
[hleg, icns] =legend('Junk', 'Euks', 'Syn', 'LowP Euk', 'HiP Euk', 'SynEuk Coincident', 'SynEuk Coincident2', 'Noise', 'Orientation','Horizontal');
set(hleg,'Position',[0.35    .2    0.5526    0.0268])
for p = 10:2:24
    icns(p).MarkerSize = 18;
end

    %this function called for the individual subplots
    function make_plot(ax, nparX, nparY, mylim)
        loglog(ax,fcsdat(:,nparX), fcsdat(:,nparY), '.', 'markersize', 1)
        ylabel(ax,fcshdr.par(nparY).name)
        xlabel(ax,fcshdr.par(nparX).name)
        axis(ax,'square')
        hold(ax, 'on')
        loglog(ax,fcsdat(class==1,nparX), fcsdat(class==1,nparY), '.g', 'markersize', 1)
        loglog(ax,fcsdat(class==2,nparX), fcsdat(class==2,nparY), 'r.', 'markersize', 1)
        loglog(ax,fcsdat(class==3,nparX), fcsdat(class==3,nparY), 'm.',  'markersize', 1)
        loglog(ax,fcsdat(class==4,nparX), fcsdat(class==4,nparY), 'k.', 'markersize', 1)
        loglog(ax,fcsdat(class==5,nparX), fcsdat(class==5,nparY), 'c.', 'markersize', 2)   
        loglog(ax,fcsdat(class==6,nparX), fcsdat(class==6,nparY), '.y', 'markersize', 1.5)
        loglog(ax,fcsdat(class==7,nparX), fcsdat(class==7,nparY), '.', 'color', [.4 .4 .4], 'markersize', 1)
        
        axis(ax,mylim)
        set(ax, 'xtick', [1e2 1e4 1e6])
        set(ax, 'ytick', [1e2 1e4 1e6])
    end

end
