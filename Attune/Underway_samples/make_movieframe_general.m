% This function is to be called by the moviemaker for each fcs file
% it generates 11 scatter plots of useful channel data. using labels
% taken directly from

% This version based off the plotting options in assign_class_EN649

%For now, quality control thresholds are still calculated from fcsdata, at
%some point we may want to make sure that that is only calculated one
%place in the pipeline


function [Fig] = make_movieframe_general(fcsdat, fcshdr, class, ch_names, QC_flag)

% Trying to streamline quality control to only have it once
% %First, Quality Control Flags
% %check flowrate data to flag plots for low quality
% t = find(fcsdat(:,12)>200 & fcsdat(:,3)>200);
% QC_flowrate(1) = (median(fcsdat(t,3)./fcsdat(t,12)));
% QC_flowrate(2) = (std(fcsdat(t,3)./fcsdat(t,12)));
% 
% QC_flag = 0; %default bad
% if (QC_flowrate(2)<2 & QC_flowrate(1)<1.5)
%     QC_flag = 1; %set to good
% end

%trying to keep figures in the background
set(0, 'DefaultFigureVisible', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');

%we're going to call subplot_tight so that margins aren't so huge
%need some parameters
xsp=0.04; ysp=0.1;

%parameter numbers for main polygons
npar_eukX = strmatch([ch_names{1}], {fcshdr.par.name}); %BL3-H
npar_eukY = strmatch([ch_names{2}], {fcshdr.par.name}); %GL2-H %changed all 19-16 & 10-7 indices down 1 compared to TN368
npar_synX = strmatch([ch_names{3}], {fcshdr.par.name}); %FSC-H
npar_synY = strmatch([ch_names{4}], {fcshdr.par.name}); %GL2-H

%preset subplot locations 
Fig = figure(1);
clf
%set(Fig, 'Position', [115 57 1e+03 482], 'Visible', 'off')
set(0,'CurrentFigure',Fig)
a1 = 3; a2 = 5;

ax1 = subplot_tight(a1,a2,1,[ysp xsp]);

ax2 = subplot_tight(a1,a2,2,[ysp xsp]);
ax3 = subplot_tight(a1,a2,3, [ysp xsp]);
ax4 = subplot_tight(a1,a2,4, [ysp xsp]);
ax5 = subplot_tight(a1,a2,5, [ysp xsp]);
ax6 = subplot_tight(a1,a2,6, [ysp xsp]);
ax7 = subplot_tight(a1,a2,7, [ysp xsp]);
ax8 = subplot_tight(a1,a2,8,[ysp xsp]);
ax9 = subplot_tight(a1,a2,9,[ysp xsp]);
ax10 = subplot_tight(a1,a2,10,[ysp xsp]);
ax11 = subplot_tight(a1,a2,11,[ysp xsp]);

%axis limits
mylim1 = [1e1 1.1e6 1e1 1.1e6];
mylim2 = [1e3 1.1e6 1e3 1.1e6];

%plot plot plot
make_plot(ax1,npar_eukX,npar_eukY,mylim1)
title(ax1, fcshdr.filename,'fontsize', 8, 'interpreter', 'none')
make_plot(ax2,npar_synX,npar_synY,mylim1)
make_plot(ax3,npar_synY,9,mylim1)
title(ax3, [fcshdr.date, ' ', fcshdr.starttime], 'interpreter', 'none', 'fontsize', 10)
make_plot(ax4,12,npar_synY,mylim1)
make_plot(ax5,16,15,mylim1) %BL3H v GL3H
if QC_flag == 0
    title(ax5, 'POOR QUALITY flag tripped', 'color', 'r', 'fontweight', 'bold')
end
make_plot(ax6,12,15,mylim1)
make_plot(ax7,11,15, mylim1)
make_plot(ax8,16,7, mylim1)
make_plot(ax9,npar_eukY-1,npar_eukY, mylim2)
make_plot(ax10,11,12, mylim1)
make_plot(ax11,12,14, mylim1)
make_legend

        Fig = getframe(Fig); 
        
    %this function called for the individual subplots
    function make_plot(ax, nparX, nparY, mylim)
        loglog(ax,fcsdat(:,nparX), fcsdat(:,nparY), '.', 'markersize', 1)
        ylabel(ax,fcshdr.par(nparY).name)
        xlabel(ax,fcshdr.par(nparX).name)
        %axis(ax,'square')
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
        ax.XLabel.FontSize = 6;  
        ax.YLabel.FontSize = 8;
       
    end



    function make_legend
            leglabels = {'Noise'; 'Euks'; 'Syn'; 'LowP Euk'; 'HiP Euk'; 'SynEuk Coincident'; 'Noise'}; 
            ind = [0 (sum(class==[1:6])==0)]; %find classes with no counts to remove from legend
            leglabels(find(ind)) = []; 
            [hleg, icns] =legend(leglabels, 'Orientation', 'Horizontal'); 
            set(hleg,'Position',[0.35    .2    0.5526    0.0268])
            p = find(isgraphics(icns, "Line")); 
        for pp = 1:length(p)
            icns(p(pp)).MarkerSize = 18; 
        end
     end
        

end
