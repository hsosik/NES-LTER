% This function is to be called by the moviemaker for each fcs file
% it generates 11 scatter plots of useful channel data. using labels
% taken directly from

% This version based off the plotting options in assign_class_EN649

%For now, quality control thresholds are still calculated from fcsdata, at
%some point we may want to make sure that that is only calculated one
%place in the pipeline


function [Fig] = make_movieframe_density(fcsdat, fcshdr, class, ch_names, QC_flag, bounds)

% Trying to streamline quality control to only have it once
%that's why it's an input

%trying to keep figures in the background
%set(0, 'DefaultFigureVisible', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');

%we're going to call subplot_tight so that margins aren't so huge
%need some parameters
xsp=0.04; ysp=0.1;


if strcmp(ch_names, 'early') 
    ch_names = {'BL3-H', 'GL1-H', 'SSC-H', 'GL1-H'}; 
    pe_chs = [10 19 28]; %PE-A, PE-H, PE-W
elseif strcmp(ch_names, 'late'); 
    ch_names = {'BL3-H', 'GL2-H', 'SSC-H', 'GL2-H'};
    pe_chs = [9 18 27]; %PE-A, PE-H, PE-W
end


%parameter numbers for main polygons
npar_eukX = strmatch([ch_names{1}], {fcshdr.par.name}); %BL3-H
npar_eukY = strmatch([ch_names{2}], {fcshdr.par.name}); %PE-H 
npar_synX = strmatch([ch_names{3}], {fcshdr.par.name}); %SSC-H
npar_synY = strmatch([ch_names{4}], {fcshdr.par.name}); %PE-H

if npar_synY ~= pe_chs(2)
    keyboard
end


if exist('bounds', 'var')
        geuk_main_gate = bounds{1};
        gsyn_main_gate = bounds{2}; 
        synGL1H2BL3Hslope = bounds{3} ;
        synGL1H2BL3Hoffset = bounds{4} ;
        synPEA2PEHmax = bounds{5} ; 
        syneukBL3H2SSCHslope = bounds{6}; 
        syneukBL3H2SSCHoffset = bounds{7} ;
        nonsynfactorA = bounds{8} ; 
        nonsynfactorB = bounds{9}; 
end


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
mylim2 = [1e1 1.1e5 1e1 1.1e6];

%plot plot plot EP: first subplot
make_plot(ax1,npar_eukX,npar_eukY,mylim1)
title(ax1, fcshdr.filename,'fontsize', 8, 'interpreter', 'none')
if exist('bounds', 'var')
    plot(ax1,geuk_main_gate(:,1), geuk_main_gate(:,2), 'g.-')
    f = @(x) 10.^(log10(x).*synGL1H2BL3Hslope+synGL1H2BL3Hoffset);
    fplot(ax1, f, xlim(ax1))
end
make_plot(ax2,npar_synX,npar_synY,mylim1)
if exist('bounds', 'var')
    plot(ax2,gsyn_main_gate(:,1)+1e-5,gsyn_main_gate(:,2), 'r.-')
end
make_plot(ax3,pe_chs(2),pe_chs(1),mylim1)

if exist('bounds', 'var')
    line(ax3,xlim(ax3),xlim(ax3)*synPEA2PEHmax)
end

title(ax3, [fcshdr.date, ' ', fcshdr.starttime], 'interpreter', 'none', 'fontsize', 10)
make_plot(ax4,19,npar_synY,mylim1)
make_plot(ax5,12,15,mylim1) 
if exist('bounds', 'var')
    f = @(x) 10.^(log10(x).*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset);
    fplot(ax5, f, xlim(ax5))
end
if QC_flag == 0
    title(ax5, 'POOR QUALITY flag tripped', 'color', 'r', 'fontweight', 'bold')
end
make_density_plot(ax6,npar_eukX,npar_eukY,mylim1)
make_plot(ax7,11,15, mylim1)
%make_plot(ax8,17,18, mylim1)
make_plot(ax8, pe_chs(2)-1 , pe_chs(2), mylim1)
if exist('bounds', 'var')
    line(ax8, xlim(ax8), xlim(ax8)*nonsynfactorA) %5
    line(ax8, xlim(ax8), xlim(ax8)*nonsynfactorB) %2.5e
end
make_plot(ax9,20, 11, mylim2) %fsc-W vs H 
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
       
        if ~isempty(mylim)
        axis(ax,mylim)
        end
        set(ax, 'xtick', [1e2 1e4 1e6])
        set(ax, 'ytick', [1e2 1e4 1e6])
        ax.XLabel.FontSize = 6;  
        ax.YLabel.FontSize = 8;
       
    end

    function make_density_plot(ax, nparX, nparY, mylim)

        %mycenters= {logspace(log10(mylim1(1)), log10(mylim1(2)), 40) logspace(log10(mylim1(3)), log10(mylim1(4)), 40)};        
        mycenters= {log10(mylim1(1)):.1:log10(mylim1(2)) log10(mylim1(3)):.1:log10(mylim1(4))};
        if sum(class == 1)>10
        hist3(ax, [log10(fcsdat(class == 1,nparX)), log10(fcsdat(class == 1,nparY))], mycenters, 'CDataMode','auto', 'EdgeColor', 'none'); 
        end;
        view(ax, 2)
        set(ax, 'ZScale', 'log')
        
        hold(ax, 'on')
        if sum(class == 2)>10
        hist3(ax, [log10(fcsdat(class == 2,nparX)), log10(fcsdat(class == 2,nparY))], mycenters, 'CDataMode','auto', 'EdgeColor', 'none'); 
        end
        if sum(class == 3)>10
        hist3(ax, [log10(fcsdat(class == 3,nparX)), log10(fcsdat(class == 3,nparY))], mycenters, 'CDataMode','auto', 'EdgeColor', 'none'); 
        end
        if sum(class == 4)>10
        hist3(ax, [log10(fcsdat(class == 4,nparX)), log10(fcsdat(class == 4,nparY))], mycenters, 'CDataMode','auto', 'EdgeColor', 'none'); 
        end
        if sum(class == 5)>10
        hist3(ax, [log10(fcsdat(class == 5,nparX)), log10(fcsdat(class == 5,nparY))], mycenters, 'CDataMode','auto', 'EdgeColor', 'none'); 
        end

        %[values, centers] = hist3([log10(fcsdat(:,nparX)), log10(fcsdat(:,nparY))], mycenters, 'CDataMode','auto')
        % imagesc(centers{:}, values')
        %set(ax, 'Ydir', 'normal')

        ylabel(ax,fcshdr.par(nparY).name)
        xlabel(ax,fcshdr.par(nparX).name)

        ax.XLabel.FontSize = 6;  
        ax.YLabel.FontSize = 8;

               
    end



    function make_legend
            leglabels = {'Noise'; 'Euks'; 'Syn'; 'LowP Euk'; 'HiP Euk'; 'SynEuk Coincident'; 'SynEuk Coincident2'; 'Noise'}; 
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
