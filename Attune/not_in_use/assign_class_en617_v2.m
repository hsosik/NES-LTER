function [ class ] = assign_class_en617( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )
%function [ class ] = assign_class_spiropa( fcsdat, fcshdr, plot_flag )
%
%   input:  data and hdr from FCS file (fcsdat, fcshdr)
%           boolean flag for plotting of cytograms
%   output: vector of integers specifying classes These gates configured
%   for SPIROPA cruise AR29, April 2018
%
%Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2019

%Initialze class vector
    class = zeros(size(fcsdat,1),1);
    
    %parameter numbers for main euk polygon
    npar_eukX = 15; %BL3-H
    npar_eukY = 19; %GL1-H
    
    %parameter numbers for main syn polygon
    npar_synX = 11; %FSC-H
    npar_synX = 12; %SSC-H
    npar_synY = 19; %GL1-H
    %simple rectangle for syn main gate
    %gsyn_11_19 = [100 300 ; 8e3 300; 8e3 4e5; 100 4e5]; %[Xmin Ymin; X     max Ymax]
    gsyn_12_19 = [1 350 ; 3e2 350; 3e2 4e5; 1 4e5]; %[Xmin Ymin; Xmax Ymax]
    
   % if strmatch('SPIROPA_RB1904_Grazer',filename)
   % end
    eukminX = 500;
    geuk_15_19 = [eukminX 199; 3000 199; 1100000 2e4; 1100000 1; eukminX 1];
    
    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 1.2; %PE to CHL SPIROPA 1.2
    synGL1H2BL3Hoffset = -.4; %PE to CHL SPIROPA .4
    syneukBL3H2SSCHslope = 1.3; %CHL to SSC
    syneukBL3H2SSCHoffset = -.8; %-.6; %PE to CHL
    synnonsynGL1H2GL2Hhi = 1.2; %SPIROPA 1.2
    synnonsynGL1H2GL2Hlo = .6; %SPIROPA 0.6
    
    %find indices of cells within the gates
    fcsdatlog = log10(fcsdat); %use log10 to make sure inpolygon corresponds to view of polygon on log-log plots
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_15_19(:,1)),log10(geuk_15_19(:,2))); 
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_12_19(:,1)),log10(gsyn_12_19(:,2))));
    minY = prctile(fcsdat(in_syn,19),10)*.4; maxY = prctile(fcsdat(in_syn,19),90)*10;
    minY = max(minY, 250);
    eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.4;
    geuk_15_19 = [eukminX 199; 1e4 199; 1100000 2e4; 1100000 1; eukminX 1];
    if (startsWith(filename, 'NESLTER_EN617_20Jul2018A') | startsWith(filename, 'NESLTER_EN617_20Jul2018B')) %before inc PE HV or Chl HV
        synnonsynGL1H2GL2Hhi = .3; %SPIROPA 1.2
        synnonsynGL1H2GL2Hlo = .2; %SPIROPA 0.6
        synGL1H2BL3Hoffset = -.1; %PE to CHL SPIROPA .4
        minY = 200;
        eukminX = 150;
        geuk_15_19 = [eukminX 199; 5000 199; 1100000 5e4; 1100000 1; eukminX 1];
    elseif (startsWith(filename, 'NESLTER_EN617_20Jul2018C') | startsWith(filename, 'NESLTER_EN617_21Jul2018A')) %after inc PE HV, before inc Chl HV
        synGL1H2BL3Hslope = 1.9; %PE to CHL SPIROPA 1.2
        synGL1H2BL3Hoffset = -1.5; %PE to CHL SPIROPA .4
        eukminX = 150;
        geuk_15_19 = [eukminX 199; 3000 199; 1100000 2e4; 1100000 1; eukminX 1];
    end;
    %gsyn_11_19 = [100 minY ; 8e3 minY; 8e3 maxY; 100 maxY]; %[Xmin Ymin; Xmax Ymax]
    gsyn_12_19 = [1 minY ; 3e2 minY; 3e2 maxY; 1 maxY]; %[Xmin Ymin; Xmax Ymax]
    
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_15_19(:,1)),log10(geuk_15_19(:,2))); 
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_12_19(:,1)),log10(gsyn_12_19(:,2))));
    in_nonsyn_lowPE = fcsdat(:,npar_synY) > minY & fcsdat(:,npar_synY) < maxY & fcsdat(:,19)./fcsdat(:,18) < synnonsynGL1H2GL2Hlo & ~(fcsdat(:,19)<1e4 & fcsdat(:,11)>1e4);
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > maxY & fcsdat(:,19)./fcsdat(:,18) > synnonsynGL1H2GL2Hlo & (fcsdat(:,19)./fcsdat(:,18) < synnonsynGL1H2GL2Hhi | fcsdat(:,18) > 1e6);
        
    %assign values in class vector
    
    class(in_nonsyn_lowPE) = 3;
    class(in_nonsyn_hiPE) = 4;
    class(in_syn) = 2; %must be done after nonsyn
    class(in_syn & (fcsdatlog(:,19)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_15_19(:,1)))) = 5; %Syn&Euk coincident
    class(in_euk) = 1; %AFTER lowPE
    class(class == 0 & fcsdat(:,19)>minY & fcsdat(:,19)<maxY/1 & fcsdatlog(:,15)>fcsdatlog(:,12)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 6; %more Syn&Euk coincident
    class(fcsdat(:,15) < 200 & fcsdat(:,19) < 200) = 7; %noise %LAST
    
    if plot_flag
        figure(1), clf
        ax1 = subplot(2,5,1);
        ax2 = subplot(2,5,2);
        ax3 = subplot(2,5,3);
        ax4 = subplot(2,5,4);
        ax5 = subplot(2,5,5);
        ax6 = subplot(2,5,6);
        ax7 = subplot(2,5,7);
        ax8 = subplot(2,5,8);
        ax9 = subplot(2,5,9);
        ax10 = subplot(2,5,10);
        
        mylim1 = [1e1 1.1e6 1e1 1.1e6];
        mylim2 = [1e3 1.1e6 1e3 1.1e6];
        
        make_plot(ax1,npar_eukX,npar_eukY,mylim1)
        %plot(ax1,euk_x_polygon, euk_y_polygon, 'g.-')
        plot(ax1,geuk_15_19(:,1), geuk_15_19(:,2), 'g.-')
        f = @(x) 10.^(log10(x).*synGL1H2BL3Hslope+synGL1H2BL3Hoffset);
        fplot(ax1, f, xlim(ax1))
        make_plot(ax2,npar_synX,npar_synY,mylim1)
        plot(ax2,gsyn_12_19(:,1),gsyn_12_19(:,2), 'r.-')
        make_plot(ax3,19,10,mylim1)
        title(ax2, filename, 'interpreter', 'none')
        line(ax3,xlim(ax3),xlim(ax3)*synGL1A2GL1Hmax)
        make_plot(ax4,12,19,mylim1) 
        make_plot(ax5,17,15,mylim1) %BL3H v GL3H
        if QC_flag == 0
            title(ax5, 'POOR QUALITY flag tripped', 'color', 'r', 'fontweight', 'bold')
        end
        %make_plot(ax6,6,10,mylim1)
        make_plot(ax6,12,15,mylim1)
        f = @(x) 10.^(log10(x).*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset);
        fplot(ax6, f, xlim(ax6))
        make_plot(ax7,11,15, mylim1)
        %make_plot2(ax7,17,8, mylim1)
        make_plot(ax8,17,8, mylim1)
        make_plot(ax9,18,19, mylim2)
        line(ax9, xlim(ax9), xlim(ax9)*synnonsynGL1H2GL2Hhi)
        line(ax9, xlim(ax9), xlim(ax9)*synnonsynGL1H2GL2Hlo)
        make_plot(ax10,11,12, mylim1)
        pause (.01)
    end
        
    function make_plot(ax, nparX, nparY, mylim)
        loglog(ax,fcsdat(:,nparX), fcsdat(:,nparY), '.', 'markersize', 1)
        ylabel(ax,fcshdr.par(nparY).name)
        xlabel(ax,fcshdr.par(nparX).name)
        axis(ax,'square')
        hold(ax, 'on')
        loglog(ax,fcsdat(class==1,nparX), fcsdat(class==1,nparY), '.g', 'markersize', 1)
        loglog(ax,fcsdat(class==2,nparX), fcsdat(class==2,nparY), 'r.', 'markersize', 1)
        loglog(ax,fcsdat(class==3,nparX), fcsdat(class==3,nparY), '^', 'color' ,[1 .5 0],  'markersize', 1)
        loglog(ax,fcsdat(class==4,nparX), fcsdat(class==4,nparY), 'm.', 'markersize', 1)
        loglog(ax,fcsdat(class==5,nparX), fcsdat(class==5,nparY), 'c.', 'markersize', 2)
        loglog(ax,fcsdat(class==7,nparX), fcsdat(class==7,nparY), '.', 'color', [.4 .4 .4], 'markersize', 1)
        loglog(ax,fcsdat(class==6,nparX), fcsdat(class==6,nparY), '.y', 'markersize', 1.5) %'color', [.7 .2 1]
        axis(ax,mylim) 
        set(ax, 'xtick', [1e2 1e4 1e6])
    end   
    function make_plot2(ax, nparX, nparY, mylim)
        %semilogy(ax,fcsdat(:,1), fcsdat(:,nparY)./fcsdat(:,nparX), '.', 'markersize', 1)
        ylabel(ax,[fcshdr.par(nparY).name ' / ' fcshdr.par(nparX).name])
        xlabel(ax,fcshdr.par(1).name)
        axis(ax,'square')
        hold(ax, 'on')
       % semilogy(ax,fcsdat(class==1,1), fcsdat(class==1,nparY)./fcsdat(class==1,nparX), '.g', 'markersize', 1)
        semilogy(ax,fcsdat(class==2,1), fcsdat(class==2,nparY)./fcsdat(class==2,nparX), 'r*', 'markersize', 1)
        semilogy(ax,fcsdat(class==3,1), fcsdat(class==3,nparY)./fcsdat(class==3,nparX), '^', 'color' ,[1 .5 0], 'markersize', 1)
        semilogy(ax,fcsdat(class==4,1), fcsdat(class==4,nparY)./fcsdat(class==4,nparX), 'm*', 'markersize', 1)
        semilogy(ax,fcsdat(class==5,1), fcsdat(class==5,nparY)./fcsdat(class==5,nparX), 'c.')
        semilogy(ax,fcsdat(class==6,1), fcsdat(class==6,nparY)./fcsdat(class==6,nparX), 'k.', 'markersize', 1)
        ylim(ax, [0 4]) 
    end

end
    
