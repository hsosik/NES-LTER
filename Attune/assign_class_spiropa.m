function [ class ] = assign_class_spiropa( fcsdat, fcshdr, plot_flag )
%function [ class ] = assign_class_spiropa( fcsdat, fcshdr, plot_flag )
%
%   input:  data and hdr from FCS file (fcsdat, fcshdr)
%           boolean flag for plotting of cytograms
%   output: vector of integers specifying classes
%   These gates configured for SPIROPA cruise AR29, April 2018
%
%Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2019

%Initialze class vector
    class = zeros(size(fcsdat,1));
    
    %parameter numbers for main euk polygon
    npar_eukX = 12; %SSC-H
    npar_eukY = 15; %BL1-H
    %polygon vertices for euks
    euk_x_polygon = [10^1.5  10^1.5 10^6 10^6   10^5    10^4.3  10^3.7  10^2.5      10^1.5 ];
    euk_y_polygon = [10^2.8  10^6 10^6 10^5.8 10^5.5  10^4.8  10^4  10^3        10^2.8 ];

    %parameter numbers for main syn polygon
    npar_synX = 3; %SSC-A
    npar_synY = 10; %GL1-A
    %simple rectangle for syn main gate
    SynXmin= 1; SynXmax= 10^3;
    SynYmin= 10^3; SynYmax= 10^5;
    syn_x_polygon = [SynXmin SynXmin SynXmax SynXmax SynXmin];
    syn_y_polygon = [SynYmin SynYmax SynYmax SynYmin SynYmin];
    synGL1A2GL1Hmax = 4; %PE area to height
    synFSCHmax = 1e4;
    synGL1A2BL3Amax = 3.5; %PE to CHL
    
    %find indices of cells within the gates
    fcsdatlog = log10(fcsdat); %use log10 to make sure inpolygon corresponds to view of polygon on log-log plots
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(euk_x_polygon),log10(euk_y_polygon)); 
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(syn_x_polygon),log10(syn_y_polygon)) & fcsdat(:,11)<synFSCHmax & fcsdat(:,10)./fcsdat(:,19) < synGL1A2GL1Hmax & fcsdat(:,10)./fcsdat(:,6) > synGL1A2BL3Amax);

    %assign values in class vector
    class(in_euk) = 1;
    class(in_syn) = 2;
    
    if plot_flag
        figure(1), clf
        ax1 = subplot(2,3,1);
        ax2 = subplot(2,3,2);
        ax3 = subplot(2,3,3);
        ax4 = subplot(2,3,4);
        ax5 = subplot(2,3,5);
        ax6 = subplot(2,3,6);
        mylim1 = [1e1 1e6 1e1 1e6];
        make_plot(ax1,12,15,mylim1)
        plot(ax1,euk_x_polygon, euk_y_polygon, 'g.-')
        make_plot(ax2,3,10,mylim1)
        plot(ax2,syn_x_polygon,syn_y_polygon, 'r.-')
        make_plot(ax3,19,10,mylim1)
        axes(ax3)
        line(xlim,xlim*synGL1A2GL1Hmax)
        %make_plot(ax4,6,10,mylim1)
        make_plot(ax4,11,19,mylim1) 
        make_plot(ax5,15,19,mylim1) 
        make_plot(ax6,6,10,mylim1)
        axes(ax6)
        line(xlim, xlim*synGL1A2BL3Amax)
        pause (.01)
     end
        
    function make_plot(ax, nparX, nparY, mylim)
        loglog(ax,fcsdat(:,nparX), fcsdat(:,nparY), '.', 'markersize', 1)
        ylabel(ax,fcshdr.par(nparY).name)
        xlabel(ax,fcshdr.par(nparX).name)
        axis(ax,'square')
        hold(ax, 'on')
        loglog(ax,fcsdat(class==1,nparX), fcsdat(class==1,nparY), '.g', 'markersize', 1)
        loglog(ax,fcsdat(class==2,nparX), fcsdat(class==2,nparY), 'r*', 'markersize', 1)
        axis(ax,mylim) 
    end   

end
    
