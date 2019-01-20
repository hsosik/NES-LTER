function [ class ] = assign_class_spiropa( fcsdatscaled, fcshdr, plot_flag )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Vector to indicate class of each event
    class = zeros(numel(fcsdatscaled(:,1)),1);
    
    npar_eukX = 12; %SSC-H
    npar_eukY = 15; %BL3-H
    %defining the polygon gate for the Small Eukaryote Signal
    x_polygon = [10^1.5  10^1.5 10^6 10^6   10^5    10^4.3  10^3.7  10^2.5      10^1.5 ];
    y_polygon = [10^2.3  10^6 10^6 10^5.8 10^5.5  10^4.8  10^4  10^2.5        10^2.3 ];

    npar_synX = 3; %SSC-A
    npar_synY = 10; %GL1-A
    %npar_synY = 19; %GL1-H
    %defining the rectangular gate for the Synechecoccus Signal
    %SynXmin= 200; SynXmax= 10^4;
    SynXmin= 1; SynXmax= 10^3;
    SynYmin= 10^3.8; SynYmax= 10^5.5;
    
    x_rect = [SynXmin SynXmin SynXmax SynXmax SynXmin];
    y_rect = [SynYmin SynYmax SynYmax SynYmin SynYmin];
    
    %counting cells within the gates
    in_euk = (inpolygon(fcsdatscaled(:,npar_eukX),fcsdatscaled(:,npar_eukY),x_polygon,y_polygon) & fcsdatscaled(:,15) > 1e3);
    in_syn = (inpolygon(fcsdatscaled(:,npar_synX),fcsdatscaled(:,npar_synY),x_rect,y_rect) & fcsdatscaled(:,11)<1e4 & fcsdatscaled(:,10)./fcsdatscaled(:,19) < 4 & fcsdatscaled(:,10)./fcsdatscaled(:,6) > 3.5);

    %defining euks and syn in the class
    class(in_euk) = 1;
    class(in_syn) = 2;
    
    if plot_flag
        figure, clf
        ax1 = subplot(2,3,1);
        ax2 = subplot(2,3,2);
        ax3 = subplot(2,3,3);
        ax4 = subplot(2,3,4);
        ax5 = subplot(2,3,5);
        ax6 = subplot(2,3,6);
        mylim1 = [1e1 1e6 1e1 1e6];
        make_plot(ax1,12,15,mylim1)
        plot(ax1,x_polygon, y_polygon, 'g.-')
        make_plot(ax2,3,10,mylim1)
     %   make_plot(ax2,3,19,mylim1)
        plot(ax2,x_rect,y_rect, 'r.-')
        make_plot(ax3,19,10,mylim1)
        axes(ax3)
        line(xlim,xlim*4)
        %make_plot(ax4,6,10,mylim1)
        make_plot(ax4,11,19,mylim1) 
        make_plot(ax5,15,19,mylim1) 
        make_plot(ax6,6,10,mylim1)
        axes(ax6)
        line(xlim, xlim*3.5)
        pause %(.01)
     end
        
    function make_plot(ax, nparX, nparY, mylim)
        loglog(ax,fcsdatscaled(:,nparX), fcsdatscaled(:,nparY), '.', 'markersize', 1)
        ylabel(ax,fcshdr.par(nparY).name)
        xlabel(ax,fcshdr.par(nparX).name)
        axis(ax,'square')
        hold(ax, 'on')
        loglog(ax,fcsdatscaled(class==1,nparX), fcsdatscaled(class==1,nparY), '.g', 'markersize', 1)
        loglog(ax,fcsdatscaled(class==2,nparX), fcsdatscaled(class==2,nparY), 'r*', 'markersize', 1)
        axis(ax,mylim) 
    end   

end
    
