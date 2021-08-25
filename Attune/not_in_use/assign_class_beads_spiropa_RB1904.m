function [ class ] = assign_class_spiropa( fcsdat, fcshdr, plot_flag, filename, QC_flag )
%function [ class ] = assign_class_spiropa( fcsdat, fcshdr, plot_flag )
%
%   input:  data and hdr from FCS file (fcsdat, fcshdr)
%           boolean flag for plotting of cytograms
%   output: vector of integers specifying classes
%   These gates configured for SPIROPA cruise AR29, April 2018
%
%Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2019

%Initialze class vector
    class = zeros(size(fcsdat,1),1);
    
    %parameter numbers for initial bead windows
    nparX = 19; %GL1-H
    nparY = 15; %BL3-H
    
    %13:19 {'BL1-H'}    {'BL2-H'}    {'BL3-H'}    {'GL4-H'}    {'GL3-H'}    {'GL2-H'}    {GL1-H'}
%     a = 13; b = 19;
%      gmin1(13:19) = [50      90      1e6 1e6 1e6 2.3e4  1.3e5 ];
%      gmax1(13:19) = [180     275     2e6 2e6 1e6 3.3e4  2e5 ];
%      gmin2(13:19) = [4e3     1.9e3   300 25  25  350    800];
%      gmax2(13:19) = [2.5e4   2.5e4   850 250 100 1100   1800];
%      gmin3(13:19) = [38 1000 250 750 550 1.4e4  400 ];
%      gmax3(13:19) = [140 2500 525 1800 1000 2.1e4 1000 ];

    a = 13; b = 18;
    gmin1(a:b) = [50      90    1e6 1e6 1e6 2.3e4];%  1.3e5 ];
    gmax1(a:b) = [180     275   2e6 2e6 2e6 3.3e4];%  2e5 ];
    gmin2(a:b) = [4e3     1.9e3 300 25  25  350];%    800];
    gmax2(a:b) = [2.5e4   2.5e4 850 250 100 1100];%   1800];
    gmin3(a:b) = [38      1000  250 750 550 1.4e4];%  400 ];
    gmax3(a:b) = [140     2500  525 1800 1000 2.1e4];% 1000 ];

    numpar = b-a+1;
    in1 = NaN(size(fcsdat,1),numpar); in2 = in1; in3 = in1;
    for cc = a:b
        in1(:,cc-a+1) = (fcsdat(:,cc)>gmin1(cc)&fcsdat(:,cc)<gmax1(cc));
        in2(:,cc-a+1) = (fcsdat(:,cc)>gmin2(cc)&fcsdat(:,cc)<gmax2(cc));
        in3(:,cc-a+1) = (fcsdat(:,cc)>gmin3(cc)&fcsdat(:,cc)<gmax3(cc));
    end
    
%      gg = load('c:\work\SPIROPA\RB1904\bead_win_ssc160');
%      g1 = gg.g1_19_15;
%      g2 = gg.g2_19_15;
%      g3 = gg.g3_19_15;
    
    %find indices of cells within the gates
        
%     fcsdatlog = log10(fcsdat); %use log10 to make sure inpolygon corresponds to view of polygon on log-log plots
%     in1 = find(isinpoly(fcsdatlog(:,19), fcsdatlog(:,15), log10(g1(:,1)), log10(g1(:,2))));
%     in2 = find(isinpoly(fcsdatlog(:,19), fcsdatlog(:,15), log10(g2(:,1)), log10(g2(:,2))));
%     in3 = find(isinpoly(fcsdatlog(:,19), fcsdatlog(:,15), log10(g3(:,1)), log10(g3(:,2))));
    
    
    %assign values in class vector
    class(sum(in1,2)==numpar) = 1;
    class(sum(in2,2)==numpar) = 2;
    class(sum(in3,2)==numpar) = 3;
        
    if plot_flag
        figure(1), clf
        ax1 = subplot(2,3,1);
        ax2 = subplot(2,3,2);
        ax3 = subplot(2,3,3);
        ax4 = subplot(2,3,4);
        ax5 = subplot(2,3,5);
        ax6 = subplot(2,3,6);
        mylim1 = [1e1 1.1e6 1e1 1.1e6];
        make_plot(ax1,nparX,nparY,mylim1)
        %plot(ax1,euk_x_polygon, euk_y_polygon, 'g.-')
        make_plot(ax2,12,nparY,mylim1)
%         plot(ax2,syn_x_polygon,syn_y_polygon, 'r.-')
         make_plot(ax3,19,10,mylim1)
         title(ax2, filename, 'interpreter', 'none')
        %make_plot(ax4,12,14,mylim1)
        make_plot(ax4,14,nparY,mylim1)
        make_histo(ax5,12,mylim1)
        make_histo(ax6,14,mylim1)
        if QC_flag == 0
            title(ax5, 'POOR QUALITY flag tripped', 'color', 'r', 'fontweight', 'bold')
        end
        pause %(.01)
    end
        
    function make_plot(ax, nparX, nparY, mylim)
        loglog(ax,fcsdat(:,nparX), fcsdat(:,nparY), '.', 'markersize', 1)
        ylabel(ax,fcshdr.par(nparY).name)
        xlabel(ax,fcshdr.par(nparX).name)
        axis(ax,'square')
        hold(ax, 'on')
        loglog(ax,fcsdat(class==1,nparX), fcsdat(class==1,nparY), '.g', 'markersize', 1)
        loglog(ax,fcsdat(class==2,nparX), fcsdat(class==2,nparY), 'r*', 'markersize', 1)
        loglog(ax,fcsdat(class==3,nparX), fcsdat(class==3,nparY), '^', 'color' ,[1 .5 0],  'markersize', 1)
        axis(ax,mylim) 
    end   

    function make_histo(ax, nparX, mylim)
        bins = logspace(2,6,128);
        h = histogram(ax, fcsdat(class==1,nparX),bins,'facecolor', 'g', 'edgecolor', 'none');
        [~,tt] = max(h.Values); l(1) = mean(h.BinEdges(tt:tt+1));
        hold(ax, 'on')
        h = histogram(ax, fcsdat(class==2,nparX),bins,'facecolor', 'r', 'edgecolor', 'none');
        [~,tt] = max(h.Values); l(2) = mean(h.BinEdges(tt:tt+1));
        title(ax,['Mode = ' num2str(round(mean(h.BinEdges(tt:tt+1))))], 'color', 'r')
        h = histogram(ax, fcsdat(class==3,nparX),bins,'facecolor', [1 .5 0], 'edgecolor', 'none');
        [~,tt] = max(h.Values); l(3) = mean(h.BinEdges(tt:tt+1));
        set(ax, 'xscale', 'log')
        xlabel(ax,fcshdr.par(nparX).name)
        legend(ax, num2str(round(l')))
    end

end
    
