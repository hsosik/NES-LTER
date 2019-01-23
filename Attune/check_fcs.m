function [ output_args ] = check_fcs( filename )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

basepath = '\\sosiknas1\Lab_data\Attune\AR29\';
fpath = [basepath filesep 'FCS_backup_AR29' filesep];
outpath = [basepath filesep 'Summary' filesep 'class' filesep];
%f = 'SFD_AR29_23Apr2018A_Group_day0_Sample(111).fcs'


    [~,fcshdr,fcsdatscaled] =fca_readfcs([fpath filename]);
    load([outpath regexprep(filename,'.fcs', '')], 'class');
        figure(1), clf
        ax1 = subplot(2,3,1);
        ax2 = subplot(2,3,2);
        ax3 = subplot(2,3,3);
        ax4 = subplot(2,3,4);
        ax5 = subplot(2,3,5);
        ax6 = subplot(2,3,6);
        mylim1 = [1e1 1e6 1e1 1e6];
        make_plot(ax1,12,15,mylim1)
     %   plot(ax1,x_polygon, y_polygon, 'g.-')
        make_plot(ax2,3,10,mylim1)
     %   plot(ax2,x_rect,y_rect, 'r.-')
        make_plot(ax3,19,10,mylim1)
        axes(ax3)
        line(xlim,xlim*4)
        %make_plot(ax4,6,10,mylim1)
        make_plot(ax4,11,19,mylim1) 
        make_plot(ax5,15,19,mylim1) 
        make_plot(ax6,6,10,mylim1)
        axes(ax6)
        line(xlim, xlim*3.5)
        pause(.01)
        
        figure(2), clf
        subplot(1,2,1)
        make_plot(gca,12,3,mylim1)
        line(xlim, xlim)
        subplot(1,2,2)
        t = find(fcsdatscaled(:,12)>500 & fcsdatscaled(:,3)>100);
        disp(median(fcsdatscaled(t,3)./fcsdatscaled(t,12)))
        histogram(fcsdatscaled(t,3)./fcsdatscaled(t,12));
        
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

