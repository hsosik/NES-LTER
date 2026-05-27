%Function formatted to be modular bead processing called by
%process_wrapper_2021
%identifies clusters of beads, saves statistics, and creates bead_plots
%which are saved to output directory

%Inputs:
% outpath is location of output directory used in wrapper

%FCSfileinfo is output of FCS_DateTimeList(fpath) or loaded from .mat file

%beadfiles2include is the string that starts bead files for the relevant
%cruise, e.g. sometimes its {'Daily bead check'} or {fcb_bead'}

%beadtype can be either 'PT' or 'FCB' -- haven't made this work for PT
%beads yet.


% output details and changes
% bead plots include dot plots for all bead files considered, with center
% values plotted for each of the 3 clusters.

% March 2026, Heidi
% Start process_beads_PT_adjust_3.m from process_beads_PT_adjust_2.m
% working on new bead clustering method for standard case with OD2 on GL1 

function [] = process_beads_PT_adjust_3(basepath, FCSfileinfo, beadfiles2include, beadtype, OD2setting)

% make output directory
outpath = [basepath filesep 'bead_calibrated' filesep];
beadfigpath = [outpath filesep 'bead_plots_2026'];
if ~exist(beadfigpath, 'dir')
    mkdir(beadfigpath)
end

%now cut down to just focus on beads in our include list
t = [];
for iii = 1:length(beadfiles2include)
    t = [t strmatch(beadfiles2include{iii}, FCSfileinfo.fcslist)];
end
beadlist = FCSfileinfo.fcslist(t);
fpath = [outpath filesep '..' filesep 'FCS' filesep];

%initialize bead table
beadstat_2026 = table;

%run through bead files in FCS directory (assumed to be one up from outpath)
for b = 1:length(beadlist)
    %if b == 20 
    %    continue
    %end
    clear temp_table
    disp(b)
    
    [fcsdat, fcshdr] = fca_readfcs([fpath beadlist{b}]);
    fcsdat = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});

    %ssca_chn = strmatch('SSC-A', {fcshdr.par.name});
    %gl1a_chn = strmatch('GL1-A', {fcshdr.par.name});

    bead_ch_names = {};
    if strcmp(OD2setting,'GL1')
        bead_ch_names{1} = 'GL1-A';
        %sidescatter_ch = 'GL1-A';
    else %OD2 filter is on SSC or not on at all
        bead_ch_names{1} = 'SSC-A';
        %sidescatter_ch = 'SSC-A';
    end
    bead_ch_names{4} = 'SSC-A'; 

    %ch = nan(1,3);
    %ch(1) = sidescatter_ch;
    %ch(2) = strmatch('BL3-H', {fcshdr.par.name});  %chlorophyl channel
    bead_ch_names{2} = 'BL3-H';
    if fcshdr.date < datetime(2019,08,05) %3rd channel changes with cruise settings
        %ch(3) = strmatch('GL2-H', {fcshdr.par.name});
        bead_ch_names{3} = 'GL2-H';
    else
        %ch(3) = strmatch('GL3-H', {fcshdr.par.name});
        bead_ch_names{3} = 'GL3-H';
    end
    
    
    temp_table = table;
    temp_table.filename = beadlist(b);
    temp_table.time = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    temp_table.parname = {fcshdr.par.name};
    temp_table.OD2_filter = {OD2setting};
    temp_table.hv = cell2mat({fcshdr.par.hv});
    

    figure(99),clf
    set(gcf,'Position', [50 50 1000 550])
            
    if strcmp(beadtype, 'FCB')
        n_clust = 3;
             
        % cluster analysis
        chl_noise_cut = 200; %BL3-H
        GL3_bead_cut = 500; %hard cutoff between 0.5 and 1 micron beads on GL3-H 
        class = zeros(size(fcsdat,1),1); % c is a vector of cluster indices, 0 indicates not bead
       
        bd_fig_flag = 0;
        %find the 0.5 micron beads (class=1), force SSC-A>0 for weird noise for EN720
        tempi = find(fcsdat.(bead_ch_names{2})>chl_noise_cut & fcsdat.(bead_ch_names{3})>GL3_bead_cut & fcsdat.(bead_ch_names{2})<1e6 & fcsdat.(bead_ch_names{4})>0); %consider particles above threshhold
        if ~isempty(tempi)
            %    keyboard
            minpts = max([10 floor(length(tempi)/100)]); %dbscan parameters
            X = fcsdat{tempi,bead_ch_names([4,2])};
            ctemp2 = dbscan(real(log10(X+1)), .001, minpts, 'distance', 'squaredeuclidean');
            if max(ctemp2)>0
                [~,cc] = max(histcounts(ctemp2,1:max(ctemp2)+1)); % cc is the cluster with highest count
                class(tempi(ctemp2==cc)) = 1;
            end
            if bd_fig_flag, figure(1),clf, gscatter(X(:,1),X(:,2),ctemp2); set(gca, 'yscale', 'log', 'xscale', 'log'), xlabel(bead_ch_names(4)), ylabel(bead_ch_names(2)), pause, end
        end

        % find the 1 micron beads (class=2)
        tempi = find(fcsdat.(bead_ch_names{2})>chl_noise_cut & fcsdat.(bead_ch_names{3})<GL3_bead_cut & fcsdat.(bead_ch_names{4})>0); %consider particles above threshhold
        minpts = max([10 floor(length(tempi)/100)]); %dbscan parameters
        X = fcsdat{tempi,bead_ch_names([4,2])};
        ctemp2 = dbscan(real(log10(X+1)), .001, minpts, 'distance', 'squaredeuclidean');
        if max(ctemp2)>0
            [~,cc] = max(histcounts(ctemp2,1:max(ctemp2)+1)); % cc is the cluster with highest count
            class(tempi(ctemp2==cc)) = 2;
        end
        if bd_fig_flag, figure(1),clf, gscatter(X(:,1),X(:,2),ctemp2); set(gca, 'yscale', 'log', 'xscale', 'log'), xlabel(bead_ch_names(4)), ylabel(bead_ch_names(2)), pause, end

        % now cluster the stuff offscale on chlorophyll to find 6 micron beads (class = 3)
        tempi = find(fcsdat.(bead_ch_names{2})>=1e6);
        if length(tempi) ~= 0
            minpts = max([10 floor(length(tempi)/20)]); %dbscan parameters
            X = fcsdat{tempi,bead_ch_names([1,3])};
            ctemp2 = dbscan(real(log10(X+1)),.002, minpts, 'distance', 'squaredeuclidean');
            if max(ctemp2)>0
                [~,cc] = max(histcounts(ctemp2,1:max(ctemp2)+1)); % cc is the cluster with highest count
                class(tempi(ctemp2==cc)) = 3;
            end
            if bd_fig_flag, figure(1),clf, gscatter(X(:,1),X(:,2),ctemp2); set(gca, 'yscale', 'log', 'xscale', 'log'),xlabel(bead_ch_names(1)), ylabel(bead_ch_names(3)), pause, end
        end

        
            %% get center values for each bead cluster on OD2 channel and not 
            %depends on OD2settings
            
            if strcmp(OD2setting,'GL1')
                figure(99)
                fh = subplot(2,4,5);
                m05_noOd2A = centend2(fcsdat, class, 1, 'SSC-A');
                m1_noOd2A = centend2(fcsdat, class, 2, 'SSC-A', 1);
                m6_noOd2A = centend2(fcsdat, class, 3, 'SSC-A');
                xlabel('SSC-A, no OD'), title('1 micron beads'), xlim([3.5e4 6.5e4]), grid
                m05_noOd2H = centend2(fcsdat, class, 1, 'SSC-H');
                m1_noOd2H = centend2(fcsdat, class, 2, 'SSC-H', 1);
                m6_noOd2H = centend2(fcsdat, class, 3, 'SSC-H');
                
                fh = subplot(2,4,6);
                m05_od2A = centend2(fcsdat, class, 1, 'GL1-A');
                m1_od2A = centend2(fcsdat, class, 2, 'GL1-A', 1);
                m6_od2A = centend2(fcsdat, class, 3, 'GL1-A');
                xlabel('GL1-A, OD2 SSC'), title('1 micron beads'), xlim([600 1400]), grid
                m05_od2H = centend2(fcsdat, class, 1, 'GL1-H');
                m1_od2H = centend2(fcsdat, class, 2, 'GL1-H', 1);
                m6_od2H = centend2(fcsdat, class, 3, 'GL1-H');

                subplot(2,4,7)
                p(1) = plot(fcsdat.('GL1-A'), fcsdat.('SSC-A'), '.k', 'markersize', 2);
                hold on
                if sum(class==1)
                    p(2) = plot(fcsdat{class==1,'GL1-A'}, fcsdat{class==1,'SSC-A'}, '.');
                else
                    p(2) = plot(NaN, NaN, '.');
                end
                if sum(class==2)
                    p(3) = plot(fcsdat{class==2,'GL1-A'}, fcsdat{class==2,'SSC-A'}, '.');
                else
                    p(3) = plot(NaN, NaN, '.');
                end
                p(4) = line(xlim,m1_noOd2A/m1_od2A*xlim, 'linewidth',2); grid
                axis([-500 2000 -1e4 12e4])
                ylabel('SSC-A, no OD'), xlabel('GL1-A, OD2')
                legend(p(2:4), '0.5 \mum beads', '1 \mum beads', '1\mum, mode,SSC/GL1', 'Location','northwest', 'fontsize',6)

                subplot(2,4,8)
                ii = fcsdat.('SSC-A')>m1_noOd2A*.5 & fcsdat.('SSC-A')<m1_noOd2A*1.5;
                histogram(fcsdat{ii,'SSC-A'}./fcsdat{ii,'GL1-A'},30:1:80) %.01:.0005:.03)
                line(m1_noOd2A/m1_od2A*[1 1], ylim, 'linewidth',2)
                hold on
                histogram(fcsdat{class==2,'SSC-A'}./fcsdat{class==2,'GL1-A'},30:1:80) %.01:.0005:.03)
                xlim([30 80])
                xlabel('SSC-A/GL1-A') %xlabel('GL1-A/SSC-A, +/-25% 1 \mum bd mode')
                legend('+/-25% bd mode', '1\mum bds', 'Location','southeast', 'fontsize',6)

                hv_od2 = fcshdr.par(strcmp('GL1-A', {fcshdr.par.name})).hv;
                hv_noOd2 = fcshdr.par(strcmp('SSC-A', {fcshdr.par.name})).hv;
                
                od2_chn = 'GL1-A';
                noOd2_chn = 'SSC-A'; 
                
            elseif strcmp(OD2setting, 'SSC')
                m05_noOd2 = NaN;
                m1_noOd2 = NaN;
                m6_noOd2 = NaN; 
                
                m05_od2 = centend2(fcsdat, class, 1, ssca_chn);
                m1_od2 = centend2(fcsdat, class, 2, ssca_chn);
                m6_od2 = centend2(fcsdat, class, 3, ssca_chn);
            
                hv_od2 = fcshdr.par(ssca_chn).hv;
                hv_noOd2 = NaN; 
                
                od2_chn = ssca_chn; 
                noOd2_chn = NaN; 
                
            else %OD2setting == 'none'
                m05_noOd2 = centend2(fcsdat, class, 1, ssca_chn);
                m1_noOd2 = centend2(fcsdat, class, 2, ssca_chn);
                m6_noOd2 = centend2(fcsdat, class, 3, ssca_chn);
                
                m05_od2 = NaN;
                m1_od2 = NaN;
                m6_od2 = NaN;                 
                   
                hv_noOd2 = fcshdr.par(ssca_chn).hv;
                hv_od2 = NaN; 
                
                od2_chn = NaN; 
                noOd2_chn = ssca_chn; 
                
            end
        
          %%  plot bead clusters for visual verification
            
            
            figure(99)
            colormap(lines(4));
            subplot(2,4,1)
            %scatter(fcsdat(:,ch(1)), fcsdat(:,ch(2)), 1, class, 'filled')
            scatter(fcsdat.('SSC-A'), fcsdat.(bead_ch_names{2}), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 1.2e6]); ylim([10 1.2e6]);
            xticks([10.^[1:6]]), yticks([10.^[1:6]]), grid
            %xlabel(bead_ch_names{1}); ylabel(bead_ch_names{2})
            xlabel('SSC-A'); ylabel(bead_ch_names{2})
            
            subplot(2,4,2)
            scatter(fcsdat.(bead_ch_names{1}), fcsdat.(bead_ch_names{3}), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 1.2e6]); ylim([10 1.2e6]);
            xticks([10.^[1:6]]), yticks([10.^[1:6]]), grid
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{3})
            
            if ~isnan(od2_chn)
                
            subplot(2,4,3); hold on;
            cmap = colormap(lines(4));
            bins = logspace(1, 6.1, 256);
            histogram(fcsdat{class==1, od2_chn}, bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
            histogram(fcsdat{class==2, od2_chn}, bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            histogram(fcsdat{class==3, od2_chn}, bins, 'FaceColor', cmap(4,:), 'EdgeColor', 'none')
            yl = ylim;
            xlim([10 1.2e6]), xticks([10.^[1:6]]), grid
            plot([m05_od2A m05_od2A], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m1_od2A m1_od2A], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m6_od2A m6_od2A], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            set(gca, 'XScale', 'log')
            legend('0.5 \mum', '1 \mum', '6 \mum')
            xlabel('OD2 SSC')
            
            end
            if ~isnan(noOd2_chn)
            
            subplot(2,4,4); hold on;
            cmap = colormap(lines(4));
            bins = logspace(1, 6.1, 256);
            histogram(fcsdat{class==1, noOd2_chn}, bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
            histogram(fcsdat{class==2, noOd2_chn}, bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            histogram(fcsdat{class==3, noOd2_chn}, bins, 'FaceColor', cmap(4,:), 'EdgeColor', 'none')
            yl = ylim;
            xlim([10 1.2e6]), xlim([10 1.2e6]), xticks([10.^[1:6]]), grid
            plot([m05_noOd2A m05_noOd2A], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m1_noOd2A m1_noOd2A], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m6_noOd2A m6_noOd2A], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            set(gca, 'XScale', 'log')
            legend('0.5 \mum', '1 \mum', '6 \mum')
            xlabel('SSC No Filter')

            end
            
            
            subplot(2,4,1), title(beadlist(b), 'interpreter', 'none', 'HorizontalAlignment', 'left','FontSize',10)
            subplot(2,4,3), title([datestr(temp_table.time) '; ' bead_ch_names{1} ' hv = ' num2str(fcshdr.par(strcmp(bead_ch_names{1}, {fcshdr.par.name})).hv) ' (' bead_ch_names{1} ')'], 'HorizontalAlignment', 'left','FontSize',10)
            print(figure(99), fullfile(beadfigpath, regexprep(beadlist{b}, '.fcs', '.png')), '-dpng')
            
               
        % HEIDI: DOES THIS STILL WORK??
        % quality control on SSC-A. Depends on OD2 filter setup
        % Without OD2 filter, we are really only capturing 1micron bead
        % cluster... 
        % Heidi: April 2026 --this does not look useful as is
        if m1_od2A<(m6_od2A/5) && m05_od2A<(m1_od2A/5)
            QC_flag = 1;
        else
            QC_flag = 0;
        end
        temp_table.QC_flag = logical(QC_flag);
        %1 is bad, 0 is good.

        
    else %need to do PT bead clustering here
        n_clust = 2 ; % 2 sizes of PT beads 
        
        c = -1*ones(size(fcsdat,1),1); % c is a vector of cluster indeces, -1 indicates noise
        tempi = find(fcsdat(:,ch(2))<1e7); %consider particles below threshold 
        minpts = floor(length(tempi)/20); %this is a parameter for clustering algorithm
        ctemp = dbscan(real(log10(fcsdat(tempi,ch)+1)), .1, minpts, 'distance', 'squaredeuclidean');
        c(tempi) = ctemp;
        nc = max(c); %this is number of cluster indexes
        if nc < 2
            minpts = floor(length(tempi)/40);
            ctemp = dbscan(real(log10(fcsdat(tempi,ch)+1)), .1, minpts, 'distance', 'squaredeuclidean');
            c(tempi) = ctemp;
            nc = max(c);
        end
        if nc<2
            keyboard
        end
        meanc = NaN(length(ch), nc); %get mean signal on each channel for each cluster
        %loglog(fcsdat(:,ssca_chn), fcsdat(:,ch(2)), '.')
        for iii = 1:nc
            meanc(:,iii) = mean(fcsdat(c==iii,ch));
            %hold on 
            %loglog(fcsdat((ctemp == iii),ssca_chn), fcsdat((ctemp==iii),ch(2)), '.')
        end
      
        [~,c1] = min(meanc(2,:)); %first cluster is smallest on flouresc.
        c2 = find(meanc(2,:)>meanc(2,c1)); 
                
        class = zeros(length(c), 1); %reset class
        class(ismember(c, c1)) = 1; 
        class(ismember(c, c2)) = 2; %class now has two groups assinged 
        % 1 is 2.4 microns, 2 is 3.2 microns
        
    %% 
        if strcmp(OD2setting, 'GL1') %then we have to do projection
            
            od2_chn = gl1a_chn;
            noOd2_chn = ssca_chn; 
                
            % central tendency of the bead clusters
            m24_noOd2 = centend2(fcsdat, class, 1, ssca_chn);
            m32_noOd2 = centend2(fcsdat, class, 2, ssca_chn);
            
            m24_od2 = centend2(fcsdat, class, 1, gl1a_chn);
            m32_od2 = centend2(fcsdat, class, 2, gl1a_chn);

             hv_od2 = fcshdr.par(gl1a_chn).hv;
             hv_noOd2 = fcshdr.par(ssca_chn).hv;
             
        elseif strcmp(OD2setting, 'SSC')        
            
             m24_noOd2 = NaN;
             m32_noOd2 = NaN;
                
            % central tendency of the bead clusters
            m24_od2 = centend2(fcsdat, class, 1, ssca_chn);
            m32_od2 = centend2(fcsdat, class, 2, ssca_chn);
            
             hv_od2 = fcshdr.par(ssca_chn).hv;
             hv_noOd2 = NaN; 
                
             od2_chn = ssca_chn; 
             noOd2_chn = NaN; 
        else %OD2setting == 'none'
                m24_noOd2 = centend2(fcsdat, class, 1, ssca_chn);
                m32_noOd2 = centend2(fcsdat, class, 2, ssca_chn);
                
                m24_od2 = NaN;
                m32_od2 = NaN;
                   
                hv_noOd2 = fcshdr.par(ssca_chn).hv;
                hv_od2 = NaN; 
                
                od2_chn = NaN; 
                noOd2_chn = ssca_chn; 
        end
            
        
        
        %plot for visual verification 
            figure(99), clf
            colormap(lines(4));
            subplot(2,3,1)
            scatter(fcsdat(:,ch(1)), fcsdat(:,ch(2)), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 1.2e6]); ylim([10 1.2e6]);
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{2})
            
            subplot(2,3,2)
            scatter(fcsdat(:,ch(1)), fcsdat(:,ch(3)), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 2e6]); ylim([10 1.2e6]);
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{3})
            
            
            if ~isnan(od2_chn)
                
            subplot(2,3,4); hold on;
            cmap = colormap(lines(4));
            bins = logspace(1, 6.1, 256);
            histogram(fcsdat(class==1, od2_chn), bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==2, od2_chn), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            yl = ylim;
            xlim([10 2e6])
            plot([m24_od2 m24_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m32_od2 m32_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            set(gca, 'XScale', 'log')
            legend('2.4 \mum', '3.2 \mum')
            xlabel('OD2 SSC')
            
            end
            if ~isnan(noOd2_chn)
            
            subplot(2,3,5); hold on;
            cmap = colormap(lines(4));
            bins = logspace(1, 6.1, 256);
            histogram(fcsdat(class==1, noOd2_chn), bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==2, noOd2_chn), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            yl = ylim;
            xlim([10 2e6])
            plot([m24_noOd2 m24_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m32_noOd2 m32_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            set(gca, 'XScale', 'log')
            legend('2.4 \mum', '3.2 \mum')
            xlabel('SSC No Filter')
            
            subplot(2,3,3)
            histogram(fcsdat(class==2, od2_chn), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            hold on
            plot([m1_od2 m1_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            
            subplot(2,3,6)
            histogram(fcsdat(class==2, noOd2_chn), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            hold on
            plot([m1_noOd2 m1_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)

            end

        if m24_noOd2>=(m32_noOd2)
            QC_flag = 1;
        else
            QC_flag = 0;
        end
        temp_table.QC_flag = logical(QC_flag);


        %1 is bad, 0 is good.
            subplot(2,2,1), title(beadlist(b), 'interpreter', 'none')
            subplot(2,2,3), title([datestr(temp_table.time) '; ' bead_ch_names{1} ' hv = ' num2str(fcshdr.par(ch(1)).hv) ' (' fcshdr.par(ch(1)).name ')'])
            print(figure(99), fullfile(beadfigpath, regexprep(beadlist{b}, '.fcs', '.png')), '-dpng')
            
           
    end
    %%
    
    for iii = 1:n_clust
        nstr = num2str(iii);
        temp_table.(['mean' nstr]) = mean(fcsdat(class==iii,:));
        temp_table.(['std' nstr]) = std(fcsdat(class==iii,:));
        temp_table.(['median' nstr]) = median(fcsdat(class==iii,:));
        temp_table.(['number' nstr]) = sum(class==iii);
        temp_table.(['centend' nstr])= centend2(fcsdat, class, iii, bead_ch_names{1});
    end
    clear iii
    
    if strcmp(beadtype, 'FCB')
        %Store center on SSC-H channel only
        temp_table.OD2centersA = [m05_od2A m1_od2A m6_od2A];
        temp_table.OD2centersH = [m05_od2H m1_od2H m6_od2H];
        temp_table.OD2_hv = hv_od2;
    
        temp_table.NoOD2centersA = [m05_noOd2A m1_noOd2A m6_noOd2A];
        temp_table.NoOD2centersH = [m05_noOd2H m1_noOd2H m6_noOd2H];
        temp_table.NoOD2_hv = hv_noOd2;
        
        X = fcsdat.("SSC-A")(class==2)./fcsdat.("GL1-A")(class==2);
        if ~isempty(X)
            temp_table.mean1micron_SSCA2GL1A = mean(X);
            t = 30:1:80; %.01:.0005:.03; %t = .01:.0005:.05; %.05 needed for first file on en668
            temp_table.mode1micron_SSCA2GL1A = t(mode(discretize(X,t)));
        else
            temp_table.mean1micron_SSCA2GL1A = NaN;
            temp_table.mode1micron_SSCA2GL1A = NaN;
        end
        
    else %PT beads, only 2 values
       %Store center on SSC-H channel only
        temp_table.OD2centers = [m24_od2 m32_od2];
        temp_table.OD2_hv = hv_od2;
    
        temp_table.NoOD2centers = [m24_noOd2 m32_noOd2];
        temp_table.NoOD2_hv = hv_noOd2;
    end
    
    
    %store beadstat statistics
    beadstat_2026(b,:) = temp_table;
    
end


%want to add a saved plot of ssc values over time compared to average
%need to separate by hv, so we'll just use mode

notes = ['beads processed ' string(datetime)];

% save beadstat
save([outpath 'beadstat_2026'],'beadstat_2026', 'notes', 'beadtype', 'bead_ch_names', 'beadfiles2include')

disp('Result file saved:')
disp([outpath 'beadstat_2026'])


figure(100)
set(gcf,'Position', [100 200 1200 350])
clf
subplot(1,3,1)
%scatter(beadstat_2026.time, beadstat_2026.OD2centers(:,2), 12, beadstat_2026.hv(:, od2_chn), 'filled')
scatter(beadstat_2026.time, beadstat_2026.OD2centersA(:,2), 12, beadstat_2026.hv(:, strcmp(beadstat_2026.parname(1,2:end), od2_chn)), 'filled')
hold on
ylabel('GL1-A, OD2 SSC')
colorbar
yl1 = mean(beadstat_2026.OD2centersA(:,2),"omitmissing")*[.9 1.1]; yl2 = ylim;
line(xlim, yl1(1)*[1 1], 'linestyle', '--'), line(xlim, yl1(2)*[1 1], 'linestyle', '--') 
ylim([min([yl1(1) yl2(1)]) max([yl1(2), yl2(2)])])
grid
subplot(1,3,2)
%scatter(beadstat_2026.time, beadstat_2026.NoOD2centers(:,2), 12, beadstat_2026.hv(:, ssca_chn), 'filled')
scatter(beadstat_2026.time, beadstat_2026.NoOD2centersA(:,2), 12, beadstat_2026.hv(:, strcmp(beadstat_2026.parname(1,2:end), noOd2_chn)), 'filled')
ylabel('SSC-A, No OD')
colorbar 
hold on
yl1 = mean(beadstat_2026.NoOD2centersA(:,2),"omitmissing")*[.9 1.1]; yl2 = ylim;
line(xlim, yl1(1)*[1 1], 'linestyle', '--'), line(xlim, yl1(2)*[1 1], 'linestyle', '--') 
ylim([min([yl1(1) yl2(1)]) max([yl1(2), yl2(2)])])
title('1 \mum beads')
grid
subplot(1,3,3)
clear p
p(1) = scatter(beadstat_2026.time, beadstat_2026.NoOD2centersA(:,2)./beadstat_2026.OD2centersA(:,2), 12, double(beadstat_2026.QC_flag), 'filled');
hold on
p(2) = plot(beadstat_2026.time, beadstat_2026.mean1micron_SSCA2GL1A, '*-');
%p(3) = plot(beadstat_2026.time, beadstat_2026.mode1micron_SSCA2GL1A,'+-');
caxis([0 1])
ylabel('SSC-A/GL1-A')
colorbar
yl1 = mean(beadstat_2026.NoOD2centersA(:,2)./beadstat_2026.OD2centersA(:,2),"omitmissing")*[.9 1.1]; yl2 = ylim;
line(xlim, yl1(1)*[1 1], 'linestyle', '--'), line(xlim, yl1(2)*[1 1], 'linestyle', '--') 
ylim([min([yl1(1) yl2(1)]) max([yl1(2), yl2(2)])])
grid
colormap("winter")
legend(p,'bd mode ratio', 'mean ratio', 'mode ratio', 'fontsize',6)
print(figure(100), fullfile(beadfigpath, 'bead_timeseries.png'), '-dpng')


end

function c = centend2(fcsdat, class_vec, class, channel, fig_flag)

y = fcsdat{class_vec==class, channel};
m = mean(y); s = std(y);
ind = y<m+s*2 & y>m-s*2;
[n, edg] = histcounts(y(ind),50);
edg = edg(1:end-1) + (edg(2) - edg(1))/2;
%s = smooth(medfilt1(n',6),7);
%s = smooth(medfilt1(n',25),7);
if m>0
    sm = smooth(medfilt1(n',round(s/m*200)),7);
else
    sm = NaN;
end
[m, ind] = max(sm);
ind = (sm==m); %find all points equal to max
c = mean(edg(ind)); %take the midpoint of equal to max
if exist("fig_flag", 'var')
    plot(edg, n, '.-'), hold on
    plot(edg, sm, 'm-', 'linewidth', 2)
    line([c c], ylim, 'linewidth', 2)
    title(['bead ' num2str(class) ' A and H'])
    %disp(s/m*100)
    %pause  
end
end
