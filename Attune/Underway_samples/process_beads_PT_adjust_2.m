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


function [] = process_beads_PT_adjust_2(basepath, FCSfileinfo, beadfiles2include, beadtype, OD2setting)

% make output directory
outpath = [basepath filesep 'bead_calibrated' filesep];
beadfigpath = [outpath filesep 'bead_plots_2021'];
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
beadstat_2021 = table;


%run through bead files in FCS directory (assumed to be one up from outpath)
for b = 1:length(beadlist)
    if b == 20 
        continue
    end
    clear temp_table
    disp(b)
    
    [fcsdat, fcshdr] = fca_readfcs([fpath beadlist{b}]);
    
    ssca_chn = strmatch('SSC-A', {fcshdr.par.name});
    gl1a_chn = strmatch('GL1-A', {fcshdr.par.name});

    bead_ch_names = {};
    if strcmp(OD2setting,'GL1')
        bead_ch_names{1} = 'GL1-A';
        sidescatter_ch = gl1a_chn;
    else %OD2 filter is on SSC or not on at all
        bead_ch_names{1} = 'SSC-A';
        sidescatter_ch = ssca_chn;
    end
    
    ch = nan(1,3);
    ch(1) = sidescatter_ch;
    ch(2) = strmatch('BL3-H', {fcshdr.par.name});  %chlorophyl channel
    bead_ch_names{2} = 'BL3-H';
    if fcshdr.date < datetime(2019,08,05) %3rd channel changes with cruise settings
        ch(3) = strmatch('GL2-H', {fcshdr.par.name});
        bead_ch_names{3} = 'GL2-H';
    else
        ch(3) = strmatch('GL3-H', {fcshdr.par.name});
        bead_ch_names{3} = 'GL3-H';
    end
    
    
    temp_table = table;
    temp_table.filename = beadlist(b);
    temp_table.time = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    temp_table.parname = {fcshdr.par.name};
    temp_table.OD2_filter = {OD2setting};
    temp_table.hv = cell2mat({fcshdr.par.hv});
    
    if strcmp(beadtype, 'FCB')
        n_clust = 3;
             
        % cluster analysis
        chl_noise_cut = 200;
        c = -1*ones(size(fcsdat,1),1); % c is a vector of cluster indeces, -1 indicates noise
        tempi = find(fcsdat(:,ch(2))>chl_noise_cut & fcsdat(:,ch(2))<1e6); %consider particles above threshhold
        minpts = floor(length(tempi)/75); %this is a parameter for clustering algorithm
        ctemp = dbscan(real(log10(fcsdat(tempi,ch)+1)), .005, minpts, 'distance', 'squaredeuclidean');
        c(tempi) = ctemp;
        nc = max(c); %this is number of cluster indexes
        meanc = NaN(length(ch), nc); %get mean signal on each channel for each cluster
        for iii = 1:nc
            meanc(:,iii) = mean(fcsdat(c==iii,ch));
        end
        clear iii
        [~,j] = sort(meanc(3,:));
        if (nc > 1 & meanc(1:2,j(1))./meanc(1:2,j(2)) > .8 & meanc(1:2,j(1))./meanc(1:2,j(2)) <1.5 & meanc(3,j(1))./meanc(3,j(2)) > 0.3 & meanc(3,j(1))./meanc(3,j(2)) < 1.5)
            c(c==j(1)) = j(2);
            for iii = 1:nc
                meanc(:,iii) = mean(fcsdat(c==iii,ch));
            end
        end
        clear iii
        [~,c2] = min(meanc(3,:)); %smallest on GL3, 1 micron
        tt = find(meanc(1,:) < meanc(1,c2));
        [~,c1] = min(meanc(2,tt)); %smallest on SSC (GL1), 0.5 micron
        c1 = tt(c1);
        
        
        % second cluster the stuff offscale on chlorophyll
        tempi = find(fcsdat(:,ch(2))>=1e6);
        if length(tempi) ~= 0
            minpts = max([10 floor(length(tempi)/15)]); %dbscan parameters
            ctemp2 = dbscan(real(log10(fcsdat(tempi,ch([1,3]))+1)), .003, minpts, 'distance', 'squaredeuclidean');
            c(tempi(ctemp2>0)) = ctemp2(ctemp2>0)+max(ctemp);
            nc = max(ctemp2);
            meanc = NaN(length(ch), nc);
            for ii = 1:nc
                meanc(:,ii) = mean(fcsdat(c==ii+max(ctemp),ch));
            end
        end
        [~,c3] = min(meanc(1,:)); %smallest on SSC (GL1)
        c3 = c3+max(ctemp);
        
        
        class = zeros(length(c), 1);
        class(c==c1) = 1; %0.5 micron
        class(c==c2) = 2; %1 micron
        class(c==c3) = 3; %6 micron
        
        
            %% get center values for each bead cluster on OD2 channel and not 
            %depends on OD2settings
            
            if strcmp(OD2setting,'GL1')
                m05_noOd2 = centend(fcsdat, class, 1, ssca_chn);
                m1_noOd2 = centend(fcsdat, class, 2, ssca_chn);
                m6_noOd2 = centend(fcsdat, class, 3, ssca_chn);
            
                m05_od2 = centend(fcsdat, class, 1, gl1a_chn);
                m1_od2 = centend(fcsdat, class, 2, gl1a_chn);
                m6_od2 = centend(fcsdat, class, 3, gl1a_chn);
                
                hv_od2 = fcshdr.par(gl1a_chn).hv;
                hv_noOd2 = fcshdr.par(ssca_chn).hv;
                
                od2_chn = gl1a_chn;
                noOd2_chn = ssca_chn; 
                
            elseif strcmp(OD2setting, 'SSC')
                m05_noOd2 = NaN;
                m1_noOd2 = NaN;
                m6_noOd2 = NaN; 
                
                m05_od2 = centend(fcsdat, class, 1, ssca_chn);
                m1_od2 = centend(fcsdat, class, 2, ssca_chn);
                m6_od2 = centend(fcsdat, class, 3, ssca_chn);
            
                hv_od2 = fcshdr.par(ssca_chn).hv;
                hv_noOd2 = NaN; 
                
                od2_chn = ssca_chn; 
                noOd2_chn = NaN; 
                
            else %OD2setting == 'none'
                m05_noOd2 = centend(fcsdat, class, 1, ssca_chn);
                m1_noOd2 = centend(fcsdat, class, 2, ssca_chn);
                m6_noOd2 = centend(fcsdat, class, 3, ssca_chn);
                
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
            clf
            set(gcf,'Position', [100 125 1100 650])
            
            
            figure(99), clf
            colormap(lines(4));
            subplot(2,2,1)
            scatter(fcsdat(:,ch(1)), fcsdat(:,ch(2)), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 1.2e6]); ylim([10 1.2e6]);
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{2})
            
            subplot(2,2,2)
            scatter(fcsdat(:,ch(1)), fcsdat(:,ch(3)), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 1.2e6]); ylim([10 1.2e6]);
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{3})
            
            
            if ~isnan(od2_chn)
                
            subplot(2,2,3); hold on;
            cmap = colormap(lines(4));
            bins = logspace(1, 6.1, 256);
            histogram(fcsdat(class==1, od2_chn), bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==2, od2_chn), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==3, od2_chn), bins, 'FaceColor', cmap(4,:), 'EdgeColor', 'none')
            yl = ylim;
            xlim([10 1.2e6])
            plot([m05_od2 m05_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m1_od2 m1_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m6_od2 m6_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            set(gca, 'XScale', 'log')
            legend('0.5 \mum', '1 \mum', '6 \mum')
            xlabel('OD2 SSC')
            
            end
            if ~isnan(noOd2_chn)
            
            subplot(2,2,4); hold on;
            cmap = colormap(lines(4));
            bins = logspace(1, 6.1, 256);
            histogram(fcsdat(class==1, noOd2_chn), bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==2, noOd2_chn), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==3, noOd2_chn), bins, 'FaceColor', cmap(4,:), 'EdgeColor', 'none')
            yl = ylim;
            xlim([10 1.2e6])
            plot([m05_noOd2 m05_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m1_noOd2 m1_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m6_noOd2 m6_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            set(gca, 'XScale', 'log')
            legend('0.5 \mum', '1 \mum', '6 \mum')
            xlabel('SSC No Filter')
            
            
            end
            
            
            subplot(2,2,1), title(beadlist(b), 'interpreter', 'none')
            subplot(2,2,3), title([datestr(temp_table.time) '; ' bead_ch_names{1} ' hv = ' num2str(fcshdr.par(ch(1)).hv) ' (' fcshdr.par(ch(1)).name ')'])
            print(figure(99), fullfile(beadfigpath, regexprep(beadlist{b}, '.fcs', '.png')), '-dpng')
            
               
        % HEIDI: DOES THIS STILL WORK??
        % quality control on SSC-A. Depends on OD2 filter setup
        % Without OD2 filter, we are really only capturing 1micron bead
        % cluster... 
        if m1_od2<(m6_od2/5) && m05_od2<(m1_od2/5)
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
            m24_noOd2 = centend(fcsdat, class, 1, ssca_chn);
            m32_noOd2 = centend(fcsdat, class, 2, ssca_chn);
            
            m24_od2 = centend(fcsdat, class, 1, gl1a_chn);
            m32_od2 = centend(fcsdat, class, 2, gl1a_chn);

             hv_od2 = fcshdr.par(gl1a_chn).hv;
             hv_noOd2 = fcshdr.par(ssca_chn).hv;
             
        elseif strcmp(OD2setting, 'SSC')        
            
             m24_noOd2 = NaN;
             m32_noOd2 = NaN;
                
            % central tendency of the bead clusters
            m24_od2 = centend(fcsdat, class, 1, ssca_chn);
            m32_od2 = centend(fcsdat, class, 2, ssca_chn);
            
             hv_od2 = fcshdr.par(ssca_chn).hv;
             hv_noOd2 = NaN; 
                
             od2_chn = ssca_chn; 
             noOd2_chn = NaN; 
        else %OD2setting == 'none'
                m24_noOd2 = centend(fcsdat, class, 1, ssca_chn);
                m32_noOd2 = centend(fcsdat, class, 2, ssca_chn);
                
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
            subplot(2,2,1)
            scatter(fcsdat(:,ch(1)), fcsdat(:,ch(2)), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 1.2e6]); ylim([10 1.2e6]);
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{2})
            
            subplot(2,2,2)
            scatter(fcsdat(:,ch(1)), fcsdat(:,ch(3)), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 2e6]); ylim([10 1.2e6]);
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{3})
            
            
            if ~isnan(od2_chn)
                
            subplot(2,2,3); hold on;
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
            
            subplot(2,2,4); hold on;
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
        temp_table.(['centend' nstr])= centend(fcsdat, class, iii, ch(1));
    end
    clear iii
    
    if strcmp(beadtype, 'FCB')
        %Store center on SSC-H channel only
        temp_table.OD2centers = [m05_od2 m1_od2 m6_od2];
        temp_table.OD2_hv = hv_od2;
    
        temp_table.NoOD2centers = [m05_noOd2 m1_noOd2 m6_noOd2];
        temp_table.NoOD2_hv = hv_noOd2;
    else %PT beads, only 2 values
       %Store center on SSC-H channel only
        temp_table.OD2centers = [m24_od2 m32_od2];
        temp_table.OD2_hv = hv_od2;
    
        temp_table.NoOD2centers = [m24_noOd2 m32_noOd2];
        temp_table.NoOD2_hv = hv_noOd2;
    end
    
    
    %store beadstat statistics
    beadstat_2021(b,:) = temp_table;
    
end


%want to add a saved plot of ssc values over time compared to average
%need to separate by hv, so we'll just use mode

notes = ['beads processed ' string(datetime)];

% save beadstat
save([outpath 'beadstat_2021'],'beadstat_2021', 'notes', 'beadtype', 'bead_ch_names', 'beadfiles2include')

disp('Result file saved:')
disp([outpath 'beadstat_2021'])


figure(100)
clf
subplot(1,2,1)
scatter(beadstat_2021.time, beadstat_2021.OD2centers(:,2), 12, beadstat_2021.hv(:, ssca_chn), 'filled')
ylabel('OD2 Signal')
colorbar
subplot(1,2,2)
scatter(beadstat_2021.time, beadstat_2021.NoOD2centers(:,2), 12, beadstat_2021.hv(:, ssca_chn), 'filled')
ylabel('No Filter Signal')
colorbar 
hold on
datetick
print(figure(100), fullfile(beadfigpath, 'bead_timeseries.png'), '-dpng')


end



function c = centend(fcsdat, class_vec, class, channel)

bins = logspace(1, 6.1, 1000);
[n, edg] = histcounts(fcsdat(class_vec==class, channel), bins);
edg = edg(1:end-1) + (edg(2) - edg(1))/2;
s = smooth(n', 5);
[~, ind] = max(s);
c = edg(ind);

end

