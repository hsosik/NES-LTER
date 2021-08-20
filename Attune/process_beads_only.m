%Function formatted to be modular bead processing called by
%process_wrapper_2021
%identifies clusters of beads, saves statistics, and creates bead_plots
%which are saved to output directory

%Inputs: 
% outpath is location of output directory used in wrapper

% bead_ch_names = names of cnannels to be used for classifying beads. May
% vary between cruises. Example: bead_ch_names = {'SSC-A', 'BL3-H', 'GL3-H'}

    %this replaced bead_channels input which 3x1 vector with indeces of channels
    % previously we had:
    %ssc_ch(ii) = [19 15 17]; %GL1-H, low sensitivity SSC, with OD2,
    %BL3-H, Chlorophyl,  then %GL3-H
    %if startdate < datenum('1-Aug-2019')   
    %    ssc_ch(ii) = [12 15 18]; %here its SSC-H, Chlorophyll, GL2-H
    % but want to allow for things to change

%FCSfileinfo is output of FCS_DateTimeList(fpath) or loaded from .mat file 

%beadfiles2include is the string that starts bead files for the relevant
%cruise, e.g. sometimes its {'Daily bead check'} or {fcb_bead'} 

%beadtype can be either 'PT' or 'FCB'? -- haven't made this work for PT
%beads yet. 


% output details and changes 

%always want to save SSC-H mean. Channel 12. 
%changed 7/8/21 so that this is true regardless of whether GL1 is used for
%clustering. We use SSC for the calibration and conversion to volume. 

% bead plots include dot plots for all bead files considered, with center
% values plotted for each of the 3 clusters. Also inlcudes bead_timeseries.png 
% which has SSC value of the 1 micron beads over the course of the
% cruise compared to the median, which is what will be used to calibrate.


function [] = process_beads_only(outpath, bead_ch_names, FCSfileinfo, beadfiles2include, beadtype)

% make output directory
beadfigpath = [outpath filesep 'bead_plots'];
if ~exist(beadfigpath, 'dir')
    mkdir(beadfigpath)
end

%now cut down to just focus on beads in our include
t = [];
for iii = 1:length(beadfiles2include)
    t = [t strmatch(beadfiles2include{iii}, FCSfileinfo.fcslist)];
end
beadlist = FCSfileinfo.fcslist(t);

%initialize some variables of interest
bead_ssch = NaN(length(beadlist),1);
bead_qc = bead_ssch;
hv_ssch = bead_ssch;
hv_ssca = bead_ssch; 
bead_time = NaT(length(beadlist), 1);

beadstat = table;
beadstat_temp = table;
fpath = [outpath filesep '..' filesep 'FCS' filesep];

%run through bead files in FCS directory (assumed to be one up from
%outpath)

for b = 1:length(beadlist)
    disp(b)
    [fcsdat, fcshdr] = fca_readfcs([fpath beadlist{b}]);

    if b == 1 %only need to do this once to check which channel numbers 
        ch = nan(1,3); 
        for i = 1:3
            ch(i) = strmatch([bead_ch_names{i}], {fcshdr.par.name});
        end       
    end
    
    
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
    
    temp_table = table;
    parname = {fcshdr.par.name};
    for iii = 1:3
        nstr = num2str(iii);
        temp_table.(['mean' nstr]) = mean(fcsdat(class==iii,:));
        temp_table.(['std' nstr]) = std(fcsdat(class==iii,:));
        temp_table.(['median' nstr]) = median(fcsdat(class==iii,:));
        temp_table.(['number' nstr]) = sum(class==iii);
    end
    clear iii
    
    % central tendency of the bead clusters
    m05 = centend(fcsdat, class, 1, ch(1));
    m1 = centend(fcsdat, class, 2, ch(1));
    m6 = centend(fcsdat, class, 3, ch(1));
    
    %Store center on SSC-H channel only 
    temp_table.SSCHcenter1 = centend(fcsdat, class, 1, 12);
    temp_table.SSCHcenter2 = centend(fcsdat, class, 2, 12);
    temp_table.SSCHcenter3 = centend(fcsdat, class, 3, 12);
    
    % HEIDI: DOES THIS STILL WORK??
    % quality control metrics
    if m1<(m6/5) && m05<(m1/5)
        QC_flag = 1;
    else
        QC_flag = 0;
    end
    QC_flag = logical(QC_flag);
    %1 is bad, 0 is good. 
    
    % plot bead clusters for visual verification
    figure(99)
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
    subplot(2,1,2); hold on;
    cmap = colormap(lines(4));
    bins = logspace(1, 6.1, 256);
    histogram(fcsdat(class==1, ch(1)), bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
    histogram(fcsdat(class==2, ch(1)), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
    histogram(fcsdat(class==3, ch(1)), bins, 'FaceColor', cmap(4,:), 'EdgeColor', 'none')
    yl = ylim;
    xlim([10 1.2e6])
    plot([m05 m05], [yl(1) yl(2)], '-k', 'LineWidth', 1)
    plot([m1 m1], [yl(1) yl(2)], '-k', 'LineWidth', 1)
    plot([m6 m6], [yl(1) yl(2)], '-k', 'LineWidth', 1)
    set(gca, 'XScale', 'log')
    legend('0.5 \mum', '1 \mum', '6 \mum')
    
    
    %store some statistics
    beadstat(b,:) = temp_table;
    beadstat_temp.hv(b,:) = {fcshdr.par.hv};
    bead_ssch(b) = centend(fcsdat, class, 2, 12); %SSCH channel always 7/8/21
    bead_ssca(b) = centend(fcsdat, class, 2, 3); %SSCA channel too 7/13/21
    bead_qc(b) = QC_flag;
    bead_time(b) = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    hv_ssch(b) = fcshdr.par(12).hv; %store ssc-h hv separately
    hv_ssca(b) = fcshdr.par(3).hv; %store ssc-a hv separately
    
    %Save figure
    subplot(2,2,1), title(beadlist(b), 'interpreter', 'none')
    subplot(2,2,2), title([datestr(bead_time(b)) '; ' bead_ch_names{1} ' hv = ' num2str(hv_ssch(b)) ' (' fcshdr.par(ch(1)).name ')'])
    print(figure(99), fullfile(beadfigpath, regexprep(beadlist{b}, '.fcs', '.png')), '-dpng')
    
    
end


for ii = 1:length(beadstat_temp.hv(:)), if length(beadstat_temp.hv{ii})==0, beadstat_temp.hv{ii} = [NaN]; end; end;
clear ii
beadstat.hv = cell2mat(beadstat_temp.hv); clear beadstat_temp
bead_qc = logical(bead_qc);
bead_med_ssch = NaN(length(unique(hv_ssch)), 2);
bead_med_ssch(:,1) = unique(hv_ssch);
bead_med_ssca = NaN(length(unique(hv_ssca)), 2);
bead_med_ssca(:,1) = unique(hv_ssca);
for ii = 1:size(bead_med_ssch,1) %heidi
    bead_med_ssch(ii,2) = nanmedian(bead_ssch(~bead_qc & hv_ssch==bead_med_ssch(ii,1)));
end
clear ii 
for ii = 1:size(bead_med_ssca,1) 
    bead_med_ssca(ii,2) = nanmedian(bead_ssca(~bead_qc & hv_ssca==bead_med_ssca(ii,1)));
end
clear ii


%want to add a saved plot of ssc values over time compared to average
%need to separate by hv, so we'll just use mode 
hval = bead_med_ssch(bead_med_ssch(:,1) == mode(hv_ssch), 2); 
aval = bead_med_ssca(bead_med_ssca(:,1) == mode(hv_ssca), 2); 
 

figure(99)
clf
subplot(1,2,1)
scatter(bead_time, bead_ssch, 12, hv_ssch,'filled')
hold on
plot(bead_time, ones(1,length(bead_time)).*hval)
ylabel('SSC-H')
datetick

subplot(1,2,2)
scatter(bead_time, bead_ssca, 12, hv_ssca, 'filled')
hold on 
plot(bead_time, ones(1,length(bead_time)).*aval) 
datetick
ylabel('SSC-A')
legend({'cluster average'; 'median over time'}, 'Location', 'southwest')
print(figure(99), fullfile(beadfigpath, 'bead_timeseries.png'), '-dpng')


% save beadstat
save([outpath 'beadstat'],'bead*', 'parname', 'hv_*')

disp(['Result file saved:'])
disp([outpath 'beadstat'])

end



function c = centend(fcsdat, class_vec, class, channel)

bins = logspace(1, 6.1, 1000);
[n, edg] = histcounts(fcsdat(class_vec==class, channel), bins);
edg = edg(1:end-1) + (edg(2) - edg(1))/2;
s = smooth(n', 5);
[~, ind] = max(s);
c = edg(ind);

end
