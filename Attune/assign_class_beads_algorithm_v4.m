function [class, m05, m1, m6, beadstat, QC_flag] = assign_class_beads_algorithm_v4(fcsdat, fcshdr, plot_flag)
%function [class, m05, m1, m6, beadstat, QC_flag] = assign_class_beads_algorithm_v4(fcsdat, fcshdr, plot_flag)

% INPUTS
% fcsdat = attune data structure from fca_readfcs.m
% fcshdr = attune hdr structure from fca_readfcs.m
% plot_flag = logical to tell function to output plots of the data

% OUTPUTS
% class = vector of class values for each observation (0=junk, 1=0.5um, 2=1um, 3=6um)
% m05 = size of 0.5-um beads
% m1 = size of 1-um beads
% m6 = size of 6-um beads
% beadstat table
% QC_flag = quality control (1=good, 0=bad)

% NOTE: laser channels were adjusted after cruise TN368, so you may need
% change the channel references depending on which cruise you are doing

% Heidi M. Sosik, Woods Hole Oceanographic Institution, June 2020: 
% modified from assign_class_beads_algorithm_v3 (Kevin's version)
%   refine bead clustering with new distance metric
%   automatically set minpts scaled from number of data points above chl threshold
%   compute and save a range of bead stats
%   export bead plots with titles

startdate = datetime([fcshdr.date, ' ', fcshdr.starttime]);
% set up channels
if startdate < datetime('1-Aug-2019')
    ch = [12 15 18];
else
    %ch1 = 19; %SSC on GL1-H; %ch2 = 15; %Chl; %ch3 = 17; %GL3-H
    ch = [19 15 17];
end

% cluster analysis
%c = dbscan(log10(fcsdat(:, [ch1 ch2])), epsilon, minpts); %Kevin's line
% first cluster only stuff not offscale on chl
chl_noise_cut = 200;
c = -1*ones(size(fcsdat,1),1);
tempi = find(fcsdat(:,ch(2))>chl_noise_cut & fcsdat(:,ch(2))<1e6);
minpts = floor(length(tempi)/75);
ctemp = dbscan(log10(fcsdat(tempi,ch)+1), .005, minpts, 'distance', 'squaredeuclidean'); 
c(tempi) = ctemp;
nc = max(c);
meanc = NaN(length(ch), nc);
for ii = 1:nc
    meanc(:,ii) = mean(fcsdat(c==ii,ch));
end
[~,c1] = min(meanc(1,:)); %smallest on SSC (GL1), 0.5 micron
[~,c2] = min(meanc(3,:)); %smallest on GL3, 1 micron

% second cluster the stuff offscale on chl
tempi = find(fcsdat(:,ch(2))>=1e6);
minpts = max([10 floor(length(tempi)/15)]);
ctemp2 = dbscan(log10(fcsdat(tempi,ch([1,3]))+1), .003, minpts, 'distance', 'squaredeuclidean');
c(tempi(ctemp2>0)) = ctemp2(ctemp2>0)+max(ctemp);
nc = max(ctemp2);
meanc = NaN(length(ch), nc);
for ii = 1:nc
    meanc(:,ii) = mean(fcsdat(c==ii+max(ctemp),ch));
end
[~,c3] = min(meanc(1,:)); %smallest on SSC (GL1)
c3 = c3+max(ctemp);

% figure(199), clf
% subplot(1,2,1), scatter(fcsdat(:,ch(1)), fcsdat(:,ch(3)), 5, c, 'filled'), set(gca, 'yscale', 'log', 'xscale', 'log')
% subplot(1,2,2), scatter(fcsdat(:,ch(1)), fcsdat(:,ch(2)), 5, c, 'filled'), set(gca, 'yscale', 'log', 'xscale', 'log')

class = zeros(length(c), 1);
class(c==c1) = 1; %0.5 micron
class(c==c2) = 2; %1 micron
class(c==c3) = 3; %6 micron

beadstat = table;
parname = {fcshdr.par.name};
for ii = 1:3
    nstr = num2str(ii);
    beadstat.(['mean' nstr]) = mean(fcsdat(class==ii,:));
    beadstat.(['std' nstr]) = std(fcsdat(class==ii,:));
    beadstat.(['median' nstr]) = median(fcsdat(class==ii,:));
end

% central tendency of the bead clusters
m05 = centend(fcsdat, class, 1, ch(1));
m1 = centend(fcsdat, class, 2, ch(1));
m6 = centend(fcsdat, class, 3, ch(1));

beadstat.SSCcenter1 = m05;
beadstat.SSCcenter2 = m1;
beadstat.SSCcenter3 = m6;

% HEIDI: DOES THIS STILL WORK??
% quality control metrics
if m1<(m6/5) && m05<(m1/5)
    QC_flag = 1;
else
    QC_flag = 0;
end
QC_flag = logical(QC_flag);

% plot bead clusters for visual verification
if plot_flag
    figure(99), clf
    set(gcf,'Position', [100 100 1100 650])

    colormap(lines(4));
    subplot(2,2,1)
    scatter(fcsdat(:,ch(1)), fcsdat(:,ch(2)), 1, class, 'filled')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlim([10 1.2e6]); ylim([10 1.2e6]);
    xlabel('SSC-H'); ylabel('680/40')
    subplot(2,2,2)
    scatter(fcsdat(:,ch(1)), fcsdat(:,ch(3)), 1, class, 'filled')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlim([10 1.2e6]); ylim([10 1.2e6]);
    xlabel('SSC-H'); ylabel('620/15')
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
end

function c = centend(fcsdat, class_vec, class, channel)

bins = logspace(1, 6.1, 256);
[n, edg] = histcounts(fcsdat(class_vec==class, channel), bins);
edg = edg(1:end-1) + (edg(2) - edg(1))/2;
s = smooth(n, 5);
[~, ind] = max(s);
c = edg(ind);

end

end

