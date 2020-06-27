function [class, m05, m1, m6, QC_flag] = assign_class_beads_algorithm_v3(fcsdat, epsilon, minpts, startdate, plot_flag)

% INPUTS
% fcsdat = attune data structure from fca_readfcs.m
% epsilon = cluster neighborhood value used in DBSCAN (log value)
% minpts = minimum number of points in a cluster
% startdate = used to determine which channels to use
% plot_flag = logical to tell function to output plots of the data

% OUTPUTS
% class = vector of class values for each observation (0=junk, 1=0.5um, 2=1um, 3=6um)
% m05 = size of 0.5-um beads
% m1 = size of 1-um beads
% m6 = size of 6-um beads
% QC_flag = quality control (1=good, 0=bad)

% NOTE: laser channels were adjusted after cruise TN368, so you may need
% change the channel references depending on which cruise you are doing

% set up channels
if datetime(startdate) < datetime('1-Aug-2019')
    ch1 = 12;
    ch2 = 15;
    ch3 = 18;
else
    ch1 = 19;
    ch2 = 15;
    ch3 = 17;
end

% cluster analysis
c = dbscan(log10(fcsdat(:, [ch1 ch2])), epsilon, minpts);

% identify bead clusters
pd = NaN(length(unique(c))-1, 1);
for i=1:length(pd)
    pd(i) = sqrt(mean(log10(fcsdat(c==i, ch1)))^2 + mean(log10(fcsdat(c==i, ch2)))^2);
end
[~,pdind] = sort(pd, 'descend');

bssc = NaN(3,1);
for j=1:3
    bssc(j) = mean(fcsdat(c==pdind(j), ch1));
end
class = zeros(length(c), 1);
class(c==pdind(bssc==min(bssc))) = 1;
class(c==pdind(bssc~=max(bssc) & bssc~=min(bssc))) = 2;
class(c==pdind(bssc==max(bssc))) = 3;

% central tendency of the bead clusters
m05 = centend(fcsdat, class, 1, ch1);
m1 = centend(fcsdat, class, 2, ch1);
m6 = centend(fcsdat, class, 3, ch1);

% quality control metrics
if m1<(m6/5) && m05<(m1/5)
    QC_flag = 1;
else
    QC_flag = 0;
end
QC_flag = logical(QC_flag);

% plot bead clusters for visual verification
if plot_flag
figure('Position', [100 100 1100 650])

colormap(lines(4));
subplot(2,2,1)
scatter(fcsdat(:,ch1), fcsdat(:,ch2), 1, class, 'filled')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlim([10 1e7]); ylim([10 1e7]);
xlabel('SSC-H'); ylabel('680/40')
subplot(2,2,2)
scatter(fcsdat(:,ch1), fcsdat(:,ch3), 1, class, 'filled')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlim([10 1e7]); ylim([10 1e7]);
xlabel('SSC-H'); ylabel('620/15')
subplot(2,1,2); hold on;
cmap = colormap(lines(4));
bins = logspace(1, 6.1, 256);
histogram(fcsdat(class==1, ch1), bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
histogram(fcsdat(class==2, ch1), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
histogram(fcsdat(class==3, ch1), bins, 'FaceColor', cmap(4,:), 'EdgeColor', 'none')
yl = ylim;
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

