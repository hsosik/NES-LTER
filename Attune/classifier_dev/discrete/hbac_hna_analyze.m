function R = hbac_hna_analyze(tag)
% function R = hbac_hna_analyze(tag)
% HNA / LNA analysis of the hbac_hna_survey events (bacteria, auto class 9).
% Per run:
%   1-D SYBR (BL1) density (ksdensity, bandwidth cfg bw): the two strongest
%   modes at least min_sep apart, the density minimum between them (valley)
%   and the dip ratio (valley density / lower mode density; 1 = no dip)
%   2-D GMM in [SSC BL1], 1 vs 2 components: delta BIC, and the 2-group
%   separation along BL1 (mean difference / pooled SD)
% Draws output/discrete/hbac_hna_<tag>_ridge.png (BL1 density per run, runs in
% time order, valleys marked) and _cyto.png (example SSC-BL1 cytograms).
%
% Nov 2026, discrete-sample classifier development (HNA/LNA split)

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 1, tag = 'v1'; end
load(fullfile(od, 'hbac_hna_survey.mat'), 'R', 'E')
bw = 0.04; min_sep = 0.15; xg = (-2.3:0.01:0)';
nr = height(R);
[R.lo_mode, R.hi_mode, R.valley, R.dip, R.hna_frac, R.dbic, R.gsep, R.g_lo, R.g_hi, R.g_frac] = deal(nan(nr, 1));
D = nan(numel(xg), nr);
opts = statset('MaxIter', 300);
for r = 1:nr
    e = E(E.run == r & E.class == 9, :);
    if height(e) < 500, continue, end
    x = double(e.BL1);
    f = ksdensity(x, xg, 'bandwidth', bw); D(:, r) = f / max(f);
    [pk, loc] = findpeaks(f, 'MinPeakProminence', 0.02 * max(f));
    if numel(pk) >= 2
        [~, o] = sort(pk, 'descend'); p1 = loc(o(1));
        cand = loc(o(2:end)); cand = cand(abs(xg(cand) - xg(p1)) >= min_sep);
        if ~isempty(cand)
            p2 = cand(1); a = min(p1, p2); b = max(p1, p2);
            [fv, k] = min(f(a:b)); v = a + k - 1;
            R.lo_mode(r) = xg(a); R.hi_mode(r) = xg(b); R.valley(r) = xg(v);
            R.dip(r) = fv / min(f(a), f(b));
            R.hna_frac(r) = mean(x > xg(v));
        end
    end
    X = [double(e.SSC) x];
    if height(X) > 8000, X = X(randperm(height(X), 8000), :); end
    try
        g1 = fitgmdist(X, 1); g2 = fitgmdist(X, 2, 'Replicates', 3, 'Options', opts, 'RegularizationValue', 1e-4);
        R.dbic(r) = g1.BIC - g2.BIC;
        [~, lo] = min(g2.mu(:, 2)); hi = 3 - lo;
        R.g_lo(r) = g2.mu(lo, 2); R.g_hi(r) = g2.mu(hi, 2);
        R.gsep(r) = diff(g2.mu([lo hi], 2)) / sqrt(mean(squeeze(g2.Sigma(2, 2, :))));
        R.g_frac(r) = g2.ComponentProportion(hi);
    catch
    end
end
R.month = month(datetime(R.time, 'convertfrom', 'datenum'));
save(fullfile(od, ['hbac_hna_' tag '.mat']), 'R', 'D', 'xg')
fprintf('%d runs with >= 500 bacteria; two modes found in %d; dip < 0.7 in %d, < 0.9 in %d\n', ...
    sum(any(isfinite(D))), sum(isfinite(R.dip)), sum(R.dip < 0.7), sum(R.dip < 0.9))
fprintf('valley BL1: median %.2f, 10-90%% %.2f to %.2f; low mode %.2f, high mode %.2f; HNA fraction median %.2f\n', ...
    median(R.valley, 'omitnan'), prctile(R.valley, [10 90]), median(R.lo_mode, 'omitnan'), median(R.hi_mode, 'omitnan'), median(R.hna_frac, 'omitnan'))
G = groupsummary(R, 'cruise', {'median' 'min' 'max'}, {'valley' 'dip'});
disp(G)

% ridge image
[~, o] = sort(R.time);
fig = figure('visible', 'off', 'position', [50 50 1200 900]);
imagesc(1:nr, xg, D(:, o)), axis xy, colormap(parula), hold on
plot(1:nr, R.valley(o), 'w.', 'markersize', 8)
plot(1:nr, R.edge_BL1(o), 'r-')
c = R.cruise(o); b = [1; find(c(2:end) ~= c(1:end - 1)) + 1];
set(gca, 'xtick', b, 'xticklabel', c(b), 'xticklabelrotation', 90)
ylabel('BL1 (SYBR), bead-relative log10'), ylim([-2.1 -0.2])
title('bacteria (class 9) SYBR density per run (peak = 1), white = valley, red = BL1 trigger edge')
exportgraphics(fig, fullfile(od, ['hbac_hna_' tag '_ridge.png']), 'resolution', 100), close(fig)

% example cytograms: 20 runs spread over time
k = o(unique(round(linspace(1, nr, 20))));
fig = figure('visible', 'off', 'position', [20 20 1600 1250]);
for m = 1:numel(k)
    r = k(m); subplot(4, 5, m)
    e = E(E.run == r, :);
    sb = e.class == 9;
    plot(e.SSC(~sb), e.BL1(~sb), '.', 'color', [.8 .8 .8], 'markersize', 1), hold on
    scatter(e.SSC(sb), e.BL1(sb), 1, ksdens2(double(e.SSC(sb)), double(e.BL1(sb))), 'filled')
    if isfinite(R.valley(r)), yline(R.valley(r), 'r-'); end
    yline(R.edge_BL1(r), 'k:');
    xlim([-3.2 -0.9]), ylim([-2.1 -0.2]), grid on
    title(sprintf('%s %s %gm\ndip %.2f HNA %.2f', R.cruise(r), datestr(R.time(r), 'mmm yy'), R.depth(r), R.dip(r), R.hna_frac(r)), 'fontsize', 8)
end
exportgraphics(fig, fullfile(od, ['hbac_hna_' tag '_cyto.png']), 'resolution', 90), close(fig)
end

function d = ksdens2(x, y)
% point density on a 2-D histogram (for colouring)
[n, ex, ey] = histcounts2(x, y, 60);
n = imgaussfilt(n, 1);
ix = discretize(x, ex); iy = discretize(y, ey);
d = n(sub2ind(size(n), ix, iy));
d = log10(d + 1);
end
