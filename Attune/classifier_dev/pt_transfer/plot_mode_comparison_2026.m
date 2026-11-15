function F = plot_mode_comparison_2026()
% function F = plot_mode_comparison_2026()
% July/Aug 2026 culture calibration: how the distributions and modes look
% when the mode is taken per linear bin (as in compile_attune_histograms /
% the 10 Sep compile_cc_histograms) or in log space.
%   Figure 1 (output/pt_transfer/modes_attune_2026.png): per culture row of
%     the 11 Sep fit, the scatter channel the fit uses (SSC-H, or GL1-H when
%     the culture's SSC mode is > 1x the 1 um bead), as a histogram of log10
%     signal; red dashed = AttuneMode (linear-bin histogram mode), black =
%     log-space KDE mode (bead_peak). Pooled rows pool their component runs.
%   Figure 2 (output/pt_transfer/modes_coulter_2026.png): pooled Coulter
%     counts per culture on a log diameter axis; grey = count per channel
%     (linear diameter bins), black = count x diameter (per log diameter);
%     red dashed = CCMode 10 Sep, black = CCMode 14 Nov (per log diameter).
%   Also refits the 11 Sep rows (log-space Attune modes, diatom kept) with
%   the 10 Sep and 14 Nov Coulter modes.
%
% Nov 2026, PT transfer / volume calibration

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
od = fullfile(dev, 'output', 'pt_transfer');
p26 = ['\' '\sosiknas1\Lab_data\Attune\size_calibration_July2026\'];
gd = fullfile(p26, 'Attune', 'Gated_FCS');
s = load(fullfile(p26, 'AttuneSummary_saved_11Sep2026.mat'), 'AttuneMode'); A = s.AttuneMode;
c0 = load(fullfile(p26, 'CCSummary_saved_10Sep2026.mat'), 'CCMode');
c1 = load(fullfile(p26, 'CCSummary_saved_14Nov2026.mat'), 'CCMode', 'ccdat', 'CCSizeCal_smoothhists');

% ---- rows and bead normalization as in size_calibration_fit_jul2026.m
ibd = contains(A.Row, 'bead');
MT = A(~ibd, :);
bt = removevars(A(contains(A.label, '1_um'), :), {'label'});
bd = groupsummary(bt, 'run_date', {'mean'});
[~, ib] = ismember(MT.run_date, bd.run_date); bm = bd(ib, :);
use = ~contains(MT.Row, {'Micromonas(2)' 'T_weissflogii_CCMP' 'Iso' '28-Jul-2026_I_gal'}) & ...
    ~(contains(MT.Row, 'Nano') & ~contains(MT.Row, 'pooled')) & ismember(MT.label, c1.CCMode.Row);
use = use | contains(MT.Row, {'Isochrysis_Isochrysis' '28-Jul-2026_I_galbana'}); % shown; averaged in the fit
rows = find(use);

% ---- figure 1: Attune distributions
f1 = figure('visible', 'off', 'position', [0 0 1600 900]);
t1 = tiledlayout(3, 5, 'TileSpacing', 'compact');
kS = nan(height(MT), 1); kG = kS;
for r = rows'
    stem = MT.Row{r};
    if contains(stem, 'pooled'), L = dir(fullfile(gd, [extractBefore(stem, '_pooled') '*.fcs'])); else, L = dir(fullfile(gd, [stem '.fcs'])); end
    L = L(~contains({L.name}, 'pooled'));
    S = []; G = [];
    for k = 1:numel(L)
        [d, h] = fca_readfcs(fullfile(gd, L(k).name)); n = {h.par.name};
        S = [S; d(:, strcmp(n, 'SSC-H'))]; G = [G; d(:, strcmp(n, 'GL1-H'))]; %#ok<AGROW>
    end
    kS(r) = 10^bead_peak(log10(S(S > 0))); kG(r) = 10^bead_peak(log10(G(G > 0)));
    useG = MT.("SSC-H")(r) / bm.("mean_SSC-H")(r) > 1;
    if useG, v = G; ch = "GL1-H"; km = kG(r); else, v = S; ch = "SSC-H"; km = kS(r); end
    nexttile; hold on
    lv = log10(v(v > 0));
    histogram(lv, linspace(prctile(lv, 0.5), prctile(lv, 99.5), 80), 'Normalization', 'pdf', ...
        'EdgeColor', 'none', 'FaceColor', [0.6 0.75 0.95])
    xline(log10(MT.(ch)(r)), 'r--', 'linewidth', 1.5)
    xline(log10(km), 'k-', 'linewidth', 1.5)
    title(sprintf('%s (%s)', strrep(MT.label(r), '_', ' '), ch), 'fontsize', 9)
    subtitle(sprintf('linear-bin %.0f | log %.0f (%+.2f log)', MT.(ch)(r), km, log10(km / MT.(ch)(r))), 'fontsize', 7)
    grid on
end
xlabel(t1, 'log10 scatter (raw)'); ylabel(t1, 'events per log10 unit')
title(t1, 'Attune culture distributions: red dashed = linear-bin histogram mode (AttuneMode), black = log-space mode')
exportgraphics(f1, fullfile(od, 'modes_attune_2026.png'), 'Resolution', 110)

% ---- figure 2: Coulter distributions
H = c1.CCSizeCal_smoothhists; fn = c1.ccdat.Properties.VariableNames;
f2 = figure('visible', 'off', 'position', [0 0 1600 900]);
t2 = tiledlayout(3, 5, 'TileSpacing', 'compact');
for i = 1:height(H)
    ind = find(contains(fn, H.filestr{i}));
    x = c1.ccdat.(fn{ind(end)}).diameter_micron; y = zeros(size(x));
    for k = ind, y = y + interp1(c1.ccdat.(fn{k}).diameter_micron, c1.ccdat.(fn{k}).count, x, 'linear', 0); end
    y(1:10) = NaN;
    nexttile; hold on
    yl = y .* x;
    plot(x, y / max(y), '.', 'color', [0.6 0.6 0.6], 'markersize', 5)
    plot(x, yl / max(yl), 'k.', 'markersize', 5)
    hs = H.histogram{i};
    plot(hs(:, 1), hs(:, 2) / max(hs(:, 2)), 'k-')
    m0 = c0.CCMode{H.Row{i}, 'diameter_mode'}; m1 = c1.CCMode{H.Row{i}, 'diameter_mode'};
    xline(m0, 'r--', 'linewidth', 1.5); xline(m1, 'k-', 'linewidth', 1.5)
    set(gca, 'xscale', 'log'); xlim([m1 / 2.2 m1 * 2.2]); ylim([0 1.1]); grid on
    title(strrep(H.Row{i}, '_', ' '), 'fontsize', 9)
    subtitle(sprintf('10 Sep %.3f | 14 Nov %.3f um (%+.3f log vol)', m0, m1, 3 * log10(m1 / m0)), 'fontsize', 7)
end
xlabel(t2, 'Coulter diameter (\mum, log axis)'); ylabel(t2, 'normalized count')
title(t2, 'Coulter pooled distributions: grey = per channel (linear bins), black = per log diameter; red dashed = 10 Sep mode, black = 14 Nov mode')
exportgraphics(f2, fullfile(od, 'modes_coulter_2026.png'), 'Resolution', 110)

% ---- refit: 11 Sep rows, log-space Attune modes, diatom kept, both Coulter mode sets
x = kS ./ bm.("mean_SSC-H"); g = kG ./ bm.("mean_GL1-H");
sw = MT.("SSC-H") ./ bm.("mean_SSC-H") > 1; x(sw) = g(sw); % same switch as the 11 Sep script
fitrows = use & ~contains(MT.Row, {'Isochrysis_Isochrysis' '28-Jul-2026_I_galbana'});
ig = contains(MT.Row, {'Isochrysis_Isochrysis' '28-Jul-2026_I_galbana'});
for c = {c0, c1}
    cc = c{1}.CCMode;
    vol = nan(height(MT), 1);
    [ok, j] = ismember(MT.label, cc.Row); vol(ok) = pi / 6 * cc.diameter_mode(j(ok)).^3;
    X = [log10(x(fitrows)); mean(log10(x(ig)))]; Y = [log10(vol(fitrows)); log10(mean(vol(ig)))];
    [p, st] = fit(X, Y, 'poly1'); ci = confint(p);
    fprintf('n=%d  p1 %.4f [%.3f %.3f]  p2 %.4f [%.3f %.3f]  rmse %.3f\n', numel(X), p.p1, ci(:, 1), p.p2, ci(:, 2), st.rmse)
end
% final calibration plot (14 Nov Coulter modes = last loop pass), as vol_cal_SSC-H_merge
lab = [MT.label(fitrows); "I galbana/Isochrysis"];
f3 = figure('visible', 'off', 'position', [0 0 760 700]); hold on
xx = linspace(-1.5, 2.2, 100)';
pb = predint(p, xx, 0.95, 'observation', 'off');
plot(X, Y, 'k.', 'markersize', 18, 'DisplayName', 'cultures (log-space modes)')
plot(xx, p(xx), 'b-', 'linewidth', 2, 'DisplayName', 'fit')
plot(xx, pb, 'b--', 'HandleVisibility', 'off'); plot(nan, nan, 'b--', 'DisplayName', 'prediction bounds 95%')
plot(xx, 1.1696 * xx + 1.1592, '-', 'color', [0.6 0.6 0.6], 'DisplayName', '11 Sep 2026 fit (linear-bin modes)')
text(X + 0.04, Y - 0.06, strrep(lab, '_', ' '), 'fontsize', 8)
text(-1.35, 3.6, sprintf('y = %.4f x + %.4f', p.p1, p.p2), 'color', 'r', 'fontsize', 12)
text(-1.35, 3.3, sprintf('r^2 = %.3f, rmse = %.3f, n = %d', st.rsquare, st.rmse, numel(X)), 'fontsize', 10)
axis([-1.5 2.2 -0.5 4]); axis square; grid on
xlabel('log10 [Attune SSC-H (merge), 1 \mum bead normalized]'); ylabel('log10 [Coulter Counter volume (\mum^3)]')
legend('location', 'southeast')
title('July/Aug 2026 SSC-H volume calibration (Attune\_SC2vol, 14 Nov 2026)')
exportgraphics(f3, fullfile(od, 'vol_cal_SSC-H_merge_14Nov2026.png'), 'Resolution', 120)
F = struct('kS', kS, 'kG', kG);
end
