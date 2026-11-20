% Checks on analyze_preservation_shift output (Nov 2026, preservation-correction analysis):
%  1. size dependence without regression to the mean: d vs the mean of preserved
%     and live size (d vs live alone is biased negative by noise in live)
%  2. storage time and season at the cruise level (cruise medians)
%  3. robust within-cruise spread
%  4. preserved vs live size slope (robust): pres = a + b*live
here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
P = readtable(fullfile(od, 'preservation_correction_pairs_disc.csv'), 'TextType', 'string');
P.avg_size = (P.pres_size + P.live_size) / 2;
R = table;
for p = ["Syn" "Euk" "crypto_dim" "crypto_bright" "Pro"]
    x = P(P.pop == p, :);
    lm = fitlm(x.avg_size, x.d_size, 'RobustOpts', 'on');
    lp = fitlm(x.live_size, x.pres_size, 'RobustOpts', 'on'); %live noise is small (live null), so OLS slope is ~unbiased
    ci = coefCI(lp); mlive = median(x.live_size);
    dfit = lp.Coefficients.Estimate(1) + (lp.Coefficients.Estimate(2) - 1) * prctile(x.live_size, [10 50 90]);
    C = groupsummary(x, 'cruise', 'median', {'d_size' 'storage_days'});
    C.season = arrayfun(@(c) mode(categorical(x.season(x.cruise == c))), C.cruise);
    [rs, ps] = corr(C.median_d_size, log10(C.median_storage_days), 'type', 'Spearman', 'rows', 'complete');
    % within cruise: storage spread is small, so pool deviations from the cruise median
    [~, ic] = ismember(x.cruise, C.cruise);
    dev = x.d_size - C.median_d_size(ic); sdev = x.storage_days - C.median_storage_days(ic);
    [rw, pw] = corr(dev, sdev, 'type', 'Spearman', 'rows', 'complete');
    ps_season = NaN;
    if numel(unique(C.season)) > 1, ps_season = kruskalwallis(C.median_d_size, C.season, 'off'); end
    R = [R; table(p, height(x), lm.Coefficients.Estimate(2), lm.Coefficients.pValue(2), rs, ps, rw, pw, ps_season, ...
        1.4826 * mad(dev, 1), lp.Coefficients.Estimate(2), ci(2, 1), ci(2, 2), mlive, dfit(1), dfit(2), dfit(3), 'VariableNames', {'pop' 'n' 'slope_vs_avgsize' 'p_avgsize' 'rho_cruise_storage' 'p_cruise_storage' ...
        'rho_within_storage' 'p_within_storage' 'p_season_cruise' 'within_cruise_robust_sd' ...
        'slope_pres_on_live' 'slope_ci_lo' 'slope_ci_hi' 'median_live_size' 'dfit_live_p10' 'dfit_live_p50' 'dfit_live_p90'})]; %#ok<AGROW>
    disp(p), disp(groupsummary(C, 'season', 'median', 'median_d_size'))
end
disp(R)
writetable(R, fullfile(od, 'preservation_correction_checks.csv'))

