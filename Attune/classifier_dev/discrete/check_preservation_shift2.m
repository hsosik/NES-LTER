% Is the preserved-vs-live size compression within cruises or only between them?
% And do Syn and Euk shifts co-vary per bottle / per cruise (a shared bead or
% calibration offset would move both)? Nov 2026, preservation-correction analysis
here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
P = readtable(fullfile(od, 'preservation_correction_pairs_disc.csv'), 'TextType', 'string');
R = table;
for p = ["Syn" "Euk"]
    x = P(P.pop == p, :);
    [g, cr] = findgroups(x.cruise);
    ml = splitapply(@median, x.live_size, g); mp = splitapply(@median, x.pres_size, g);
    lw = fitlm(x.live_size - ml(g), x.pres_size - mp(g), 'RobustOpts', 'on'); ciw = coefCI(lw);
    lb = fitlm(ml, mp, 'RobustOpts', 'on'); cib = coefCI(lb);
    R = [R; table(p, lw.Coefficients.Estimate(2), ciw(2,1), ciw(2,2), std(x.live_size - ml(g)), ...
        lb.Coefficients.Estimate(2), cib(2,1), cib(2,2), std(ml), 'VariableNames', ...
        {'pop' 'within_slope' 'w_lo' 'w_hi' 'within_live_sd' 'between_slope' 'b_lo' 'b_hi' 'between_live_sd'})]; %#ok<AGROW>
end
disp(R)
writetable(R, fullfile(od, 'preservation_correction_within_between.csv'))
% Syn vs Euk shift, same bottle
s = P(P.pop == "Syn", {'cruise' 'cast' 'niskin' 'd_size' 'live_size'}); e = P(P.pop == "Euk", {'cruise' 'cast' 'niskin' 'd_size' 'live_size'});
J = innerjoin(s, e, 'Keys', {'cruise' 'cast' 'niskin'});
[r, pv] = corr(J.d_size_s, J.d_size_e, 'type', 'Spearman');
fprintf('Syn vs Euk shift, same bottle: rho %.2f (p %.2g, n %d)\n', r, pv, height(J))
C = groupsummary(J, 'cruise', 'median', {'d_size_s' 'd_size_e' 'live_size_s' 'live_size_e'});
[r, pv] = corr(C.median_d_size_s, C.median_d_size_e, 'type', 'Spearman');
fprintf('Syn vs Euk shift, cruise medians: rho %.2f (p %.2g, n %d)\n', r, pv, height(C))
[r, pv] = corr(C.median_live_size_s, C.median_live_size_e, 'type', 'Spearman');
fprintf('Syn vs Euk live size, cruise medians: rho %.2f (p %.2g)\n', r, pv)
disp(sortrows(C, 'median_d_size_s'))

% Correction forms: residual bottle-level error (preserved corrected - live)
%   none, constant (cruise-median offset), linear in preserved size (live ~ pres)
% and within-population compression (p10-p90 spread, preserved / live)
p1 = 1.1962;
F = table; col = lines(5);
fig = figure('visible', 'off', 'position', [50 50 1300 420]);
tiledlayout(1, 3, 'TileSpacing', 'compact')
for p = ["Syn" "Euk" "Pro"]
    x = P(P.pop == p, :);
    C = groupsummary(x, 'cruise', 'median', 'd_size');
    off = median(C.median_d_size);
    li = fitlm(x.pres_size, x.live_size, 'RobustOpts', 'on'); %prediction of live from preserved
    r0 = x.d_size; r1 = x.d_size - off; r2 = -(predict(li, x.pres_size) - x.live_size);
    rsd = @(r) 1.4826 * mad(r, 1);
    spread = median((x.pres_size_p90 - x.pres_size_p10) ./ (x.live_size_p90 - x.live_size_p10));
    F = [F; table(p, height(x), off, p1 * off, 10^(-p1 * off), median(r0), rsd(r0), rsd(r1), rsd(r2), ...
        li.Coefficients.Estimate(1), li.Coefficients.Estimate(2), spread, 'VariableNames', ...
        {'pop' 'n' 'offset_size' 'offset_logvol' 'vol_correction_factor' 'bias_none' 'rsd_none' 'rsd_constant' 'rsd_linear' ...
        'live_on_pres_a' 'live_on_pres_b' 'within_pop_spread_ratio'})]; %#ok<AGROW>
    nexttile, hold on
    [g, cr] = findgroups(x.cruise);
    scatter(x.live_size, x.pres_size, 14, g, 'filled', 'MarkerFaceAlpha', 0.7)
    colormap(gca, turbo(numel(cr)))
    ax = [min([x.live_size; x.pres_size]) max([x.live_size; x.pres_size])] + [-0.05 0.05];
    plot(ax, ax, 'k:', 'displayname', '1:1')
    plot(ax, ax + off, 'k--')
    lp = fitlm(x.live_size, x.pres_size, 'RobustOpts', 'on');
    plot(ax, predict(lp, ax'), 'r-')
    axis([ax ax]), axis square, grid on
    xlabel('live median size (log_{10} rel. 1 \mum bead)'), ylabel('preserved median size')
    title(sprintf('%s: n %d bottles, %d cruises; offset %.3f, slope %.2f', p, height(x), numel(cr), off, lp.Coefficients.Estimate(2)), 'fontsize', 9)
end
exportgraphics(fig, fullfile(od, 'preservation_correction_pres_vs_live.png'), 'resolution', 90)
close(fig)
disp(F)
writetable(F, fullfile(od, 'preservation_correction_forms.csv'))
