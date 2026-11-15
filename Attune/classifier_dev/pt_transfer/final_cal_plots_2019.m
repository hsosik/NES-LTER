function final_cal_plots_2019()
% function final_cal_plots_2019()
% Final volume calibrations for the pre-Aug 2019 cruises (one 488 nm SSC
% channel), drawn as they would be applied:
%   left:  OD2 on SSC (EN608, AR28B, AR29, EN617, AR31A, TN368): the linear
%          fit log10 vol = 1.4066 log10(SSC/bead) + 1.2109 (sizecal_2019)
%   right: no OD2 (AR34B, RB1904 at SSC 160 V; EN627 200 V): the same linear
%          fit applied to decompressed SSC (ssc_decompress). Cultures (200 V) at
%          their measured SSC (open) and decompressed (filled); the effective
%          curves vol(measured SSC) at 160 and 200 V end where the correction is
%          flagged (larger cells get no volume).
% Draws ../output/pt_transfer/final_cal_2019.png.
%
% Nov 2026, PT transfer / volume calibration

od = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'output', 'pt_transfer');
s = load(fullfile(od, 'sizecal_2019.mat'), 'R'); R = s.R;
pO = R.Properties.UserData.od2; st = R.Properties.UserData.od2_stats;
lin = @(x) pO.p1 * x + pO.p2;
xx = linspace(-1.6, 2.2, 200)';
pb = predint(pO, xx, 0.95, 'observation', 'off');
lab = strrep(R.culture, '_', ' ');

f = figure('visible', 'off', 'position', [0 0 1400 680]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');

nexttile; hold on
plot(xx, lin(xx), 'b-', 'linewidth', 2, 'DisplayName', sprintf('fit: %.4f x + %.4f', pO.p1, pO.p2))
plot(xx, pb, 'b--', 'HandleVisibility', 'off'); plot(nan, nan, 'b--', 'DisplayName', '95% prediction bounds')
plot(R.od2_norm, R.lvol, 'k.', 'markersize', 20, 'DisplayName', 'cultures, OD2 on SSC (260 V)')
text(R.od2_norm + 0.04, R.lvol - 0.07, lab, 'fontsize', 8)
text(-1.5, 3.75, sprintf('rmse %.3f, r^2 %.3f, n %d', st.rmse, st.rsquare, height(R)), 'fontsize', 10)
xlabel('log10 [SSC-H / 1 \mum bead]  (OD2 on SSC)'); ylabel('log10 Coulter volume (\mum^3)')
title('OD2 on SSC: EN608, AR28B, AR29, EN617, AR31A, TN368')
legend('location', 'southeast'); grid on; axis([-1.6 2.2 -0.5 4]); axis square

nexttile; hold on
[tc, fc] = ssc_decompress(R.neu_norm, 200);
plot(xx, lin(xx), '-', 'color', [0.6 0.6 1], 'linewidth', 1, 'DisplayName', 'linear fit vs decompressed SSC')
cols = [0.85 0.33 0.1; 0 0 0];
hvs = [160 200];
for k = 1:2
    sm = linspace(-1.6, 1.2, 600)';
    [t, fl] = ssc_decompress(sm, hvs(k));
    plot(sm(~fl), lin(t(~fl)), '-', 'color', cols(k, :), 'linewidth', 2.5, ...
        'DisplayName', sprintf('applied: vol vs measured SSC, %d V', hvs(k)))
    last = find(~fl, 1, 'last');
    plot(sm(last), lin(t(last)), 'v', 'color', cols(k, :), 'markerfacecolor', cols(k, :), 'markersize', 9, 'HandleVisibility', 'off')
    xline(sm(last), ':', sprintf('flag above (%d V)', hvs(k)), 'color', cols(k, :), 'HandleVisibility', 'off', ...
        'LabelVerticalAlignment', 'bottom', 'fontsize', 8)
end
plot(R.neu_norm, R.lvol, 'ko', 'markersize', 7, 'DisplayName', 'cultures, measured (no OD2, 200 V)')
q = ~fc;
plot(tc(q), R.lvol(q), 'k.', 'markersize', 20, 'DisplayName', 'same, decompressed')
plot([R.neu_norm(q) tc(q)]', [R.lvol(q) R.lvol(q)]', 'k-', 'HandleVisibility', 'off')
text(R.neu_norm + 0.04, R.lvol - 0.07, lab, 'fontsize', 8)
xlabel('log10 [SSC-H / 1 \mum bead]  (no OD2; measured, or decompressed for filled points)'); ylabel('log10 Coulter volume (\mum^3)')
title('No OD2: AR34B, RB1904 (160 V), EN627 (200 V, flagged)')
legend('location', 'southeast', 'fontsize', 8); grid on; axis([-1.6 2.2 -0.5 4]); axis square
title(tl, 'Pre-Aug 2019 cruises: one linear 488 nm calibration; decompression first where there is no OD2')
exportgraphics(f, fullfile(od, 'final_cal_2019.png'), 'Resolution', 110)
fprintf('flagged cultures at 200 V: %s\n', strjoin(R.culture(fc), ', '))
end
