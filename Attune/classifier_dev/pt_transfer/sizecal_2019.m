function R = sizecal_2019()
% function R = sizecal_2019()
% Volume calibrations for the cruises before the Aug 2019 layout change (one
% scatter channel, SSC 488/10), from the March/April 2019 culture runs:
%   OD2 on SSC (SSC 260 V): EN608, AR28B, AR29, EN617, AR31A, TN368
%   no OD2     (SSC 200 V): AR34B, RB1904 (EN627: flagged)
% As for 2026 (Attune_SC2vol, 14 Nov 2026), all modes are in log space:
%   - culture: KDE peak (bead_peak) of log10 SSC-H, events gated on
%     fluorescence only, identically in both runs (chl BL3 > 10^2.5; Syn also
%     PE (old-layout GL1-H) > 10^3.5) - see culture_od2_fluorgate
%   - 1 um FCB bead: find_beads on the same day's bead runs, KDE peak of SSC-H
%   - Coulter: mode of count x diameter (per log diameter), smoothed as in
%     compile_cc_histograms, above CCsize_threshold
% Syn = strain 1820i (9 Apr 2019; the Coulter "Syn" of the 2019 fit), with
% 9 Apr beads; the other cultures with 20 Mar beads.
% Fits log10 vol = p1 log10(SSC / 1 um bead) + p2 for each configuration
% (all 12 cultures, as the composite best fit); for no OD2 also only the
% cultures below the compression of SSC without OD2 (raw SSC-H < 10^4.8 at
% 200 V, culture_od2_check). Draws output/pt_transfer/sizecal_2019.png.
%
% Nov 2026, PT transfer / volume calibration

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
od = fullfile(dev, 'output', 'pt_transfer');
d0 = ['\' '\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\FCS\'];
s = load(['\' '\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\FitResults_saved_SSC-H_01Apr2026.mat'], ...
    'CCsizecal', 'CCsize_threshold', 'ModeTable');
CC = s.CCsizecal; thr = s.CCsize_threshold; MT19 = s.ModeTable;

% ---- 1 um bead peaks (SSC-H), per configuration and day
sig = bead_signature();
bf = struct('od2_20Mar', 'size_cal_OD2_Mar2019_Group_day0_Bead Solution.fcs', ...
    'od2_9Apr', 'size_cal_OD2_Mar2019_Group_day0_Bead_Solution2.fcs', ...
    'neu_20Mar', 'Size_cal_nuetral_Mar2019_Group_day0_Bead Solution.fcs', ...
    'neu_9Apr', 'Size_cal_nuetral_Mar2019_Group_day0_Bead Solution2.fcs');
bead = struct();
for k = fieldnames(bf)'
    [d, h] = fca_readfcs([d0 bf.(k{1})]);
    if startsWith(k{1}, 'od2'), conf = attune_config('TN368'); else, conf = attune_config('AR34B'); end
    [cl, pops] = find_beads(double(d), h, sig, conf);
    v = d(cl == 2, strcmp({h.par.name}, 'SSC-H'));
    if pops(2).found && numel(v) >= 100
        bead.(k{1}) = 10^bead_peak(log10(v(v > 0)));
    else
        bead.(k{1}) = NaN;
    end
    fprintf('1 um bead %-9s n=%4d QC=%d SSC-H %.0f (10^%.3f)\n', k{1}, numel(v), pops(2).QC, bead.(k{1}), log10(bead.(k{1})))
end

% no-OD2 bead run of 20 Mar not usable (5915 events, 1 um not found): use 9 Apr,
% as the 2019 analysis did (beads_1um_mode_NoOD)
if isnan(bead.neu_20Mar), bead.neu_20Mar = bead.neu_9Apr; fprintf('neu_20Mar: using the 9 Apr no-OD2 bead\n'), end
fprintf('OD2 bead 20 Mar vs 9 Apr: %+.3f log (day-to-day check)\n', log10(bead.od2_20Mar / bead.od2_9Apr))

% ---- cultures
cult = {'Syn' '1820i' '9Apr'; 'Micromonas' 'Micromonas' '20Mar'; 'Nano' 'Nano' '20Mar'; 'Pavlova' 'Pavlova' '20Mar';
    'I_galbana' 'I_galbana' '20Mar'; 'Prymn_pav' 'PRYMN_PAV' '20Mar'; 'Dun' 'Dun' '20Mar'; 'Crypto_spmc' 'CRYPTO_SPMC' '20Mar';
    'Htriq' 'HTriq' '20Mar'; 'H_akashiwo' 'H_Akashiwo_CCMP3374' '20Mar'; 'Ambopo' 'AMBOPO' '20Mar'; 'Pptp02' 'PPTP02' '20Mar'};
rows = {};
for i = 1:size(cult, 1)
    r = struct('culture', string(cult{i, 1}));
    for cfgname = ["neu" "od2"]
        if cfgname == "neu", pat = ['Size_cal_nuetral_Mar2019_Group_day0_' cult{i, 2} '*.fcs'];
        else, pat = ['size_cal_OD2_Mar2019_Group_day0_' cult{i, 2} '*.fcs']; end
        F = dir([d0 pat]);
        if strcmp(cult{i, 2}, '1820i') && cfgname == "neu", F = F(contains({F.name}, '_2')); end
        [d, h] = fca_readfcs([d0 F(1).name]);
        L = @(n) log10(max(d(:, strcmp({h.par.name}, n)), 1));
        g = L('BL3-H') > 2.5;
        if ismember(cult{i, 1}, {'Syn'}), g = g & L('GL1-H') > 3.5; end
        v = L('SSC-H'); v = v(g & v > 0);
        r.(cfgname + "_n") = numel(v);
        r.(cfgname + "_raw") = bead_peak(v);
        r.(cfgname + "_norm") = r.(cfgname + "_raw") - log10(bead.(cfgname + "_" + cult{i, 3}));
        r.(cfgname + "_hv") = h.par(strcmp({h.par.name}, 'SSC-H')).hv;
    end
    % Coulter: log-diameter mode (count x diameter), smoothed as compile_cc_histograms
    x = CC.(cult{i, 1}); dd = x(:, 1); y = x(:, 2) .* dd;
    if ismember(cult{i, 1}, {'Syn' 'Micromonas' 'Nano'}), ysm = smooth(medfilt1(y, 3), 8); else, ysm = smooth(medfilt1(y, 7), 20); end
    ysm(dd <= thr(i)) = NaN; ysm(1:2) = NaN;
    [~, im] = max(ysm);
    r.cc_d = dd(im); r.vol = pi / 6 * dd(im)^3;
    r.vol_old = MT19.cc_mode_vol(strcmp(MT19.files, cult{i, 1}));
    rows{end+1, 1} = r; %#ok<AGROW>
end
R = struct2table([rows{:}]');
fc = figure('visible', 'off', 'position', [0 0 1500 800]); tc = tiledlayout(3, 4, 'TileSpacing', 'compact');
for i = 1:size(cult, 1)
    x = CC.(cult{i, 1}); dd = x(:, 1); y = x(:, 2); yl = y .* dd;
    nexttile; hold on
    plot(dd, y / max(y(dd > thr(i))), '.', 'color', [0.6 0.6 0.6], 'markersize', 5)
    plot(dd, yl / max(yl(dd > thr(i))), 'k.', 'markersize', 5)
    xline((6 / pi * R.vol_old(i))^(1/3), 'r--', 'linewidth', 1.5); xline(R.cc_d(i), 'k-', 'linewidth', 1.5)
    xline(max(thr(i), dd(1)), 'b:')
    set(gca, 'xscale', 'log'); xlim([R.cc_d(i) / 2.5 R.cc_d(i) * 2.5]); ylim([0 1.2]); grid on
    title(strrep(cult{i, 1}, '_', ' ')); subtitle(sprintf('2019 fit %.3g | per log d %.3g um^3', R.vol_old(i), R.vol(i)), 'fontsize', 7)
end
title(tc, '2019 Coulter: grey per channel, black per log diameter; red dashed = 2019 ModeTable mode, black = log-diameter mode; blue = CCsize\_threshold')
exportgraphics(fc, fullfile(od, 'modes_coulter_2019.png'), 'Resolution', 110)
R.lvol = log10(R.vol);
disp(R(:, {'culture', 'neu_n', 'od2_n', 'neu_raw', 'od2_raw', 'neu_norm', 'od2_norm', 'vol_old', 'vol'}))

% ---- fits
lin = R.neu_raw < 4.8;
[pO, sO] = fit(R.od2_norm, R.lvol, 'poly1');
[pN, sN] = fit(R.neu_norm, R.lvol, 'poly1');
[pL, sL] = fit(R.neu_norm(lin), R.lvol(lin), 'poly1');
pr = @(nm, p, st) fprintf('%-36s p1 %.4f [%.3f %.3f]  p2 %.4f [%.3f %.3f]  rmse %.3f  r2 %.3f  n %d\n', nm, ...
    p.p1, confint(p) * [1 0]' .* [1; 1], p.p2, confint(p) * [0 1]', st.rmse, st.rsquare, st.dfe + 2);
pr('OD2 on SSC (260 V), all', pO, sO);
pr('no OD2 (200 V), all', pN, sN);
pr('no OD2 (200 V), raw SSC < 10^4.8', pL, sL);
R.Properties.UserData = struct('bead', bead, 'od2', pO, 'od2_stats', sO, 'neu_all', pN, 'neu_lin', pL, 'neu_lin_stats', sL);
save(fullfile(od, 'sizecal_2019.mat'), 'R')

% ---- plot
f = figure('visible', 'off', 'position', [0 0 1400 650]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
xx = linspace(-1.6, 2.2, 100)';
nexttile; hold on
plot(R.od2_norm, R.lvol, 'k.', 'markersize', 18, 'DisplayName', 'cultures, OD2 on SSC (260 V)')
plot(xx, pO(xx), 'b-', 'linewidth', 2, 'DisplayName', sprintf('fit: %.4f x + %.4f', pO.p1, pO.p2))
pb = predint(pO, xx, 0.95, 'observation', 'off'); plot(xx, pb, 'b--', 'HandleVisibility', 'off')
plot(xx, 1.1962 * xx + 1.0826, '-', 'color', [0.6 0.6 0.6], 'DisplayName', '2026 SSC/GL1 merge (Attune\_SC2vol)')
text(R.od2_norm + 0.04, R.lvol - 0.07, strrep(R.culture, '_', ' '), 'fontsize', 8)
title(sprintf('OD2 on SSC: rmse %.3f, r^2 %.3f, n %d', sO.rmse, sO.rsquare, height(R)))
xlabel('log10 [SSC-H / 1 \mum bead] (OD2 on SSC)'); ylabel('log10 Coulter volume (\mum^3)')
legend('location', 'northwest'); grid on; axis([-1.6 2.2 -0.5 4]); axis square
nexttile; hold on
plot(R.neu_norm(lin), R.lvol(lin), 'k.', 'markersize', 18, 'DisplayName', 'cultures, raw SSC < 10^{4.8}')
plot(R.neu_norm(~lin), R.lvol(~lin), 'ko', 'DisplayName', 'compressed (raw SSC \geq 10^{4.8})')
plot(xx, pL(xx), 'b-', 'linewidth', 2, 'DisplayName', sprintf('fit (linear range): %.4f x + %.4f', pL.p1, pL.p2))
pb = predint(pL, xx, 0.95, 'observation', 'off'); plot(xx, pb, 'b--', 'HandleVisibility', 'off')
plot(xx, pN(xx), 'b:', 'DisplayName', sprintf('fit (all): %.4f x + %.4f', pN.p1, pN.p2))
plot(xx, 1.1962 * xx + 1.0826, '-', 'color', [0.6 0.6 0.6], 'DisplayName', '2026 SSC/GL1 merge (Attune\_SC2vol)')
text(R.neu_norm + 0.04, R.lvol - 0.07, strrep(R.culture, '_', ' '), 'fontsize', 8)
xline(4.8 - log10(bead.neu_20Mar), ':', 'compression', 'HandleVisibility', 'off')
title(sprintf('No OD2: rmse %.3f (linear range), n %d', sL.rmse, sum(lin)))
xlabel('log10 [SSC-H / 1 \mum bead] (no OD2, 200 V)'); ylabel('log10 Coulter volume (\mum^3)')
legend('location', 'northwest'); grid on; axis([-1.6 2.2 -0.5 4]); axis square
title(tl, 'March/April 2019 cultures: volume calibrations for the pre-Aug 2019 cruises (log-space modes)')
exportgraphics(f, fullfile(od, 'sizecal_2019.png'), 'Resolution', 110)
end
