function C = ssc_compression_fit()
% function C = ssc_compression_fit()
% Compression of SSC-H without OD2 for large particles, from the March/April
% 2019 cultures run without OD2 (SSC 200 V) and with OD2 on SSC (260 V; linear
% for these signals). Each culture is gated on fluorescence only, identically
% in both runs, and its SSC-H quantiles (10-90%) are paired across the two
% runs (quantile-quantile); each run is normalized to its own 1 um FCB bead
% (find_beads; no-OD2 bead of 9 Apr, as sizecal_2019), so
%   t = log10 true SSC / bead   (OD2 run)
%   s = log10 measured SSC / bead (no-OD2 run)
% Pairs with raw OD2 SSC-H < 10^2.4 are left out (low-signal excess with OD2).
% Model, linear below a knee t0 and rolling off above it:
%   s = t - c/b * log10(1 + 10^(b (t - t0)))
% (slope 1 below t0, 1 - c above; b sets the sharpness). Fitted by least
% squares. Check: FCB beads of the 21 Jan 2020 filter swap (SSC 220 V, the
% at-sea voltage), each run normalized to its own 1 um bead.
% Saves ../output/pt_transfer/ssc_compression.mat (C.par = [t0 b c]) and
% ssc_compression.png. Use ssc_decompress to invert.
%
% Nov 2026, PT transfer / volume calibration

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
od = fullfile(dev, 'output', 'pt_transfer');
s19 = load(fullfile(od, 'sizecal_2019.mat'), 'R'); bead = s19.R.Properties.UserData.bead;
d0 = ['\' '\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\FCS\'];
cult = {'Syn' '1820i' '9Apr'; 'Micromonas' 'Micromonas' '20Mar'; 'Nano' 'Nano' '20Mar'; 'Pavlova' 'Pavlova' '20Mar';
    'I_galbana' 'I_galbana' '20Mar'; 'Prymn_pav' 'PRYMN_PAV' '20Mar'; 'Dun' 'Dun' '20Mar'; 'Crypto_spmc' 'CRYPTO_SPMC' '20Mar';
    'Htriq' 'HTriq' '20Mar'; 'H_akashiwo' 'H_Akashiwo_CCMP3374' '20Mar'; 'Ambopo' 'AMBOPO' '20Mar'; 'Pptp02' 'PPTP02' '20Mar'; 'SA2' 'SA2' '20Mar'};
q = 10:5:90;
rows = {};
for i = 1:size(cult, 1)
    v = struct();
    for cf = ["neu" "od2"]
        if cf == "neu", F = dir([d0 'Size_cal_nuetral_Mar2019_Group_day0_' cult{i, 2} '*.fcs']);
        else, F = dir([d0 'size_cal_OD2_Mar2019_Group_day0_' cult{i, 2} '*.fcs']); end
        if strcmp(cult{i, 2}, '1820i') && cf == "neu", F = F(contains({F.name}, '_2')); end
        [d, h] = fca_readfcs([d0 F(1).name]);
        L = @(n) log10(max(d(:, strcmp({h.par.name}, n)), 1));
        g = L('BL3-H') > 2.5;
        if strcmp(cult{i, 1}, 'Syn'), g = g & L('GL1-H') > 3.5; end
        x = L('SSC-H'); v.(cf) = x(g & x > 0);
    end
    qn = prctile(v.neu, q); qo = prctile(v.od2, q);
    for j = 1:numel(q)
        rows{end+1, 1} = struct('culture', string(cult{i, 1}), 'q', q(j), 'raw_neu', qn(j), 'raw_od2', qo(j), ...
            's', qn(j) - log10(bead.(['neu_' cult{i, 3}])), 't', qo(j) - log10(bead.(['od2_' cult{i, 3}]))); %#ok<AGROW>
    end
end
P = struct2table([rows{:}]');
use = P.raw_od2 >= 2.4;
model = @(p, t) t - p(3) / p(2) * log10(1 + 10.^(p(2) * (t - p(1))));
obj = @(p) sum((P.s(use) - model(p, P.t(use))).^2);
p = fminsearch(obj, [0.6 4 0.8], optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4));
P.resid = P.s - model(p, P.t);
rmse = sqrt(mean(P.resid(use).^2));
fprintf('compression fit: knee t0 %.3f (x%.2f bead), sharpness b %.2f, slope above knee %.3f, rmse %.3f log (n=%d)\n', ...
    p(1), 10^p(1), p(2), 1 - p(3), rmse, sum(use))
G = groupsummary(P(use, :), 'culture', 'median', {'t', 'resid'}); disp(sortrows(G, 'median_t'))

% check: 21 Jan 2020 filter-swap beads at SSC 220 V (culture_od2_check style peaks, find_beads labels)
fsw = ['\' '\sosiknas1\lab_data\Attune\ATX_FSC_misc_experiments\bead_hv_experiment2\filter_swap_experiment\'];
[a, ha] = fca_readfcs([fsw 'filter experiment neutral.fcs']); [b2, hb] = fca_readfcs([fsw 'filter experiment od2.fcs']);
sig = bead_signature();
[cla] = find_beads(double(a), ha, sig, 'GL1');
pk = @(d, h, c) bead_peak(log10(max(d(c, strcmp({h.par.name}, 'SSC-H')), 1)));
% OD2 run: label by fluorescence clusters (find_beads gates on SSC, which OD2 moves); use the
% populations' fluorescence signatures: 1 um bright BL1, 6 um saturating BL3 with dim BL1
Lb = @(n) log10(max(b2(:, strcmp({hb.par.name}, n)), 1));
one = Lb('BL1-H') > 4 & Lb('BL1-H') < 4.6 & Lb('BL3-H') < 3.8 & Lb('FSC-H') < 4.6;
six = Lb('BL3-H') > 5.8 & Lb('BL1-H') < 3.2;
s6 = pk(a, ha, cla == 3) - pk(a, ha, cla == 2);
t6 = bead_peak(Lb('SSC-H')) * 0 + bead_peak(Lb('SSC-H')); %#ok<NASGU>
t6 = bead_peak(log10(max(b2(six, strcmp({hb.par.name}, 'SSC-H')), 1))) - bead_peak(log10(max(b2(one, strcmp({hb.par.name}, 'SSC-H')), 1)));
fprintf('2020 filter swap, SSC 220 V, 6 um bead: true %.3f, measured %.3f, model at 200 V %.3f\n', t6, s6, model(p, t6))

C = struct('par', p, 'rmse', rmse, 'pairs', P, 'model', 'log10 meas = t - c/b log10(1 + 10^(b (t - t0)))', ...
    'check_6um_220V', [t6 s6 model(p, t6)]);
save(fullfile(od, 'ssc_compression.mat'), 'C')

f = figure('visible', 'off', 'position', [0 0 1200 560]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on
cs = unique(P.culture); col = [lines(7); 0.5 0.5 0.5; 0 0 0; 0.6 0.4 0.2; 0.9 0.5 0.8; 0.3 0.7 0.3; 0.7 0.7 0.2];
for i = 1:numel(cs)
    k = P.culture == cs(i) & use;
    plot(P.t(k), P.s(k), '.-', 'color', col(i, :), 'DisplayName', strrep(cs(i), '_', ' '))
end
tt = linspace(-1.5, 2.2, 200);
plot(tt, model(p, tt), 'k-', 'linewidth', 2.5, 'DisplayName', 'fit')
plot(tt, tt, 'k:', 'DisplayName', '1:1')
plot(t6, s6, 'kp', 'markersize', 14, 'markerfacecolor', 'y', 'DisplayName', '6 \mum bead, 2020, 220 V')
xlabel('true SSC (with OD2), log10 bead units'); ylabel('SSC without OD2 (200 V), log10 bead units')
legend('location', 'northwest', 'fontsize', 7); grid on; axis square; axis([-1.5 2.2 -1.5 2.2])
title(sprintf('knee %.2f (x%.1f bead), slope above %.2f, rmse %.3f', p(1), 10^p(1), 1 - p(3), rmse))
nexttile; hold on
for i = 1:numel(cs)
    k = P.culture == cs(i) & use;
    plot(P.t(k), P.resid(k), '.-', 'color', col(i, :))
end
yline(0, 'k'); grid on; axis square
xlabel('true SSC, log10 bead units'); ylabel('residual (log10)')
title('residuals by culture (quantiles 10-90%)')
title(tl, 'SSC-H compression without OD2: same cultures with and without OD2 (March/April 2019)')
exportgraphics(f, fullfile(od, 'ssc_compression.png'), 'Resolution', 110)
end
