function T = culture_pigment_ratios()
% function T = culture_pigment_ratios()
% Can pigment groups be told apart in the 2026 cultures (standard layout) by
%   (1) chl fluorescence excited at 532 vs 488 nm: log10(GL4-H / BL3-H)
%       (both 680/40; a two-point excitation spectrum), each brought to the
%       reference voltages (BL3 340, GL4 260 V; cfg.k)
%   (2) side scatter at 532 vs 488 nm: log10(GL1/bead) - log10(SSC/bead), with
%       SSC decompressed (ssc_decompress at the file's SSC voltage; events
%       beyond the invertible range dropped)
% per culture: medians over the gated events. Groups: cyanobacteria (Syn),
% green, haptophyte, cryptophyte, diatom, dinoflagellate. Mie (classifier_dev/
% mie_scatter) predicts scatter ratio green/haptophyte > diatom/dino >
% cryptophyte for large cells, ~0 for small cells.
% Draws ../output/pt_transfer/culture_pigment_ratios.png.
%
% Nov 2026, PT transfer / pigment ratios

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
od = fullfile(dev, 'output', 'pt_transfer');
cfg = classifier_config();
p26 = ['\' '\sosiknas1\Lab_data\Attune\size_calibration_July2026\'];
gd = fullfile(p26, 'Attune', 'Gated_FCS');
s = load(fullfile(p26, 'AttuneSummary_saved_11Sep2026.mat'), 'AttuneMode'); A = s.AttuneMode; A.label = string(A.label);
s = load(fullfile(p26, 'CCSummary_saved_14Nov2026.mat'), 'CCMode'); CC = s.CCMode;
grp = struct('Syn7335', "cyanobacteria", 'Syn8113', "cyanobacteria", 'Syn_WH8018', "cyanobacteria", ...
    'Micromonas', "green", 'C2Micromonas', "green", 'Nano', "green", 'Dun', "green", ...
    'I_galbana', "haptophyte", 'Isochrysis', "haptophyte", 'Pavlova', "haptophyte", ...
    'Crypto', "cryptophyte", 'Rhodomonas', "cryptophyte", 'T_weissflogii', "diatom", ...
    'Heterocapsa', "dinoflagellate", 'SE62', "coccolithophore (SE62)");
rows = {};
for i = 1:height(A)
    lab = A.label(i); stem = A.Properties.RowNames{i};
    if ~isfield(grp, lab) || contains(stem, 'pooled'), continue, end
    fn = fullfile(gd, [stem '.fcs']);
    if ~isfile(fn), continue, end
    b = A(A.run_date == A.run_date(i) & A.label == "1_um", :);
    bS = median(b.("SSC-H")); bG = median(b.("GL1-H"));
    [d, h] = fca_readfcs(fn); n = {h.par.name};
    col = @(c) find(strcmp(n, c), 1);
    hv = @(c) h.par(col(c)).hv;
    v = @(c) d(:, col(c));
    ok = v('BL3-H') > 10 & v('GL4-H') > 10 & v('SSC-H') > 10 & v('GL1-H') > 10;
    bl3 = log10(v('BL3-H')) - cfg.k.BL3 * log10(hv('BL3-H') / 340);
    gl4 = log10(v('GL4-H')) - cfg.k.GL4 * log10(hv('GL4-H') / 260);
    sm = log10(v('SSC-H') / bS);
    [st, fl] = ssc_decompress(sm, hv('SSC-H'));
    g1 = log10(v('GL1-H') / bG);
    okS = ok & ~fl;
    r = struct('culture', lab, 'file', string(stem), 'group', grp.(lab), 'n', sum(ok), ...
        'frac_ssc_flagged', mean(fl(ok)), 'gl1_norm', median(g1(ok)), ...
        'chl_ex532_488', median(gl4(ok) - bl3(ok)), ...
        'scat_532_488', median(g1(okS) - st(okS)), 'scat_532_488_raw', median(g1(ok) - sm(ok)), ...
        'vol', NaN);
    if ismember(lab, CC.Properties.RowNames), r.vol = pi / 6 * CC{lab, 'diameter_mode'}^3; end
    rows{end+1, 1} = r; %#ok<AGROW>
end
T = struct2table([rows{:}]');
T = sortrows(T, {'group', 'vol'});
disp(T(:, {'culture', 'group', 'n', 'vol', 'gl1_norm', 'chl_ex532_488', 'scat_532_488', 'scat_532_488_raw', 'frac_ssc_flagged'}))
save(fullfile(od, 'culture_pigment_ratios.mat'), 'T')

gs = ["cyanobacteria" "green" "haptophyte" "cryptophyte" "diatom" "dinoflagellate" "coccolithophore (SE62)"];
col = [0.5 0.5 0.5; 0.2 0.65 0.2; 0.9 0.6 0; 0.8 0.1 0.5; 0.1 0.4 0.8; 0.6 0.3 0.1; 0 0 0];
f = figure('visible', 'off', 'position', [0 0 1500 500]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
lv = log10(T.vol); lv(isnan(lv)) = -0.2;
for pnl = 1:3
    nexttile; hold on
    for k = 1:numel(gs)
        q = T.group == gs(k);
        if ~any(q), continue, end
        switch pnl
            case 1, x = lv(q); y = T.chl_ex532_488(q);
            case 2, x = lv(q); y = T.scat_532_488(q);
            case 3, x = T.chl_ex532_488(q); y = T.scat_532_488(q);
        end
        plot(x, y, 'o', 'markersize', 8, 'color', col(k, :), 'markerfacecolor', col(k, :), 'DisplayName', gs(k))
        if pnl == 3, text(x + 0.01, y, strrep(T.culture(q), '_', ' '), 'fontsize', 7), else, text(x + 0.04, y, strrep(T.culture(q), '_', ' '), 'fontsize', 7), end
    end
    grid on
    switch pnl
        case 1, xlabel('log10 Coulter volume (\mum^3)'); ylabel('log10 chl fluorescence, ex 532 / ex 488 (GL4/BL3)'); title('chl excitation ratio vs size')
        case 2, xlabel('log10 Coulter volume (\mum^3)'); ylabel('log10 scatter 532 / 488 (GL1 / decompressed SSC), bead units'); title('scatter ratio vs size'); yline(0, 'k:', 'HandleVisibility', 'off')
        case 3, xlabel('log10 chl fluorescence, ex 532 / ex 488'); ylabel('log10 scatter 532 / 488'); title('scatter ratio vs excitation ratio'); legend('location', 'eastoutside')
    end
end
title(tl, 'July/Aug 2026 cultures: pigment signals at 532 vs 488 nm (medians of gated events)')
exportgraphics(f, fullfile(od, 'culture_pigment_ratios.png'), 'Resolution', 110)
end
