function M = field_pigment_map(cruises, nfiles)
% function M = field_pigment_map(cruises, nfiles)
% Chl fluorescence excitation ratio (532 vs 488 nm: GL4 - BL3 in log10 bead
% units, attune_features) against cell volume for eukaryotes in field samples,
% with the July/Aug 2026 culture medians on the same scale as reference.
% Bead units remove drift in the two lasers' relative power between cruises.
%   - volume: merged SSC/GL1 scatter (cfg.gl1_merge_bead window, as
%     assign_class_v2) -> Attune_SC2vol('SSC-H')
%   - classes from assign_class_v2: 1 = euk without PE, 3/4 = cryptophytes
%   - cultures: same features, each culture day's ref from that day's gated
%     FCB bead files (1 um for BL3/SSC/GL1, 0.5 um for GL4, as cfg.bead_pop)
%   - March/April 2019 cultures (old layout, OD2 on SSC; squares): green chl
%     680/40 was file GL3-H, mapped to physical GL4 by attune_config through
%     attune_features; bead ref from the same day's FCB bead run (find_beads);
%     volume from SSC-H / 1 um bead with Attune_SC2vol('SSC-H_488_2019');
%     events gated on chl (BL3-H > 10^2.5; Syn also PE > 10^3.5) as sizecal_2019
% Current-layout cruises only (GL1 = side scatter through OD2). Files named
% Dockwater / config_change are skipped (EN644: dock samples of 6-15 Aug 2019,
% before the layout change, read with the wrong channels before Nov 2026).
% Draws ../output/pt_transfer/field_pigment_map.png.
% field_pigment_map('replot') redraws from the saved field_pigment_map.mat.
%
% Nov 2026, PT transfer / pigment ratios

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
warning('off', 'all')
od = fullfile(dev, 'output', 'pt_transfer');
cfg = classifier_config();
if nargin > 0 && ischar(cruises) && strcmp(cruises, 'replot')
    s = load(fullfile(od, 'field_pigment_map.mat'), 'M'); M = s.M;
    plot_map(M, od)
    return
end
root = ['\' '\sosiknas1\lab_data\Attune\cruise_data\'];
if nargin < 1 || isempty(cruises)
    cruises = {'EN644' 'EN649' 'EN655' 'EN661' 'EN668' 'AT46' 'EN688' 'EN695' 'EN706' 'EN712' 'EN720' 'EN727' 'AR92' 'AR95' 'AR104'};
end
if nargin < 2 || isempty(nfiles), nfiles = 8; end
D = dir(root); D = {D([D.isdir]).name};

% ---- cultures
p26 = ['\' '\sosiknas1\Lab_data\Attune\size_calibration_July2026\'];
gd = fullfile(p26, 'Attune', 'Gated_FCS');
s = load(fullfile(od, 'culture_pigment_ratios.mat'), 'T'); CT = s.T;
CT.group(CT.culture == "SE62") = "coccolithophore (SE62)"; % haptophyte; saved table may carry the old 'unknown' label
C = {};
for i = 1:height(CT)
    fn = fullfile(gd, CT.file(i) + ".fcs");
    day = extractBefore(CT.file(i), 12);
    ref = day_ref(gd, day, cfg);
    if isempty(ref), continue, end
    [d, h] = fca_readfcs(char(fn));
    F = attune_features(d, h, ref);
    sz = merged_size(F, cfg);
    ok = F.BL3 > -2 & F.GL4 > -3;
    C{end+1, 1} = struct('culture', CT.culture(i), 'group', CT.group(i), ...
        'lvol', median(log10(Attune_SC2vol(10.^sz(ok), 'SSC-H', 0))), 'pig', median(F.GL4(ok) - F.BL3(ok))); %#ok<AGROW>
end
C = struct2table([C{:}]');
C = groupsummary(C, {'culture', 'group'}, 'median', {'lvol', 'pig'});
C.year = repmat(2026, height(C), 1);
C19 = cultures_2019(cfg);
C = [C(:, {'culture' 'group' 'median_lvol' 'median_pig' 'year'}); C19];

% ---- field
E = {};
for ic = 1:numel(cruises)
    m = D(endsWith(D, ['_' cruises{ic}]));
    if isempty(m), continue, end
    bp = fullfile(root, m{1}, 'bead_calibrated');
    try
        bs = load_beadstat(bp); S = load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo'); Fi = S.FCSfileinfo;
    catch
        continue
    end
    Fi = Fi(Fi.QC_flag == 1 & ~contains(Fi.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete' 'Dockwater' 'config_change'}), :);
    if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
    pick = unique(round(linspace(1, height(Fi), min(nfiles, height(Fi)))));
    for k = pick
        try
            [d, h] = fca_readfcs(fullfile(root, m{1}, 'FCS', Fi.fcslist{k}));
            ref = bead_reference(bs, Fi.matdate_start(k), cfg);
            cl = assign_class_v2(d, h, ref, cfg, Fi.pro_measured(k));
        catch
            continue
        end
        F = attune_features(d, h, ref);
        sz = merged_size(F, cfg);
        q = ismember(cl, [1 3 4]) & F.GL4 > -3;
        mo = month(datetime(Fi.matdate_start(k), 'ConvertFrom', 'datenum'));
        E{end+1, 1} = table(repmat(string(cruises{ic}), sum(q), 1), repmat(mo, sum(q), 1), cl(q), ...
            log10(Attune_SC2vol(10.^sz(q), 'SSC-H', 0)), F.GL4(q) - F.BL3(q), ...
            'VariableNames', {'cruise', 'month', 'class', 'lvol', 'pig'}); %#ok<AGROW>
    end
    fprintf('%s done\n', cruises{ic})
end
E = vertcat(E{:});
M = struct('cultures', C, 'events', E);
save(fullfile(od, 'field_pigment_map.mat'), 'M', '-v7.3')
plot_map(M, od)
end

function plot_map(M, od)
C = M.cultures; E = M.events;

% ---- plot: winter (Dec-Apr) and summer (Jun-Oct), euks and cryptophytes
gs = ["cyanobacteria" "green" "haptophyte" "cryptophyte" "diatom" "dinoflagellate" "raphidophyte" "coccolithophore (SE62)" "unconfirmed (2019)"];
gc = [0.5 0.5 0.5; 0.2 0.65 0.2; 0.9 0.6 0; 0.8 0.1 0.5; 0.1 0.4 0.8; 0.6 0.3 0.1; 0.5 0.8 0.8; 0 0 0; 1 1 1];
seas = {'winter-spring (Dec-May)', @(m) ismember(m, [12 1 2 3 4 5]); 'summer-fall (Jun-Nov)', @(m) ismember(m, 6:11)};
f = figure('visible', 'off', 'position', [0 0 1450 620]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
xe = -0.5:0.05:4; ye = -1.2:0.02:0.6;
for sp = 1:2
    nexttile; hold on
    q = seas{sp, 2}(E.month);
    H = histcounts2(E.lvol(q & E.class == 1), E.pig(q & E.class == 1), xe, ye)';
    imagesc(xe(1:end-1) + 0.025, ye(1:end-1) + 0.01, log10(H + 1)); set(gca, 'ydir', 'normal')
    colormap(flipud(gray));
    Hc = histcounts2(E.lvol(q & ismember(E.class, [3 4])), E.pig(q & ismember(E.class, [3 4])), xe, ye)';
    if any(Hc(:)), contour(xe(1:end-1) + 0.025, ye(1:end-1) + 0.01, log10(Hc + 1), 4, 'm', 'linewidth', 1, 'DisplayName', 'cryptophytes (classes 3/4)'), end
    for yr = [2026 2019]
        mk = 'o'; if yr == 2019, mk = 's'; end
        for k = 1:numel(gs)
            r = C.group == gs(k) & C.year == yr;
            if ~any(r), continue, end
            plot(C.median_lvol(r), C.median_pig(r), mk, 'markersize', 9, 'markeredgecolor', 'k', ...
                'markerfacecolor', gc(k, :), 'DisplayName', gs(k) + " " + yr)
            lab = strrep(C.culture(r), '_', ' '); if yr == 2019, lab = lab + " '19"; end
            text(C.median_lvol(r) + 0.05, C.median_pig(r), lab, 'fontsize', 7)
        end
    end
    xlabel('log10 volume (\mum^3): field and 2026 merged SSC/GL1; 2019 SSC 488'); ylabel('chl ex 532 - ex 488 (GL4 - BL3, log10 bead units)')
    title(sprintf('%s: %d euk events, %s', seas{sp, 1}, sum(q & E.class == 1), strjoin(unique(E.cruise(q))', ' ')), 'fontsize', 8)
    axis([-0.5 4 -1.2 0.6]); grid on; box on
    if sp == 2, legend('location', 'eastoutside', 'fontsize', 8), end
end
title(tl, 'Field eukaryotes (density, grey) and cryptophytes (magenta contours) vs culture medians (circles 2026, squares 2019)')
exportgraphics(f, fullfile(od, 'field_pigment_map.png'), 'Resolution', 110)
end

function sz = merged_size(F, cfg)
sz = F.SSC;
hi = F.GL1 > log10(cfg.gl1_merge_bead(2));
mid = F.GL1 > log10(cfg.gl1_merge_bead(1)) & ~hi;
sz(hi) = F.GL1(hi);
sz(mid) = log10((10.^F.SSC(mid) + 10.^F.GL1(mid)) / 2);
end

function ref = day_ref(gd, day, cfg)
% bead reference from that day's gated FCB bead files (KDE peaks)
ref = [];
f1 = dir(fullfile(gd, [char(day) '_FCB_bead_check*_1_um_beads.fcs']));
f05 = dir(fullfile(gd, [char(day) '_FCB_bead_check*_0.5_um_beads.fcs']));
f6 = dir(fullfile(gd, [char(day) '_FCB_bead_check*_McLane.fcs']));
if isempty(f1) || isempty(f05), return, end
fs = {f05(1), f1(1)}; if ~isempty(f6), fs{3} = f6(1); end
ref.chans = cfg.chans; ref.floor_raw = cfg.floor_raw; ref.k = cfg.k;
for c = cfg.chans
    ch = c{1}; pidx = min(cfg.bead_pop.(ch), numel(fs));
    [d, h] = fca_readfcs(fullfile(gd, fs{pidx}.name));
    j = strcmp({h.par.name}, [ch '-H']);
    ref.bead_log.(ch) = bead_peak(log10(max(d(:, j), 1)));
    ref.bead_hv.(ch) = h.par(j).hv;
end
end

function C = cultures_2019(cfg)
% March/April 2019 culture runs with OD2 on SSC (old layout), on the field scale
d0 = ['\' '\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\FCS\'];
conf = attune_config('TN368', datenum(2019, 3, 20)); %#ok<DATNM> % OD2 on SSC, v2019 layout
sig = bead_signature();
bf = struct('d20Mar', 'size_cal_OD2_Mar2019_Group_day0_Bead Solution.fcs', ...
    'd9Apr', 'size_cal_OD2_Mar2019_Group_day0_Bead_Solution2.fcs');
refs = struct();
for k = fieldnames(bf)'
    [d, h] = fca_readfcs([d0 bf.(k{1})]);
    cl = find_beads(double(d), h, sig, conf);
    ref = struct('chans', {cfg.chans}, 'floor_raw', cfg.floor_raw, 'k', cfg.k, 'chan', conf.chan, 'layout', conf.layout);
    for c = cfg.chans
        ch = c{1}; j = [];
        if ~isempty(conf.chan.(ch)), j = find(strcmp({h.par.name}, conf.chan.(ch)), 1); end
        v = []; if ~isempty(j), v = d(cl == cfg.bead_pop.(ch) & d(:, j) > 0, j); end
        if numel(v) >= 50
            ref.bead_log.(ch) = bead_peak(log10(v)); ref.bead_hv.(ch) = h.par(j).hv;
        else
            ref.bead_log.(ch) = NaN; ref.bead_hv.(ch) = 1;
        end
    end
    refs.(k{1}) = ref;
end
cult = {'Syn' '1820i' 'd9Apr' 'cyanobacteria'; 'Micromonas' 'Micromonas' 'd20Mar' 'green'; 'Nano' 'Nano' 'd20Mar' 'green';
    'Pavlova' 'Pavlova' 'd20Mar' 'haptophyte'; 'I_galbana' 'I_galbana' 'd20Mar' 'haptophyte';
    'Prymn_pav' 'PRYMN_PAV' 'd20Mar' 'haptophyte'; 'Dun' 'Dun' 'd20Mar' 'green'; 'Crypto_spmc' 'CRYPTO_SPMC' 'd20Mar' 'cryptophyte';
    'Htriq' 'HTriq' 'd20Mar' 'dinoflagellate'; 'H_akashiwo' 'H_Akashiwo_CCMP3374' 'd20Mar' 'raphidophyte';
    'Ambopo' 'AMBOPO' 'd20Mar' 'dinoflagellate'; 'Pptp02' 'PPTP02' 'd20Mar' 'unconfirmed (2019)'};
C = table('Size', [0 5], 'VariableTypes', {'string' 'string' 'double' 'double' 'double'}, ...
    'VariableNames', {'culture' 'group' 'median_lvol' 'median_pig' 'year'});
for i = 1:size(cult, 1)
    F0 = dir([d0 'size_cal_OD2_Mar2019_Group_day0_' cult{i, 2} '*.fcs']);
    if isempty(F0), continue, end
    [d, h] = fca_readfcs([d0 F0(1).name]);
    ref = refs.(cult{i, 3});
    F = attune_features(d, h, ref);
    L = @(n) log10(max(d(:, strcmp({h.par.name}, n)), 1));
    g = L('BL3-H') > 2.5 & F.GL4 > -3 & isfinite(F.SSC);
    if strcmp(cult{i, 1}, 'Syn'), g = g & L('GL1-H') > 3.5; end
    C(end+1, :) = {string(cult{i, 1}), string(cult{i, 4}), ...
        median(log10(Attune_SC2vol(10.^F.SSC(g), 'SSC-H_488_2019', 0))), median(F.GL4(g) - F.BL3(g)), 2019}; %#ok<AGROW>
end
disp(C)
end
