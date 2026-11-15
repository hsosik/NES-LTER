function [C, P] = syn_spectral_types(nPerCruise, outdir)
% function [C, P] = syn_spectral_types(nPerCruise, outdir)
% Are there discrete Syn pigment types that recur across samples and cruises?
% For a sample of files from every cruise with beadstat.mat, classify with
% assign_class_v2 and describe the Syn by bead-normalized spectral ratios:
%   ex  = BL2 - GL2   PE fluorescence excited at 488 vs 532 nm: a proxy for the
%                     PUB:PEB excitation ratio that separates types 3a/3b/3c
%   em  = GL2 - GL3   green-excited PE, two emission bands (PE-I/PE-II or PC)
%   pc  = GL4 - GL2   green-excited far-red (PC / allophycocyanin side) vs PE
%   pchl = GL2 - BL3  PE per chl
% Per file, a 1-3 component Gaussian mixture on [ex em] (BIC) gives the Syn
% spectral groups; C has one row per group (centre, spread, share), P holds a
% per-file subsample of Syn events for pooled density plots.
% Ratios are bead-normalized per channel (0.5 um bead for BL2/GL2/GL4, 6 um for GL3),
% so they compare across cruises if the bead lots' own spectra are stable.
%
% Nov 2026, classifier redesign exploration

if ~exist('nPerCruise', 'var') || isempty(nPerCruise), nPerCruise = 8; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'))
if ~exist('outdir', 'var') || isempty(outdir), outdir = fullfile(here, 'output'); end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
exclude = {'bead' 'Bead' 'Cast' 'cast' 'test' 'Dockwater' 'discrete' 'Rinses' 'Filter config' 'Cultures' 'Grazer' 'Dilution'};
D = dir([root '20*_*']); D = D([D.isdir]);
D = D(arrayfun(@(d) exist(fullfile(root, d.name, 'bead_calibrated', 'beadstat.mat'), 'file') == 2 && ...
    exist(fullfile(root, d.name, 'bead_calibrated', 'FCSfileinfo.mat'), 'file') == 2, D));
fprintf('%d cruises with beadstat.mat\n', numel(D))
Cc = cell(numel(D), 1); Pc = cell(numel(D), 1);
parfor i = 1:numel(D)
    try
        [Cc{i}, Pc{i}] = do_cruise(fullfile(root, D(i).name), D(i).name, nPerCruise, exclude);
        fprintf('%s: %d files, %d groups\n', D(i).name, numel(unique(Cc{i}.file)), height(Cc{i}));
    catch ME
        fprintf('%s: FAILED %s\n', D(i).name, ME.message);
    end
end
C = vertcat(Cc{:}); P = vertcat(Pc{:});
save(fullfile(outdir, 'syn_spectral_types.mat'), 'C', 'P')
end

function [C, P] = do_cruise(cpath, cdir, nPerCruise, exclude)
cfg = classifier_config();
bp = fullfile(cpath, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
F = FCSfileinfo;
keep = ~contains(F.fcslist, exclude);
if ismember('QC_flag', F.Properties.VariableNames), keep = keep & F.QC_flag == 1; end
F = sortrows(F(keep, :), 'matdate_start');
if ~ismember('pro_measured', F.Properties.VariableNames), F.pro_measured = zeros(height(F), 1); end
ii = find(F.pro_measured == 0); if isempty(ii), ii = (1:height(F))'; end
pick = ii(unique(round(linspace(1, numel(ii), min(nPerCruise, numel(ii))))));
uw = load_uw(bp);
C = table; P = table;
for k = pick'
    fn = F.fcslist{k};
    try
        [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
        ref = bead_reference(bs, F.matdate_start(k), cfg);
        [new, info] = assign_class_v2(d, h, ref, cfg, F.pro_measured(k));
    catch
        continue
    end
    if ~isfield(info, 'features'), continue, end
    X = info.features; pn = {h.par.name};
    sat = false(size(new));
    for ch = {'BL2' 'GL2' 'GL3' 'GL4'}
        j = strcmp(pn, [ch{1} '-H']); if any(j), sat = sat | d(:, j) >= cfg.sig.sat_frac * max(d(:, j)); end
    end
    ok = new == 2 & ~sat & X.BL2 > info.floor.BL2 + 1 & X.GL2 > info.floor.GL2 + 1 & X.GL3 > info.floor.GL3 + 1;
    if sum(ok) < 200, continue, end
    ex = X.BL2(ok) - X.GL2(ok); em = X.GL2(ok) - X.GL3(ok); pc = X.GL4(ok) - X.GL2(ok); pchl = X.GL2(ok) - X.BL3(ok);
    % per-file mixture on the two spectral ratios
    Z = [ex em];
    rng(1, 'twister')
    best = []; bic = Inf;
    for K = 1:3
        try
            g = fitgmdist(Z(randperm(size(Z, 1), min(20000, size(Z, 1))), :), K, 'RegularizationValue', 1e-4, 'Replicates', 3);
            if g.BIC < bic, best = g; bic = g.BIC; end
        catch
        end
    end
    if isempty(best), continue, end
    lab = cluster(best, Z);
    lat = NaN; lon = NaN; temp = NaN;
    if ~isempty(uw), u = find(strcmp(uw.file, fn), 1); if ~isempty(u), lat = uw.lat(u); lon = uw.lon(u); temp = uw.temp(u); end, end
    for c = 1:best.NumComponents
        m = lab == c;
        if sum(m) < 50, continue, end
        C = [C; table(string(cdir), string(fn), F.matdate_start(k), lat, lon, temp, best.NumComponents, c, sum(m), sum(m) / numel(m), ...
            median(ex(m)), median(em(m)), median(pc(m)), median(pchl(m)), 1.4826 * mad(ex(m), 1), 1.4826 * mad(em(m), 1), median(X.size(ok)), ...
            'VariableNames', {'cruise' 'file' 'time' 'lat' 'lon' 'temp' 'K' 'comp' 'n' 'share' 'ex' 'em' 'pc' 'pchl' 'ex_sd' 'em_sd' 'syn_size'})]; %#ok<AGROW>
    end
    q = randperm(numel(ex), min(1000, numel(ex)));
    P = [P; table(repmat(string(cdir), numel(q), 1), repmat(string(fn), numel(q), 1), repmat(temp, numel(q), 1), ex(q), em(q), pc(q), pchl(q), ...
        'VariableNames', {'cruise' 'file' 'temp' 'ex' 'em' 'pc' 'pchl'})]; %#ok<AGROW>
end
end

function uw = load_uw(bp)
uw = [];
fn = fullfile(bp, 'AttuneTable_uw_match.mat');
if ~exist(fn, 'file'), return, end
s = load(fn, 'Attune_uw_match');
if ~isfield(s, 'Attune_uw_match'), return, end
t = s.Attune_uw_match; v = t.Properties.VariableNames;
uw.file = cellstr(t.Filename);
uw.lat = first_good(t, v(contains(v, 'lat', 'IgnoreCase', true) & ~contains(v, 'flag', 'IgnoreCase', true)));
uw.lon = first_good(t, v(contains(v, 'lon', 'IgnoreCase', true) & ~contains(v, 'flag', 'IgnoreCase', true)));
uw.temp = first_good(t, v(contains(v, {'sbe48t' 'tsg1_temperature' 'tsg1_sst' 'sbe45t' 'sst' 'temperature'}, 'IgnoreCase', true) & ~contains(v, 'air', 'IgnoreCase', true)));
end

function x = first_good(t, names)
x = nan(height(t), 1);
for i = 1:numel(names)
    y = t.(names{i});
    if isnumeric(y) && any(isfinite(y)), x = double(y); return, end
end
end
