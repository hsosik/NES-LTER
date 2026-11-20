function S = preservation_shift(cruisedir, maxdepth, maxmin, outdir)
% function S = preservation_shift(cruisedir, maxdepth, maxmin, outdir)
% How do preserved cells sit on the bead-normalized scale compared with live
% cells from the same water? For each surface Niskin bottle with preserved
% runs, and the underway (live) files taken within maxmin minutes of it:
%   - preserved runs: bead-relative features (attune_features via
%     assign_class_v2, lab bead runs from lab_beadstat), population medians
%     for the hand gates (Gated_Table class: 1 Euk, 2 Syn, 4 Pro, 5 low-PE euk,
%     6 high-PE euk) and for v2 run unchanged on the preserved run; lab bead
%     runs from the same day only (discrete_bead_reference: runs without one
%     are skipped; a population that failed QC that day comes from the
%     nearest good run, and days_to_beads records how far away)
%   - underway runs: the same medians for v2 classes (cruise bead table); for
%     bottles with a preserved pro run, also the nearest Pro-setting underway
%     file within 3 h
% e.g. S = preservation_shift('20240906_EN720')
%
% cruisedir - cruise folder under cruise_data, e.g. '20240906_EN720'
% maxdepth  - bottles shallower than this (m) count as surface (default 7)
% maxmin    - underway files starting within this many minutes (default 30)
% outdir    - default classifier_dev/output/discrete
% S         - long table: one row per file x scheme x class, with n and
%             medians of the features; saved as presshift_<cruise>.mat
%
% Nov 2026, discrete-sample classifier development

if nargin < 2 || isempty(maxdepth), maxdepth = 7; end
if nargin < 3 || isempty(maxmin), maxmin = 30; end
maxmin_pro = 180; %Pro-setting underway files are sparser (alternating settings)
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
if nargin < 4 || isempty(outdir), outdir = fullfile(fileparts(here), 'output', 'discrete'); end
if ~exist(outdir, 'dir'), mkdir(outdir), end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cruise = extractAfter(cruisedir, '_');
cfg = classifier_config();
cfg.cruise = cruise;
cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
lab = load(fullfile(fileparts(here), 'output', 'discrete', 'lab_beadstat.mat'));
labbt = datenum(lab.beadstat.time);

% preserved runs at surface bottles
g = load(fullfile(root, cruisedir, 'preserved', 'outputs', 'Gated_Table.mat'));
G = g.gated_table;
G.Properties.VariableNames = lower(G.Properties.VariableNames);
if ~all(ismember({'date_sampled' 'depth_m'}, G.Properties.VariableNames))
    % no step 3 yet: bottle time and depth from the NES-LTER API
    bt = readtable(['https://nes-lter-api.whoi.edu/api/ctd/bottles/' lower(cruise) '.csv'], 'delimiter', ',');
    bt = renamevars(bt, {'date' 'depsm'}, {'date_sampled' 'depth_m'});
    G = outerjoin(G, bt, 'keys', {'cast' 'niskin'}, 'RightVariables', {'date_sampled' 'depth_m'}, 'Type', 'left');
    fprintf('%s: bottle metadata from the NES-LTER API\n', cruise)
end
ts = todatenum(G.date_sampled);
surf = G.depth_m <= maxdepth & isfinite(ts) & ~cellfun(@isempty, G.class) & ~contains(G.fcslist, 'hbac');
G = G(surf, :); ts = ts(surf);
fprintf('%s: %d preserved phyto runs at %d surface bottles\n', cruise, height(G), numel(unique(ts)))

% underway runs
bp = fullfile(root, cruisedir, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & ~contains(FCSfileinfo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end

S = table;
[ub, ~, ib] = unique([G.cast G.niskin], 'rows');
for b = 1:size(ub, 1)
    rows = find(ib == b);
    t0 = ts(rows(1));
    bottle = {string(cruise), ub(b, 1), ub(b, 2), G.depth_m(rows(1)), t0};
    % preserved
    for r = rows'
        fn = G.fcslist{r};
        rt = runtype(fn);
        try
            [fcsdat, fcshdr] = fca_readfcs(fullfile(root, cruisedir, 'preserved', 'FCS', fn));
        catch err
            warning('%s: %s', fn, err.message), continue
        end
        trun = datenum([fcshdr.date ' ' fcshdr.starttime]);
        ref = discrete_bead_reference(lab.beadstat, trun, cfg); %same-day lab beads; failed populations from the nearest good run
        if isempty(ref), warning('%s: no lab bead run on %s, skipped', fn, datestr(trun, 1)), continue, end
        fb = cell2mat(struct2cell(ref.fallback));
        [v2, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, strcmp(rt, 'pro'));
        hand = double(G.class{r}(:));
        if numel(hand) ~= numel(v2), warning('%s: class length mismatch', fn), continue, end
        extra = {"preserved", rt, string(fn), trun, max(fb)}; %days_to_beads: 0 same day, else the largest fallback
        S = [S; summarize(info.features, hand, "hand", bottle, extra); summarize(info.features, v2, "v2", bottle, extra)]; %#ok<AGROW>
    end
    % live (underway) within maxmin minutes
    near = find(abs(Fi.matdate_start - t0) * 1440 <= maxmin);
    [~, o] = sort(abs(Fi.matdate_start(near) - t0)); near = near(o(1:min(2, end)));
    % Pro-setting underway file for the preserved pro run: nearest within maxmin_pro
    if any(runtype(string(G.fcslist(rows))) == "pro") && ~any(Fi.pro_measured(near))
        pn = find(Fi.pro_measured & abs(Fi.matdate_start - t0) * 1440 <= maxmin_pro);
        [~, o] = min(abs(Fi.matdate_start(pn) - t0)); near = [near; pn(o)];
    end
    for k = near'
        fn = Fi.fcslist{k};
        try
            [fcsdat, fcshdr] = fca_readfcs(fullfile(root, cruisedir, 'FCS', fn));
        catch err
            warning('%s: %s', fn, err.message), continue
        end
        ref = bead_reference(bs, Fi.matdate_start(k), cfg);
        [v2, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, Fi.pro_measured(k));
        rt = "underway"; if Fi.pro_measured(k), rt = "underway_pro"; end
        extra = {"live", rt, string(fn), Fi.matdate_start(k), NaN};
        S = [S; summarize(info.features, v2, "v2", bottle, extra)]; %#ok<AGROW>
    end
    fprintf('bottle C%dN%d: %d preserved, %d underway files\n', ub(b, 1), ub(b, 2), numel(rows), numel(near))
end
save(fullfile(outdir, ['presshift_' cruise '.mat']), 'S', 'maxdepth', 'maxmin')
end

function T = summarize(F, cls, scheme, bottle, extra)
feat = {'size' 'SSC' 'FSC' 'GL1' 'BL1' 'BL2' 'BL3' 'GL2' 'GL3' 'GL4' 'PEchl' 'GL2GL3' 'chlSSC'};
F.BL2GL2 = F.BL2 - F.GL2;
feat = [feat {'BL2GL2'}];
T = table;
for c = setdiff(unique(cls(:))', 0)
    ix = cls == c;
    m = varfun(@median, F(ix, feat)); m.Properties.VariableNames = feat;
    sd = (prctile(F.size(ix), 84.1) - prctile(F.size(ix), 15.9)) / 2;
    row = [cell2table([bottle extra {scheme, c, sum(ix), sd}], 'VariableNames', ...
        {'cruise' 'cast' 'niskin' 'depth_m' 'time_sampled' 'kind' 'runtype' 'file' 'time_run' 'days_to_beads' 'scheme' 'class' 'n' 'size_sd'}) m];
    T = [T; row]; %#ok<AGROW>
end
end

function rt = runtype(fn)
if contains(fn, 'pro', 'IgnoreCase', true), rt = "pro";
elseif contains(fn, 'phyto_PE'), rt = "PE";
elseif contains(fn, 'phyto_CHL'), rt = "CHL";
else, rt = "other";
end
end

function t = todatenum(d)
if isdatetime(d), t = datenum(d); return, end
if isnumeric(d), t = d; return, end
d = string(d);
t = nan(size(d));
for i = 1:numel(d)
    if ismissing(d(i)) || strlength(d(i)) == 0, continue, end
    s = extractBefore(d(i) + "+", "+"); s = replace(s, "T", " ");
    try
        t(i) = datenum(datetime(s, 'InputFormat', 'yyyy-MM-dd HH:mm:ss'));
    catch
        try, t(i) = datenum(datetime(s)); catch, end
    end
end
end
