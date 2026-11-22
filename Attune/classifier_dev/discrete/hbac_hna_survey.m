function [R, E] = hbac_hna_survey(cruisedirs, nper, seed, outname)
% function [R, E] = hbac_hna_survey(cruisedirs, nper, seed, outname)
% HNA / LNA survey of SYBR-stained bacteria runs (hbac_SYBR_SSC, preserved).
% For each cruise, nper hbac runs drawn at random (seeded) from the Gated_Table
% list; each run is classified with assign_class_hbac and its small events
% (SSC < -0.8) are kept, subsampled to 30000, with the auto class, so that
% bacteria (9) can be examined in SYBR (BL1) vs size (SSC).
% Events before the hand time window (first 20% of the run) are dropped
% where the run was hand gated.
%
% R - one row per run: cruise, file, date, depth, cast, niskin, counts, BL1 edge
% E - events: run index, BL1, SSC, BL3, GL2, class
% Saves output/discrete/<outname>.mat (default hbac_hna_survey)
%
% Nov 2026, discrete-sample classifier development (HNA/LNA split)

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 2 || isempty(nper), nper = 8; end
if nargin < 3, seed = 1; end
if nargin < 4, outname = 'hbac_hna_survey'; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = discrete_config();
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
R = table; E = table;
for ci = 1:numel(cruisedirs)
    cd_ = cruisedirs{ci}; cruise = extractAfter(cd_, '_'); cfg.cruise = cruise;
    S = load(fullfile(root, cd_, 'preserved', 'outputs', 'Gated_Table.mat'));
    fn = fieldnames(S); G = S.(fn{find(cellfun(@(q) istable(S.(q)), fn), 1)});
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    v = G.Properties.VariableNames;
    ix = find(contains(G.fcslist, 'hbac') & ~contains(G.fcslist, {'skip' 'test' 'badly' 'virus'}));
    rng(seed + ci);
    ix = ix(sort(randperm(numel(ix), min(nper, numel(ix)))));
    dv = v(startsWith(v, 'depth_m'));
    for i = ix'
        try, [dat, hdr] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', G.fcslist{i})); catch, continue, end
        t = datenum([hdr.date ' ' hdr.starttime]);
        ref = discrete_bead_reference(lab, t, cfg);
        if isempty(ref), fprintf('%s no beads\n', G.fcslist{i}), continue, end
        try, [cl, info] = assign_class_hbac(dat, hdr, ref, cfg); catch me, fprintf('%s %s\n', G.fcslist{i}, me.message), continue, end
        F = info.features;
        keep = true(height(F), 1);
        if ismember('class', v) && ~isempty(G.class{i}) && numel(G.class{i}) == height(F) && any(G.class{i})
            h = double(G.class{i}(:)); keep = dat(:, 1) >= min(dat(h ~= 0, 1));
        end
        dep = NaN; if ~isempty(dv), d = G.(dv{1})(i); if iscell(d), d = str2double(d{1}); end, dep = double(d); end
        cs = NaN; nk = NaN;
        if ismember('cast', v), cs = double(G.cast(i)); end
        if ismember('niskin', v), nk = double(G.niskin(i)); end
        R = [R; table(string(cruise), string(G.fcslist{i}), t, dep, cs, nk, sum(keep), sum(cl == 9 & keep), sum(cl == 7 & keep), ...
            sum(cl == 2 & keep), info.edges.BL1, info.edges.SSC, ...
            'VariableNames', {'cruise' 'file' 'time' 'depth' 'cast' 'niskin' 'n' 'bact' 'pro' 'syn' 'edge_BL1' 'edge_SSC'})]; %#ok<AGROW>
        j = find(keep & F.SSC < -0.8 & isfinite(F.BL1) & isfinite(F.SSC));
        if numel(j) > 30000, j = j(sort(randperm(numel(j), 30000))); end
        E = [E; table(repmat(height(R), numel(j), 1), single(F.BL1(j)), single(F.SSC(j)), single(F.BL3(j)), single(F.GL2(j)), int8(cl(j)), ...
            'VariableNames', {'run' 'BL1' 'SSC' 'BL3' 'GL2' 'class'})]; %#ok<AGROW>
        fprintf('%s %s depth %g bact %d\n', cruise, G.fcslist{i}, dep, R.bact(end))
    end
end
save(fullfile(od, [outname '.mat']), 'R', 'E')
end
