function T = compare_synpart_edge(basedir)
% function T = compare_synpart_edge(basedir)
% Syn on particle must have chl above the chl trigger edge (Nov 2026). Compares
% assign_class_v2 with a copy of the code before that rule (assign_class_v2_presynpart
% in basedir) on: the 96 non-Pro files of rescore_euk_pe_rule, 8 Pro-setting files
% each from AR92, AR95 and HRS2609, and the EN617 files of validate_early_layout.
% Per file: events changed, Syn on particle before/after, agreement with the old
% classes before/after (where class files exist).
% Saves output/early_layout_robust/synpart_edge.mat
dev = fileparts(mfilename('fullpath')); addpath(fullfile(dev, '..')); addpath(basedir)
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
out = fullfile(dev, 'output', 'early_layout_robust');
dirs = dir(root); dirs = {dirs([dirs.isdir]).name};
P = load(fullfile(out, 'euk_pe_rule.mat')).T;
L = table(string(cellfun(@(c) dirs{endsWith(dirs, ['_' char(c)])}, cellstr(P.cruise), 'UniformOutput', false)), P.file, false(height(P), 1), ...
    'VariableNames', {'cruisedir' 'file' 'pro'});
for cd = {'20250814_AR92' '20251018_AR95' '20260814_HRS2609'}
    bp = fullfile(root, cd{1}, 'bead_calibrated'); load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
    Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & FCSfileinfo.pro_measured == 1 & ~contains(FCSfileinfo.fcslist, {'test' 'bead' 'Bead'}), :);
    Fi = Fi(cellfun(@(c) isfile(fullfile(bp, 'class', regexprep(c, '\.fcs$', '.mat'))), Fi.fcslist), :); Fi = sortrows(Fi, 'matdate_start');
    f = string(Fi.fcslist(unique(round(linspace(1, height(Fi), 8)))));
    L = [L; table(repmat(string(cd{1}), numel(f), 1), f, true(numel(f), 1), 'VariableNames', {'cruisedir' 'file' 'pro'})]; %#ok<AGROW>
end
V = load(fullfile(out, 'agreement_EN617_layout.mat')).V;
L = [L; table(repmat("20180720_EN617", height(V), 1), V.file, false(height(V), 1), 'VariableNames', {'cruisedir' 'file' 'pro'})];
T = table; lastcd = "";
for k = 1:height(L)
    cd = char(L.cruisedir(k)); cruise = extractAfter(cd, '_'); fn = char(L.file(k));
    bp = fullfile(root, cd, 'bead_calibrated');
    if ~strcmp(lastcd, cd)
        bs = load_beadstat(bp); FI = [];
        if isfile(fullfile(bp, 'FCSfileinfo.mat')), FI = load(fullfile(bp, 'FCSfileinfo.mat')).FCSfileinfo; end
        cfg = classifier_config(); cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
        lastcd = string(cd);
    end
    [d, h] = fca_readfcs(fullfile(root, cd, 'FCS', fn));
    if isempty(FI), t = datenum([h.date ' ' h.starttime]); else, t = FI.matdate_start(strcmp(FI.fcslist, fn)); end %#ok<DATNM>
    ref = bead_reference(bs, t, cfg);
    a = assign_class_v2(d, h, ref, cfg, double(L.pro(k)));
    b = assign_class_v2_presynpart(d, h, ref, cfg, double(L.pro(k)));
    cf = fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat'));
    agb = NaN; aga = NaN;
    if isfile(cf)
        o = load(cf).class(:, 1); o(~ismember(o, [0 1 2 3 4 7])) = 0;
        sc = @(x) sum(o == x & o ~= 0) / sum(o ~= 0 | x ~= 0); agb = sc(b); aga = sc(a);
    end
    T = [T; table(string(cruise), string(fn), L.pro(k), sum(a ~= b), sum(b == 8), sum(a == 8), sum(b == 2), sum(a == 2), agb, aga, ...
        'VariableNames', {'cruise' 'file' 'pro' 'n_changed' 'synpart_before' 'synpart_after' 'syn_before' 'syn_after' 'ag_before' 'ag_after'})]; %#ok<AGROW>
    fprintf('%s %s pro %d | changed %d | Syn on particle %d -> %d | Syn %d -> %d | ag %.3f -> %.3f\n', cruise, extractBefore(fn, '_Group'), L.pro(k), T.n_changed(end), ...
        T.synpart_before(end), T.synpart_after(end), T.syn_before(end), T.syn_after(end), agb, aga);
end
save(fullfile(out, 'synpart_edge.mat'), 'T')
end
