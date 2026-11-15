function T = rescore_euk_pe_rule()
% function T = rescore_euk_pe_rule()
% After removing the Pro-tail cut and adding "Euks carry no Syn-like PE" (Nov 2026):
% re-score the 96 non-Pro files of assess_pro_tail_removal against the old classes
% and compare with the agreement before (Pro-tail cut on) and with the cut simply off.
% Saves output/early_layout_robust/euk_pe_rule.mat
dev = fileparts(mfilename('fullpath')); addpath(fullfile(dev, '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
out = fullfile(dev, 'output', 'early_layout_robust');
P = load(fullfile(out, 'pro_tail_removal.mat')).T;
[~, ia] = unique(P.file, 'stable'); P = P(ia, :);
dirs = dir(root); dirs = {dirs([dirs.isdir]).name};
T = table;
for c = unique(P.cruise, 'stable')'
    cd = dirs{endsWith(dirs, ['_' char(c)])};
    bp = fullfile(root, cd, 'bead_calibrated');
    bs = load_beadstat(bp); load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
    cfg = classifier_config(); cfg.cruise = char(c); cfg.gl1_is_od2_scatter = strcmpi(attune_config(char(c)).OD2, 'GL1');
    for k = find(P.cruise == c)'
        fn = char(P.file(k)); r = strcmp(FCSfileinfo.fcslist, fn);
        [d, h] = fca_readfcs(fullfile(root, cd, 'FCS', fn));
        ref = bead_reference(bs, FCSfileinfo.matdate_start(r), cfg);
        [a, info] = assign_class_v2(d, h, ref, cfg, 0);
        o = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat'))).class(:, 1);
        o(~ismember(o, [0 1 2 3 4 7])) = 0;
        ag = sum(o == a & o ~= 0) / sum(o ~= 0 | a ~= 0);
        T = [T; table(c, P.file(k), P.ag_with_cut(k), P.ag_without(k), ag, info.n_euk_synpe, sum(a == 1), sum(a == 2), sum(o == 1), sum(o == 2), ...
            P.euk_with_cut(k), 'VariableNames', {'cruise' 'file' 'ag_before' 'ag_cut_off' 'ag_now' 'n_euk_synpe' 'euk_now' 'syn_now' 'old_euk' 'old_syn' 'euk_before'})]; %#ok<AGROW>
        fprintf('%s %s ag %.3f -> %.3f (cut off only %.3f) | Euk %d -> %d (old %d), synpe removed %d\n', c, extractBefore(fn, '_Group'), T.ag_before(end), ag, T.ag_cut_off(end), ...
            T.euk_before(end), T.euk_now(end), T.old_euk(end), T.n_euk_synpe(end));
    end
end
save(fullfile(out, 'euk_pe_rule.mat'), 'T')
end
