function T = scan_edge_euks()
% function T = scan_edge_euks()
% With the Pro-tail cut gone (Nov 2026): per file, Euks that sit where the cut used
% to act (no PE, chl within edge_pro_band of the chl trigger edge, no larger than
% Syn + 0.1), whether the Euk chl-gap rule found a gap, and what the old classes
% called those events. Files: the 96 non-Pro files of rescore_euk_pe_rule.
% Saves output/early_layout_robust/edge_euks.mat
dev = fileparts(mfilename('fullpath')); addpath(fullfile(dev, '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
out = fullfile(dev, 'output', 'early_layout_robust');
P = load(fullfile(out, 'euk_pe_rule.mat')).T;
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
        if ~isfield(info.trigger_edges, 'BL3') || info.ncomp == 0, continue, end
        o = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat'))).class(:, 1);
        F = info.features; e = info.trigger_edges.BL3;
        hp = F.GL2 > info.rules.pe_min & F.PEchl > cfg.sig.pe_min_PEchl;
        ee = a == 1 & ~hp & F.BL3 < e + cfg.sig.edge_pro_band & F.size < info.syn_ref.size + 0.1;
        T = [T; table(c, string(cd), string(fn), sum(ee), sum(ee & o == 1), sum(ee & o == 0), sum(a == 1), isfinite(info.euk_chl_gap), ...
            'VariableNames', {'cruise' 'cruisedir' 'file' 'n_edge_euk' 'old_euk' 'old_none' 'euk' 'gap_found'})]; %#ok<AGROW>
        fprintf('%s %s edge Euks %d (old Euk %d, old none %d) of %d Euks, gap %d\n', c, extractBefore(fn, '_Group'), T.n_edge_euk(end), T.old_euk(end), T.old_none(end), T.euk(end), T.gap_found(end));
    end
end
save(fullfile(out, 'edge_euks.mat'), 'T')
end
