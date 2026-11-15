function T = assess_pro_tail_removal(nfiles)
% (ran against the code before the Pro-tail cut was removed, Nov 2026; the edge_pro_tail switch no longer exists)
% function T = assess_pro_tail_removal(nfiles)
% Can the Pro-tail cut go now that the Euk chl-gap rule exists? For non-Pro files
% of early and post-2019 cruises, classify with the cut on (default) and off
% (cfg.sig.edge_pro_tail = false; the capped noise-pile rule stays on) and score
% both against the old classes. Also counts the events the cut removes that the
% gap rule would catch anyway.
% Saves output/early_layout_robust/pro_tail_removal.mat
% Nov 2026
dev = fileparts(mfilename('fullpath')); addpath(fullfile(dev, '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
out = fullfile(dev, 'output', 'early_layout_robust');
C = {'20190705_TN368' '20190414_AR34B' '20190512_RB1904' '20190201_EN627' '20250124_EN727' '20250814_AR92' '20251018_AR95' '20260814_HRS2609'};
T = table;
for i = 1:numel(C)
    cruise = extractAfter(C{i}, '_'); bp = fullfile(root, C{i}, 'bead_calibrated');
    bs = load_beadstat(bp); load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
    if ~ismember('pro_measured', FCSfileinfo.Properties.VariableNames), FCSfileinfo.pro_measured = zeros(height(FCSfileinfo), 1); end
    f = fullfile(out, ['agreement_' cruise '_layout.mat']);
    if isfile(f)
        files = load(f).V.file;
    else
        Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & FCSfileinfo.pro_measured == 0 & ~contains(FCSfileinfo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
        Fi = sortrows(Fi, 'matdate_start');
        files = string(Fi.fcslist(unique(round(linspace(1, height(Fi), nfiles)))));
        if strcmp(cruise, 'AR95'), files = unique([files; string(FCSfileinfo.fcslist(contains(FCSfileinfo.fcslist, '19Oct2025_AltPro(19)') & contains(FCSfileinfo.fcslist, 'Euk(20)')))]); end
        if strcmp(cruise, 'HRS2609'), files = unique([files; string(FCSfileinfo.fcslist(contains(FCSfileinfo.fcslist, '17Aug') & contains(FCSfileinfo.fcslist, 'Euk(8)')))]); end
    end
    cfg = classifier_config(); cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
    for k = 1:numel(files)
        r = strcmp(FCSfileinfo.fcslist, files(k));
        if FCSfileinfo.pro_measured(r), continue, end
        [d, h] = fca_readfcs(fullfile(root, C{i}, 'FCS', char(files(k))));
        ref = bead_reference(bs, FCSfileinfo.matdate_start(r), cfg);
        [a, ia] = assign_class_v2(d, h, ref, cfg, 0);
        c2 = cfg; c2.sig.edge_pro_tail = false;
        [b, ib] = assign_class_v2(d, h, ref, c2, 0);
        o = load(fullfile(bp, 'class', regexprep(char(files(k)), '\.fcs$', '.mat'))).class(:, 1);
        o(~ismember(o, [0 1 2 3 4 7])) = 0;
        sc = @(x) sum(o == x & o ~= 0) / sum(o ~= 0 | x ~= 0);
        ch = a ~= b; %events the cut changes
        T = [T; table(string(cruise), files(k), ia.edge_pile_capped, sc(a), sc(b), sum(ch), sum(ch & o == 1), sum(ch & o == 0), ...
            sum(a == 1), sum(b == 1), sum(o == 1), sum(ch & b ~= 1), ib.euk_chl_gap - ia.trigger_edges.BL3, ...
            'VariableNames', {'cruise' 'file' 'capped' 'ag_with_cut' 'ag_without' 'n_changed' 'changed_old_euk' 'changed_old_none' ...
            'euk_with_cut' 'euk_without' 'old_euk' 'changed_not_euk_without' 'gap_above_edge_without'})]; %#ok<AGROW>
        fprintf('%s %s ag %.3f -> %.3f without | changed %d (old Euk %d, old none %d) | Euk %d / %d / old %d | gap %.2f\n', cruise, extractBefore(files(k), '_Group'), ...
            T.ag_with_cut(end), T.ag_without(end), T.n_changed(end), T.changed_old_euk(end), T.changed_old_none(end), T.euk_with_cut(end), T.euk_without(end), T.old_euk(end), T.gap_above_edge_without(end));
    end
end
save(fullfile(out, 'pro_tail_removal.mat'), 'T')
end
