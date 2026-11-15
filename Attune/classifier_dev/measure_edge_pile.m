function P = measure_edge_pile(nfiles)
% function P = measure_edge_pile(nfiles)
% Fraction of chlorophyll-bearing events (assign_class_v2 has_chl) that sit in the
% pile at the chl trigger edge: no PE, BL3 within cfg.sig.edge_pro_band of the
% edge. Early cruises (files scored in output/early_layout) and post-2019
% cruises (files spread through each). Saves output/early_layout/edge_pile.mat
% Nov 2026, noise-pile robustness
dev = fileparts(mfilename('fullpath')); addpath(fullfile(dev, '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
out = fullfile(dev, 'output', 'early_layout');
C = {'20190705_TN368' '20190414_AR34B' '20190512_RB1904' '20190201_EN627' '20250124_EN727' '20250814_AR92' '20251018_AR95'};
P = table;
for i = 1:numel(C)
    cruise = extractAfter(C{i}, '_'); bp = fullfile(root, C{i}, 'bead_calibrated');
    bs = load_beadstat(bp); load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
    f = fullfile(out, ['agreement_' cruise '_layout.mat']);
    if isfile(f)
        V = load(f).V; files = V.file; ag = V.v2_agreement;
    else
        Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & ~contains(FCSfileinfo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
        Fi = sortrows(Fi, 'matdate_start');
        files = string(Fi.fcslist(unique(round(linspace(1, height(Fi), nfiles))))); ag = nan(size(files));
    end
    cfg = classifier_config(); cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
    for k = 1:numel(files)
        t = FCSfileinfo.matdate_start(strcmp(FCSfileinfo.fcslist, files(k)));
        [d, h] = fca_readfcs(fullfile(root, C{i}, 'FCS', char(files(k))));
        ref = bead_reference(bs, t, cfg);
        [F, fi] = attune_features(d, h, ref);
        e = NaN; eg = NaN;
        tr = [fi.trigger fi.trigger2]; lv = [fi.trigger_level fi.trigger2_level];
        for j = 1:2
            ch = char(regexp(tr(j), '(?<=_)\w+$', 'match', 'once'));
            if strcmp(ch, 'BL3'), e = fi.floor.BL3 + log10(lv(j)) - log10(cfg.floor_raw); end
            if strcmp(ch, 'GL2'), eg = fi.floor.GL2 + log10(lv(j)) - log10(cfg.floor_raw); end
        end
        pe_min = max(cfg.sig.pe_min_GL2, eg);
        has_pe = F.GL2 > pe_min & F.GL2 - F.BL3 > cfg.sig.pe_min_PEchl & F.BL3 > cfg.sig.chl_min;
        has_chl = (F.BL3 > max(cfg.sig.chl_min, e) | has_pe) & F.SSC > fi.floor.SSC + 0.1;
        pile = has_chl & ~has_pe & F.BL3 < e + cfg.sig.edge_pro_band;
        P = [P; table(string(cruise), files(k), ag(k), height(F), sum(has_chl), sum(pile), sum(pile) / sum(has_chl), ...
            sum(has_chl & ~pile), lv(1), 'VariableNames', {'cruise' 'file' 'agreement' 'n' 'n_has_chl' 'n_pile' 'pile_frac' 'n_nonpile' 'trig1_level'})]; %#ok<AGROW>
        fprintf('%s %s ag %.2f n %d chl %d pile %.3f nonpile %d\n', cruise, files(k), ag(k), height(F), sum(has_chl), P.pile_frac(end), P.n_nonpile(end));
    end
end
save(fullfile(out, 'edge_pile.mat'), 'P')
end
