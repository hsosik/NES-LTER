function T = measure_pro_tail(nfiles)
% (ran against the code before the Pro-tail cut was removed, Nov 2026; edge_pro_size_above_syn no longer exists)
% function T = measure_pro_tail(nfiles)
% How the Pro-tail cut (non-Pro files: no-PE events within edge_pro_band of the
% chl trigger edge, no larger than Syn + edge_pro_size_above_syn; flag 3) acts per
% file, against the old classes: events cut, how many the old classes called Euk,
% and the shape of the chl distribution of the small no-PE events just above the
% edge (does it fall from the edge = a pile, or rise to a cluster peak above it?).
% Early cruises (output/early_layout_robust file lists) and post-2019 cruises.
% Saves output/early_layout_robust/pro_tail.mat
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
    end
    cfg = classifier_config(); cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
    for k = 1:numel(files)
        r = strcmp(FCSfileinfo.fcslist, files(k));
        if FCSfileinfo.pro_measured(r), continue, end
        [d, h] = fca_readfcs(fullfile(root, C{i}, 'FCS', char(files(k))));
        ref = bead_reference(bs, FCSfileinfo.matdate_start(r), cfg);
        [c, info] = assign_class_v2(d, h, ref, cfg, 0);
        if ~isfield(info.trigger_edges, 'BL3') || info.ncomp == 0, continue, end
        o = load(fullfile(bp, 'class', regexprep(char(files(k)), '\.fcs$', '.mat'))).class(:, 1);
        F = info.features; e = info.trigger_edges.BL3; ss = info.syn_ref.size;
        pt = info.flag == 3;
        % small no-PE events from the edge up to 1 log10 above it
        nope = F.GL2 < info.rules.pe_min | F.PEchl < cfg.sig.pe_min_PEchl;
        cand = nope & F.size < ss + cfg.sig.edge_pro_size_above_syn & F.BL3 > e & F.BL3 < e + 1 & ismember(c, [0 1]) & F.SSC > info.floor.SSC + 0.1;
        pk = NaN; dr = NaN;
        if sum(cand) >= 50
            [fk, gk] = ksdensity(F.BL3(cand), 'Bandwidth', 0.05, 'NumPoints', 200);
            in = gk >= e & gk <= e + 1; fk = fk(in); gk = gk(in);
            [fm, ip] = max(fk); pk = gk(ip) - e;
            dr = fk(find(gk >= e + 0.05, 1)) / fm; %density just above the edge / peak
        end
        T = [T; table(string(cruise), files(k), info.edge_pile_frac, e, ss, sum(pt), sum(pt & o == 1), sum(o == 1), sum(c == 1), sum(cand), pk, dr, ...
            'VariableNames', {'cruise' 'file' 'pile_frac' 'edge' 'syn_size' 'n_cut' 'n_cut_old_euk' 'old_euk' 'v2_euk' 'n_cand' 'peak_above_edge' 'edge_to_peak'})]; %#ok<AGROW>
        fprintf('%s %s cut %d (old Euk %d) | old Euk %d v2 %d | cand %d peak %.2f above edge, edge/peak %.2f\n', cruise, extractBefore(files(k), '_Group'), T.n_cut(end), T.n_cut_old_euk(end), T.old_euk(end), T.v2_euk(end), T.n_cand(end), pk, dr);
    end
end
save(fullfile(out, 'pro_tail.mat'), 'T')
end
