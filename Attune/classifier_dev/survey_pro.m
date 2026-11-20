function survey_pro(cruises, nper, outname)
% function survey_pro(cruises, nper, outname)
% Survey mixture components in the Pro window (dim, no PE, absolute size limits)
% across Pro-setting files of many cruises: Pro-named files that are Pro-capable by
% attune_pro_measured (voltage rule), nper spread through each cruise, classified
% with pro_measured = 1. One row per component (C: size, chl, chl:scatter and size
% relative to the file's Syn components, height above the chl trigger edge, SDs,
% old-class Pro inside it) and one per file (Fsum). Saves output/<outname>.mat
% e.g. survey_pro({'20231108_AR79' '20240401_AR82A' '20261019_AR104'}, 12, 'pro_survey')
%
% Nov 2026, classifier redesign prototype
sp = fullfile(fileparts(mfilename('fullpath')), 'output');
addpath(fileparts(mfilename('fullpath')), fullfile(fileparts(mfilename('fullpath')), '..'))
diary(fullfile(sp, [outname '.log'])); diary on
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
C = table; Fsum = table;
for ci = 1:numel(cruises)
    cd0 = cruises{ci};
    cpath = fullfile(root, cd0); bp = fullfile(cpath, 'bead_calibrated');
    cfg = classifier_config();
    cfg.cruise = extractAfter(cd0, '_');
    cfg.gl1_is_od2_scatter = strcmpi(attune_config(cfg.cruise).OD2, 'GL1');
    sig = cfg.sig;
    try
        bs = load_beadstat(bp);
    catch me
        fprintf('%s: no beadstat (%s)\n', cd0, me.message); continue
    end
    s = load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo'); Fi = s.FCSfileinfo;
    Fi = Fi(~contains(Fi.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete' 'PT'}), :);
    if ismember('QC_flag', Fi.Properties.VariableNames), Fi = Fi(Fi.QC_flag == 1, :); end
    Fi = Fi(contains(Fi.fcslist, {'Pro' 'pro'}), :);
    % voltage-based Pro capability (wrapper rule)
    Fi = attune_pro_measured(Fi, fullfile(cpath, 'FCS'));
    Fi = sortrows(Fi(Fi.pro_measured == 1, :), 'matdate_start');
    fprintf('%s: %d Pro-capable Pro-named files\n', cd0, height(Fi));
    if isempty(Fi), continue, end
    pick = unique(round(linspace(1, height(Fi), min(nper, height(Fi)))));
    for k = pick
        fn = Fi.fcslist{k};
        try
            [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
            ref = bead_reference(bs, Fi.matdate_start(k), cfg);
            [new, info] = assign_class_v2(d, h, ref, cfg, 1);
        catch me
            fprintf('  %s: error %s\n', fn, me.message); continue
        end
        cf = fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat'));
        nold = NaN;
        old7 = false(size(new));
        if exist(cf, 'file'), c = load(cf, 'class'); old7 = c.class(:, 1) == 7; nold = sum(old7); end
        F = info.features;
        e = NaN; if isfield(info.trigger_edges, 'BL3'), e = info.trigger_edges.BL3; end
        m = info.comp;
        if isempty(m) || ~ismember('label', m.Properties.VariableNames), continue, end
        syn = m.label == 2;
        sync = NaN; if any(syn), sync = sum(m.chlSSC(syn) .* m.weight(syn)) / sum(m.weight(syn)); end
        synsize = NaN; if any(syn), synsize = sum(m.SSC(syn) .* m.weight(syn)) / sum(m.weight(syn)); end
        synbl3 = NaN; if any(syn), synbl3 = sum(m.BL3(syn) .* m.weight(syn)) / sum(m.weight(syn)); end
        pe_min = sig.pe_min_GL2;
        win = m.SSC < sig.pro_max_SSC & m.BL3 < sig.pro_max_BL3 & ~(m.PEchl > sig.pe_min_PEchl & m.GL2 > pe_min) & m.chlSSC >= sig.detritus_max_chlSSC;
        hc = info.trace.has_chl;
        X = [F.size(hc) F.BL3(hc) F.PEchl(hc) F.GL2GL3(hc)];
        comp = zeros(size(new)); comp(hc) = cluster(info.gmm, X);
        for j = find(win)'
            S = info.gmm.Sigma(:, :, j);
            inj = comp == j;
            C = [C; table(string(cd0), string(fn), Fi.matdate_start(k), j, m.weight(j), sum(inj), sum(new == 7 & inj), sum(new == 7), nold, ...
                m.SSC(j), m.BL3(j), m.GL2(j), m.chlSSC(j), sqrt(S(1, 1)), sqrt(S(2, 2)), m.label(j), sync, synsize, synbl3, e, info.n_has_chl, ...
                sum(inj & old7), ...
                'VariableNames', {'cruise' 'file' 'time' 'comp' 'weight' 'n_comp' 'n_pro_comp' 'n_pro_file' 'n_old_pro' ...
                'size' 'BL3' 'GL2' 'chlSSC' 'sd_size' 'sd_BL3' 'label' 'syn_chlSSC' 'syn_size' 'syn_BL3' 'edge' 'n_has_chl' 'n_oldpro_in_comp'})]; %#ok<AGROW>
        end
        Fsum = [Fsum; table(string(cd0), string(fn), Fi.matdate_start(k), sum(new == 7), nold, sum(win), sync, synsize, synbl3, e, ...
            'VariableNames', {'cruise' 'file' 'time' 'n_pro' 'n_old_pro' 'n_wincomp' 'syn_chlSSC' 'syn_size' 'syn_BL3' 'edge'})]; %#ok<AGROW>
        fprintf('  %s: v2 Pro %d, old %d, window comps %d\n', fn, sum(new == 7), nold, sum(win));
    end
    save(fullfile(sp, [outname '.mat']), 'C', 'Fsum')
end
save(fullfile(sp, [outname '.mat']), 'C', 'Fsum')
diary off
end
