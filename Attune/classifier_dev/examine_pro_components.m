function P = examine_pro_components(tags, outdir)
% function P = examine_pro_components(tags, outdir)
% For every Pro-setting file scored by score_agreement for the cruises in tags,
% describe each mixture component named Pro: weight, events it gives class 7,
% spread (SD) in size and chl, and how many of the file's old-class Pro there are.
% Used to set a rule for when a file has a real Pro cluster.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
P = table;
for t = tags
    a = load(fullfile(outdir, ['agreement_' t{1} '.mat']), 'A'); A = a.A(a.A.pro_measured == 1, :);
    for i = 1:height(A)
        cpath = fullfile(root, char(A.cruise(i))); fn = char(A.file(i));
        bs = load_beadstat(fullfile(cpath, 'bead_calibrated'));
        [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
        ref = bead_reference(bs, A.time(i), cfg);
        [new, info] = assign_class_v2(d, h, ref, cfg, 1);
        c = load(fullfile(cpath, 'bead_calibrated', 'class', regexprep(fn, '\.fcs$', '.mat')), 'class');
        g = info.gmm; m = info.comp;
        % events given class 7 by each component (after all checks)
        hc = info.trace.has_chl; comp = zeros(size(new)); comp(hc) = cluster(g, [info.features.size(hc) info.features.BL3(hc) info.features.PEchl(hc) info.features.GL2GL3(hc)]);
        for j = find(m.label == 7)'
            S = g.Sigma(:, :, j);
            P = [P; table(string(t{1}), string(fn), j, m.weight(j), sum(new == 7 & comp == j), sum(new == 7), sum(c.class(:, 1) == 7), ...
                m.SSC(j), m.BL3(j), m.GL2(j), sqrt(S(1, 1)), sqrt(S(2, 2)), S(1, 2) / sqrt(S(1, 1) * S(2, 2)), info.n_has_chl, ...
                'VariableNames', {'cruise' 'file' 'comp' 'weight' 'n_pro_from_comp' 'n_pro_file' 'n_pro_old' 'size' 'BL3' 'GL2' 'sd_size' 'sd_BL3' 'r_size_BL3' 'n_has_chl'})]; %#ok<AGROW>
        end
        if ~any(m.label == 7)
            fprintf('%s: no Pro component (old Pro %d)\n', fn, sum(c.class(:, 1) == 7));
        end
    end
end
disp(P(:, 3:end))
save(fullfile(outdir, 'pro_components.mat'), 'P')
end
