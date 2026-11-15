function T = examine_peeuk(tags, outdir)
% function T = examine_peeuk(tags, outdir)
% Do PE-containing Euks fall into two groups (e.g. high-PE cryptophyte-type and
% low-PE with a Syn-like spectrum)? Pools v2 PE euks (class 3) and particle-
% associated Syn (class 8) from the files scored by score_agreement and plots
% PE:chl vs GL2:GL3, PE vs chl, and PE (GL3) vs size per cruise, flagging GL2-
% saturated events (whose GL2:GL3 is underestimated).
% Writes outdir/peeuk_<cruise>.png and returns per-event table T.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
T = table;
for t = tags
    a = load(fullfile(outdir, ['agreement_' t{1} '.mat']), 'A'); A = a.A;
    for i = 1:height(A)
        cpath = fullfile(root, char(A.cruise(i))); fn = char(A.file(i));
        bs = load_beadstat(fullfile(cpath, 'bead_calibrated'));
        [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
        ref = bead_reference(bs, A.time(i), cfg);
        [new, info] = assign_class_v2(d, h, ref, cfg, A.pro_measured(i));
        F = info.features;
        raw2 = d(:, strcmp({h.par.name}, 'GL2-H'));
        sat = raw2 >= max(raw2) * 0.98;
        k = find(ismember(new, [3 8]));
        if isempty(k), continue, end
        T = [T; table(repmat(string(t{1}), numel(k), 1), repmat(i, numel(k), 1), new(k), F.size(k), F.BL3(k), F.GL2(k), F.GL3(k), ...
            F.PEchl(k), F.GL2GL3(k), sat(k), F.GL3(k) > info.floor.GL3 + cfg.sig.ratio_min_above_floor, repmat(info.syn_ref.size, numel(k), 1), ...
            repmat(info.syn_ref.GL2GL3, numel(k), 1), 'VariableNames', {'cruise' 'filei' 'class' 'size' 'BL3' 'GL2' 'GL3' 'PEchl' 'GL2GL3' 'GL2sat' 'gl3ok' 'syn_size' 'syn_GL2GL3'})]; %#ok<AGROW>
    end
    C = T(T.cruise == t{1}, :);
    p3 = C.class == 3; ok = p3 & ~C.GL2sat & C.gl3ok;
    fprintf('%s: PE euk %d (GL2 saturated %d = %.0f%%), class 8 %d; Syn GL2:GL3 median %.2f\n', t{1}, sum(p3), sum(p3 & C.GL2sat), 100 * mean(C.GL2sat(p3)), sum(C.class == 8), median(C.syn_GL2GL3));
    if sum(ok) > 50
        X = [C.PEchl(ok) C.GL2GL3(ok)];
        g1 = fitgmdist(X, 1); g2 = fitgmdist(X, 2, 'Replicates', 5, 'RegularizationValue', 1e-4);
        [~, ix] = sort(g2.mu(:, 1));
        fprintf('   unsaturated PE euks %d: BIC 1 comp %.0f, 2 comp %.0f; 2-comp means PE:chl/GL2:GL3 %s, weights %s\n', sum(ok), g1.BIC, g2.BIC, ...
            mat2str(g2.mu(ix, :), 2), mat2str(g2.ComponentProportion(ix), 2));
        fprintf('   GL2:GL3 percentiles 10/50/90 %s; PE:chl %s; GL2 %s\n', mat2str(prctile(C.GL2GL3(ok), [10 50 90]), 2), mat2str(prctile(C.PEchl(ok), [10 50 90]), 2), mat2str(prctile(C.GL2(ok), [10 50 90]), 2));
    end
    f = figure('Visible', 'off', 'Position', [20 20 1500 450]); tiledlayout(1, 3, 'TileSpacing', 'compact');
    cl = {C.class == 3 & ~C.GL2sat, C.class == 3 & C.GL2sat, C.class == 8}; cc = [.85 0 .75; 0 0 0; .6 .1 .1];
    ln = {'PE euk' 'PE euk, GL2 saturated' 'Syn on particle'};
    pv = {'PEchl' 'GL2GL3' 'PE:chl' 'GL2:GL3'; 'BL3' 'GL2' 'chl (BL3)' 'PE (GL2)'; 'size' 'GL3' 'size (SSC/GL1)' 'GL3'};
    for v = 1:3
        nexttile; hold on
        for j = 1:3
            scatter(C.(pv{v, 1})(cl{j}), C.(pv{v, 2})(cl{j}), 4, cc(j, :), 'filled', 'MarkerFaceAlpha', .35)
        end
        if v == 1, yline(median(C.syn_GL2GL3), 'k--', 'Syn median'); end
        xlabel(pv{v, 3}); ylabel(pv{v, 4}); grid on; box on
        if v == 3, legend(ln, 'Location', 'southeast', 'FontSize', 7); end
    end
    sgtitle(sprintf('%s: v2 PE euks and Syn-on-particle, %d files', t{1}, height(A)))
    exportgraphics(f, fullfile(outdir, ['peeuk_' t{1} '.png']), 'Resolution', 90); close(f)
end
end
