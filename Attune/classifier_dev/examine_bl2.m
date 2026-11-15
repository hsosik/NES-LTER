function S = examine_bl2(tags, outdir)
% function S = examine_bl2(tags, outdir)
% Is BL2 (blue-laser ~574 nm, PE excited less efficiently than by the green laser)
% a usable PE measure where GL2 saturates? For v2 Syn (2), PE euks (3) and
% Syn-on-particle (8) pooled per cruise: how often GL2 and BL2 saturate, how
% tightly BL2 tracks GL2 where both are on scale, and whether the BL2:GL2 ratio
% differs between Syn and PE euks. Writes outdir/bl2_<cruise>.png.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
S = table;
for t = tags
    a = load(fullfile(outdir, ['agreement_' t{1} '.mat']), 'A'); A = a.A;
    C = table;
    for i = 1:height(A)
        cpath = fullfile(root, char(A.cruise(i))); fn = char(A.file(i));
        bs = load_beadstat(fullfile(cpath, 'bead_calibrated'));
        [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
        ref = bead_reference(bs, A.time(i), cfg);
        [new, info] = assign_class_v2(d, h, ref, cfg, A.pro_measured(i));
        F = info.features; pn = {h.par.name};
        r2 = d(:, strcmp(pn, 'GL2-H')); rb = d(:, strcmp(pn, 'BL2-H'));
        satG = r2 >= max(r2) * 0.98; satB = rb >= max(rb) * 0.98;
        k = find(ismember(new, [2 3 8]));
        C = [C; table(repmat(i, numel(k), 1), new(k), F.size(k), F.BL3(k), F.GL2(k), F.GL3(k), F.BL2(k), satG(k), satB(k), ...
            F.GL3(k) > info.floor.GL3 + cfg.sig.ratio_min_above_floor, F.BL2(k) > info.floor.BL2 + 1, ...
            'VariableNames', {'filei' 'class' 'size' 'BL3' 'GL2' 'GL3' 'BL2' 'GL2sat' 'BL2sat' 'gl3ok' 'bl2ok'})]; %#ok<AGROW>
    end
    C.BL2GL2 = C.BL2 - C.GL2;
    pe = C.class == 3; sy = C.class == 2;
    both = ~C.GL2sat & ~C.BL2sat & C.bl2ok;
    fprintf('\n%s: PE euks %d, GL2 saturated %d; of those BL2 saturated %d, BL2 > floor+1 %d\n', t{1}, sum(pe), sum(pe & C.GL2sat), sum(pe & C.GL2sat & C.BL2sat), sum(pe & C.GL2sat & C.bl2ok))
    for cl = [2 3 8]
        g = C.class == cl & both;
        if sum(g) < 20, continue, end
        p = polyfit(C.GL2(g), C.BL2(g), 1); res = C.BL2(g) - polyval(p, C.GL2(g));
        fprintf('  class %d on scale in both (%d): BL2-GL2 median %.2f (p10 %.2f, p90 %.2f); BL2 on GL2 slope %.2f, residual SD %.2f\n', ...
            cl, sum(g), median(C.BL2GL2(g)), prctile(C.BL2GL2(g), 10), prctile(C.BL2GL2(g), 90), p(1), std(res))
    end
    % spectral ratios side by side for unsaturated PE euks: BL2:GL2 vs GL2:GL3
    g = pe & both & C.gl3ok;
    if sum(g) > 20
        fprintf('  PE euks unsaturated: corr(BL2-GL2, GL2-GL3) %.2f\n', corr(C.BL2GL2(g), C.GL2(g) - C.GL3(g)))
    end
    S = [S; table(string(t{1}), sum(pe), sum(pe & C.GL2sat), sum(pe & C.GL2sat & C.BL2sat), 'VariableNames', {'cruise' 'n_peeuk' 'n_GL2sat' 'n_both_sat'})]; %#ok<AGROW>

    f = figure('Visible', 'off', 'Position', [20 20 1500 450]); tiledlayout(1, 3, 'TileSpacing', 'compact');
    grp = {sy & ~C.GL2sat, C.class == 8, pe & ~C.GL2sat, pe & C.GL2sat}; gc = [1 .55 0; .6 .1 .1; .85 0 .75; 0 0 0];
    gn = {'Syn' 'Syn on particle' 'PE euk' 'PE euk, GL2 saturated'};
    pv = {'GL2' 'BL2' 'PE (GL2)' 'BL2 (PE, blue laser)'; 'BL3' 'BL2' 'chl (BL3)' 'BL2 (PE, blue laser)'; 'GL3' 'BL2' 'GL3' 'BL2 (PE, blue laser)'};
    for v = 1:3
        nexttile; hold on
        for j = 1:4
            x = find(grp{j}); x = x(1:max(1, ceil(numel(x) / 20000)):end);
            scatter(C.(pv{v, 1})(x), C.(pv{v, 2})(x), 3, gc(j, :), 'filled', 'MarkerFaceAlpha', .3)
        end
        xlabel(pv{v, 3}); ylabel(pv{v, 4}); grid on; box on
        if v == 1, legend(gn, 'Location', 'northwest', 'FontSize', 7); end
    end
    sgtitle(sprintf('%s: BL2 as a PE measure, %d files', t{1}, height(A)))
    exportgraphics(f, fullfile(outdir, ['bl2_' t{1} '.png']), 'Resolution', 90); close(f)
end
disp(S)
end
