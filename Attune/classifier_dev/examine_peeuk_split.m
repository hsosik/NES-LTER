function T = examine_peeuk_split(tags, outdir)
% function T = examine_peeuk_split(tags, outdir)
% Pool v2 PE euks (class 3) with their old class (3 low PE / 4 high PE) and the
% candidate split variables, all relative to the file's Syn:
%   dBL2GL3  BL2:GL3 minus the Syn median (blue- vs green-excited PE; both channels
%            on scale for nearly all PE euks; cryptophyte-type PE sits below Syn)
%   dGL2GL3  GL2:GL3 minus the Syn median (only where GL2 is not saturated)
%   PEchl    PE:chl from the best on-scale PE channel (GL2, else BL2 scaled to GL2)
% Writes outdir/peeuk_split_<cruise>.png and returns the event table.
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
        c = load(fullfile(cpath, 'bead_calibrated', 'class', regexprep(fn, '\.fcs$', '.mat')), 'class'); old = c.class(:, 1);
        F = info.features; pn = {h.par.name};
        satG = d(:, strcmp(pn, 'GL2-H')) >= 0.98 * max(d(:, strcmp(pn, 'GL2-H')));
        satB = d(:, strcmp(pn, 'BL2-H')) >= 0.98 * max(d(:, strcmp(pn, 'BL2-H')));
        satG3 = d(:, strcmp(pn, 'GL3-H')) >= 0.98 * max(d(:, strcmp(pn, 'GL3-H')));
        BL2GL3 = F.BL2 - F.GL3;
        sy = new == 2 & ~satG & ~satB & F.GL3 > info.floor.GL3 + 1 & F.BL2 > info.floor.BL2 + 1;
        if sum(sy) < 50, continue, end
        synB = median(BL2GL3(sy)); synG = median(F.GL2GL3(sy));
        synBG = median(F.BL2(sy) - F.GL2(sy));
        k = find(new == 3);
        pe = F.GL2(k); pe(satG(k)) = F.BL2(k(satG(k))) - synBG;  %BL2 stands in for saturated GL2
        T = [T; table(repmat(string(t{1}), numel(k), 1), repmat(i, numel(k), 1), old(k), F.size(k), F.BL3(k), pe, pe - F.BL3(k), ...
            BL2GL3(k) - synB, F.GL2GL3(k) - synG, satG(k), satB(k) | satG3(k), ...
            'VariableNames', {'cruise' 'filei' 'old' 'size' 'BL3' 'PE' 'PEchl' 'dBL2GL3' 'dGL2GL3' 'GL2sat' 'BL2orGL3sat'})]; %#ok<AGROW>
    end
    C = T(T.cruise == t{1} & ~T.BL2orGL3sat, :);
    fprintf('\n%s: PE euks %d (old 3: %d, old 4: %d, old other: %d)\n', t{1}, height(C), sum(C.old == 3), sum(C.old == 4), sum(~ismember(C.old, [3 4])))
    for o = [3 4]
        g = C.old == o;
        fprintf('  old %d: dBL2GL3 p10/50/90 %s  dGL2GL3(unsat) %s  PEchl %s  PE %s\n', o, mat2str(prctile(C.dBL2GL3(g), [10 50 90]), 2), ...
            mat2str(prctile(C.dGL2GL3(g & ~C.GL2sat), [10 50 90]), 2), mat2str(prctile(C.PEchl(g), [10 50 90]), 2), mat2str(prctile(C.PE(g), [10 50 90]), 2))
    end
    if height(C) > 50
        X = [C.dBL2GL3 C.PEchl];
        g2 = fitgmdist(X, 2, 'Replicates', 5, 'RegularizationValue', 1e-4); g1 = fitgmdist(X, 1);
        [~, ix] = sort(g2.mu(:, 1));
        fprintf('  [dBL2GL3 PEchl]: BIC 1 comp %.0f, 2 comp %.0f; means %s weights %s\n', g1.BIC, g2.BIC, mat2str(g2.mu(ix, :), 2), mat2str(g2.ComponentProportion(ix), 2))
    end
    f = figure('Visible', 'off', 'Position', [20 20 1500 450]); tiledlayout(1, 3, 'TileSpacing', 'compact');
    grp = {C.old == 3, C.old == 4, ~ismember(C.old, [3 4])}; gc = [.2 .6 1; .85 0 .2; .6 .6 .6]; gn = {'old 3 (low PE)' 'old 4 (high PE)' 'old other'};
    pv = {'dBL2GL3' 'PEchl' 'BL2:GL3 - Syn' 'PE:chl (GL2, BL2 if saturated)'; 'BL3' 'PE' 'chl (BL3)' 'PE (GL2, BL2 if saturated)'; 'dBL2GL3' 'dGL2GL3' 'BL2:GL3 - Syn' 'GL2:GL3 - Syn'};
    for v = 1:3
        nexttile; hold on
        for j = [3 1 2]
            g = grp{j}; if v == 3, g = g & ~C.GL2sat; end
            scatter(C.(pv{v, 1})(g), C.(pv{v, 2})(g), 5, gc(j, :), 'filled', 'MarkerFaceAlpha', .4)
        end
        xlabel(pv{v, 3}); ylabel(pv{v, 4}); grid on; box on
        if v == 1, legend(gn([3 1 2]), 'Location', 'southwest', 'FontSize', 7); end
    end
    sgtitle(sprintf('%s: v2 PE euks coloured by old class', t{1}))
    exportgraphics(f, fullfile(outdir, ['peeuk_split_' t{1} '.png']), 'Resolution', 90); close(f)
end
end
