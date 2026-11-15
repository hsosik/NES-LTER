function M = discriminate_noise(cruisedir, nfiles, outdir)
% function M = discriminate_noise(cruisedir, nfiles, outdir)
% Can noise be separated from small Euks using all Attune parameters?
% Trains classifiers on events both old and v2 agree on (Euk vs 0), within
% the region where they overlap (small, dim), using
%   (a) all parameters
%   (b) only parameters neither classifier uses (BL1, BL2, GL4, GL1,
%       pulse area/height ratios, FSC-W) -- independent evidence, since the
%       labels were made from SSC/FSC/BL3/GL2/GL3
% and reports cross-validated AUC, which parameters matter, and how the
% disputed events (old 0, v2 Euk) score.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
cpath = fullfile(root, cruisedir);
bp = fullfile(cpath, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & ~contains(FCSfileinfo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
Fi = Fi(cellfun(@(c) exist(fullfile(bp, 'class', regexprep(c, '\.fcs$', '.mat')), 'file') == 2, Fi.fcslist), :);
Fi = sortrows(Fi, 'matdate_start');
pick = unique(round(linspace(1, height(Fi), min(nfiles, height(Fi)))));

chans = cfg.chans;
X = []; y = []; fileid = [];
for k = pick
    fn = Fi.fcslist{k};
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, Fi.matdate_start(k), cfg);
    [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, 0);
    c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class');
    old = c.class(:, 1);
    F = info.features;
    pn = {fcshdr.par.name};
    T = table;
    for ch = chans
        T.([ch{1} '_H']) = F.(ch{1});
        a = fcsdat(:, strcmp(pn, [ch{1} '-A'])); h = fcsdat(:, strcmp(pn, [ch{1} '-H']));
        T.([ch{1} '_AH']) = log10(max(a, 1) ./ max(h, 1));   %pulse shape
    end
    T.FSC_W = log10(max(fcsdat(:, strcmp(pn, 'FSC-W')), 1));
    % group: 1 agreed euk, 0 agreed non-phyto, 2 disputed (old 0, v2 euk)
    lab = nan(size(old));
    lab(old == 1 & new == 1) = 1; lab(old == 0 & new == 0) = 0; lab(old == 0 & new == 1) = 2;
    % restrict to the overlap region: small and dim, where the decision is hard
    region = F.size < 0.5 & F.BL3 > -1 & F.BL3 < 1.2;
    keep = find(region & ~isnan(lab));
    keep = keep(1:max(1, ceil(numel(keep) / 3000)):end);
    X = [X; T(keep, :)]; y = [y; lab(keep)]; fileid = [fileid; k * ones(numel(keep), 1)]; %#ok<AGROW>
end
fprintf('overlap region events: agreed euk %d, agreed non-phyto %d, disputed %d\n', sum(y == 1), sum(y == 0), sum(y == 2))

used = {'SSC_H' 'FSC_H' 'BL3_H' 'GL2_H' 'GL3_H'};
vars = X.Properties.VariableNames;
sets = {vars, setdiff(vars, used, 'stable')};
setname = {'all parameters' 'independent parameters only'};
tr = y ~= 2;
M = struct();
f = figure('Visible', 'off', 'Position', [40 40 1600 500]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
for s2 = 1:2
    Xs = X{:, sets{s2}};
    Xs(~isfinite(Xs)) = 0;
    % cross-validate by file so neighbouring events in one file don't leak
    files = unique(fileid);
    foldOf = containers.Map(num2cell(files), num2cell(mod(0:numel(files) - 1, 5) + 1));
    fo = cellfun(@(v) foldOf(v), num2cell(fileid));
    score = nan(size(y));
    for kf = 1:5
        trn = tr & fo ~= kf; tst = fo == kf;
        mdl = fitcensemble(Xs(trn, :), y(trn), 'Method', 'Bag', 'NumLearningCycles', 100);
        [~, sc] = predict(mdl, Xs(tst, :));
        score(tst) = sc(:, 2);
    end
    [~, ~, ~, auc] = perfcurve(y(tr), score(tr), 1);
    mdl = fitcensemble(Xs(tr, :), y(tr), 'Method', 'Bag', 'NumLearningCycles', 100);
    imp = oobPermutedPredictorImportance(mdl);
    [~, o] = sort(imp, 'descend');
    M(s2).set = setname{s2}; M(s2).auc = auc;
    M(s2).top = sets{s2}(o(1:min(8, end)));
    M(s2).disputed_score = score(y == 2);
    fprintf('\n%s: cross-validated AUC (agreed euk vs agreed non-phyto) = %.3f\n', setname{s2}, auc)
    fprintf('  most informative: %s\n', strjoin(M(s2).top, ', '))
    d = score(y == 2);
    fprintf('  disputed events: %.0f%% score euk-like (>0.5), %.0f%% clearly euk (>0.8), %.0f%% clearly noise (<0.2)\n', ...
        100 * mean(d > .5), 100 * mean(d > .8), 100 * mean(d < .2))
    nexttile; hold on
    edges = 0:.025:1;
    histogram(score(y == 1), edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0 .6 0], 'LineWidth', 1.5)
    histogram(score(y == 0), edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [.5 .5 .5], 'LineWidth', 1.5)
    histogram(d, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [.85 0 0], 'LineWidth', 1.5)
    xlabel('P(euk), cross-validated'); title(sprintf('%s (AUC %.2f)', setname{s2}, auc)); grid on
    legend({'agreed euk' 'agreed non-phyto' 'disputed'}, 'Location', 'north')
end
% where do disputed events fall on the strongest independent parameters?
nexttile; hold on
top = M(2).top(1:2);
g = {y == 0, y == 1, y == 2}; cc = [.6 .6 .6; 0 .6 0; .85 0 0];
for j = 1:3
    scatter(X.(top{1})(g{j}), X.(top{2})(g{j}), 3, cc(j, :), 'filled', 'MarkerFaceAlpha', .3)
end
xlabel(top{1}, 'Interpreter', 'none'); ylabel(top{2}, 'Interpreter', 'none'); grid on
title('top two independent parameters')
title(tl, sprintf('%s: separating small Euks from noise in the overlap region (%d files)', cruisedir, numel(unique(fileid))), 'Interpreter', 'none')
exportgraphics(f, fullfile(outdir, ['noise_discrimination_' extractAfter(cruisedir, '_') '.png']), 'Resolution', 90)
save(fullfile(outdir, ['noise_discrimination_' extractAfter(cruisedir, '_') '.mat']), 'X', 'y', 'fileid', 'M')
end
