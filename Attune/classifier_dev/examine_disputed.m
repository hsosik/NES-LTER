function D = examine_disputed(cruisedir, nfiles, outdir)
% function D = examine_disputed(cruisedir, nfiles, outdir)
% Events the existing classifier left as 0 but assign_class_v2 calls Euk:
% are they real small eukaryotes or noise? Compares them with agreed Euks,
% agreed Syn and agreed non-phytoplankton on cluster shape, pigment
% spectrum, pulse shape, and how their numbers vary through the cruise.
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

grp = {'agreed euk' 'disputed (old 0, v2 euk)' 'agreed syn' 'agreed 0'};
P = cell(1, 4);   %pooled events per group (subsampled)
D = table;
for k = pick
    fn = Fi.fcslist{k};
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, Fi.matdate_start(k), cfg);
    [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, 0);
    c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class');
    old = c.class(:, 1);
    F = info.features;
    pn = {fcshdr.par.name};
    F.BL3AH = log10(max(fcsdat(:, strcmp(pn, 'BL3-A')), 1) ./ max(fcsdat(:, strcmp(pn, 'BL3-H')), 1));
    F.FSCW = fcsdat(:, strcmp(pn, 'FSC-W'));
    F.BL1BL3 = F.BL1 - F.BL3;
    g = {old == 1 & new == 1, old == 0 & new == 1, old == 2 & new == 2, old == 0 & new == 0};
    for j = 1:4
        idx = find(g{j});
        idx = idx(1:max(1, ceil(numel(idx) / 1500)):end);
        P{j} = [P{j}; F(idx, {'size' 'SSC' 'BL1' 'BL3' 'GL2' 'chlSSC' 'BL1BL3' 'BL3AH' 'FSCW'})];
    end
    D = [D; table(string(fn), Fi.matdate_start(k), Fi.vol_analyzed(k), sum(g{1}), sum(g{2}), sum(g{3}), sum(g{4}), ...
        'VariableNames', {'file' 'time' 'vol' 'n_agreed_euk' 'n_disputed' 'n_agreed_syn' 'n_agreed_0'})]; %#ok<AGROW>
end

% summary stats
fprintf('\n%-28s %8s %8s %8s %8s %8s %8s\n', 'group', 'size', 'BL3', 'chl:size', 'BL1-BL3', 'A/H BL3', 'FSC-W')
for j = 1:4
    q = P{j};
    fprintf('%-28s %8.2f %8.2f %8.2f %8.2f %8.3f %8.1f\n', grp{j}, median(q.size), median(q.BL3), median(q.chlSSC), median(q.BL1BL3), median(q.BL3AH), median(q.FSCW))
end
r1 = corr(D.n_disputed ./ D.vol, D.n_agreed_euk ./ D.vol, 'Type', 'Spearman');
r2 = corr(D.n_disputed ./ D.vol, D.n_agreed_0 ./ D.vol, 'Type', 'Spearman');
r3 = corr(D.n_disputed ./ D.vol, D.n_agreed_syn ./ D.vol, 'Type', 'Spearman');
fprintf('\ndisputed concentration vs agreed euk r=%.2f, vs agreed syn r=%.2f, vs agreed non-phyto r=%.2f (Spearman, %d files)\n', r1, r3, r2, height(D))
fprintf('disputed as fraction of v2 euks: median %.2f, range %.2f-%.2f\n', median(D.n_disputed ./ (D.n_disputed + D.n_agreed_euk)), min(D.n_disputed ./ (D.n_disputed + D.n_agreed_euk)), max(D.n_disputed ./ (D.n_disputed + D.n_agreed_euk)))

% figure: projections and distributions
cols = [0 .6 0; .85 0 0; .9 .5 0; .6 .6 .6];
f = figure('Visible', 'off', 'Position', [40 40 1700 950]);
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
pairs = {'size' 'BL3'; 'BL3' 'GL2'; 'BL3' 'BL1'; 'size' 'FSCW'};
for p = 1:4
    nexttile; hold on
    for j = [4 3 1 2]
        scatter(P{j}.(pairs{p, 1}), P{j}.(pairs{p, 2}), 3, cols(j, :), 'filled', 'MarkerFaceAlpha', .35)
    end
    xlabel(pairs{p, 1}); ylabel(pairs{p, 2}); grid on
    if p == 1, legend(grp([4 3 1 2]), 'Location', 'northwest', 'FontSize', 7); end
end
hv = {'chlSSC' 'chl per size (log10, bead units)'; 'BL1BL3' 'green / chl fluorescence (BL1-BL3)'; 'BL3AH' 'BL3 pulse area/height (log10)'};
for h = 1:3
    nexttile; hold on
    for j = [1 2 4]
        histogram(P{j}.(hv{h, 1}), 60, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', cols(j, :), 'LineWidth', 1.5)
    end
    xlabel(hv{h, 2}); grid on
end
nexttile; hold on
tm = datetime(D.time, 'ConvertFrom', 'datenum');
plot(tm, D.n_agreed_euk ./ D.vol, 'o-', 'Color', cols(1, :))
plot(tm, D.n_disputed ./ D.vol, 'o-', 'Color', cols(2, :))
plot(tm, D.n_agreed_0 ./ D.vol, 'o-', 'Color', cols(4, :))
set(gca, 'YScale', 'log'); ylabel('events per ml'); grid on
legend({'agreed euk' 'disputed' 'agreed non-phyto'}, 'Location', 'best', 'FontSize', 7)
title(sprintf('disputed vs euk r=%.2f, vs non-phyto r=%.2f', r1, r2))
title(tl, sprintf('%s: events old classifier = 0 but v2 = Euk (red), %d files', cruisedir, height(D)), 'Interpreter', 'none')
exportgraphics(f, fullfile(outdir, ['disputed_' extractAfter(cruisedir, '_') '.png']), 'Resolution', 90)
end
