function [D, P] = examine_syn_smear(cruisedir, nfiles, outdir)
% function [D, P] = examine_syn_smear(cruisedir, nfiles, outdir)
% Events assign_class_v2 calls Syn but the existing classifier left as 0
% (the winter "smear" toward larger scatter): apart from SSC, do they match
% the file's core Syn cluster? Compares, relative to each file's agreed Syn:
%  - every fluorescence channel and the pigment ratios (PE:chl, GL2:GL3, ...)
%  - scatter in other channels (FSC, GL1) and pulse shape (A/H, FSC-W)
%  - how pigment changes with scatter along the smear (doublet: both rise
%    ~0.3; Syn on a particle: scatter rises, pigment doesn't)
%  - how the smear's concentration varies with Syn and non-phyto particles
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

grp = {'core syn (both)' 'smear (old 0, v2 syn)' 'old-only syn (old syn, v2 0)' 'agreed 0, Syn size range'};
vars = {'size' 'SSC' 'FSC' 'GL1' 'BL1' 'BL2' 'BL3' 'GL2' 'GL3' 'GL4' 'PEchl' 'GL2GL3' 'chlSSC' 'PEsize' ...
    'GL4BL3' 'BL2BL3' 'BL1BL3' 'GL1SSC' 'FSCSSC' 'SSCAH' 'BL3AH' 'GL2AH' 'FSCW'};
P = cell(1, 4);
D = table;
H = []; Hn = []; nsmall = [];
for k = pick
    fn = Fi.fcslist{k};
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, Fi.matdate_start(k), cfg);
    [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, 0);
    c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class');
    old = c.class(:, 1);
    F = info.features;
    pn = {fcshdr.par.name};
    ah = @(ch) log10(max(fcsdat(:, strcmp(pn, [ch '-A'])), 1) ./ max(fcsdat(:, strcmp(pn, [ch '-H'])), 1));
    F.SSCAH = ah('SSC'); F.BL3AH = ah('BL3'); F.GL2AH = ah('GL2');
    F.FSCW = log10(max(fcsdat(:, strcmp(pn, 'FSC-W')), 1));
    F.PEsize = F.GL2 - F.size;
    F.GL4BL3 = F.GL4 - F.BL3; F.BL2BL3 = F.BL2 - F.BL3; F.BL1BL3 = F.BL1 - F.BL3;
    F.GL1SSC = F.GL1 - F.SSC; F.FSCSSC = F.FSC - F.SSC;
    core = old == 2 & new == 2;
    if sum(core) < 30, continue, end
    % express everything relative to this file's core Syn so files can be pooled
    c0 = varfun(@median, F(core, vars)); c0.Properties.VariableNames = vars;
    R = F(:, vars);
    for v = vars, R.(v{1}) = R.(v{1}) - c0.(v{1}); end
    synrange = R.size > -0.5 & R.size < 1.2 & F.BL3 > info.trigger_edges.BL3;
    % events with Syn's pigment spectrum and amount, at any scatter: does the smear stop at the v2 size cut?
    synspec = abs(R.PEchl) < 0.3 & abs(R.GL2GL3) < 0.1 & abs(R.BL3) < 0.5;
    sedges = -1:0.25:3;
    H(end + 1, :) = histcounts(R.size(synspec), sedges); %#ok<AGROW>
    Hn(end + 1, :) = histcounts(R.size(synspec & new == 2), sedges); %#ok<AGROW>
    % dim small events in the smear (scatter well below Syn)
    nsmall(end + 1) = sum(old == 0 & new == 2 & R.size < -0.5); %#ok<AGROW>
    g = {core, old == 0 & new == 2, old == 2 & new == 0, old == 0 & new == 0 & synrange};
    for j = 1:4
        idx = find(g{j});
        idx = idx(1:max(1, ceil(numel(idx) / 1500)):end);
        P{j} = [P{j}; R(idx, :)];
    end
    D = [D; table(string(fn), Fi.matdate_start(k), Fi.vol_analyzed(k), sum(g{1}), sum(g{2}), sum(g{3}), sum(g{4}), height(F), ...
        'VariableNames', {'file' 'time' 'vol' 'n_core' 'n_smear' 'n_oldonly' 'n_agreed0' 'n_events'})]; %#ok<AGROW>
end

% medians relative to core Syn (log10)
fprintf('\nmedian offset from the file''s core Syn (log10, bead units)\n%-10s', '')
fprintf(' %12s', 'core', 'smear', 'old-only', 'agreed 0'); fprintf('\n')
for v = vars
    fprintf('%-10s', v{1})
    for j = 1:4, fprintf(' %12.2f', median(P{j}.(v{1}), 'omitnan')); end
    fprintf('\n')
end
% spread: is the smear as tight as the core in each parameter?
fprintf('\nIQR (log10)\n')
for v = vars
    fprintf('%-10s', v{1})
    for j = 1:2, fprintf(' %12.2f', iqr(P{j}.(v{1}))); end
    fprintf('\n')
end
% along the smear: pigment vs scatter offset
fprintf('\nsmear binned by scatter offset from core Syn: median offsets\n%-10s %6s %7s %7s %7s %7s %7s %7s %7s\n', 'dSize', 'n', 'BL3', 'GL2', 'GL4', 'PEchl', 'GL2GL3', 'FSC', 'FSCW')
q = P{2}; edges = -0.2:0.2:1.2;
for b = 1:numel(edges) - 1
    i = q.size >= edges(b) & q.size < edges(b + 1);
    fprintf('%4.1f-%4.1f  %6d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n', edges(b), edges(b + 1), sum(i), ...
        median(q.BL3(i)), median(q.GL2(i)), median(q.GL4(i)), median(q.PEchl(i)), median(q.GL2GL3(i)), median(q.FSC(i)), median(q.FSCW(i)))
end
% fit pigment vs scatter slope within the smear
bb = [ones(height(q), 1) q.size] \ q.BL3;
bp2 = [ones(height(q), 1) q.size] \ q.GL2;
fprintf('\nsmear slope vs scatter offset: BL3 %.2f, GL2 %.2f (1 = pigment scales with scatter, 0 = constant pigment)\n', bb(2), bp2(2))

fprintf('\nevents with Syn pigment spectrum and amount, by scatter offset from core Syn (all files)\n%-10s %8s %8s\n', 'dSize', 'all', 'v2 Syn')
sedges = -1:0.25:3;
for b = 1:numel(sedges) - 1
    fprintf('%5.2f-%5.2f %8d %8d\n', sedges(b), sedges(b + 1), sum(H(:, b)), sum(Hn(:, b)))
end
q = P{2}(P{2}.size < -0.5, :);
fprintf('\nsmall dim smear events (scatter >0.5 below core Syn): %d total (%.1f%% of smear); median offsets BL3 %.2f, GL2 %.2f, PEchl %.2f, GL2GL3 %.2f, BL1BL3 %.2f, BL2BL3 %.2f, FSCW %.2f\n', ...
    sum(nsmall), 100 * sum(nsmall) / sum(D.n_smear), median(q.BL3), median(q.GL2), median(q.PEchl), median(q.GL2GL3), median(q.BL1BL3), median(q.BL2BL3), median(q.FSCW))

% concentrations
cs =D.n_core ./ D.vol; sm = D.n_smear ./ D.vol; a0 = D.n_agreed0 ./ D.vol;
fprintf('\nsmear fraction of v2 Syn: median %.2f, range %.2f-%.2f (%d files)\n', median(D.n_smear ./ (D.n_smear + D.n_core)), ...
    min(D.n_smear ./ (D.n_smear + D.n_core)), max(D.n_smear ./ (D.n_smear + D.n_core)), height(D))
fprintf('Spearman r: smear vs core Syn %.2f, vs non-phyto in Syn size range %.2f, vs core*non-phyto %.2f, vs total event rate %.2f\n', ...
    corr(sm, cs, 'Type', 'Spearman'), corr(sm, a0, 'Type', 'Spearman'), corr(sm, cs .* a0, 'Type', 'Spearman'), corr(sm, D.n_events ./ D.vol, 'Type', 'Spearman'))
fprintf('Spearman r: smear fraction vs core Syn %.2f, vs non-phyto %.2f, vs total event rate %.2f\n', ...
    corr(sm ./ (sm + cs), cs, 'Type', 'Spearman'), corr(sm ./ (sm + cs), a0, 'Type', 'Spearman'), corr(sm ./ (sm + cs), D.n_events ./ D.vol, 'Type', 'Spearman'))

% figure
cols = [.9 .5 0; .85 0 0; 0 .4 1; .6 .6 .6];
f = figure('Visible', 'off', 'Position', [40 40 1800 950]);
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
pairs = {'size' 'BL3'; 'size' 'GL2'; 'BL3' 'GL2'; 'GL2GL3' 'PEchl'; 'size' 'FSC'; 'size' 'FSCW'; 'BL3' 'BL2BL3'};
for p = 1:size(pairs, 1)
    nexttile; hold on
    for j = [4 1 3 2]
        scatter(P{j}.(pairs{p, 1}), P{j}.(pairs{p, 2}), 3, cols(j, :), 'filled', 'MarkerFaceAlpha', .3)
    end
    if p <= 2, plot([-0.5 1.5], [-0.5 1.5], 'k:'); end
    xlabel([pairs{p, 1} ' − core Syn']); ylabel([pairs{p, 2} ' − core Syn']); grid on
    if p == 1, legend(grp([4 1 3 2]), 'Location', 'northwest', 'FontSize', 7); end
end
nexttile; hold on
tm = datetime(D.time, 'ConvertFrom', 'datenum');
plot(tm, cs, 'o-', 'Color', cols(1, :)); plot(tm, sm, 'o-', 'Color', cols(2, :)); plot(tm, a0, 'o-', 'Color', cols(4, :))
set(gca, 'YScale', 'log'); ylabel('events per ml'); grid on
legend({'core Syn' 'smear' 'non-phyto, Syn size'}, 'Location', 'best', 'FontSize', 7)
title(tl, sprintf('%s: Syn smear vs core Syn, offsets from each file''s core Syn, %d files (dotted: 1:1)', cruisedir, height(D)), 'Interpreter', 'none')
exportgraphics(f, fullfile(outdir, ['syn_smear_' extractAfter(cruisedir, '_') '.png']), 'Resolution', 90)
end
