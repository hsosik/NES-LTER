function S = examine_synpart(cruisedir, nfiles, outdir)
% function S = examine_synpart(cruisedir, nfiles, outdir)
% Is class 8 (particle-associated Syn) a separate group, or just the large
% end of the file's Syn cluster? Two tests, relative to each file's Syn:
%  1. pigment vs scatter: within Syn, bigger cells carry more pigment. If class 8
%     is the Syn tail, its pigment continues that trend; if it is one Syn on a
%     larger particle, its pigment stays near one Syn while scatter rises.
%  2. scatter distribution: compare the number of Syn-spectrum events above the
%     class-8 size cut with what the Syn cluster's own (Gaussian) shape predicts.
% Also compares FSC, GL1, FSC-W and pulse shape.
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
Fi = sortrows(Fi, 'matdate_start');
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
pick = unique(round(linspace(1, height(Fi), min(nfiles, height(Fi)))));

vars = {'size' 'FSC' 'GL1' 'BL3' 'GL2' 'GL4' 'PEchl' 'GL2GL3' 'FSCW' 'SSCAH' 'BL3AH'};
P = cell(1, 2); S = table;
sedges = -1:0.1:3;
Hs = zeros(1, numel(sedges) - 1); H8 = Hs;
for k = pick
    fn = Fi.fcslist{k};
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, Fi.matdate_start(k), cfg);
    [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, Fi.pro_measured(k));
    syn = new == 2; sp = new == 8;
    if sum(syn) < 100, continue, end
    F = info.features;
    pn = {fcshdr.par.name};
    ah = @(ch) log10(max(fcsdat(:, strcmp(pn, [ch '-A'])), 1) ./ max(fcsdat(:, strcmp(pn, [ch '-H'])), 1));
    F.SSCAH = ah('SSC'); F.BL3AH = ah('BL3');
    F.FSCW = log10(max(fcsdat(:, strcmp(pn, 'FSC-W')), 1));
    c0 = varfun(@median, F(syn, vars)); c0.Properties.VariableNames = vars;
    R = F(:, vars);
    for v = vars, R.(v{1}) = R.(v{1}) - c0.(v{1}); end
    % test 1: slope of pigment on scatter within Syn (robust: binned medians)
    b = robustfit(R.size(syn), R.BL3(syn)); slope = b(2);
    % pigment class 8 would have if it continued the Syn trend
    pred8 = slope * median(R.size(sp));
    % test 2: Gaussian tail of Syn scatter above the class-8 cut
    sd = 1.4826 * mad(R.size(syn), 1);
    expect = sum(syn) * (1 - normcdf(max(cfg.sig.synpart_min_size_above, cfg.sig.synpart_size_sd * sd), 0, sd));
    Hs = Hs + histcounts(R.size(syn), sedges); H8 = H8 + histcounts(R.size(sp), sedges);
    S = [S; table(string(fn), Fi.matdate_start(k), sum(syn), sum(sp), sd, slope, median(R.size(sp)), median(R.BL3(sp)), pred8, expect, ...
        'VariableNames', {'file' 'time' 'n_syn' 'n_synpart' 'syn_size_sd' 'BL3_slope_in_syn' 'synpart_size' 'synpart_BL3' 'synpart_BL3_if_tail' 'tail_expected'})]; %#ok<AGROW>
    grp = {syn, sp};
    for j = 1:2
        g = find(grp{j});
        g = g(1:max(1, ceil(numel(g) / 1500)):end);
        P{j} = [P{j}; R(g, :)];
    end
end
fprintf('\n%s, %d files\n', cruisedir, height(S))
fprintf('class 8 as fraction of Syn: median %.3f (range %.3f-%.3f)\n', median(S.n_synpart ./ S.n_syn), min(S.n_synpart ./ S.n_syn), max(S.n_synpart ./ S.n_syn))
[~, im] = max(S.n_synpart ./ S.n_syn); fprintf('highest: %s (%d class 8, %d Syn)\n', S.file(im), S.n_synpart(im), S.n_syn(im))
fprintf('Syn scatter SD (robust) median %.2f log10\n', median(S.syn_size_sd))
fprintf('test 1: chl vs scatter slope within Syn median %.2f; class 8 median offsets: size %.2f, BL3 %.2f, BL3 if it continued the Syn trend %.2f\n', ...
    median(S.BL3_slope_in_syn), median(S.synpart_size), median(S.synpart_BL3), median(S.synpart_BL3_if_tail))
fprintf('test 2: class 8 observed %d vs Gaussian Syn tail beyond the cut %.0f (ratio %.1f)\n', sum(S.n_synpart), sum(S.tail_expected), sum(S.n_synpart) / sum(S.tail_expected))
fprintf('\nmedian offset from the file''s Syn\n%-8s %8s %8s\n', '', 'syn', 'class 8')
for v = vars
    fprintf('%-8s %8.2f %8.2f\n', v{1}, median(P{1}.(v{1})), median(P{2}.(v{1})))
end
fprintf('\nclass 8 binned by scatter offset: median BL3, FSC offsets\n')
q = P{2};
for e = 0.3:0.2:1.9
    i = q.size >= e & q.size < e + 0.2;
    if sum(i) > 10, fprintf('%4.1f-%4.1f n=%5d BL3 %5.2f FSC %5.2f\n', e, e + 0.2, sum(i), median(q.BL3(i)), median(q.FSC(i))); end
end

f = figure('Visible', 'off', 'Position', [40 40 1600 420]);
tiledlayout(1, 4, 'TileSpacing', 'compact')
nexttile; hold on
scatter(P{1}.size, P{1}.BL3, 3, [.9 .5 0], 'filled', 'MarkerFaceAlpha', .2)
scatter(P{2}.size, P{2}.BL3, 4, [.8 0 0], 'filled', 'MarkerFaceAlpha', .4)
x = [-1 2]; plot(x, median(S.BL3_slope_in_syn) * x, 'k--')
xlabel('size − file Syn'); ylabel('BL3 − file Syn'); grid on; legend({'Syn' 'class 8' 'Syn chl–size trend'}, 'Location', 'northwest', 'FontSize', 7)
nexttile; hold on
scatter(P{1}.size, P{1}.FSC, 3, [.9 .5 0], 'filled', 'MarkerFaceAlpha', .2)
scatter(P{2}.size, P{2}.FSC, 4, [.8 0 0], 'filled', 'MarkerFaceAlpha', .4)
xlabel('size − file Syn'); ylabel('FSC − file Syn'); grid on
nexttile; hold on
ctr = sedges(1:end - 1) + 0.05;
stairs(ctr, Hs + 0.5, 'Color', [.9 .5 0]); stairs(ctr, H8 + 0.5, 'Color', [.8 0 0])
plot(ctr, sum(Hs) * 0.1 * normpdf(ctr, 0, median(S.syn_size_sd)), 'k--')
set(gca, 'YScale', 'log'); ylim([0.5 Inf]); xlabel('size − file Syn'); ylabel('events'); grid on
legend({'Syn' 'class 8' 'Gaussian fit to Syn'}, 'FontSize', 7)
nexttile; hold on
scatter(P{1}.size, P{1}.FSCW, 3, [.9 .5 0], 'filled', 'MarkerFaceAlpha', .2)
scatter(P{2}.size, P{2}.FSCW, 4, [.8 0 0], 'filled', 'MarkerFaceAlpha', .4)
xlabel('size − file Syn'); ylabel('FSC-W − file Syn'); grid on
sgtitle(sprintf('%s: class 8 vs the file''s Syn, %d files', cruisedir, height(S)), 'Interpreter', 'none')
exportgraphics(f, fullfile(outdir, ['synpart_' extractAfter(cruisedir, '_') '.png']), 'Resolution', 90)
end
