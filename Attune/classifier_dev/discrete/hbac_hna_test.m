function R = hbac_hna_test(survey, tag, exclude)
% function R = hbac_hna_test(survey, tag, exclude)
% Test hbac_hna_split on the bacteria (class 9) of a hbac_hna_survey file.
% Per run: cut and method of the per-run rule, HNA fraction with that cut and
% with the fixed default cut, and the HNA fraction a 2-component GMM in
% [SSC BL1] gives (independent check). Runs listed in the survey file
% 'exclude' (e.g. the tuning set) are dropped.
% Draws output/discrete/hbac_hna_<tag>_test.png (fractions) and
% _examples.png (cytograms of fallback runs and the largest rule-vs-fixed
% differences, cut marked).
%
% Nov 2026, discrete-sample classifier development (HNA/LNA split)

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
load(fullfile(od, [survey '.mat']), 'R', 'E')
if nargin >= 3 && ~isempty(exclude)
    X = load(fullfile(od, [exclude '.mat']), 'R');
    drop = ismember(R.file, X.R.file);
    fprintf('dropping %d runs also in %s\n', sum(drop), exclude)
else
    drop = false(height(R), 1);
end
nr = height(R); R.run = (1:nr)';
[R.cut, R.dip, R.hna, R.hna_fixed, R.hna_gmm, R.nb] = deal(nan(nr, 1)); R.method = strings(nr, 1);
opts = statset('MaxIter', 500);
for r = find(~drop)'
    e = E(E.run == r, :); b = e.class == 9;
    R.nb(r) = sum(b);
    if sum(b) < 200, continue, end
    [bt, s] = hbac_hna_split(double(e.BL1), double(e.SSC), b);
    R.cut(r) = s.cut; R.dip(r) = s.dip; R.method(r) = s.method; R.hna(r) = s.hna_frac;
    [~, s0] = hbac_hna_split(double(e.BL1), double(e.SSC), b, struct('slope', 0, 'win', [0 0], 'flank', 0.35, 'bw', 0.04, 'min_n', Inf, 'max_dip', 0, 'default_cut', -1.27));
    R.hna_fixed(r) = s0.hna_frac;
    X = double([e.SSC(b) e.BL1(b)]); if height(X) > 8000, rng(r); X = X(randperm(height(X), 8000), :); end
    try
        w = warning('off', 'all');
        g = fitgmdist(X, 2, 'Replicates', 3, 'Options', opts, 'RegularizationValue', 1e-4); warning(w)
        [~, hi] = max(g.mu(:, 2)); R.hna_gmm(r) = g.ComponentProportion(hi);
    catch
    end
end
R = R(~drop & isfinite(R.cut), :);
v = R.method == "valley";
fprintf('%d runs: valley cut in %d (%.0f%%), fixed fallback in %d\n', height(R), sum(v), 100 * mean(v), sum(~v))
fprintf('valley cuts: median %.3f, 5-95%% %.2f to %.2f\n', median(R.cut(v)), prctile(R.cut(v), [5 95]))
fprintf('HNA fraction: median %.2f (5-95%% %.2f-%.2f)\n', median(R.hna), prctile(R.hna, [5 95]))
fprintf('|rule - fixed| HNA fraction: median %.3f, 90th %.3f, max %.3f\n', median(abs(R.hna - R.hna_fixed)), prctile(abs(R.hna - R.hna_fixed), 90), max(abs(R.hna - R.hna_fixed)))
fprintf('rule vs GMM HNA fraction: r = %.2f, median |diff| %.3f\n', corr(R.hna, R.hna_gmm, 'rows', 'complete'), median(abs(R.hna - R.hna_gmm), 'omitnan'))
G = groupsummary(R, 'cruise', 'median', {'cut' 'hna' 'dip'}); G.n_fixed = groupsummary(R, 'cruise', @(m) sum(m ~= "valley"), 'method').fun1_method;
disp(G)
save(fullfile(od, ['hbac_hna_' tag '_test.mat']), 'R')

fig = figure('visible', 'off', 'position', [50 50 1300 400]);
subplot(1, 3, 1), plot(R.hna_fixed(v), R.hna(v), 'b.', R.hna_fixed(~v), R.hna(~v), 'ro'), refline(1, 0), grid on
xlabel('HNA fraction, flat fixed cut -1.27'), ylabel('HNA fraction, per-run rule'), legend('valley', 'fixed fallback', 'location', 'northwest')
subplot(1, 3, 2), plot(R.hna_gmm, R.hna, 'k.'), refline(1, 0), grid on, xlabel('HNA fraction, 2-group GMM'), ylabel('per-run rule')
subplot(1, 3, 3), histogram(R.cut(v), -1.46:0.02:-1.08), xline(-1.27, 'r'), xlabel('valley cut (BL1)'), ylabel('runs'), grid on
exportgraphics(fig, fullfile(od, ['hbac_hna_' tag '_test.png']), 'resolution', 90), close(fig)

% examples: all fallbacks (up to 8), then the largest rule-vs-fixed differences, then random valley runs
[~, o] = sort(abs(R.hna - R.hna_fixed), 'descend');
k = [find(~v); o(v(o))]; k = unique(k(1:min(12, end)), 'stable');
rng(7); rest = setdiff(find(v), k); k = [k; rest(randperm(numel(rest), min(20 - numel(k), numel(rest))))];
fig = figure('visible', 'off', 'position', [20 20 1600 1250]);
for m = 1:numel(k)
    r = k(m); subplot(4, 5, m)
    e = E(E.run == R.run(r), :); b = e.class == 9;
    plot(e.SSC(~b), e.BL1(~b), '.', 'color', [.8 .8 .8], 'markersize', 1), hold on
    lo = b & e.BL1 - 0.35 * (e.SSC + 2) <= R.cut(r);
    plot(e.SSC(lo), e.BL1(lo), '.', 'color', [.2 .4 .9], 'markersize', 1)
    plot(e.SSC(b & ~lo), e.BL1(b & ~lo), '.', 'color', [.9 .35 .1], 'markersize', 1)
    xs = [-3.2 -0.9]; plot(xs, R.cut(r) + 0.35 * (xs + 2), 'k-'); yline(-1.27, 'k:');
    xlim([-3.2 -0.9]), ylim([-2.1 -0.2]), grid on
    title(sprintf('%s %s %gm\n%s dip %.2f HNA %.2f (fixed %.2f)', R.cruise(r), extractBetween(R.file(r), 'SSC_', '.fcs'), round(R.depth(r)), R.method(r), R.dip(r), R.hna(r), R.hna_fixed(r)), 'fontsize', 8, 'interpreter', 'none')
end
annotation('textbox', [0 0.95 1 0.04], 'string', 'bacteria (class 9): blue LNA, orange HNA; solid = cut used (tilted), dotted = flat fixed -1.27; gray other events. x SSC, y BL1 (bead-relative log10)', 'edgecolor', 'none', 'horizontalalignment', 'center')
exportgraphics(fig, fullfile(od, ['hbac_hna_' tag '_examples.png']), 'resolution', 90), close(fig)
end

