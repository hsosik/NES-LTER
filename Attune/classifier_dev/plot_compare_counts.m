function plot_compare_counts(tags, outdir)
% function plot_compare_counts(tags, outdir)
% Per-file counts by product group, existing classifier vs assign_class_v2,
% from compare_v2 output (e.g., tags = {'AR92' 'EN727'}).
%
% Nov 2026, classifier redesign prototype

groups = {'syn' 'euk' 'pe_euk' 'pro'};
f = figure('Visible', 'off', 'Position', [50 50 1600 400 * numel(tags)]);
tl = tiledlayout(numel(tags), numel(groups), 'TileSpacing', 'compact');
for t = 1:numel(tags)
    load(fullfile(outdir, ['compare_v2_' tags{t} '.mat']), 'R', 'C')
    fprintf('\n%s (%d files)\n', tags{t}, height(R)); disp(C)
    for g = 1:numel(groups)
        nexttile; hold on
        o = R.(['old_' groups{g}]); n = R.(['new_' groups{g}]);
        tm = datetime(R.time, 'ConvertFrom', 'datenum');
        plot(tm, o, 'ok', 'MarkerSize', 5)
        plot(tm, n, '.r', 'MarkerSize', 12)
        set(gca, 'YScale', 'log'); grid on
        tot = sum(o); totn = sum(n);
        title(sprintf('%s %s: old %d, v2 %d (%+.0f%%)', tags{t}, groups{g}, tot, totn, 100 * (totn - tot) / max(tot, 1)), 'Interpreter', 'none')
        if g == 1, ylabel('events per file'); legend({'existing' 'v2'}, 'Location', 'best'); end
        r = (n + 1) ./ (o + 1);
        fprintf('%-7s files within 10%%: %d of %d; median v2/old %.3f\n', groups{g}, sum(abs(r - 1) < .1 & o > 20), sum(o > 20), median(r(o > 20)));
    end
end
exportgraphics(f, fullfile(outdir, 'compare_counts.png'), 'Resolution', 90)
end
