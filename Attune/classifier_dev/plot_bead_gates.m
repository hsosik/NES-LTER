function plot_bead_gates(fcsdat, fcshdr, class, pops, ttl, w)
% function plot_bead_gates(fcsdat, fcshdr, class, pops, ttl, w)
% One panel per bead population, on its gating channels: all events grey,
% gated beads coloured, search window dashed, expected position (+) and
% found centre (o). Title line per panel gives n, robust SD and QC.
% w - search window half-width (default 0.35)
%
% Nov 2026, bead stats redesign prototype

if nargin < 6, w = 0.35; end
names = {fcshdr.par.name};
col = lines(3);
clf
for ip = 1:numel(pops)
    p = pops(ip);
    subplot(1, numel(pops), ip)
    if numel(p.gatech) < 2, continue, end
    c1 = strcmp(names, p.gatech{1}); c2 = strcmp(names, p.gatech{2});
    X = log10(max(fcsdat(:, [find(c1) find(c2)]), 1));
    plot(X(:, 1), X(:, 2), '.', 'color', [.6 .6 .6], 'markersize', 2), hold on
    for j = 1:numel(pops)
        if any(class == j)
            plot(X(class == j, 1), X(class == j, 2), '.', 'color', col(j, :), 'markersize', 3)
        end
    end
    e = p.expected;
    if all(isfinite(e))
        plot(e(1) + w * [-1 1 1 -1 -1], e(2) + w * [-1 -1 1 1 -1], 'k--')
        plot(e(1), e(2), 'k+', 'markersize', 12, 'linewidth', 1.5)
    end
    if all(isfinite(p.center)), plot(p.center(1), p.center(2), 'ko', 'markersize', 10, 'linewidth', 1.5), end
    axis([0 6.1 0 6.1]), grid on, axis square
    xlabel([p.gatech{1} ' (log_{10})']), ylabel([p.gatech{2} ' (log_{10})'])
    title(sprintf('%s n=%d rsd=%.3f/%.3f QC=%d\n%s', p.name, p.n, p.rsd, p.QC, p.why), ...
        'fontsize', 8, 'interpreter', 'none')
end
if nargin > 4 && ~isempty(ttl), sgtitle(ttl, 'interpreter', 'none', 'fontsize', 9), end
end
