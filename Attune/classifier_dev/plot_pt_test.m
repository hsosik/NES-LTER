function plot_pt_test(fcsdat, fcshdr, P, conf, ttl)
% function plot_pt_test(fcsdat, fcshdr, P, conf, ttl)
% One Attune performance test (PT beads), for the bead review of cruises
% without FCB bead runs. Panels: FSC-H x SSC-H with the singlet gate
% (pt_test_peaks: within 0.15 log10 of the medians), then pairs of
% fluorescence channels for the singlets, showing the PT intensity levels;
% the events each channel's peak comes from (brightest level) are coloured
% and the peak is marked.
% P    - pt_test_peaks output
% conf - attune_config for the test's date
%
% Nov 2026, bead review for PT-only cruises

names = {fcshdr.par.name};
L = @(c) log10(max(fcsdat(:, c), 1));
pairs = {{'FSC' 'SSC'} {'BL1' 'BL3'} {'BL2' 'GL2'} {'GL3' 'GL4'}};
pairs = pairs(cellfun(@(p) ~isempty(P.ch.(p{1}).col) && ~isempty(P.ch.(p{2}).col), pairs));
g = P.sing;
clf
for k = 1:numel(pairs)
    a = P.ch.(pairs{k}{1}); b = P.ch.(pairs{k}{2});
    x = L(a.col); y = L(b.col);
    subplot(1, numel(pairs), k)
    if k == 1
        plot(x, y, '.', 'color', [.6 .6 .6], 'markersize', 2), hold on
        plot(x(g), y(g), '.', 'color', [1 .5 .05], 'markersize', 3)
        mx = median(x); my = median(y);
        plot(mx + 0.15 * [-1 1 1 -1 -1], my + 0.15 * [-1 -1 1 1 -1], 'k--')
        ttl1 = sprintf('%d events, %d singlets', numel(x), sum(g));
    else
        plot(x(g), y(g), '.', 'color', [.6 .6 .6], 'markersize', 2), hold on
        plot(x(a.use), y(a.use), '.', 'color', [.12 .47 .71], 'markersize', 3)
        plot(x(b.use), y(b.use), '.', 'color', [.89 .1 .11], 'markersize', 3)
        ttl1 = 'singlets; blue/red: brightest level (peak)';
    end
    xline(a.pk, 'k-'); yline(b.pk, 'k-');
    axis([0 6.1 0 6.1]), grid on, axis square
    xlabel(sprintf('%s = %s (log_{10}) %g V', pairs{k}{1}, names{a.col}, a.hv), 'interpreter', 'none')
    ylabel(sprintf('%s = %s (log_{10}) %g V', pairs{k}{2}, names{b.col}, b.hv), 'interpreter', 'none')
    title(sprintf('%s\npeak %.3f / %.3f, CV %.3f / %.3f', ttl1, a.pk, b.pk, a.rcv, b.rcv), 'fontsize', 8)
end
if nargin > 4 && ~isempty(ttl), sgtitle(sprintf('%s (layout %s)', ttl, conf.layout), 'interpreter', 'none', 'fontsize', 9), end
end
