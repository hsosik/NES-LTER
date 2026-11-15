function plot_gl1_lowsignal(G)
% function plot_gl1_lowsignal(G)
% Figure for gl1_lowsignal_check: bead-normalized GL1 (OD2) - SSC (no OD2)
% against the raw GL1 signal, for all particles (binned per file) and for
% Syn (one point per file). Saved as output/pt_transfer/gl1_lowsignal.png.
%
% Nov 2026, PT transfer / OD2 low-signal check

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'pt_transfer');
if nargin < 1, s = load(fullfile(od, 'gl1_lowsignal.mat'), 'G'); G = s.G; end
B = G.bins; S = G.syn;
cr = unique(B.cruise); col = turbo(numel(cr));
f = figure('visible', 'off', 'position', [0 0 1500 550]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');

nexttile; hold on
for i = 1:numel(cr)
    q = B.cruise == cr(i);
    plot(B.raw_gl1(q), B.dGL1(q), '.', 'color', col(i, :), 'markersize', 4, 'HandleVisibility', 'off')
end
e = 1:0.1:5; [~, ~, bi] = histcounts(B.raw_gl1, e);
m = accumarray(bi(bi > 0), B.dGL1(bi > 0), [numel(e) - 1 1], @median, NaN);
plot(e(1:end-1) + 0.05, m, 'k-', 'linewidth', 2, 'DisplayName', 'median, all cruises')
plot(S.raw_gl1, S.d, 'kp', 'markerfacecolor', 'y', 'DisplayName', 'Syn (per file)')
xlabel('raw GL1-H (log10 counts)'); ylabel('GL1_{norm} - SSC_{norm} (log10)')
title('All particles (bins of SSC) and Syn'); grid on; ylim([-0.6 1]); legend('location', 'northeast')

nexttile; hold on
for i = 1:numel(cr)
    q = B.cruise == cr(i);
    plot(B.ssc_norm(q), B.dGL1(q), '.', 'color', col(i, :), 'markersize', 4, 'DisplayName', cr(i))
end
xlabel('SSC_{norm} (log10, 1 \mum bead = 0)'); ylabel('GL1_{norm} - SSC_{norm}')
title('Same, against particle size (SSC)'); grid on; ylim([-0.6 1])
legend('location', 'eastoutside', 'fontsize', 7)

nexttile; hold on
for i = 1:numel(cr)
    q = S.cruise == cr(i);
    plot(S.raw_gl1(q), S.d(q), 'o', 'color', col(i, :), 'markerfacecolor', col(i, :), 'markersize', 4)
end
xlabel('Syn raw GL1-H (log10)'); ylabel('Syn GL1_{norm} - SSC_{norm}')
title('Syn per file, by cruise'); grid on
title(tl, 'Side scatter through OD2 (GL1, 532 nm) vs without OD2 (SSC, 488 nm), each normalized to the 1 \mum bead on its own channel')
exportgraphics(f, fullfile(od, 'gl1_lowsignal.png'), 'Resolution', 110)
end
