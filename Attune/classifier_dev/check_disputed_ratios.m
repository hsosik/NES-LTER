% check_disputed_ratios
% Brightness-independent check of the disputed events (old 0, v2 Euk):
% GL4 and BL3 both see chlorophyll fluorescence (green vs blue excitation),
% so chlorophyll-containing cells should share a GL4:BL3 ratio at any
% brightness. Uses output of discriminate_noise.
%
% Nov 2026, classifier redesign prototype

load(fullfile('output', 'noise_discrimination_EN727.mat'), 'X', 'y')
X.GL4BL3 = X.GL4_H - X.BL3_H;
X.BL2BL3 = X.BL2_H - X.BL3_H;
X.GL1SSC = X.GL1_H - X.SSC_H;
nm = {'agreed non-phyto' 'agreed euk' 'disputed'};
code = [0 1 2];
fprintf('%-18s %8s %8s %8s %8s %8s %8s\n', 'group', 'BL3', 'GL4', 'GL4-BL3', 'BL2-BL3', 'BL3 A/H', 'FSC-W')
for j = 1:3
    q = X(y == code(j), :);
    fprintf('%-18s %8.2f %8.2f %8.2f %8.2f %8.3f %8.2f\n', nm{j}, median(q.BL3_H), median(q.GL4_H), median(q.GL4BL3), median(q.BL2BL3), median(q.BL3_AH), median(q.FSC_W))
end
% euks binned by BL3 brightness: does GL4-BL3 depend on brightness?
e = X(y == 1, :); d = X(y == 2, :);
edges = -0.4:0.2:1.2;
fprintf('\nGL4-BL3 by BL3 brightness bin (median, n)\n%-10s %-16s %-16s %-16s\n', 'BL3', 'agreed euk', 'disputed', 'agreed non-phyto')
nz = X(y == 0, :);
for b = 1:numel(edges) - 1
    ie = e.BL3_H >= edges(b) & e.BL3_H < edges(b + 1);
    id = d.BL3_H >= edges(b) & d.BL3_H < edges(b + 1);
    in = nz.BL3_H >= edges(b) & nz.BL3_H < edges(b + 1);
    fprintf('%4.1f-%4.1f  %6.2f (%5d)   %6.2f (%5d)   %6.2f (%5d)\n', edges(b), edges(b + 1), median(e.GL4BL3(ie)), sum(ie), median(d.GL4BL3(id)), sum(id), median(nz.GL4BL3(in)), sum(in))
end

f = figure('Visible', 'off', 'Position', [40 40 1500 480]);
tiledlayout(1, 3, 'TileSpacing', 'compact')
cc = [.6 .6 .6; 0 .6 0; .85 0 0];
pairs = {'BL3_H' 'GL4_H'; 'BL3_H' 'GL4BL3'; 'GL4BL3' 'BL3_AH'};
for p = 1:3
    nexttile; hold on
    for j = 1:3
        q = X(y == code(j), :);
        scatter(q.(pairs{p, 1}), q.(pairs{p, 2}), 3, cc(j, :), 'filled', 'MarkerFaceAlpha', .3)
    end
    xlabel(pairs{p, 1}, 'Interpreter', 'none'); ylabel(pairs{p, 2}, 'Interpreter', 'none'); grid on
    if p == 1, legend(nm, 'Location', 'northwest'); end
end
exportgraphics(f, fullfile('output', 'disputed_chl_ratio_EN727.png'), 'Resolution', 90)
