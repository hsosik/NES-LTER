function plot_perf_transfer(V, C)
% function plot_perf_transfer(V, C)
% Figures for the PT performance-test transfer (saved in output/pt_transfer):
%   1. validation: residual (observed FCB peak - prediction from the nearest
%      performance test) vs time, per channel, for the bead population used
%      by the classifier on that channel (cfg.bead_pop)
%   2. early cruises: Syn SSC normalized with the new references (and with
%      the current PT-conversion references where they exist) by month,
%      against the FCB cruises in output/syn_ssc_check.mat
%
% Nov 2026, PT transfer for early cruises

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev)
od = fullfile(dev, 'output', 'pt_transfer');
if nargin < 1 || isempty(V), s = load(fullfile(dev, 'output', 'perf_transfer_validation.mat'), 'V'); V = s.V; end
if nargin < 2 || isempty(C), s = load(fullfile(od, 'compare_perf_refs.mat'), 'C'); C = s.C; end
cfg = classifier_config();
P = {'05' '1' '6'};
chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL2' 'GL3' 'GL4'};

% 1. validation residuals
f = figure('visible', 'off', 'position', [0 0 1400 900]);
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
cr = unique(V.cruise); col = lines(numel(cr));
for k = 1:numel(chans)
    nexttile; hold on
    ch = chans{k}; pop = P{cfg.bead_pop.(ch)};
    for i = 1:numel(cr)
        g = V.chan == ch & V.pop == pop & V.cruise == cr(i);
        plot(V.time(g), V.resid(g), '.', 'color', col(i, :), 'DisplayName', cr(i) + " (" + V.state(find(g, 1)) + ")")
    end
    g = V.chan == ch & V.pop == pop;
    title(sprintf('%s (%s um): rsd %.3f', ch, strrep(pop, '05', '0.5'), 1.4826 * mad(V.resid(g), 1)))
    ylim([-0.12 0.12]); grid on; yline(0, 'k')
    if k == 1, legend('location', 'southwest', 'fontsize', 7), end
end
ylabel(tl, 'FCB peak - performance-test prediction (log10)')
title(tl, 'Performance tests as a transfer standard: FCB bead runs vs nearest performance test (offset per channel and optical state)')
exportgraphics(f, fullfile(od, 'validation_residuals.png'), 'Resolution', 110)

% 2. Syn SSC normalized, by month, vs the FCB cruises
s = load(fullfile(dev, 'output', 'syn_ssc_check.mat'), 'S'); S = s.S;
f = figure('visible', 'off', 'position', [0 0 1200 600]); hold on
fcb = S.source == "FCB";
tn = S.cruise == "TN368";
plot(S.month(fcb & ~tn) + 0.3 * (rand(sum(fcb & ~tn), 1) - 0.5), S.norm_ssc(fcb & ~tn), '.', 'color', [0.7 0.7 0.7], 'DisplayName', 'FCB cruises, SSC without OD2')
plot(S.month(tn) + 0.3 * (rand(sum(tn), 1) - 0.5), S.norm_ssc(tn), 'kp', 'markerfacecolor', 'k', 'DisplayName', 'TN368: FCB beads, OD2 on SSC')
ce = unique(C.cruise); col = lines(numel(ce));
for i = 1:numel(ce)
    g = C.cruise == ce(i) & isfinite(C.norm_new);
    plot(C.month(g) + 0.3 * (rand(sum(g), 1) - 0.5), C.norm_new(g), 'o', 'color', col(i, :), 'DisplayName', ce(i) + " perf transfer")
    if any(isfinite(C.norm_old(g)))
        plot(C.month(g) + 0.3 * (rand(sum(g), 1) - 0.5), C.norm_old(g), 'x', 'color', col(i, :), 'DisplayName', ce(i) + " PT conversion")
    end
end
xlim([0.5 12.5]); xticks(1:12); ylim([-1.6 -0.5]); grid on
text(0.6, -0.55, 'EN608, AR28B, AR29, EN617, AR31A, TN368: OD2 on SSC; EN627: no OD2 (degraded scatter, Feb 2019). Off-scale Syn gate failures not shown.', 'fontsize', 8)
xlabel('month'); ylabel('Syn SSC - 1 um bead SSC (log10)')
legend('location', 'eastoutside', 'fontsize', 8)
title('Bead-normalized Syn SSC: early cruises (performance-test transfer vs PT conversion) against FCB-bead cruises')
exportgraphics(f, fullfile(od, 'syn_ssc_early_cruises.png'), 'Resolution', 110)
end
