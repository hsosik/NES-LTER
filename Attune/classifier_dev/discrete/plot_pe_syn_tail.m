function plot_pe_syn_tail(outfile)
% function plot_pe_syn_tail(outfile)
% Hand gates only (no classifier): how the dim, low-chl Syn tail at the PE trigger
% was gated in PE-triggered preserved runs. Per run, the chl-vs-size trend of hand
% Syn clear of the trigger (GL2 > edge + 0.6) is fitted; events within 0.6 of the
% trigger and 0.5-1.5 below the trend, Syn-sized, with Syn-like PE:chl, are the tail.
%   red    hand Syn (core)
%   blue   hand Syn in the tail region (gated as Syn)
%   orange tail-region events the hand gate left unassigned
%   grey   other events
% Tuning cruises (top) vs EN617 / TN368 (bottom). The GL2 trigger edge is taken as the
% 0.5th percentile of GL2 (PE-triggered runs). Writes output/discrete/examples/pe_syn_tail.png
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 1, outfile = fullfile(od, 'examples', 'pe_syn_tail.png'); end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
runs = {'20220729_EN687' 'C10N18'; '20220729_EN687' 'C10N13'; '20230807_EN706' 'C13N18'; '20240906_EN720' 'C15N16'; ...
    '20180720_EN617' 'C13N18(1)'; '20180720_EN617' 'C21N20'; '20190705_TN368' '15N24(1)'; '20190705_TN368' 'C61N24'};
fig = figure('visible', 'off', 'position', [20 20 1500 1900]);
tl = tiledlayout(fig, size(runs, 1), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:size(runs, 1)
    cdir = runs{k, 1}; cr = extractAfter(cdir, '_');
    G = load(fullfile(root, cdir, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table; G.Properties.VariableNames = lower(G.Properties.VariableNames);
    gi = find(contains(G.fcslist, 'phyto_PE_SSC') & endsWith(G.fcslist, ['_' runs{k, 2} '.fcs']) & ~cellfun(@isempty, G.class), 1);
    fn = G.fcslist{gi};
    cfg = discrete_config(); cfg.cruise = cr; cfg.gl1_is_od2_scatter = true;
    [d, h] = fca_readfcs(fullfile(root, cdir, 'preserved', 'FCS', fn));
    ref = discrete_bead_reference(lab, datenum([h.date ' ' h.starttime]), cfg);
    F = attune_features(d, h, ref); F.PEchl = F.GL2 - F.BL3;
    % size as in assign_class_v2: SSC, merged with GL1 (OD2 scatter) for large particles
    F.size = F.SSC;
    if any(isfinite(F.GL1))
        hi = F.GL1 > log10(cfg.gl1_merge_bead(2)); mid = F.GL1 > log10(cfg.gl1_merge_bead(1)) & ~hi;
        F.size(hi) = F.GL1(hi); F.size(mid) = log10((10.^F.SSC(mid) + 10.^F.GL1(mid)) / 2);
    end
    hc = double(G.class{gi}(:)); late = d(:, 1) >= min(d(hc ~= 0, 1));
    F = F(late, :); hc = hc(late);
    e = prctile(F.GL2, 0.5);
    syn = hc == 2;
    core = syn & F.GL2 > e + 0.6;
    b = robustfit(F.size(core), F.BL3(core));
    sz = median(F.size(core)); pe = median(F.PEchl(core));
    r = F.BL3 - (b(1) + b(2) * F.size);
    region = F.GL2 <= e + 0.6 & r < -0.5 & r > -1.5 & abs(F.size - sz) < 0.5 & F.PEchl > pe - 0.5 & F.PEchl < pe + 0.9;
    tail_in = syn & region; tail_out = ~syn & region; core_s = syn & ~region; oth = ~syn & ~region;
    pairs = {'size' 'BL3'; 'BL3' 'GL2'; 'size' 'PEchl'};
    lims = struct('size', [-2.4 0.5], 'BL3', [-2 1.5], 'GL2', [-1.8 1], 'PEchl', [-1.5 1.2]);
    for p = 1:3
        ax = nexttile(tl); hold(ax, 'on')
        x = pairs{p, 1}; y = pairs{p, 2};
        plot(ax, F.(x)(oth), F.(y)(oth), '.', 'color', [.8 .8 .8], 'markersize', 2)
        plot(ax, F.(x)(core_s), F.(y)(core_s), '.', 'color', [.85 .1 .1], 'markersize', 2)
        plot(ax, F.(x)(tail_out), F.(y)(tail_out), '.', 'color', [1 .55 0], 'markersize', 3)
        plot(ax, F.(x)(tail_in), F.(y)(tail_in), '.', 'color', [0 .35 1], 'markersize', 3)
        if strcmp(y, 'GL2') || strcmp(x, 'GL2'), yline(ax, e, ':k'); end
        xlim(ax, lims.(x)); ylim(ax, lims.(y)); grid(ax, 'on'); box(ax, 'on')
        xlabel(ax, x); ylabel(ax, y); set(ax, 'fontsize', 8)
        if p == 1
            grp = 'tuning cruise'; if k > 4, grp = 'unseen cruise'; end
            title(ax, sprintf('%s %s (%s): tail gated as Syn %d, left out %d', cr, runs{k, 2}, grp, sum(tail_in), sum(tail_out)), 'fontsize', 9, 'interpreter', 'none')
        end
    end
end
title(tl, {'PE-triggered runs, hand gates only: red = hand Syn, blue = dim low-chl tail gated as Syn,' 'orange = same tail region left unassigned, grey = other; dotted line = PE (GL2) trigger edge'}, 'fontsize', 11)
if ~exist(fileparts(outfile), 'dir'), mkdir(fileparts(outfile)), end
exportgraphics(fig, outfile, 'resolution', 90); close(fig)
fprintf('%s\n', outfile);
end
