function show_pro_tail_examples(list, outfile)
% function show_pro_tail_examples(list, outfile)
% Events that become Euk when the Pro-tail cut is switched off
% (cfg.sig.edge_pro_tail = false), one row per file:
%   1 size vs chl    v2 Euks (cut on) green, events added without the cut red,
%                    other events grey; chl trigger edge (dashed), edge + band
%                    (dotted), Syn size + edge_pro_size_above_syn (dash-dot)
%   2 chl vs PE      same colours (added events have no PE)
%   3 chl histogram  small no-PE events (size < Syn + 0.1) above the edge:
%                    Euks kept (green), added without the cut (red), old-class
%                    Euks among the added (black outline); Euk chl gap if found
% list - n x 2 cell: cruise folder, FCS file name
% Nov 2026, Pro-tail assessment
addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
nf = size(list, 1);
f = figure('Visible', 'off', 'Position', [50 50 1500 380 * nf]);
tl = tiledlayout(nf, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nf
    cruise = extractAfter(list{i, 1}, '_'); bp = fullfile(root, list{i, 1}, 'bead_calibrated');
    bs = load_beadstat(bp); load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
    r = strcmp(FCSfileinfo.fcslist, list{i, 2});
    cfg = classifier_config(); cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
    [d, h] = fca_readfcs(fullfile(root, list{i, 1}, 'FCS', list{i, 2}));
    ref = bead_reference(bs, FCSfileinfo.matdate_start(r), cfg);
    [a, ia] = assign_class_v2(d, h, ref, cfg, 0);
    c2 = cfg; c2.sig.edge_pro_tail = false;
    b = assign_class_v2(d, h, ref, c2, 0);
    o = load(fullfile(bp, 'class', regexprep(list{i, 2}, '\.fcs$', '.mat'))).class(:, 1);
    F = ia.features; e = ia.trigger_edges.BL3; ss = ia.syn_ref.size; sg = cfg.sig;
    add = a ~= 1 & b == 1; keep = a == 1;
    oth = ~add & ~keep & F.SSC > ia.floor.SSC + 0.1 & F.BL3 > e - 0.6;
    % 1 size vs chl
    nexttile; hold on
    scatter(F.size(oth), F.BL3(oth), 2, [.75 .75 .75], 'filled')
    scatter(F.size(keep), F.BL3(keep), 3, [0 .6 0], 'filled')
    scatter(F.size(add), F.BL3(add), 8, [.9 0 0], 'filled')
    yline(e, 'k--'); yline(e + sg.edge_pro_band, 'k:'); xline(ss + sg.edge_pro_size_above_syn, 'k-.');
    xlim([-2.5 1.5]); ylim([e - 0.6 e + 2.5]); grid on; box on
    xlabel('size (log10 rel. 1 um bead)'); ylabel('chl (BL3)')
    title(sprintf('%s  %s', cruise, strrep(extractBefore(list{i, 2}, '_Group'), '_', '\_')), 'FontWeight', 'normal')
    text(0.02, 0.97, sprintf('Euk %d (cut on) -> %d (cut off)\nadded %d: old Euk %d, old none %d, old other %d', sum(a == 1), sum(b == 1), sum(add), ...
        sum(add & o == 1), sum(add & o == 0), sum(add & o ~= 0 & o ~= 1)), 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8)
    % 2 chl vs PE
    nexttile; hold on
    scatter(F.BL3(oth), F.GL2(oth), 2, [.75 .75 .75], 'filled')
    scatter(F.BL3(keep), F.GL2(keep), 3, [0 .6 0], 'filled')
    scatter(F.BL3(add), F.GL2(add), 8, [.9 0 0], 'filled')
    xline(e, 'k--'); xline(e + sg.edge_pro_band, 'k:'); yline(ia.rules.pe_min, 'k--');
    xlim([e - 0.6 e + 2.5]); grid on; box on
    xlabel('chl (BL3)'); ylabel('PE (GL2)')
    % 3 chl histogram of small no-PE events
    nexttile; hold on
    sm = F.size < ss + sg.edge_pro_size_above_syn & F.BL3 > e - 0.05 & F.BL3 < e + 1.5;
    ed = e - 0.05:0.025:e + 1.5;
    hk = histcounts(F.BL3(sm & keep), ed); ha = histcounts(F.BL3(sm & add), ed); ho = histcounts(F.BL3(sm & ~keep & ~add), ed);
    bar(ed(1:end - 1) + 0.0125, [hk; ha; ho]', 1, 'stacked', 'EdgeColor', 'none');
    colororder(gca, [0 .6 0; .9 0 0; .75 .75 .75])
    stairs(ed, [histcounts(F.BL3(sm & add & o == 1), ed) 0], 'k', 'LineWidth', 1)
    xline(e, 'k--'); xline(e + sg.edge_pro_band, 'k:');
    if isfinite(ia.euk_chl_gap), xline(ia.euk_chl_gap, 'b-', 'Euk gap'); end
    grid on; box on; xlabel('chl (BL3), size < Syn + 0.1'); ylabel('events')
    legend({'Euk (cut on)' 'added (cut off)' 'other' 'added, old Euk'}, 'Location', 'northeast', 'FontSize', 7)
end
title(tl, 'Events that become Euk without the Pro-tail cut (red); dashed = chl trigger edge, dotted = edge + 0.4, dash-dot = Syn size + 0.1')
exportgraphics(f, outfile, 'Resolution', 110); close(f)
end
