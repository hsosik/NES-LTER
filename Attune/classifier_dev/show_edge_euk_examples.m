function show_edge_euk_examples(list, outfile)
% function show_edge_euk_examples(list, outfile)
% Euks without PE at the chl trigger edge, now that the Pro-tail cut is gone
% (Nov 2026): Euks with no PE, chl within edge_pro_band of the chl edge and no
% larger than Syn + 0.1 ("edge Euks"). One row per file:
%   1 size vs chl    other Euks green, edge Euks red (old class none) or black
%                    (old Euk), other events grey; chl edge dashed, edge + band
%                    dotted, Syn size + 0.1 dash-dot
%   2 chl vs PE      same colours
%   3 chl histogram  small no-PE events (size < Syn + 0.1) from the edge up: other
%                    Euks, edge Euks, unassigned; the Euk chl gap where one was found
% list - n x 2 cell: cruise folder, FCS file name
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
    o = load(fullfile(bp, 'class', regexprep(list{i, 2}, '\.fcs$', '.mat'))).class(:, 1);
    F = ia.features; e = ia.trigger_edges.BL3; ss = ia.syn_ref.size; band = cfg.sig.edge_pro_band;
    hp = F.GL2 > ia.rules.pe_min & F.PEchl > cfg.sig.pe_min_PEchl;
    ee = a == 1 & ~hp & F.BL3 < e + band & F.size < ss + 0.1;
    eu = a == 1 & ~ee;
    oth = a ~= 1 & F.SSC > ia.floor.SSC + 0.1 & F.BL3 > e - 0.6;
    eo = ee & o == 1; en = ee & o ~= 1;
    % 1 size vs chl
    nexttile; hold on
    scatter(F.size(oth), F.BL3(oth), 2, [.75 .75 .75], 'filled')
    scatter(F.size(eu), F.BL3(eu), 3, [0 .6 0], 'filled')
    scatter(F.size(en), F.BL3(en), 9, [.9 0 0], 'filled')
    scatter(F.size(eo), F.BL3(eo), 9, [0 0 0], 'filled')
    yline(e, 'k--'); yline(e + band, 'k:'); xline(ss + 0.1, 'k-.');
    xlim([-2.5 1.5]); ylim([e - 0.6 e + 2.5]); grid on; box on
    xlabel('size (log10 rel. 1 um bead)'); ylabel('chl (BL3)')
    title(sprintf('%s  %s', cruise, strrep(extractBefore(list{i, 2}, '_Group'), '_', '\_')), 'FontWeight', 'normal')
    gtxt = 'no Euk chl gap found'; if isfinite(ia.euk_chl_gap), gtxt = sprintf('Euk chl gap at %.2f (edge %.2f)', ia.euk_chl_gap, e); end
    text(0.02, 0.97, sprintf('Euk %d; edge Euks %d: old none %d, old Euk %d, old other %d\n%s', sum(a == 1), sum(ee), sum(ee & o == 0), sum(eo), sum(ee & o ~= 0 & o ~= 1), gtxt), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8)
    % 2 chl vs PE
    nexttile; hold on
    scatter(F.BL3(oth), F.GL2(oth), 2, [.75 .75 .75], 'filled')
    scatter(F.BL3(eu), F.GL2(eu), 3, [0 .6 0], 'filled')
    scatter(F.BL3(en), F.GL2(en), 9, [.9 0 0], 'filled')
    scatter(F.BL3(eo), F.GL2(eo), 9, [0 0 0], 'filled')
    xline(e, 'k--'); xline(e + band, 'k:'); yline(ia.rules.pe_min, 'k--');
    xlim([e - 0.6 e + 2.5]); grid on; box on
    xlabel('chl (BL3)'); ylabel('PE (GL2)')
    % 3 chl histogram, small no-PE events
    nexttile; hold on
    sm = ~hp & F.size < ss + 0.1 & F.BL3 > e - 0.05 & F.BL3 < e + 1.5 & F.SSC > ia.floor.SSC + 0.1;
    ed = e - 0.05:0.025:e + 1.5;
    hb = [histcounts(F.BL3(sm & eu), ed); histcounts(F.BL3(sm & ee), ed); histcounts(F.BL3(sm & a ~= 1), ed)];
    b = bar(ed(1:end - 1) + 0.0125, hb', 1, 'stacked', 'EdgeColor', 'none');
    b(1).FaceColor = [0 .6 0]; b(2).FaceColor = [.9 0 0]; b(3).FaceColor = [.75 .75 .75];
    xline(e, 'k--'); xline(e + band, 'k:');
    if isfinite(ia.euk_chl_gap), xline(ia.euk_chl_gap, 'b-', 'Euk gap'); end
    grid on; box on; xlabel('chl (BL3), no PE, size < Syn + 0.1'); ylabel('events')
    legend(b, {'other Euks' 'edge Euks' 'not Euk'}, 'Location', 'northeast', 'FontSize', 7)
end
title(tl, 'Euks with no PE at the chl trigger edge (red: old classes left them out; black: old Euk); dashed = chl trigger edge, dotted = edge + 0.4, dash-dot = Syn size + 0.1')
exportgraphics(f, outfile, 'Resolution', 110); close(f)
end
