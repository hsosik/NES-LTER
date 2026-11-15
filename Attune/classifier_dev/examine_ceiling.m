function T = examine_ceiling(tags, nper, outdir)
% function T = examine_ceiling(tags, nper, outdir)
% Why do large events pile up at size ~ +1 even with the GL1 (OD2 side scatter)
% substitute for saturated SSC? For the nper lowest-agreement files of each
% cruise in tags (from score_agreement output), report how SSC and GL1 behave
% at the top of their ranges and how the merge offset is estimated.
% Writes outdir/ceiling_<file>.png and returns a summary table.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
T = table;
for t = tags
    a = load(fullfile(outdir, ['agreement_' t{1} '.mat']), 'A'); A = sortrows(a.A, 'agreement');
    for i = 1:min(nper, height(A))
        cpath = fullfile(root, char(A.cruise(i))); fn = char(A.file(i));
        bs = load_beadstat(fullfile(cpath, 'bead_calibrated'));
        [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
        ref = bead_reference(bs, A.time(i), cfg);
        [~, info] = assign_class_v2(d, h, ref, cfg, A.pro_measured(i));
        F = info.features;
        [~, finfo] = attune_features(d, h, ref); info.floor = finfo.floor;
        % GL1 on the 1 um bead scale (your batch_mergedSSC2volume convention: SSC and GL1
        % each relative to their own 1 um bead peak are taken as equivalent)
        cfg1 = cfg; cfg1.bead_pop.GL1 = 2; ref1 = bead_reference(bs, A.time(i), cfg1);
        GL1b = F.GL1 + (ref.bead_log.GL1 - ref1.bead_log.GL1) + ref.k.GL1 * log10(ref1.bead_hv.GL1 / ref.bead_hv.GL1);
        pn = {h.par.name};
        rawSSC = d(:, strcmp(pn, 'SSC-H')); rawGL1 = d(:, strcmp(pn, 'GL1-H'));
        rawSSCA = d(:, strcmp(pn, 'SSC-A'));
        sat = max(F.SSC); satv = ~(F.SSC < sat - 0.05);
        gsat = max(F.GL1); gsatv = ~(F.GL1 < gsat - 0.05);
        ov = ~satv & F.SSC > -0.54 & F.GL1 > info.floor.GL1 + 1;
        dlt = F.SSC(ov) - F.GL1(ov);
        % slope of SSC on GL1 in the overlap (1 if a constant offset is right)
        p = polyfit(F.GL1(ov), F.SSC(ov), 1);
        % offset in the upper part of the overlap only
        up = ov & F.SSC > sat - 0.8;
        dlt1 = F.SSC(ov) - GL1b(ov);
        r = table(string(A.cruise(i)), string(fn), A.agreement(i), height(F), sum(satv), max(rawSSC), sum(gsatv), max(rawGL1), ...
            sum(satv & gsatv), sum(ov), median(dlt), iqr(dlt), p(1), median(F.SSC(up) - F.GL1(up)), sum(up), ...
            median(F.GL1(satv) - info.floor.GL1), sum(satv & F.GL1 < info.floor.GL1 + 1), sat, prctile(F.size(satv), [5 50 95]), ...
            h.par(strcmp(pn, 'SSC-H')).hv, h.par(strcmp(pn, 'GL1-H')).hv, median(dlt1), iqr(dlt1), ...
            median(F.SSC(up) - GL1b(up)), prctile(GL1b(satv), [5 50 95]), ...
            'VariableNames', {'cruise' 'file' 'agreement' 'n' 'n_SSCsat' 'raw_SSC_max' 'n_GL1sat' 'raw_GL1_max' 'n_both_sat' 'n_overlap' ...
            'offset' 'offset_iqr' 'slope' 'offset_upper' 'n_upper' 'GL1_sat_above_floor' 'n_sat_GL1_near_floor' 'SSC_sat' 'size_sat_p5_50_95' 'hv_SSC' 'hv_GL1' 'offset_1um' 'offset_1um_iqr' 'offset_1um_upper' 'GL1_1um_sat_p5_50_95'});
        T = [T; r]; %#ok<AGROW>
        disp(r)

        f = figure('Visible', 'off', 'Position', [20 20 1500 450]);
        tiledlayout(1, 3, 'TileSpacing', 'compact');
        nexttile; hold on
        k = 1:max(1, ceil(height(F) / 60000)):height(F);
        scatter(F.GL1(k), F.SSC(k), 2, [.6 .6 .6], 'filled', 'MarkerFaceAlpha', .3)
        scatter(F.GL1(ov), F.SSC(ov), 2, [0 .4 1], 'filled', 'MarkerFaceAlpha', .3)
        scatter(F.GL1(satv), F.SSC(satv), 4, [.85 0 0], 'filled')
        plot(xlim, xlim + median(dlt), 'k-'); plot(xlim, polyval(p, xlim), 'k--')
        xline(info.floor.GL1 + 1, ':'); yline(sat - 0.05, ':')
        xlabel('GL1 (bead rel.)'); ylabel('SSC (bead rel.)'); grid on
        title(sprintf('blue = overlap (n %d), red = SSC saturated (n %d)', sum(ov), sum(satv)), 'FontSize', 8)
        legend({'all' 'overlap' 'SSC sat' 'constant offset' 'fit'}, 'Location', 'northwest', 'FontSize', 7)
        nexttile; hold on
        edges = -3:0.02:3;
        histogram(F.SSC, edges, 'DisplayStyle', 'stairs'); histogram(F.size, edges, 'DisplayStyle', 'stairs');
        histogram(F.GL1 + median(dlt), edges, 'DisplayStyle', 'stairs');
        histogram(GL1b, edges, 'DisplayStyle', 'stairs');
        set(gca, 'YScale', 'log'); xlabel('bead rel. (1 um)'); legend({'SSC' 'size (merged)' 'GL1 + fitted offset' 'GL1 on 1 um bead'}, 'FontSize', 7); grid on
        nexttile; hold on
        scatter(F.size(k), F.BL3(k), 2, [.6 .6 .6], 'filled', 'MarkerFaceAlpha', .3)
        scatter(F.size(satv), F.BL3(satv), 4, [.85 0 0], 'filled')
        xlabel('size'); ylabel('chl (BL3)'); grid on; title('red = SSC saturated', 'FontSize', 8)
        sgtitle(fn, 'Interpreter', 'none', 'FontSize', 9)
        exportgraphics(f, fullfile(outdir, ['ceiling_' regexprep(fn, '\.fcs$', '.png')]), 'Resolution', 80)
        close(f)
    end
end
end
