function P = check_pro(cruisedirs, nPerCruise, outdir)
% function P = check_pro(cruisedirs, nPerCruise, outdir)
% Is Pro separated robustly from noise and non-phytoplankton particles?
% For Pro-setting files (pro_measured = 1):
%  - Pro concentration (existing and v2) vs water temperature; Pro should
%    be essentially absent in cold shelf water, so counts there are false
%    positives
%  - how close the v2 Pro cluster sits to the BL3 trigger edge (truncation)
%  - pigment / pulse-shape signature of v2 Pro vs unassigned events of
%    similar size, using channels the classifiers don't use
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
P = table; E = table;
for c = 1:numel(cruisedirs)
    cpath = fullfile(root, cruisedirs{c});
    bp = fullfile(cpath, 'bead_calibrated');
    bs = load_beadstat(bp);
    load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
    Fi = FCSfileinfo(FCSfileinfo.pro_measured == 1 & FCSfileinfo.QC_flag == 1, :);
    Fi = Fi(cellfun(@(x) exist(fullfile(bp, 'class', regexprep(x, '\.fcs$', '.mat')), 'file') == 2, Fi.fcslist), :);
    Fi = sortrows(Fi, 'matdate_start');
    u = load(fullfile(bp, 'AttuneTable_uw_match.mat'), 'Attune_uw_match'); u = u.Attune_uw_match;
    v = u.Properties.VariableNames;
    tn = v(contains(v, {'sbe48t' 'tsg1_temperature' 'tsg1_sst' 'sbe45t' 'sst' 'temperature'}, 'IgnoreCase', true) & ~contains(v, 'air', 'IgnoreCase', true));
    temp = nan(height(u), 1);
    for t = tn
        x = u.(t{1}); if isnumeric(x) && mean(isfinite(x)) > .5, temp = double(x); break, end
    end
    pick = unique(round(linspace(1, height(Fi), min(nPerCruise, height(Fi)))));
    for k = pick
        fn = Fi.fcslist{k};
        [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
        ref = bead_reference(bs, Fi.matdate_start(k), cfg);
        [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, 1);
        o = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class'); old = o.class(:, 1);
        F = info.features;
        edgeBL3 = NaN; if isfield(info.trigger_edges, 'BL3'), edgeBL3 = info.trigger_edges.BL3; end
        ui = find(strcmp(cellstr(u.Filename), fn), 1);
        tt = NaN; if ~isempty(ui), tt = temp(ui); end
        pro = new == 7;
        r = table(string(cruisedirs{c}), string(fn), Fi.matdate_start(k), tt, Fi.vol_analyzed(k), sum(old == 7), sum(pro), ...
            median(F.BL3(pro)) - edgeBL3, mean(F.BL3(pro) < edgeBL3 + 0.2), string(fcshdr.tr1_par), fcshdr.tr1_level, ...
            'VariableNames', {'cruise' 'file' 'time' 'temp' 'vol' 'old_pro' 'new_pro' 'pro_BL3_above_edge' 'frac_pro_near_edge' 'trigger' 'trigger_level'});
        P = [P; r]; %#ok<AGROW>
        % events of Pro size: v2 Pro vs unassigned, for signature comparison
        small = F.size < -1.24 & F.size > -2.44;
        pn = {fcshdr.par.name};
        ah = @(ch) log10(max(fcsdat(:, strcmp(pn, [ch '-A'])), 1) ./ max(fcsdat(:, strcmp(pn, [ch '-H'])), 1));
        T = table(F.size, F.BL3, F.BL1, F.BL2, F.GL4, F.GL2, ah('BL3'), ah('SSC'), log10(max(fcsdat(:, strcmp(pn, 'FSC-W')), 1)), ...
            'VariableNames', {'size' 'BL3' 'BL1' 'BL2' 'GL4' 'GL2' 'BL3_AH' 'SSC_AH' 'FSC_W'});
        T.edge = repmat(edgeBL3, height(T), 1);
        T.group = repmat("other", height(T), 1);
        T.group(pro) = "v2 pro"; T.group(small & new == 0) = "unassigned, Pro-size";
        T.cold = repmat(tt < 18, height(T), 1);
        keep = find(T.group ~= "other"); keep = keep(1:max(1, ceil(numel(keep) / 2000)):end);
        E = [E; T(keep, :)]; %#ok<AGROW>
        fprintf('%s %s T=%.1f old %d new %d, Pro median BL3 %.2f above edge, %.0f%% within 0.2\n', cruisedirs{c}, fn, tt, sum(old == 7), sum(pro), r.pro_BL3_above_edge, 100 * r.frac_pro_near_edge);
    end
end
save(fullfile(outdir, 'check_pro.mat'), 'P', 'E')

f = figure('Visible', 'off', 'Position', [40 40 1700 900]);
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
nexttile; hold on
cr = unique(P.cruise); mk = 'os^d';
for i = 1:numel(cr)
    q = P(P.cruise == cr(i), :);
    scatter(q.temp, (q.new_pro + 1) ./ q.vol, 30, 'r', mk(i), 'filled')
    scatter(q.temp, (q.old_pro + 1) ./ q.vol, 30, 'k', mk(i))
end
set(gca, 'YScale', 'log'); xlabel('temperature (°C)'); ylabel('Pro events per ml (+1)'); grid on
title('Pro vs temperature (red v2, black existing)')
nexttile
scatter(P.temp, P.frac_pro_near_edge, 30, 'filled'); xlabel('temperature (°C)'); ylabel('fraction of v2 Pro within 0.2 of BL3 trigger edge'); grid on
pairs = {'size' 'BL3'; 'BL3' 'BL2'; 'BL3' 'GL4'; 'size' 'FSC_W'; 'BL3' 'BL3_AH'; 'BL3' 'BL1'};
for p = 1:numel(pairs) / 2
    nexttile; hold on
    g = E.group == "unassigned, Pro-size";
    scatter(E.(pairs{p, 1})(g), E.(pairs{p, 2})(g), 3, [.6 .6 .6], 'filled', 'MarkerFaceAlpha', .3)
    g = E.group == "v2 pro" & ~E.cold;
    scatter(E.(pairs{p, 1})(g), E.(pairs{p, 2})(g), 3, [0 .3 1], 'filled', 'MarkerFaceAlpha', .3)
    g = E.group == "v2 pro" & E.cold;
    scatter(E.(pairs{p, 1})(g), E.(pairs{p, 2})(g), 5, [.9 0 0], 'filled', 'MarkerFaceAlpha', .5)
    xlabel(pairs{p, 1}, 'Interpreter', 'none'); ylabel(pairs{p, 2}, 'Interpreter', 'none'); grid on
    if p == 1, legend({'unassigned, Pro-size' 'v2 Pro, warm' 'v2 Pro, <18 °C'}, 'Location', 'northwest', 'FontSize', 7); end
end
title(tl, 'Pro in Pro-setting files: temperature check and signature vs unassigned events of Pro size')
exportgraphics(f, fullfile(outdir, 'check_pro.png'), 'Resolution', 90)
end
