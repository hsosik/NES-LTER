function T = bead_ratio_history(outdir)
% function T = bead_ratio_history(outdir)
% Bead size ratios per FCB bead run across all cruises, to look for steps in
% scatter geometry (blocker bar / laser alignment) at the Thermo service dates.
% For each run: log10(1 um / 0.5 um) and log10(6 um / 1 um) peak ratios in
% FSC-H, SSC-H and (current layout only) GL1-H = 532 nm scatter through OD2.
% Both beads in a ratio come from the same run and hv, so gain drops out.
% Pre-Aug 2019 layout: GL1-H is PE fluorescence, so GL1 ratios are NaN there.
%
% outdir - figure/mat folder (default classifier_dev/output/pt_transfer)
% T - one row per run: cruise, time, hv SSC/FSC/GL1, n05/n1/n6, QC, ratios
%
% Nov 2026, scatter geometry check against service history

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
if nargin < 1, outdir = fullfile(dev, 'output', 'pt_transfer'); end
root = ['\' '\sosiknas1\lab_data\Attune\cruise_data\'];
d = dir(fullfile(root, '20*_*'));
d = d([d.isdir]);
T = table();
for i = 1:numel(d)
    try
        [b, m] = load_beadstat(fullfile(root, d(i).name));
    catch
        continue
    end
    if height(b) == 0, continue, end
    n = height(b);
    cruise = extractAfter(string(d(i).name), '_');
    src = "FCB";
    if isfield(m, 'source'), src = string(m.source); end
    newlay = b.time >= datetime(2019, 8, 19);
    r = table(repmat(cruise, n, 1), b.time, repmat(src, n, 1), b.hv.FSC, b.hv.SSC, b.hv.GL1, ...
        b.n05, b.n1, b.n6, b.QC05 & b.QC1, b.QC6 & b.QC1, newlay, ...
        'VariableNames', {'cruise' 'time' 'source' 'hvFSC' 'hvSSC' 'hvGL1' 'n05' 'n1' 'n6' 'QC_1_05' 'QC_6_1' 'newlayout'});
    for c = ["FSC_H" "SSC_H" "GL1_H"]
        r.("r1_05_" + c) = log10(b.peak1.(c) ./ b.peak05.(c));
        r.("r6_1_" + c) = log10(b.peak6.(c) ./ b.peak1.(c));
        if c == "GL1_H"
            r.("r1_05_" + c)(~newlay) = NaN; r.("r6_1_" + c)(~newlay) = NaN;
        end
    end
    T = [T; r]; %#ok<AGROW>
end
T = sortrows(T, 'time');
save(fullfile(outdir, 'bead_ratio_history.mat'), 'T')

% service events (Thermo reports; shared notes attune_optics_changes.md)
ev = {datetime(2018,7,15) 'blocker bar'; datetime(2019,8,6) 'blocker bar x2 + blue align';
    datetime(2019,10,23) 'PM'; datetime(2021,2,22) 'blocker bar'; datetime(2022,5,13) 'PM';
    datetime(2022,9,22) 'blue laser replaced'; datetime(2022,11,29) 'blue realign';
    datetime(2024,6,21) 'blue realign'; datetime(2024,7,15) 'PZT + flow cell (log)';
    datetime(2025,2,14) 'PZT power (log)'; datetime(2026,2,12) 'PM'};

f = figure('Position', [50 50 1300 900], 'Color', 'w', 'Visible', 'off');
pan = {"r1_05_SSC_H" "SSC 1/0.5 \mum" "QC_1_05"; "r1_05_FSC_H" "FSC 1/0.5 \mum" "QC_1_05";
    "r1_05_GL1_H" "GL1 (532 OD2) 1/0.5 \mum" "QC_1_05"; "r6_1_SSC_H" "SSC 6/1 \mum" "QC_6_1";
    "r6_1_FSC_H" "FSC 6/1 \mum" "QC_6_1"; "r6_1_GL1_H" "GL1 (532 OD2) 6/1 \mum" "QC_6_1"};
for p = 1:size(pan, 1)
    ax = subplot(3, 2, (mod(p-1,3))*2 + 1 + (p > 3)); hold(ax, 'on')
    y = T.(pan{p,1}); g = T.(pan{p,3});
    plot(ax, T.time(~g), y(~g), '.', 'Color', [.75 .75 .75])
    plot(ax, T.time(g), y(g), 'k.', 'MarkerSize', 9)
    % cruise medians of good runs
    [G, cr] = findgroups(T.cruise(g)); tt = T.time(g); yy = y(g);
    if ~isempty(G)
        mt = splitapply(@(x) median(x), datenum(tt), G); my = splitapply(@(x) median(x, 'omitnan'), yy, G);
        plot(ax, datetime(mt, 'ConvertFrom', 'datenum'), my, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
    end
    ok = y(g); ok = ok(isfinite(ok));
    if ~isempty(ok), ylim(ax, prctile(ok, [1 99]) + [-0.05 0.05]), end
    yl = ylim(ax);
    for e = 1:size(ev, 1)
        xline(ax, ev{e,1}, ':', 'Color', [0 0.4 0.8]);
    end
    xline(ax, datetime(2019,8,19), '-', 'Color', [0.8 0.4 0], 'LineWidth', 1);
    title(ax, pan{p,2}); ylabel(ax, 'log10 ratio'); grid(ax, 'on'); ylim(ax, yl)
end
sgtitle(f, {'FCB bead size ratios per run (black: QC good; red: cruise median)', ...
    'blue dotted: service events; orange: filter layout change 19 Aug 2019'})
exportgraphics(f, fullfile(outdir, 'bead_ratio_history.png'), 'Resolution', 110)
close(f)
end
