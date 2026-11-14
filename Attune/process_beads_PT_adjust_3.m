function beadstat = process_beads_PT_adjust_3(basepath, FCSfileinfo, beadfiles2include, beadtype, OD2setting, makeplots)
% function beadstat = process_beads_PT_adjust_3(basepath, FCSfileinfo, beadfiles2include, beadtype, OD2setting, makeplots)
% FCB bead statistics for one cruise, called by process_wrapper_2024b (step 4).
% Finds the 0.5, 1 and 6 um FCB bead populations in every bead-check file,
% measures their peak on every channel and QCs each population separately.
%
% basepath          - cruise folder (contains FCS\ and bead_calibrated\)
% FCSfileinfo       - table from FCS_DateTimeList; [] to list FCS\ directly
% beadfiles2include - cellstr; files whose names contain any of these are bead runs
% beadtype          - 'FCB' (PT beads: use process_beads_PT_adjust_2)
% OD2setting        - kept for the wrapper; the OD2 placement now comes from
%                     attune_config (configuration history), and a mismatch
%                     is reported
% makeplots         - optional, default true: gating plot per run + time series
%
% Writes bead_calibrated\beadstat.mat with
%   beadstat - one row per bead run (see measure_beads for the columns):
%              filename, time, hv, and per population P = 05, 1, 6:
%              nP, QCP, whyP, peakP, medP, rsdP
%   meta     - cruise configuration, bead signature with this cruise's learned
%              offsets, files considered, processing date
% and plots in bead_calibrated\beadstat_plots\.
% Read it with load_beadstat. The earlier beadstat_2026.mat / beadstat_2021.mat
% files are left in place.
%
% Method (classifier_dev: find_beads, bead_signature, bead_peak, attune_config,
% measure_beads):
%   - each population is searched for where the bead signature puts it at
%     the file's voltages, on two channels where it is bright (0.5 um SSC x
%     GL2, 1 um SSC x BL1, 6 um GL1 x GL3, with fallbacks when a channel is
%     saturated, too dim or not installed); a first pass learns the cruise's
%     own offsets, the second gates around them
%   - densest point, then a robust (MCD) ellipse; an upper limit on one more
%     channel drops bead aggregates (1 um with 0.5 um on GL3, 6 um with 1 um
%     on BL1)
%   - peak = KDE mode of log10 signal, so skewed tails and minor modes do not
%     move it (replaces centend2; medians kept as a check)
%   - QC per population: >= 100 events, not too broad, near its expected
%     position on the gating channels and on its other bright channels
%
% March 2026, Heidi: started from process_beads_PT_adjust_2
% Nov 2026: rewritten (bead stats redesign); see classifier_dev/BEADSTATS_HANDOFF.md

if nargin < 6 || isempty(makeplots), makeplots = true; end
if ~strcmp(beadtype, 'FCB')
    error('process_beads_PT_adjust_3 handles FCB beads only; use process_beads_PT_adjust_2 for PT beads')
end
addpath(fullfile(fileparts(mfilename('fullpath')), 'classifier_dev'))

[~, cdir] = fileparts(strip(basepath, 'right', filesep));
cruise = extractAfter(cdir, '_');
conf = attune_config(cruise);
if nargin > 4 && ~isempty(OD2setting) && ~strcmpi(OD2setting, conf.OD2)
    warning('OD2setting %s from the wrapper differs from the configuration history (%s); using %s', ...
        OD2setting, conf.OD2, conf.OD2)
end
outpath = fullfile(basepath, 'bead_calibrated');
fpath = fullfile(basepath, 'FCS');

% bead files
if isempty(FCSfileinfo)
    F = dir(fullfile(fpath, '*.fcs')); fcslist = {F.name}';
else
    fcslist = FCSfileinfo.fcslist;
end
files = fcslist(contains(fcslist, beadfiles2include));
if isempty(files) || ~conf.fcb_beads
    warning('%s: no FCB bead runs', cruise)
    beadstat = table;
    return
end

% read once
dat = cell(numel(files), 1); hdr = cell(numel(files), 1);
for f = 1:numel(files)
    try
        [d, h] = fca_readfcs(fullfile(fpath, files{f}));
        dat{f} = single(d); hdr{f} = h;
    catch err
        warning('could not read %s: %s', files{f}, err.message)
    end
end

[beadstat, meta, classes, popsall] = measure_beads(dat, hdr, files, conf);
beadstat = sortrows(beadstat, 'time');
meta.cruise = string(cruise);
meta.beadfiles2include = beadfiles2include;
meta.processed = datetime('now');
meta.code = "process_beads_PT_adjust_3 (Nov 2026 redesign)";
save(fullfile(outpath, 'beadstat.mat'), 'beadstat', 'meta')
fprintf('%s: %d bead runs; QC good 0.5 um %d, 1 um %d, 6 um %d\nsaved %s\n', cruise, height(beadstat), ...
    sum(beadstat.QC05), sum(beadstat.QC1), sum(beadstat.QC6), fullfile(outpath, 'beadstat.mat'))

if makeplots
    plotpath = fullfile(outpath, 'beadstat_plots');
    if ~exist(plotpath, 'dir'), mkdir(plotpath), end
    fig = figure('visible', 'off', 'position', [50 50 1300 470]);
    for f = 1:numel(files)
        if isempty(hdr{f}), continue, end
        plot_bead_gates(double(dat{f}), hdr{f}, double(classes{f}), popsall{f}, files{f}, meta.signature.window)
        exportgraphics(fig, fullfile(plotpath, regexprep(files{f}, '\.fcs$', '.png')), 'resolution', 70)
    end
    plot_timeseries(beadstat, cruise, plotpath)
    close(fig)
end
end

function plot_timeseries(B, cruise, plotpath)
% peak of each population on the channels the classifier normalizes by
pairs = {'05' 'BL2_H'; '05' 'GL2_H'; '05' 'GL4_H'; '1' 'FSC_H'; '1' 'SSC_H'; '1' 'GL1_H'; ...
    '1' 'BL1_H'; '1' 'BL3_H'; '6' 'GL3_H'; '1' 'SSC_A'; '1' 'GL1_A'};
fig = figure('visible', 'off', 'position', [50 50 1400 800]);
for i = 1:size(pairs, 1)
    subplot(3, 4, i)
    v = B.(['peak' pairs{i, 1}]).(pairs{i, 2}); q = B.(['QC' pairs{i, 1}]);
    semilogy(B.time(q), v(q), 'o', B.time(~q), v(~q), 'rx')
    title(sprintf('%s um, %s', strrep(pairs{i, 1}, '05', '0.5'), strrep(pairs{i, 2}, '_', '-')), 'fontsize', 8)
    grid on
end
subplot(3, 4, 12), axis off
text(0, 0.6, {cruise, 'o QC good', 'x QC failed', 'peak (KDE mode), log scale'}, 'fontsize', 9)
exportgraphics(fig, fullfile(plotpath, 'bead_timeseries.png'), 'resolution', 70)
close(fig)
end
