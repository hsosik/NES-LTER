function beadstat = generate_between_beadstats(makeplots)
% function beadstat = generate_between_beadstats(makeplots)
% Bead statistics for the FCB bead checks run between cruises (shore lab,
% discrete/preserved samples), in the same format as the cruise beadstat.mat
% files (measure_beads; read with load_beadstat).
%
% Runs: every .fcs under cruise_data\beads\FCB_bead_mix_experiment_settings
% (top level and between_cruises\) whose acquisition time is not already in
% a cruise's beadstat.mat. Excluded: the temp\ and between_cruises\skip\
% folders, size-calibration runs (size_cal_*) and duplicate copies (same
% acquisition time).
%
% Runs are processed in groups (measure_beads learns each group's bead
% offsets in its first pass): calendar year, split at the filter-layout change
% (19 Aug 2019) and at the new FCB bead mix (11 Dec 2024). OD2: none before
% Aug 2019 (AR34B/RB1904 settings), on GL1 after.
%
% Writes ...\FCB_bead_mix_experiment_settings\between_cruises\beadstat_2026\
%   beadstat.mat (beadstat, meta) and plots\ (one gating plot per run,
%   bead_timeseries.png). The older between_cruises\outputs\beadstat.mat
%   (2021 format, read by Process_Preserved_AllSteps*) is left in place.
% Extra columns: folder, group (string), bead_lot ("pre-Dec2024" or
% "new mix Dec2024").
%
% Nov 2026, bead stats redesign

if nargin < 1 || isempty(makeplots), makeplots = true; end
addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
warning('off', 'stats:robustcov:SingularCovariance')
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
base = fullfile(root, 'beads', 'FCB_bead_mix_experiment_settings');
outdir = fullfile(base, 'between_cruises', 'beadstat_2026');
newlot = datetime(2024, 12, 11);
layoutchange = datetime(2019, 8, 19);

% runs already in cruise beadstat files
D = dir([root '20*_*']); known = datetime.empty(0, 1);
for i = 1:numel(D)
    try
        B = load_beadstat(fullfile(root, D(i).name, 'bead_calibrated'), false);
        known = [known; B.time]; %#ok<AGROW>
    catch
    end
end

% candidate files
F = dir(fullfile(base, '**', '*.fcs'));
rel = strrep({F.folder}, base, '');
keep = ~contains(rel, {'\temp' '\skip' 'backup'}, 'IgnoreCase', true) & ~startsWith({F.name}, 'size_cal', 'IgnoreCase', true);
F = F(keep);
files = {}; folders = {}; dat = {}; hdr = {}; t = datetime.empty(0, 1);
for i = 1:numel(F)
    try
        [d, h] = fca_readfcs(fullfile(F(i).folder, F(i).name));
    catch
        continue
    end
    ti = datetime([h.date ' ' h.starttime]);
    if any(abs(known - ti) < minutes(1)) || any(abs(t - ti) < seconds(2)), continue, end
    files{end+1, 1} = F(i).name; folders{end+1, 1} = strrep(F(i).folder, base, ''); %#ok<AGROW>
    dat{end+1, 1} = single(d); hdr{end+1, 1} = h; t(end+1, 1) = ti; %#ok<AGROW>
end
fprintf('%d between-cruise bead runs\n', numel(files))

% groups
grp = strings(numel(t), 1);
for i = 1:numel(t)
    if t(i) < layoutchange, grp(i) = "2019 pre-Aug (v2019 layout)";
    elseif t(i) >= newlot, grp(i) = string(year(t(i))) + " new mix";
    else, grp(i) = string(year(t(i)));
    end
end
beadstat = table; metas = struct('group', {}, 'meta', {});
classes = cell(numel(t), 1); pops = cell(numel(t), 1);
for g = unique(grp, 'stable')'
    k = find(grp == g);
    if t(k(1)) < layoutchange
        conf = attune_config('', datenum(t(k(1)))); conf.OD2 = 'None';
    else
        conf = attune_config('', datenum(t(k(1)))); conf.OD2 = 'GL1';
    end
    [B, M, cl, pp] = measure_beads(dat(k), hdr(k), files(k), conf);
    B.folder = string(folders(k));
    B.group = repmat(g, height(B), 1);
    B.bead_lot = repmat(ifelse(t(k(1)) >= newlot, "new mix Dec2024", "pre-Dec2024"), height(B), 1);
    beadstat = [beadstat; B]; %#ok<AGROW>
    metas(end+1) = struct('group', g, 'meta', M); %#ok<AGROW>
    classes(k) = cl; pops(k) = pp;
    fprintf('%s: %d runs; QC good 0.5 um %d, 1 um %d, 6 um %d\n', g, height(B), sum(B.QC05), sum(B.QC1), sum(B.QC6))
end
[beadstat, o] = sortrows(beadstat, 'time');
meta.source = "between-cruise FCB bead checks";
meta.groups = metas;
meta.new_bead_mix = newlot;
meta.processed = datetime('now');
meta.code = "generate_between_beadstats (Nov 2026)";
meta.notes = ["Same format as cruise beadstat.mat (measure_beads); extra columns folder, group, bead_lot";
    "Runs already in a cruise beadstat.mat, temp\, between_cruises\skip\ and size_cal_* runs are excluded"];
if ~exist(outdir, 'dir'), mkdir(outdir), end
save(fullfile(outdir, 'beadstat.mat'), 'beadstat', 'meta')
fprintf('saved %s\n', fullfile(outdir, 'beadstat.mat'))

if makeplots
    pp = fullfile(outdir, 'plots'); if ~exist(pp, 'dir'), mkdir(pp), end
    fig = figure('visible', 'off', 'position', [50 50 1300 470]);
    for i = 1:numel(files)
        plot_bead_gates(double(dat{i}), hdr{i}, double(classes{i}), pops{i}, ...
            sprintf('%s  %s  (%s)', datestr(t(i), 'yyyy-mm-dd HH:MM'), files{i}, grp(i)), 0.25)
        exportgraphics(fig, fullfile(pp, sprintf('%s_%s.png', datestr(t(i), 'yyyymmdd_HHMM'), regexprep(files{i}, '\.fcs$', ''))), 'resolution', 60)
    end
    close(fig)
    timeseries_plot(beadstat, fullfile(pp, 'bead_timeseries.png'))
end
end

function timeseries_plot(B, fn)
pairs = {'05' 'GL2_H'; '1' 'SSC_H'; '1' 'FSC_H'; '1' 'BL1_H'; '1' 'BL3_H'; '1' 'GL1_H'; '6' 'GL3_H'; '05' 'BL2_H'};
fig = figure('visible', 'off', 'position', [50 50 1500 800]);
for i = 1:size(pairs, 1)
    subplot(2, 4, i)
    v = B.(['peak' pairs{i, 1}]).(pairs{i, 2}); q = B.(['QC' pairs{i, 1}]);
    semilogy(B.time(q), v(q), 'o', B.time(~q), v(~q), 'rx', 'markersize', 4)
    xline(datetime(2024, 12, 11), 'k--')
    title(sprintf('%s um, %s', strrep(pairs{i, 1}, '05', '0.5'), strrep(pairs{i, 2}, '_', '-')), 'fontsize', 8), grid on
end
sgtitle('Between-cruise FCB bead checks: peak (o QC good, x failed); dashed = new bead mix', 'fontsize', 9)
exportgraphics(fig, fn, 'resolution', 70)
close(fig)
end

function y = ifelse(c, a, b)
if c, y = a; else, y = b; end
end
