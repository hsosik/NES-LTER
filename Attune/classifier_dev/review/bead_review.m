function D = bead_review(cruise, opts)
% function D = bead_review(cruise, opts)
% Cruise-level bead review: write a static HTML viewer (index.html + bead_review_data.js
% + frames\<run>.png) for the cruise's FCB bead runs, from bead_calibrated\beadstat.mat
% (process_beads_PT_adjust_3). Open index.html in a browser (works from the share).
%
% The page shows
%   - each classifier channel's bead peak through the cruise, brought to a common
%     voltage, with the reference every sample file is normalized against
%     (bead_reference: median of the runs within cfg.bead_window_days) and the
%     sample times, so a bad run and the files it affects are visible together
%   - per run: the gating plot (plot_bead_gates), log10 histograms of the 0.5, 1 and
%     6 um populations on every channel with their peak (KDE mode) and median, and
%     the QC and peak table
%   - unusual runs: failed population QC, peak off the cruise median after the
%     voltage correction, broad or skewed populations, voltage changes
%   - review decisions per run (status and note), saved by the page to
%     bead_review_decisions.csv beside beadstat.mat (Chrome or Edge), as in the
%     classification review. Runs marked exclude are skipped by classification and
%     the volume step (load_beadstat -> bead_reference, batch_mergedSSC2volume) the
%     next time they run; rebuild this page to see the references without them
%
% Cruises without FCB bead runs (2018 - early 2019): beadstat.mat is the PT-transfer
% table (make_perf_beadstat, meta.source "PT performance test transfer"), one row per
% Attune performance test with FCB-equivalent peaks (PT peak + FCB-PT offset for the
% channel, bead size and optical state). The page then shows each test instead of a
% bead run: FSC x SSC with the singlet gate and the PT intensity levels (plot_pt_test),
% histograms of all events, singlets and the brightest level, the PT peak, offset and
% CV per channel, the optical state and filter set; flags for broad peaks, a state or
% filter set other than the cruise's usual, and tests before or after the sample files.
% Exclusions work the same way. Tests are read from performance_test_fcs beside
% cruise_data (opts.ptdir).
%
% Preserved (discrete) samples: beaddir = <cruise>\preservedead_calibrated, whose
% beadstat.mat is the cruise's subset of the lab bead checks (preserved_bead_step,
% meta.source "lab FCB bead checks"). Runs are read from the lab bead folder
% (meta.lab_dir), samples are the preserved runs in FCSfileinfo (runtype CHL, PE,
% pro, hbac), and each sample's reference is discrete_bead_reference (same-day runs;
% a channel without a good same-day population uses the nearest good run: orange).
% The page also lists each run day: samples, same-day runs, channels needing a
% fallback, and days without a lab bead run. Exclusions apply to this cruise.
%
% opts (all optional):
%   .root       cruise_data root passed to attune_cruise_settings
%   .beaddir    folder with beadstat.mat and FCSfileinfo.mat (default bead_calibrated);
%               without FCSfileinfo.mat the sample times are read from the FCS headers
%   .outdir     viewer folder (default <beaddir>\bead_review)
%   .ptdir      performance-test FCS folder (PT-transfer tables; default
%               <cruise_data>\..\performance_test_fcs)
% e.g. bead_review('EN727'), bead_review('AR29')
%
% Nov 2026, bead review; PT-transfer tables

if nargin < 2, opts = struct(); end
def = struct('root', '', 'beaddir', '', 'outdir', '', 'ptdir', '');
for f = fieldnames(def)', if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end, end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fullfile(fileparts(here), 'pt_transfer')) % classifier_dev, pt_test_peaks
s = attune_cruise_settings(cruise, opts.root);
bp = opts.beaddir; if isempty(bp), bp = fullfile(s.basepath, 'bead_calibrated'); end
outdir = opts.outdir; if isempty(outdir), outdir = fullfile(bp, 'bead_review'); end
decfile = fullfile(bp, 'bead_review_decisions.csv'); % where load_beadstat reads exclusions
if ~exist(fullfile(outdir, 'frames'), 'dir'), mkdir(fullfile(outdir, 'frames')); end

[bs, meta] = load_beadstat(bp, false); % as measured: the page shows QC before exclusions
src = ""; if isfield(meta, 'source'), src = string(meta.source); end
pt = src == "PT performance test transfer";
lab = src == "lab FCB bead checks"; % preserved samples: lab bead runs
if lab, addpath(fullfile(fileparts(here), 'discrete')), end
if ~pt && ~isfield(meta, 'signature')
    error('bead_review:source', ['%s: beadstat.mat is neither an FCB bead table nor a PT-transfer table ' ...
        '(meta.source "%s"); for a cruise without FCB beads rebuild it with make_perf_beadstat'], s.cruise, src)
end
if isfile(fullfile(bp, 'FCSfileinfo.mat'))
    F = load(fullfile(bp, 'FCSfileinfo.mat')); fn = fieldnames(F); Fi = F.(fn{1});
elseif lab
    error('bead_review:files', 'preserved samples: run the files step first (no FCSfileinfo.mat in %s)', bp)
else % step 1 not run yet: sample times from the FCS headers
    Fi = list_fcs(fullfile(s.basepath, 'FCS'));
end
cfg = classifier_config();
cfg.cruise = s.cruise;
conf = attune_config(s.cruise);
fcsdir = fullfile(s.basepath, 'FCS'); % where the bead runs are
if lab
    cfg.cruise = 'LAB'; conf = attune_config('LAB'); fcsdir = char(meta.lab_dir);
end
P = {'05' '1' '6'}; PN = {'0.5' '1' '6'};
chans = bs.peak1.Properties.VariableNames; % FSC_H ... GL4_H, SSC_A, GL1_A
bt = datenum(bs.time);
nr = height(bs);
if pt
    ptdir = opts.ptdir;
    if isempty(ptdir), ptdir = fullfile(fileparts(fileparts(regexprep(s.basepath, '[\\/]+$', ''))), 'performance_test_fcs'); end
end

% common voltage per channel: the most frequent among the runs (performance
% tests set their own voltages: the median)
hvc = struct();
for c = cfg.chans
    h = bs.hv.(c{1}); h = h(isfinite(h) & h > 0);
    if isempty(h), hvc.(c{1}) = NaN; elseif pt, hvc.(c{1}) = round(median(h)); else, hvc.(c{1}) = mode(h); end
end

% runs: re-gate with the cruise's signature for the plots and histograms
edges = 0:0.05:6.1;
R = cell(nr, 1);
if pt
    for i = 1:nr, R{i} = pt_run(bs, i, chans, ptdir, outdir, edges); end
    nfcb = 0;
else
    sig = meta.signature;
    nfcb = nr;
end
for i = 1:nfcb
    r = struct('file', char(bs.filename(i)), 'time', char(datetime(bs.time(i), 'Format', 'yyyy-MM-dd''T''HH:mm:ss')));
    [~, nm] = fileparts(r.file);
    r.frame = ['frames/' nm '.png'];
    r.hv = table2struct(bs.hv(i, :));
    r.pops = cell(1, 3);
    for ip = 1:3
        b = P{ip};
        q = struct('name', PN{ip}, 'n', bs.(['n' b])(i), 'QC', logical(bs.(['QC' b])(i)), 'why', char(bs.(['why' b])(i)));
        q.peak = cell2struct(num2cell(log10(table2array(bs.(['peak' b])(i, :)))), chans, 2);
        q.med = cell2struct(num2cell(log10(table2array(bs.(['med' b])(i, :)))), chans, 2);
        q.rsd = table2struct(bs.(['rsd' b])(i, :));
        q.hist = struct(); q.nnow = NaN;
        r.pops{ip} = q;
    end
    r.all = struct(); r.read = '';
    try
        [d, h] = fca_readfcs(fullfile(fcsdir, r.file));
        d = double(d); names = {h.par.name};
        rng(0) % robustcov subsets
        [cl, pops] = find_beads(d, h, sig, conf);
        fig = figure('visible', 'off', 'position', [50 50 1300 470]);
        plot_bead_gates(d, h, cl, pops, r.file, sig.window)
        exportgraphics(fig, fullfile(outdir, r.frame), 'resolution', 70)
        close(fig)
        for k = 1:numel(chans)
            j = strcmp(names, file_param(chans{k}, conf));
            if ~any(j), continue, end
            x = log10(max(d(:, j), 1));
            r.all.(chans{k}) = histcounts(x, edges);
            for ip = 1:3
                r.pops{ip}.hist.(chans{k}) = histcounts(x(cl == ip), edges);
            end
        end
        for ip = 1:3, r.pops{ip}.nnow = sum(cl == ip); end
    catch ME
        r.read = ME.message;
        warning('bead_review:read', '%s: %s', r.file, ME.message)
    end
    R{i} = r;
end

% sample files and the reference each one is normalized against (bead_reference)
isbead = contains(Fi.fcslist, {'FCB_bead' 'fcb_bead'});
if lab % preserved samples, as classify selects them
    isbead = strcmp(Fi.runtype, 'other') | contains(Fi.fcslist, s.exclude_preserved) | ~isfinite(Fi.matdate_start);
end
t = Fi.matdate_start(~isbead); names = Fi.fcslist(~isbead);
[t, o] = sort(t); names = names(o);
ns = numel(t);
S = struct('file', {names}, 'time', {cellstr(datetime(t, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'))});
S.ref = struct(); S.nruns = struct(); S.runs = struct(); S.nearest = struct();
% exactly what classification uses: bead_reference on the table with the bead
% review's exclusions applied (as decided when this page was built)
bsx = load_beadstat(bp);
warning('off', 'bead_reference:noQC')
for c = cfg.chans
    S.ref.(c{1}) = nan(ns, 1); S.nruns.(c{1}) = zeros(ns, 1); S.nearest.(c{1}) = false(ns, 1); S.runs.(c{1}) = cell(ns, 1);
end
for k = 1:ns
    if lab
        ref = discrete_bead_reference(bsx, t(k), cfg);
        if isempty(ref) % no lab bead run that day: classify skips the file
            for c = cfg.chans, S.runs.(c{1}){k} = []; end
            continue
        end
    else
        ref = bead_reference(bsx, t(k), cfg);
    end
    for c = cfg.chans
        ch = c{1};
        u = find(ismember(bsx.filename, ref.runs.(ch)));
        S.runs.(ch){k} = u(:)' - 1; % 0-based for the page
        S.nearest.(ch)(k) = ~isempty(u) && all(abs(bt(u) - t(k)) > cfg.bead_window_days);
        if lab, S.nearest.(ch)(k) = ref.fallback.(ch) > 0; end % nearest good run on another day
        % reference at the cruise's common voltage
        S.ref.(ch)(k) = ref.bead_log.(ch) + cfg.k.(ch) * log10(hvc.(ch) / ref.bead_hv.(ch));
        S.nruns.(ch)(k) = ref.ntimes.(ch);
    end
end
warning('on', 'bead_reference:noQC')

D = struct('cruise', s.cruise, 'made', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
D.basepath = s.basepath;
D.chans = cfg.chans; D.allchans = chans;
D.bead_pop = cfg.bead_pop; D.k = cfg.k; D.hvc = hvc; D.window_days = cfg.bead_window_days;
D.edges = [edges(1) edges(2) - edges(1) numel(edges) - 1];
D.mode = 'fcb';
if pt
    D.mode = 'pt'; D.qc = []; D.nshift = struct();
    D.source = char(src); D.notes = cellstr(meta.notes); D.ptdir = ptdir;
    D.state_usual = char(mode(categorical(bs.state))); D.filter_usual = char(mode(categorical(bs.filterset)));
else
    D.qc = sig.qc; D.nshift = table2struct(meta.nshift);
end
D.lab = lab;
if lab
    D.source = char(src); D.notes = cellstr(meta.notes);
    dd = meta.days; dd.day = cellstr(string(dd.day, 'yyyy-MM-dd')); dd.fallback_channels = cellstr(dd.fallback_channels);
    D.lab_days = num2cell(table2struct(dd));
    D.window_days = 0; % same day only
end
D.code = char(meta.code); D.processed = char(datetime(meta.processed, 'Format', 'yyyy-MM-dd HH:mm'));
D.runs = R;
D.excluded = cellstr(meta.excluded); % marked exclude when the page was built (references above skip them)
D.beadfile_decisions = char(meta.decisions);
D.samples = S;
D.build = D.processed;
[dd, dn, de] = fileparts(decfile); [~, ddn] = fileparts(dd);
D.decisions_path = decfile; D.decisions_name = [dn de]; D.decisions_dir = dd; D.decisions_dir_name = ddn;
D.decisions = num2cell(table2struct(review_decisions(decfile, true)));
fid = fopen(fullfile(outdir, 'bead_review_data.js'), 'w', 'n', 'UTF-8');
fprintf(fid, 'const BEADS = %s;\n', jsonencode(D, 'ConvertInfAndNaN', true));
fclose(fid);
html = fileread(fullfile(here, 'bead_review_template.html'));
html = strrep(html, '__STAMP__', char(datetime('now', 'Format', 'yyyyMMddHHmmss')));
fid = fopen(fullfile(outdir, 'index.html'), 'w', 'n', 'UTF-8'); fwrite(fid, unicode2native(html, 'UTF-8')); fclose(fid);
what = 'runs'; if pt, what = 'performance tests'; end
fprintf('%s bead review: %d %s, %d sample files\n  %s\n', s.cruise, nr, what, ns, fullfile(outdir, 'index.html'))
end

function r = pt_run(bs, i, chans, ptdir, outdir, edges)
% one performance test of a PT-transfer table: peaks from the table, PT peaks,
% CVs and histograms re-measured from the test's FCS file (pt_test_peaks)
P = {'05' '1' '6'}; PN = {'0.5' '1' '6'};
r = struct('file', char(bs.filename(i)), 'time', char(datetime(bs.time(i), 'Format', 'yyyy-MM-dd''T''HH:mm:ss')));
[~, nm] = fileparts(r.file);
r.frame = ['frames/' regexprep(nm, '[^\w\-]', '_') '.png'];
r.hv = table2struct(bs.hv(i, :));
r.pops = cell(1, 3);
for ip = 1:3
    b = P{ip};
    q = struct('name', PN{ip}, 'n', bs.(['n' b])(i), 'QC', logical(bs.(['QC' b])(i)), 'why', char(bs.(['why' b])(i)));
    q.peak = cell2struct(num2cell(log10(table2array(bs.(['peak' b])(i, :)))), chans, 2);
    q.med = cell2struct(num2cell(nan(1, numel(chans))), chans, 2);
    q.rsd = q.med; q.hist = struct(); q.nnow = NaN;
    r.pops{ip} = q;
end
phys = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
r.pt = struct('state', char(bs.state(i)), 'filterset', char(bs.filterset(i)), 'n', NaN, 'nsing', NaN, ...
    'pk', struct(), 'rcv', struct(), 'hist', struct(), 'top', struct());
for c = phys, r.pt.pk.(c{1}) = NaN; r.pt.rcv.(c{1}) = NaN; end
r.all = struct(); r.read = '';
try
    [d, h] = fca_readfcs(fullfile(ptdir, r.file));
    d = double(d); names = {h.par.name};
    tconf = attune_config('', datenum(bs.time(i))); % filter layout of the test's date
    Pk = pt_test_peaks(d, h, tconf, phys);
    r.pt.n = size(d, 1); r.pt.nsing = sum(Pk.sing);
    for c = phys, r.pt.pk.(c{1}) = Pk.ch.(c{1}).pk; r.pt.rcv.(c{1}) = Pk.ch.(c{1}).rcv; end
    fig = figure('visible', 'off', 'position', [50 50 1300 380]);
    plot_pt_test(d, h, Pk, tconf, r.file)
    exportgraphics(fig, fullfile(outdir, r.frame), 'resolution', 70)
    close(fig)
    for k = 1:numel(chans)
        j = strcmp(names, file_param(chans{k}, tconf));
        if ~any(j), continue, end
        x = log10(max(d(:, j), 1));
        r.all.(chans{k}) = histcounts(x, edges);
        r.pt.hist.(chans{k}) = histcounts(x(Pk.sing), edges);
        ch = extractBefore(chans{k}, '_');
        if endsWith(chans{k}, '_H') && ~ismember(ch, {'FSC' 'SSC'})
            r.pt.top.(chans{k}) = histcounts(x(Pk.ch.(ch).use), edges);
        end
    end
catch ME
    r.read = ME.message;
    warning('bead_review:read', '%s: %s', r.file, ME.message)
end
end

function T = list_fcs(fcsdir)
% sample-file names and start times from the FCS headers ($DATE, $BTIM), when
% the cruise has no FCSfileinfo.mat yet; bead and test files left out
L = dir(fullfile(fcsdir, '*.fcs'));
L = L(~contains({L.name}, {'bead' 'Bead' 'test' 'Test'}));
t = nan(numel(L), 1);
for i = 1:numel(L)
    fid = fopen(fullfile(fcsdir, L(i).name));
    txt = char(fread(fid, 20000, 'uint8')'); fclose(fid);
    d = regexp(txt, '\$DATE/([^/]*)/', 'tokens', 'once');
    b = regexp(txt, '\$BTIM/([^/]*)/', 'tokens', 'once');
    if isempty(d) || isempty(b), continue, end
    try, t(i) = datenum([d{1} ' ' b{1}], 'dd-mmm-yyyy HH:MM:SS'); catch, end
end
ok = isfinite(t);
T = table({L(ok).name}', t(ok), 'VariableNames', {'fcslist' 'matdate_start'});
fprintf('bead_review: no FCSfileinfo.mat; %d sample files listed from the FCS headers\n', height(T))
end

function nm = file_param(c, conf)
% physical channel name (e.g. 'GL3_H', 'SSC_A') -> this cruise's file parameter (as measure_beads)
parts = split(c, '_');
nm = conf.chan.(parts{1});
if isempty(nm), return, end
nm = [extractBefore(nm, '-') '-' parts{2}];
end
