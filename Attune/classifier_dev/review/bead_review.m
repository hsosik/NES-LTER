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
% opts (all optional):
%   .root       cruise_data root passed to attune_cruise_settings
%   .beaddir    folder with beadstat.mat and FCSfileinfo.mat (default bead_calibrated)
%   .outdir     viewer folder (default <beaddir>\bead_review)
% e.g. bead_review('EN727')
%
% Nov 2026, bead review

if nargin < 2, opts = struct(); end
def = struct('root', '', 'beaddir', '', 'outdir', '');
for f = fieldnames(def)', if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end, end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)) % classifier_dev
s = attune_cruise_settings(cruise, opts.root);
bp = opts.beaddir; if isempty(bp), bp = fullfile(s.basepath, 'bead_calibrated'); end
outdir = opts.outdir; if isempty(outdir), outdir = fullfile(bp, 'bead_review'); end
decfile = fullfile(bp, 'bead_review_decisions.csv'); % where load_beadstat reads exclusions
if ~exist(fullfile(outdir, 'frames'), 'dir'), mkdir(fullfile(outdir, 'frames')); end

[bs, meta] = load_beadstat(bp, false); % as measured: the page shows QC before exclusions
F = load(fullfile(bp, 'FCSfileinfo.mat')); fn = fieldnames(F); Fi = F.(fn{1});
cfg = classifier_config();
cfg.cruise = s.cruise;
conf = attune_config(s.cruise);
P = {'05' '1' '6'}; PN = {'0.5' '1' '6'};
chans = bs.peak1.Properties.VariableNames; % FSC_H ... GL4_H, SSC_A, GL1_A
bt = datenum(bs.time);
nr = height(bs);

% common voltage per channel: the most frequent among the runs
hvc = struct();
for c = cfg.chans
    h = bs.hv.(c{1}); h = h(isfinite(h) & h > 0);
    if isempty(h), hvc.(c{1}) = NaN; else, hvc.(c{1}) = mode(h); end
end

% runs: re-gate with the cruise's signature for the plots and histograms
edges = 0:0.05:6.1;
sig = meta.signature;
R = cell(nr, 1);
for i = 1:nr
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
        [d, h] = fca_readfcs(fullfile(s.basepath, 'FCS', r.file));
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
    ref = bead_reference(bsx, t(k), cfg);
    for c = cfg.chans
        ch = c{1};
        u = find(ismember(bsx.filename, ref.runs.(ch)));
        S.runs.(ch){k} = u(:)' - 1; % 0-based for the page
        S.nearest.(ch)(k) = ~isempty(u) && all(abs(bt(u) - t(k)) > cfg.bead_window_days);
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
D.qc = sig.qc; D.nshift = table2struct(meta.nshift);
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
fprintf('%s bead review: %d runs, %d sample files\n  %s\n', s.cruise, nr, ns, fullfile(outdir, 'index.html'))
end

function nm = file_param(c, conf)
% physical channel name (e.g. 'GL3_H', 'SSC_A') -> this cruise's file parameter (as measure_beads)
parts = split(c, '_');
nm = conf.chan.(parts{1});
if isempty(nm), return, end
nm = [extractBefore(nm, '-') '-' parts{2}];
end
