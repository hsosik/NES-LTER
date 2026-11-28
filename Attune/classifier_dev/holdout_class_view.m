function holdout_class_view(kind, outdir, nmax)
% function holdout_class_view(kind, outdir, nmax)
% Class view pages for the held-out A/B reviews: the same runs, in the same order,
% shown with the cruise review's panels, classes, Syn types and per-file QC
% (review_sample, review_viewer), classified with the current classifier. Each A/B
% page links to its class view for a run once that run has a verdict (or with blind
% off), so the blind step is kept.
%   kind 'preserved'  runs of discrete/holdout_review: <outdir>/class_view/ (CHL, PE
%                     and Pro runs) and <outdir>/class_view/hbac/ (SYBR runs);
%                     classified as holdout_review does (assign_class_discrete or
%                     assign_class_hbac, lab beads); outdir default
%                     output/discrete/holdout/review
%        'live'       files of live_holdout_review, one class view per page:
%                     <outdir>/review/class_view/ and <outdir>/review_peeuk/class_view/;
%                     assign_class_v2 with the cruise's beads; old classes go to the QC
%                     (agreement); outdir default output/live_holdout
% Run after the A/B page is built or refreshed, and again whenever the classes change.
% nmax (optional): only the first nmax runs of each page (a quick check).
% The A/B pages' index.html is rewritten from the current template, so their
% class-view links show without rebuilding their frames.
% Decisions made on a class view save to class_view_decisions.csv in its folder,
% apart from the A/B verdicts.
%
% Nov 2026, classifier redesign (held-out check)

here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here), fullfile(here, 'review'), fullfile(here, 'discrete'))
if nargin < 1 || isempty(kind), kind = 'preserved'; end
if nargin < 3 || isempty(nmax), nmax = Inf; end
[~, commit] = system(['git -C "' here '" rev-parse --short HEAD']);
build = [char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')) ' (' strtrim(commit) ')'];
switch lower(kind)
    case 'preserved'
        if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'output', 'discrete', 'holdout', 'review'); end
        preserved_view(here, outdir, build, nmax);
        ab_index(fullfile(here, 'discrete', 'holdout_review_template.html'), outdir);
    case 'live'
        if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'output', 'live_holdout'); end
        pages = {'review' 'review_peeuk'};
        for p = 1:numel(pages)
            pd = fullfile(outdir, pages{p});
            if exist(fullfile(pd, 'review_data.js'), 'file')
                live_view(outdir, pd, build, nmax);
                ab_index(fullfile(here, 'live_holdout_review_template.html'), pd);
            end
        end
    otherwise
        error('holdout_class_view:kind', 'kind must be ''preserved'' or ''live'', not %s', kind)
end
end


function preserved_view(here, outdir, build, nmax)
od = fullfile(here, 'output', 'discrete');
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
runs = ab_runs(outdir); runs = runs(1:min(end, nmax));
L = readtable(fullfile(od, 'holdout', 'holdout_lists.csv'), 'TextType', 'string');
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
cv = fullfile(outdir, 'class_view'); cvh = fullfile(cv, 'hbac');
R = review_config(); R.keep_order = true;
R.link = struct('href', '../index.html', 'text', 'A/B held-out review');
Rh = review_config_hbac(); Rh.keep_order = true;
Rh.link = struct('href', '../index.html', 'text', 'CHL, PE and Pro runs');
nh = sum(strcmp({runs.runtype}, 'hbac'));
if nh > 0, R.link(2) = struct('href', 'hbac/index.html', 'text', sprintf('SYBR runs (%d)', nh)); end
mk(cv); if nh > 0, mk(cvh); end
Q = {}; Qh = {};
warning('off', 'stats:statrobustfit:IterationLimit')
for i = 1:numel(runs)
    fn = runs(i).file;
    j = find(L.file == fn, 1);
    if isempty(j), fprintf('%s: not in holdout_lists.csv\n', fn); continue, end
    cdir = char(L.cruisedir(j)); r = char(L.runtype(j));
    cr = extractAfter(cdir, '_'); if isempty(cr), cr = cdir; end
    try
        [dat, hdr] = fca_readfcs(fullfile(root, cdir, 'preserved', 'FCS', fn));
        t = datenum([hdr.date ' ' hdr.starttime]);
        cfg = discrete_config(); cfg.cruise = cr; cfg.gl1_is_od2_scatter = true;
        ref = discrete_bead_reference(lab, t, cfg);
        if strcmp(r, 'hbac')
            [cls, info] = assign_class_hbac(dat, hdr, ref, cfg);
            Qh{end + 1} = review_sample_hbac(cls, info, cfg, cvh, fn, t, hdr, Rh); %#ok<AGROW>
        else
            [cls, info] = assign_class_discrete(dat, hdr, ref, cfg, r);
            Q{end + 1} = review_sample(dat, hdr, cls, info, [], cv, fn, t, double(strcmp(r, 'pro')), R); %#ok<AGROW>
        end
    catch err
        fprintf('%s: %s\n', fn, err.message); continue
    end
    fprintf('%d/%d %s\n', i, numel(runs), fn);
end
review_viewer(cv, Q, R, 'Held-out preserved (class view)', fullfile(cv, 'class_view_decisions.csv'), build);
if ~isempty(Qh)
    review_viewer(cvh, Qh, Rh, 'Held-out preserved (class view)', fullfile(cvh, 'class_view_decisions.csv'), build);
end
fprintf('class view: %d runs, %s\n', numel(Q) + numel(Qh), fullfile(cv, 'index.html'))
end


function live_view(outdir, pd, build, nmax)
runs = ab_runs(pd); runs = runs(1:min(end, nmax));
L = readtable(fullfile(outdir, 'holdout_lists.csv'), 'TextType', 'string');
cv = fullfile(pd, 'class_view'); mk(cv);
R = review_config(); R.keep_order = true;
R.link = struct('href', '../index.html', 'text', 'A/B held-out review');
Q = {};
warning('off', 'stats:statrobustfit:IterationLimit')
for i = 1:numel(runs)
    fn = runs(i).file;
    j = find(L.file == fn & L.set == runs(i).set, 1);
    if isempty(j), fprintf('%s: not in holdout_lists.csv\n', fn); continue, end
    cd0 = char(L.cruisedir(j));
    try
        P = cruise_paths(cd0);
        cfg = classifier_config(); cfg.cruise = extractAfter(cd0, '_');
        cfg.gl1_is_od2_scatter = strcmpi(attune_config(cfg.cruise).OD2, 'GL1');
        [dat, hdr] = fca_readfcs(fullfile(P.basepath, 'FCS', fn));
        ref = bead_reference(load_beadstat(P.bead_dir), L.time(j), cfg);
        [cls, info] = assign_class_v2(dat, hdr, ref, cfg, L.pro_measured(j));
        old = [];
        try, old = load_old_class(P, fn); catch, end %#ok<CTCH> agreement QC only
        if numel(old) ~= numel(cls), old = []; end
        Q{end + 1} = review_sample(dat, hdr, cls, info, old, cv, fn, L.time(j), L.pro_measured(j), R); %#ok<AGROW>
    catch err
        fprintf('%s: %s\n', fn, err.message); continue
    end
    fprintf('%d/%d %s\n', i, numel(runs), fn);
end
[~, pname] = fileparts(pd); ttl = 'Live held-out (class view)';
if strcmp(pname, 'review_peeuk'), ttl = 'Live held-out, PE euk only (class view)'; end
review_viewer(cv, Q, R, ttl, fullfile(cv, 'class_view_decisions.csv'), build);
fprintf('class view: %d files, %s\n', numel(Q), fullfile(cv, 'index.html'))
end


function runs = ab_runs(pd)
% the runs of an A/B page, in its order
js = fileread(fullfile(pd, 'review_data.js'));
D = jsondecode(regexprep(strtrim(js), '^window\.HOLDOUT\s*=\s*|;$', ''));
runs = D.runs; if iscell(runs), runs = [runs{:}]; end
end


function ab_index(tplfile, pd)
% the A/B page from the current template (as holdout_review / live_holdout_review write it)
tpl = fileread(tplfile);
tpl = strrep(tpl, 'review_data.js"', ['review_data.js?v=' char(datetime('now', 'Format', 'yyyyMMddHHmmss')) '"']);
fid = fopen(fullfile(pd, 'index.html'), 'w', 'n', 'UTF-8'); fwrite(fid, tpl, 'char'); fclose(fid);
end


function mk(d)
if ~exist(fullfile(d, 'frames'), 'dir'), mkdir(fullfile(d, 'frames')); end
end
