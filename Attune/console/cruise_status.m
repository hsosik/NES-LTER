function S = cruise_status(cruise, opt)
% function S = cruise_status(cruise, Sample="underway", Root="", OutPath="")
% Processing state of one cruise (underway or preserved samples) for the cruise
% console: one state and reason per step, from the outputs in bead_calibrated and
% the step records process_cruise writes (process_state.mat). Reads only small
% files and folder listings (no class files), so it is cheap enough to scan every
% cruise (cruise_status_all).
%
% S.steps  struct array, one per column of the console, in this order:
%   files        FCSfileinfo vs the FCS folder (new files not listed: partial)
%   beads        beadstat.mat; underway: bead runs acquired after the table: stale
%   bead_review  bead_review_status (page button "mark review complete")
%   classify     class_log vs the sample files; stale when made with another
%                classifier_version tag (or before tags), before the bead table, or
%                before a bead review run decision; failures (e.g. no lab beads)
%   class_review review frames in sync with the class files; files with a decision
%   volume       mergeTable.mat rows; files reclassified after it have no volume
%   table        AttuneTable / DiscreteTable newer than classify and volume;
%                underway: class review decisions newer: apply_review_excludes
%   metadata     underway match / BottleMetadata
% each with fields step, state, reason, date (of the output), detail.
% States: done, partial, stale, missing, failed, todo (needs a person: a review),
%         unknown (made before step records), na (not used for this sample type)
%
% S.cruise, S.sample, S.outpath, S.exists, S.checked, S.needs_work, S.archive,
% S.links (output folder, review pages), S.counts, S.bead_review (bead_review_status),
% S.last_steps (latest process_state record per step)
%
% Nov 2026, cruise console (phase 1: status)

arguments
    cruise {mustBeTextScalar}
    opt.Sample (1,1) string {mustBeMember(opt.Sample, ["underway" "preserved"])} = "underway"
    opt.Root (1,1) string = ""
    opt.OutPath (1,1) string = "" % output folder instead of bead_calibrated (testing, as process_cruise)
end
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(here, 'classifier_dev'), fullfile(here, 'classifier_dev', 'review'))

s = attune_cruise_settings(cruise, char(opt.Root));
pres = opt.Sample == "preserved";
if pres
    base = s.preserved_path; excl = s.exclude_preserved;
else
    base = s.basepath; excl = s.exclude;
end
op = fullfile(base, 'bead_calibrated');
if opt.OutPath ~= "", op = char(opt.OutPath); end
S = struct('cruise', string(s.cruise), 'sample', opt.Sample, 'outpath', string(op), 'exists', isfolder(op), ...
    'checked', datetime('now'), 'needs_work', false, 'archive', strings(0, 1), 'links', struct(), ...
    'counts', struct(), 'bead_review', struct(), 'last_steps', struct(), 'steps', []);
S.links = struct('outpath', string(op), 'bead_review', string(fullfile(op, 'bead_review', 'index.html')), ...
    'class_review', string(fullfile(op, 'review', 'index.html')), 'hbac_review', string(fullfile(op, 'review', 'hbac', 'index.html')));
a = dir(fullfile(base, 'bead_calibrated_archive_*'));
S.archive = string({a([a.isdir]).name})';

names = ["files" "beads" "bead_review" "classify" "class_review" "volume" "table" "metadata"];
st = struct('step', cellstr(names), 'state', 'missing', 'reason', '', 'date', NaT, 'detail', '');

% step records
last = struct();
if isfile(fullfile(op, 'process_state.mat'))
    try
        o = load(fullfile(op, 'process_state.mat'), 'steps');
        for r = o.steps(:)'
            last.(char(r.step)) = r; % later records replace earlier ones
        end
    catch
    end
end
S.last_steps = last;

% ---- files
fcs = dir(fullfile(base, 'FCS', '*.fcs'));
fcsnames = string({fcs.name})';
if isempty(fcsnames)
    % no FCS files of this sample type (e.g. a cruise with preserved samples only)
    for k = 1:numel(st), st(k).state = 'na'; st(k).reason = 'no FCS files'; end
    S.steps = st;
    S.counts = struct('n_fcs', 0, 'n_sample', 0, 'n_classified', 0, 'n_failed', 0, 'n_volume', NaN, 'n_reviewed', 0, 'n_excluded', 0);
    return
end
fi = fullfile(op, 'FCSfileinfo.mat');
Fi = table();
k = 1;
if isfile(fi)
    o = load(fi, 'FCSfileinfo'); Fi = o.FCSfileinfo;
    if ~istable(Fi), Fi = table('Size', [0 2], 'VariableTypes', {'cell' 'double'}, 'VariableNames', {'fcslist' 'matdate_start'}); end
    st(k).date = fdate(fi);
    new = setdiff(fcsnames, string(Fi.fcslist));
    if ~isempty(new)
        st(k) = set1(st(k), 'partial', sprintf('%d new FCS files not listed (run files)', numel(new)));
    else
        st(k) = set1(st(k), 'done', sprintf('%d FCS files', height(Fi)));
    end
else
    st(k).reason = 'no FCSfileinfo.mat';
end
st(k) = ledger_fail(st(k), last, 'files');
% sample files (as classify chooses them)
smp = strings(0, 1); smpclass = strings(0, 1);
if height(Fi)
    keep = ~contains(Fi.fcslist, excl);
    if pres && ismember('runtype', Fi.Properties.VariableNames), keep = keep & ~strcmp(Fi.runtype, 'other'); end
    smp = string(Fi.fcslist(keep));
    if pres && ismember('runtype', Fi.Properties.VariableNames)
        smpclass = repmat("assign_class_discrete", numel(smp), 1);
        smpclass(strcmp(Fi.runtype(keep), 'hbac')) = "assign_class_hbac";
    else
        smpclass = repmat("assign_class_v2", numel(smp), 1);
    end
end

% ---- beads
k = 2;
bf = fullfile(op, 'beadstat.mat');
bdate = fdate(bf);
if isfile(bf)
    st(k).date = bdate;
    src = '';
    try
        m = load(bf, 'meta'); if isfield(m.meta, 'source'), src = char(string(m.meta.source)); end
    catch
    end
    st(k) = set1(st(k), 'done', ['bead table' iif(isempty(src), '', [' (' src ')'])]);
    if ~pres && height(Fi)
        isb = contains(Fi.fcslist, {'FCB_bead' 'fcb_bead'});
        tb = datetime(Fi.matdate_start(isb), 'ConvertFrom', 'datenum');
        nnew = sum(tb > bdate);
        if nnew, st(k) = set1(st(k), 'stale', sprintf('%d bead runs acquired after the bead table', nnew)); end
    end
else
    st(k).reason = 'no beadstat.mat';
end
st(k) = ledger_fail(st(k), last, 'beads');

% ---- bead review
k = 3;
try
    br = bead_review_status(op);
catch ME
    br = struct('complete', false, 'state', 'error', 'reviewer', '', 'saved', NaT, 'message', ME.message);
end
S.bead_review = br;
switch br.state
    case 'complete'
        st(k) = set1(st(k), 'done', sprintf('complete (%s)', br.reviewer));
        st(k).date = br.saved;
    case 'no bead table'
        st(k) = set1(st(k), 'missing', 'no bead table');
    otherwise
        st(k) = set1(st(k), 'todo', br.state);
end
st(k).detail = br.message;
% last time a run decision changed (exclusions used by classify and volume)
bdec = decisions(fullfile(op, 'bead_review_decisions.csv'), true);

% ---- classify
k = 4;
lf = fullfile(op, 'class_log.mat');
cdir = dir(fullfile(op, 'class', '*.mat'));
L = table();
cnt = struct('n_fcs', numel(fcsnames), 'n_sample', numel(smp), 'n_classified', 0, 'n_failed', 0, ...
    'n_volume', NaN, 'n_reviewed', 0, 'n_excluded', 0);
if isfile(lf)
    o = load(lf, 'L'); L = o.L;
    st(k).date = fdate(lf);
    L = L(ismember(string(L.file), smp), :);
    good = string(L.error) == "" | ismissing(string(L.error));
    cnt.n_classified = sum(good); cnt.n_failed = sum(~good);
    if ismember('classified', L.Properties.VariableNames)
        ct = L.classified;
    else
        ct = repmat(st(k).date, height(L), 1); % before per-file times: the log's date
    end
    R = {};
    nmiss = numel(smp) - height(L);
    if nmiss > 0, R{end+1} = sprintf('%d of %d sample files not classified', nmiss, numel(smp)); end
    if cnt.n_failed
        nob = sum(contains(string(L.error(~good)), 'nobeads'));
        R{end+1} = sprintf('%d failed%s', cnt.n_failed, iif(nob > 0, sprintf(' (%d no lab beads that day)', nob), ''));
    end
    state = iif(isempty(R), 'done', 'partial');
    % stale: classifier version, bead table, bead review decisions
    stl = {};
    if ismember('classifier_tag', L.Properties.VariableNames)
        [~, ia] = ismember(string(L.file), smp);
        cur = strings(numel(ia), 1);
        for u = unique(smpclass(ia))', cur(smpclass(ia) == u) = classifier_version(u); end
        nold = sum(good & string(L.classifier_tag) ~= cur);
        if nold, stl{end+1} = sprintf('%d files from an older classifier version', nold); end
    elseif cnt.n_classified
        stl{end+1} = 'classified before classifier version tags';
    end
    if ~isnat(bdate)
        n = sum(good & ct < bdate);
        if n, stl{end+1} = sprintf('%d files classified before the bead table was rebuilt', n); end
    end
    if ~isnat(bdec)
        n = sum(good & ct < bdec);
        if n, stl{end+1} = sprintf('%d files classified before a bead review decision changed', n); end
    end
    if ~isempty(stl), state = 'stale'; R = [stl R]; end
    st(k) = set1(st(k), state, iif(isempty(R), sprintf('%d files', cnt.n_classified), strjoin(R, '; ')));
    st(k).detail = sprintf('%d of %d sample files classified, %d failed', cnt.n_classified, numel(smp), cnt.n_failed);
elseif ~isempty(cdir)
    st(k) = set1(st(k), 'stale', sprintf('%d class files but no class_log (old pipeline)', numel(cdir)));
    st(k).date = max(datetime([cdir.datenum], 'ConvertFrom', 'datenum'));
else
    st(k).reason = 'not classified';
end
st(k) = ledger_fail(st(k), last, 'classify');
cdate = st(k).date;

% ---- class review
k = 5;
rq = fullfile(op, 'review', 'review_qc.mat');
if isfile(rq)
    st(k).date = fdate(rq);
    [~, rdec, T] = decisions(fullfile(op, 'review_decisions.csv'), false);
    if height(T)
        T = T(ismember(string(T.file), smp), :);
        cnt.n_reviewed = height(T); cnt.n_excluded = sum(string(T.status) == "exclude");
    end
    if ~isnat(cdate) && st(k).date < cdate - minutes(1)
        st(k) = set1(st(k), 'stale', 'review frames older than the class files');
    elseif cnt.n_reviewed == 0
        st(k) = set1(st(k), 'todo', 'not reviewed');
    elseif cnt.n_reviewed < cnt.n_classified
        st(k) = set1(st(k), 'todo', sprintf('%d of %d files reviewed, %d excluded', cnt.n_reviewed, cnt.n_classified, cnt.n_excluded));
    else
        st(k) = set1(st(k), 'done', sprintf('all %d files reviewed, %d excluded', cnt.n_reviewed, cnt.n_excluded));
    end
    if ~isnat(rdec), st(k).detail = sprintf('last decision %s', char(rdec)); end
elseif isfile(lf) || ~isempty(cdir)
    st(k).reason = 'no classification review';
else
    st(k) = set1(st(k), 'missing', 'not classified');
end
if pres && isfile(fullfile(op, 'review', 'hbac', 'review_qc.mat'))
    st(k).detail = strtrim([st(k).detail ' (SYBR runs: review\hbac)']);
end

% ---- volume
k = 6;
mf = fullfile(op, 'mergeTable.mat');
if isfile(mf)
    st(k).date = fdate(mf);
    o = load(mf, 'mergeT'); M = o.mergeT;
    hasv = string(M.filename(isfinite(M.SSC_hv)));
    if height(L)
        good = string(L.error) == "" | ismissing(string(L.error));
        f = string(L.file(good));
        if ismember('classified', L.Properties.VariableNames)
            redone = L.classified(good) > st(k).date + seconds(5); % reclassified without volume (file dates are to the second)
        else
            redone = repmat(cdate > st(k).date + minutes(1), numel(f), 1);
        end
        withv = ismember(f, hasv) & ~redone;
        cnt.n_volume = sum(withv);
        if any(redone)
            st(k) = set1(st(k), 'stale', sprintf('%d files reclassified after the volume step', sum(redone)));
        elseif cnt.n_volume < numel(f)
            st(k) = set1(st(k), 'partial', sprintf('%d of %d classified files without a volume', numel(f) - cnt.n_volume, numel(f)));
        else
            st(k) = set1(st(k), 'done', sprintf('%d files', cnt.n_volume));
        end
    else
        st(k) = set1(st(k), 'done', sprintf('%d files', numel(hasv)));
    end
elseif isfield(last, 'volume') && last.volume.ok
    st(k) = set1(st(k), iif(~isnat(cdate) && cdate > last.volume.finished + minutes(1), 'stale', 'done'), 'step record only (no mergeTable)');
    st(k).date = last.volume.finished;
elseif pres && (isfile(lf) || ~isempty(cdir))
    st(k) = set1(st(k), 'unknown', 'volume step before step records (no mergeTable for preserved runs then)');
else
    st(k).reason = 'no volumes';
end
st(k) = ledger_fail(st(k), last, 'volume');

% ---- table
k = 7;
tf = fullfile(op, iif(pres, 'DiscreteTable.mat', 'AttuneTable.mat'));
if isfile(tf)
    st(k).date = fdate(tf);
    up = max([cdate st(6).date]);
    if ~isnat(up) && st(k).date < up - minutes(1)
        st(k) = set1(st(k), 'stale', 'older than the class files or volumes');
    else
        st(k) = set1(st(k), 'done', '');
    end
    if ~pres
        [~, rdec] = decisions(fullfile(op, 'review_decisions.csv'), false);
        if ~isnat(rdec) && rdec > st(k).date
            st(k) = set1(st(k), 'stale', strjoin({st(k).reason 'class review decisions newer (apply_review_excludes)'}, '; '));
        end
    end
else
    st(k).reason = ['no ' iif(pres, 'DiscreteTable.mat', 'AttuneTable.mat')];
end
st(k) = ledger_fail(st(k), last, 'table');

% ---- metadata
k = 8;
mdf = fullfile(op, iif(pres, 'BottleMetadata.mat', 'AttuneTable_uw_match.mat'));
if isfile(mdf)
    st(k).date = fdate(mdf);
    if ~pres && ~isnat(st(7).date) && st(k).date < st(7).date - minutes(1)
        st(k) = set1(st(k), 'stale', 'older than AttuneTable');
    else
        st(k) = set1(st(k), 'done', '');
    end
else
    st(k).reason = ['no ' iif(pres, 'BottleMetadata.mat', 'AttuneTable_uw_match.mat')];
end
st(k) = ledger_fail(st(k), last, 'metadata');

for k = 1:numel(st), st(k).reason = strip(regexprep(st(k).reason, '^; |; $', '')); end
S.steps = st;
S.counts = cnt;
S.needs_work = any(~ismember({st.state}, {'done' 'na'}));
end


function t = fdate(fn)
d = dir(fn);
if isempty(d), t = NaT; else, t = datetime(d(1).datenum, 'ConvertFrom', 'datenum'); end
end


function x = set1(x, state, reason)
x.state = state; x.reason = reason;
end


function x = ledger_fail(x, last, step)
% the latest record of this step failed: failed, unless the output is newer
if isfield(last, step) && ~last.(step).ok
    r = last.(step);
    if isnat(x.date) || x.date < r.finished
        msg = ''; if isfield(r.info, 'message'), msg = char(r.info.message); end
        x.state = 'failed';
        x.reason = sprintf('last run failed %s: %s', char(r.finished, 'dd-MMM HH:mm'), msg);
    end
end
end


function [tlast, tlast_all, T] = decisions(fn, runs_only)
% latest 'saved' time in a review decisions csv (local time); runs_only: leave out
% the bead review's completion row. T: the decisions table
tlast = NaT; tlast_all = NaT; T = table();
if ~isfile(fn), return, end
try
    o = detectImportOptions(fn, 'Delimiter', ',', 'TextType', 'string');
    o = setvartype(o, intersect(o.VariableNames, {'file' 'status' 'note' 'reviewer' 'saved' 'build'}), 'string');
    T = readtable(fn, o);
catch
    return
end
if ~height(T) || ~ismember('saved', T.Properties.VariableNames), return, end
t = parse_utc(T.saved);
if runs_only, keep = string(T.file) ~= "__review_complete__"; else, keep = true(height(T), 1); end
if runs_only, T = T(keep, :); end
if any(keep & ~isnat(t)), tlast = max(t(keep)); end
tlast_all = tlast;
end


function t = parse_utc(s)
s = string(s);
t = NaT(size(s), 'TimeZone', 'UTC');
for f = ["yyyy-MM-dd'T'HH:mm:ss.SSSXXX" "yyyy-MM-dd'T'HH:mm:ssXXX"]
    j = isnat(t) & s ~= "" & ~ismissing(s);
    if ~any(j), break, end
    try
        t(j) = datetime(s(j), 'InputFormat', f, 'TimeZone', 'UTC');
    catch
        for i = find(j)'
            try, t(i) = datetime(s(i), 'InputFormat', f, 'TimeZone', 'UTC'); catch, end
        end
    end
end
t.TimeZone = 'local'; t.TimeZone = ''; % local clock, as file dates
end


function v = iif(c, a, b)
if c, v = a; else, v = b; end
end
