function [T, S] = cruise_status_all(opt)
% function [T, S] = cruise_status_all(Name=Value)
% Processing state of every cruise (cruise_status), cached for the cruise console.
%
% T  overview table: one row per cruise and sample type (underway, preserved where
%    the cruise has a preserved folder), one column per step with its state
%    (categorical), needs_work, and counts
% S  struct array of cruise_status results (reasons, dates, links)
%
% Options
%   Refresh  true (default): scan the cruises; false: read the cache
%   Cruises  cruises to scan (default: every cruise in attune_cruise_settings.csv);
%            with a cache, only these rows are replaced
%   Samples  ["underway" "preserved"] (default both)
%   Root     folder of cruise folders (default \\sosiknas1\Lab_data\Attune\cruise_data\)
%   Show     true (default): print the overview
%
% Cache: <Root>\processing_status\status.mat (T, S, made), with status.json, status.js
% and index.html (console_publish): the cruise console page, read-only in a browser.
%
% e.g. cruise_status_all                     % scan all, print, cache
%      T = cruise_status_all(Refresh=false); % cached overview
%      cruise_status_all(Cruises="EN727")    % refresh one cruise in the cache
%
% Nov 2026, cruise console (phase 1: status)

arguments
    opt.Refresh (1,1) logical = true
    opt.Cruises string = string.empty
    opt.Samples string = ["underway" "preserved"]
    opt.Root (1,1) string = ""
    opt.Show (1,1) logical = true
end
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fileparts(mfilename('fullpath')))
root = opt.Root;
if root == "", root = "\\sosiknas1\Lab_data\Attune\cruise_data\"; end
cdir = fullfile(root, 'processing_status');
cf = fullfile(cdir, 'status.mat');

S = [];
if isfile(cf)
    o = load(cf, 'S'); S = o.S;
end
if opt.Refresh
    cr = opt.Cruises;
    if isempty(cr)
        C = readtable(fullfile(here, 'attune_cruise_settings.csv'), 'TextType', 'string', 'Delimiter', ',');
        cr = C.cruise(C.cruise ~= "");
        S = []; % full scan: start fresh
    end
    jobs = strings(0, 2);
    for c = cr(:)'
        for smp = opt.Samples(:)'
            jobs(end+1, :) = [c smp]; %#ok<AGROW>
        end
    end
    t0 = tic; nnew = 0;
    for j = 1:size(jobs, 1)
        c = jobs(j, 1); smp = jobs(j, 2);
        try
            if smp == "preserved"
                s = attune_cruise_settings(c, char(opt.Root));
                if ~isfolder(s.preserved_path), continue, end
            end
            r = cruise_status(c, Sample=smp, Root=opt.Root);
        catch ME
            warning('cruise_status_all:cruise', '%s %s: %s', c, smp, ME.message)
            continue
        end
        if ~isempty(S)
            S = S(~([S.cruise] == r.cruise & [S.sample] == r.sample));
            S = [S(:); r];
        else
            S = r;
        end
        nnew = nnew + 1;
        if opt.Show && ~rem(nnew, 10)
            fprintf('status: %d of up to %d checked (%s elapsed)\n', nnew, size(jobs, 1), char(duration(0, 0, toc(t0), 'Format', 'mm:ss')))
        end
    end
    % order: cruise start date (folder name yyyymmdd), underway before preserved
    S = sort_status(S, root);
    T = overview(S);
    made = datetime('now');
    if ~isfolder(cdir), mkdir(cdir); end
    save(cf, 'T', 'S', 'made')
    console_publish(cdir, S, made); % status.json, status.js and the page (index.html) for browsers
    if opt.Show, fprintf('status: %d rows; cache %s\n', numel(S), cf); end
else
    if isempty(S), error('cruise_status_all:nocache', 'no status cache in %s; run with Refresh=true', cdir); end
    T = overview(S);
end
if opt.Show, show(T); end
end


function S = sort_status(S, root)
if isempty(S), return, end
key = strings(numel(S), 1);
for k = 1:numel(S)
    f = regexp(char(S(k).outpath), '(\d{8})_', 'tokens', 'once');
    if isempty(f), f = {'99999999'}; end
    key(k) = string(f{1}) + S(k).cruise + iif(S(k).sample == "preserved", "2", "1");
end
[~, i] = sort(key);
S = S(i);
end


function T = overview(S)
steps = string({S(1).steps.step});
states = categorical(strings(0, 1), ["done" "partial" "stale" "todo" "failed" "missing" "unknown" "na"]);
n = numel(S);
T = table(string([S.cruise])', string([S.sample])', 'VariableNames', {'cruise' 'sample'});
for k = 1:numel(steps)
    v = strings(n, 1);
    for i = 1:n, v(i) = string(S(i).steps(k).state); end
    T.(steps(k)) = categorical(v, categories(states));
end
T.needs_work = [S.needs_work]';
cnt = [S.counts];
T.n_sample = [cnt.n_sample]'; T.n_classified = [cnt.n_classified]'; T.n_failed = [cnt.n_failed]';
T.archived = arrayfun(@(s) ~isempty(s.archive), S(:));
T.checked = [S.checked]';
end


function show(T)
% compact print: one line per cruise and sample type, a short code per step
code = containers.Map({'done' 'partial' 'stale' 'todo' 'failed' 'missing' 'unknown' 'na'}, ...
    {'ok' 'part' 'STALE' 'todo' 'FAIL' '-' '?' ''});
steps = T.Properties.VariableNames(3:10);
fprintf('\n%-8s %-9s', 'cruise', 'sample');
for s = steps, fprintf(' %-12s', s{1}); end
fprintf('\n');
for i = 1:height(T)
    fprintf('%-8s %-9s', T.cruise(i), T.sample(i));
    for s = steps, fprintf(' %-12s', code(char(T.(s{1})(i)))); end
    fprintf('\n');
end
fprintf('%d of %d need work (ok done, part partial, STALE, todo review, FAIL, - missing, ? unknown)\n\n', sum(T.needs_work), height(T));
end


function v = iif(c, a, b)
if c, v = a; else, v = b; end
end
