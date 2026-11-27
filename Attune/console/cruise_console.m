function fig = cruise_console(opt)
% function fig = cruise_console(Name=Value)
% Cruise console: processing state of every cruise (underway and preserved) in one
% window, and process_cruise runs launched from it. The page
% (cruise_console_page.html) runs inside a MATLAB uihtml component; status comes
% from the cruise_status_all cache. In the window:
%   - overview grid: one row per cruise and sample type, a coloured state per step
%     (files, beads, bead review, classify, class review, volume, table, metadata);
%     hover a cell for the reason; filters for sample type and state; rows with a
%     queued or running run are marked
%   - details of the selected cruise: per step state, reason, output date and the
%     latest step record (process_state.mat), counts, archive
%   - buttons: refresh this cruise or all cruises (cruise_status_all, updates the
%     cache and the browser copy), open the output folder and the bead, class and
%     SYBR review pages
%   - run form: steps (or "steps that need work") and options for the selected row
%     or the ticked rows; "Show command" (dry run: nothing starts) or "Run". Each
%     run is a job (console_job) in its own MATLAB process (console_run), logged to
%     <outpath>\logs\. Runs start one at a time (MaxRuns); a cruise whose lock is
%     held by a live run waits. Force, several cruises, or a cruise that already has
%     a run queued ask to confirm.
%   - run list: queued, running and recent runs, progress, log tail, cancel
%   - NAS drive: the letter this computer maps \\sosiknas1\Lab_data to (default Z,
%     kept for later sessions: nas_drive). Paths in the page and the open buttons
%     use that drive: review pages save decisions only when opened from it.
% Runs keep going when the console is closed; queued runs start when a console is
% open (it checks every few seconds). Use one console at a time.
% The same page, read-only, is written next to the cache for browsers
% (<Root>\processing_status\index.html, console_publish).
%
% Options
%   Refresh  true: scan all cruises before opening (default false: the cache; scans
%            when there is no cache)
%   Root     folder of cruise folders (default \\sosiknas1\Lab_data\Attune\cruise_data\)
%   JobRoot  folder whose processing_status\jobs holds the job files (default Root;
%            testing)
%   MaxRuns  runs at a time (default 1: a run's parallel pool takes every core)
%   Period   seconds between checks of the runs (default 4)
%   TestOutPath  testing only: every run writes to this folder instead of the cruise's
%            bead_calibrated (process_cruise OutPath)
%
% Nov 2026, cruise console (phases 2 and 3)

arguments
    opt.Refresh (1,1) logical = false
    opt.Root (1,1) string = ""
    opt.JobRoot (1,1) string = ""
    opt.MaxRuns (1,1) double {mustBePositive} = 1
    opt.Period (1,1) double {mustBePositive} = 4
    opt.TestOutPath (1,1) string = ""
end
here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here))
root = opt.Root;
if root == "", root = "\\sosiknas1\Lab_data\Attune\cruise_data\"; end
jroot = opt.JobRoot; if jroot == "", jroot = root; end
cf = fullfile(root, 'processing_status', 'status.mat');
if opt.Refresh || ~isfile(cf)
    cruise_status_all(Root=opt.Root, Show=false);
end
o = load(cf, 'S', 'made');
st = struct('S', o.S, 'made', o.made, 'root', opt.Root, 'cache', cf, 'jroot', jroot, 'maxruns', opt.MaxRuns, ...
    'jobs', {console_job('list', jroot)}, 'msg', "", 'last', "", 'testout', opt.TestOutPath, 'drive', nas_drive());

fig = uifigure('Name', 'Attune cruise console', 'Position', [80 80 1500 900]);
g = uigridlayout(fig, [1 1], 'Padding', [0 0 0 0]);
h = uihtml(g, 'HTMLSource', fullfile(here, 'cruise_console_page.html'));
h.UserData = st;
h.HTMLEventReceivedFcn = @(src, evt) console_event(src, evt);
push(h);
t = timer('ExecutionMode', 'fixedSpacing', 'Period', opt.Period, 'BusyMode', 'drop', ...
    'TimerFcn', @(~, ~) poll(h), 'Name', 'cruise_console');
fig.CloseRequestFcn = @(f, ~) close_console(f, t);
fig.UserData = struct('timer', t, 'html', h);
start(t);
end


function push(h, msg)
% send the page its data (only when something changed)
st = h.UserData;
if nargin > 1, st.msg = string(msg); end
txt = string(status_json(st.S, st.made, 'matlab', st.msg, st.jobs, st.drive));
if txt ~= st.last
    st.last = txt;
    h.UserData = st;
    h.Data = char(txt);
else
    h.UserData = st;
end
end


function poll(h)
% runs: update states, start queued runs, reload the status cache when a run
% refreshed it
if ~isvalid(h), return, end
st = h.UserData;
try
    st.jobs = console_job('list', st.jroot);
    nrun = sum(cellfun(@(j) ismember(j.state, {'starting' 'running'}), st.jobs));
    q = find(cellfun(@(j) strcmp(j.state, 'queued'), st.jobs));
    for k = flip(q(:))' % oldest first
        if nrun >= st.maxruns, break, end
        J = console_job('launch', st.jobs{k});
        if strcmp(J.state, 'starting'), nrun = nrun + 1; end
    end
    if ~isempty(q), st.jobs = console_job('list', st.jroot); end
    d = dir(st.cache);
    if ~isempty(d) && datetime(d.datenum, 'ConvertFrom', 'datenum') > st.made + seconds(1)
        o = load(st.cache, 'S', 'made'); st.S = o.S; st.made = o.made;
    end
catch ME
    st.msg = "run check failed: " + ME.message;
end
h.UserData = st;
push(h);
end


function console_event(h, evt)
% events from the page: 'refresh' {cruise, sample} ('' = all), 'open' path,
% 'run' {targets, steps, options, dryrun}, 'cancel' job id, 'drive' NAS drive letter
st = h.UserData;
fig = ancestor(h, 'figure');
switch evt.HTMLEventName
    case 'open'
        msg = console_open(char(evt.HTMLEventData));
        push(h, msg);
    case 'refresh'
        d = evt.HTMLEventData;
        cr = string(d.cruise); smp = string(d.sample);
        if cr == "", what = "all cruises (about 2 minutes)"; else, what = cr + " " + smp; end
        push(h, "refreshing " + what + " ...");
        drawnow
        try
            st = h.UserData;
            st = console_refresh(st, cr, smp);
            h.UserData = st;
            push(h, "");
        catch ME
            push(h, "refresh failed: " + ME.message);
        end
    case 'run'
        push(h, run_event(h, fig, evt.HTMLEventData));
        poll(h); % start it now if nothing else is running
    case 'drive'
        try
            st.drive = nas_drive('set', evt.HTMLEventData);
            h.UserData = st;
            push(h, "NAS drive set to " + st.drive + ": (kept for later sessions)");
        catch ME
            push(h, ME.message);
        end
    case 'jserror'
        warning('cruise_console:page', 'console page error: %s', char(evt.HTMLEventData));
    case 'cancel'
        id = char(evt.HTMLEventData);
        k = find(cellfun(@(j) strcmp(j.id, id), st.jobs), 1);
        if isempty(k), push(h, "run not found: " + id); return, end
        J = st.jobs{k};
        if ismember(J.state, {'starting' 'running'}) && ~confirm(fig, sprintf('Cancel the running %s %s run (%s)? Its processes are killed; files it was writing may be left half done.', ...
                J.cruise, J.sample, strjoin(J.steps, ' ')), 'Cancel run')
            return
        end
        J = console_job('cancel', J);
        push(h, sprintf('%s %s run %s', J.cruise, J.sample, J.state));
        poll(h);
end
end


function msg = run_event(h, fig, d)
st = h.UserData;
tg = d.targets; if ~iscell(tg), tg = num2cell(tg); end
steps = string(d.steps);
if isempty(steps), msg = "tick at least one step"; return, end
o = d.options;
spec0 = struct('steps', {cellstr(steps)}, 'options', struct('NewOnly', fo(o, 'NewOnly', false), 'Force', fo(o, 'Force', false), ...
    'Review', fo(o, 'Review', true), 'Parallel', fo(o, 'Parallel', []), 'Files', fo(o, 'Files', []), 'OutPath', char(st.testout), 'Root', char(st.root)));
if d.dryrun
    lines = strings(0, 1);
    for k = 1:numel(tg)
        J = console_job('preview', setfield(setfield(spec0, 'cruise', tg{k}.cruise), 'sample', tg{k}.sample), st.jroot); %#ok<SFLD>
        lines(end+1) = sprintf('%s %s: process_cruise("%s", %s)', J.cruise, J.sample, J.cruise, args_text(J)); %#ok<AGROW>
    end
    msg = "dry run, nothing started. " + strjoin(lines, " | ");
    return
end
% confirmations
names = strings(numel(tg), 1);
for k = 1:numel(tg), names(k) = string(tg{k}.cruise) + " " + string(tg{k}.sample); end
busy = names(arrayfun(@(k) any(cellfun(@(j) strcmp(j.cruise, tg{k}.cruise) && strcmp(j.sample, tg{k}.sample) && ...
    ismember(j.state, {'queued' 'starting' 'running'}), st.jobs)), 1:numel(tg)));
q = strings(0, 1);
if numel(tg) > 1, q(end+1) = sprintf('Queue %d runs (%s), one at a time?', numel(tg), strjoin(names, ', ')); end
if spec0.options.Force, q(end+1) = "Force is on: it rebuilds FCSfileinfo and replaces bead tables not made from FCB beads."; end
if ~isempty(busy), q(end+1) = "Already queued or running: " + strjoin(busy, ', ') + ". Queue another run after it?"; end
if ~isempty(q) && ~confirm(fig, strjoin(q, [newline newline]), 'Run process_cruise'), msg = "not started"; return, end
skipped = strings(0, 1);
for k = 1:numel(tg)
    spec = spec0; spec.cruise = tg{k}.cruise; spec.sample = tg{k}.sample;
    if strcmp(spec.sample, 'preserved') % movie and voldists are underway steps
        spec.steps = setdiff(spec.steps, {'movie' 'voldists'}, 'stable');
        if isempty(spec.steps), skipped(end+1) = names(k); continue, end %#ok<AGROW>
    end
    console_job('new', spec, st.jroot);
end
names = setdiff(names, skipped, 'stable');
if isempty(names), msg = "not started: movie and voldists are not used for preserved samples"; return, end
msg = sprintf('queued: %s (steps %s)', strjoin(names, ', '), strjoin(steps, ' '));
end


function s = args_text(J)
o = J.options;
s = sprintf('Steps=[%s], Sample="%s"', strjoin("""" + string(J.steps) + """", ' '), J.sample);
if o.NewOnly, s = s + ", NewOnly=true"; end
if o.Force, s = s + ", Force=true"; end
if ~o.Review, s = s + ", Review=false"; end
if ~isempty(o.Parallel), s = s + ", Parallel=" + string(logical(o.Parallel)); end
if ~isempty(o.Files), s = s + ", Files=" + string(o.Files); end
if ~isempty(o.OutPath), s = s + ", OutPath=""" + string(o.OutPath) + """"; end
s = s + "; runs as: " + string(J.command);
end


function ok = confirm(fig, msg, title)
if isempty(fig) || ~isvalid(fig) || strcmp(getenv('CRUISE_CONSOLE_TEST'), 'yes')
    ok = true; return % tests: no dialog
end
ok = strcmp(uiconfirm(fig, msg, title, 'Options', {'OK' 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2), 'OK');
end


function v = fo(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end


function close_console(f, t)
try, stop(t); delete(t); catch, end
delete(f);
end
