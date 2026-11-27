function out = console_job(action, varargin)
% function out = console_job(action, ...)
% process_cruise runs launched from the cruise console. Each run is a job: a JSON
% file in <Root>\processing_status\jobs\<id>.json, run in its own MATLAB process
% (matlab -batch console_run(jobfile)) with its output in <outpath>\logs\<id>.log.
%
%   J    = console_job('new', spec, root)   job file, state 'queued'. spec: cruise,
%                                           sample, steps, options (NewOnly, Force,
%                                           Review, Parallel, Files, OutPath, Root)
%   J    = console_job('preview', spec, root) the job 'new' would make, not written
%   cmd  = console_job('command', J)        the command line that runs the job
%   J    = console_job('launch', J)         starts it (refused while a live run
%                                           holds the cruise's lock: J.message)
%   J    = console_job('launch', J, true)   dry run: J.command only, nothing starts
%   J    = console_job('cancel', J)         kills the run's processes, clears its
%                                           lock, state 'cancelled'
%   jobs = console_job('list', root, days)  jobs from the last days (default 7),
%                                           newest first; running jobs whose
%                                           process is gone become 'died' (lock
%                                           cleared); running jobs get progress and
%                                           a log tail
%   J    = console_job('read', file) / console_job('write', J)
%
% States: queued, starting (process launched), running (console_run took it),
%         done, failed, cancelled, died (process gone without finishing)
%
% Nov 2026, cruise console (phase 3: runs)

switch action
    case 'new',     out = job_new(varargin{:});
    case 'preview', out = job_new(varargin{1}, varargin{2}, false);
    case 'command', out = job_command(varargin{1});
    case 'launch',  out = job_launch(varargin{:});
    case 'cancel',  out = job_cancel(varargin{1});
    case 'list',    out = job_list(varargin{:});
    case 'read',    out = job_read(varargin{1});
    case 'write',   job_write(varargin{1}); out = [];
    otherwise, error('console_job:action', 'unknown action %s', action)
end
end


function J = job_new(spec, root, write)
if nargin < 2 || root == "", root = default_root(); end
if nargin < 3, write = true; end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here))
o = struct('NewOnly', false, 'Force', false, 'Review', true, 'Parallel', [], 'Files', [], 'OutPath', '', 'Root', '');
if isfield(spec, 'options') && isstruct(spec.options)
    for f = fieldnames(spec.options)', o.(f{1}) = spec.options.(f{1}); end
end
cruise = upper(char(spec.cruise)); sample = char(spec.sample);
s = attune_cruise_settings(cruise, char(o.Root));
if ~isempty(o.OutPath)
    op = char(o.OutPath);
elseif strcmp(sample, 'preserved')
    op = fullfile(s.preserved_path, 'bead_calibrated');
else
    op = fullfile(s.basepath, 'bead_calibrated');
end
id = sprintf('%s-%s-%s', char(datetime('now'), 'yyyyMMdd-HHmmss-SSS'), cruise, sample);
jd = fullfile(char(root), 'processing_status', 'jobs');
if write && ~isfolder(jd), mkdir(jd); end
J = struct('id', id, 'cruise', cruise, 'sample', sample, 'steps', {cellstr(string(spec.steps))}, 'options', o, ...
    'state', 'queued', 'created', now_str(), 'started', '', 'finished', '', 'host', '', 'pid', [], ...
    'outpath', op, 'log', fullfile(op, 'logs', [id '.log']), 'file', fullfile(jd, [id '.json']), ...
    'command', '', 'message', '');
J.command = job_command(J);
if write, job_write(J); end
end


function cmd = job_command(J)
here = fileparts(mfilename('fullpath'));
exe = fullfile(matlabroot, 'bin', 'matlab.exe');
cmd = sprintf('"%s" -batch "addpath(''%s''); console_run(''%s'')"', exe, here, J.file);
end


function J = job_launch(J, dryrun)
if nargin < 2, dryrun = false; end
J = job_read(J.file);
J.command = job_command(J);
if dryrun, return, end
if ~strcmp(J.state, 'queued')
    J.message = sprintf('not started: state is %s', J.state); return
end
L = run_lock('check', J.outpath);
if L.exists && L.alive
    J.message = ['waiting: ' L.message];
    job_write(J);
    return
end
J.state = 'starting'; J.message = '';
job_write(J);
exe = fullfile(matlabroot, 'bin', 'matlab.exe');
args = extractAfter(J.command, strlength(sprintf('"%s" ', exe)));
psi = System.Diagnostics.ProcessStartInfo(exe, args);
psi.UseShellExecute = false; psi.CreateNoWindow = true;
pr = System.Diagnostics.Process.Start(psi);
% the launcher pid lives in its own file: once started, console_run owns the job file
fid = fopen(strrep(J.file, '.json', '.launcher'), 'w'); fprintf(fid, '%d', pr.Id); fclose(fid);
end


function J = job_cancel(J)
J = job_read(J.file);
if ismember(J.state, {'done' 'failed' 'cancelled' 'died'}), return, end
pids = [launcher_pid(J) J.pid];
for p = pids(isfinite(pids) & pids > 0)
    system(sprintf('taskkill /T /F /PID %d >NUL 2>&1', p));
end
if ~isempty(J.pid), run_lock('clear', J.outpath, J.pid); end
J = job_read(J.file); % console_run may have written meanwhile
if ~ismember(J.state, {'done' 'failed'})
    J.state = 'cancelled'; J.finished = now_str(); J.message = 'cancelled from the console';
    job_write(J);
end
end


function jobs = job_list(root, days)
if nargin < 1 || root == "", root = default_root(); end
if nargin < 2, days = 7; end
jd = fullfile(char(root), 'processing_status', 'jobs');
d = dir(fullfile(jd, '*.json'));
d = d([d.datenum] > now - days);
jobs = cell(numel(d), 1);
for k = 1:numel(d)
    try
        J = job_read(fullfile(d(k).folder, d(k).name));
    catch
        continue
    end
    J.progress = ''; J.tail = {};
    if ismember(J.state, {'starting' 'running'})
        lp = launcher_pid(J);
        alive = (~isempty(J.pid) && pid_alive(J.pid)) || (isfinite(lp) && pid_alive(lp));
        if ~alive
            J2 = job_read(J.file); % it may just have finished
            if ismember(J2.state, {'starting' 'running'})
                J2.state = 'died'; J2.finished = now_str(); J2.message = 'the run''s MATLAB process ended without finishing (crash or killed); see the log';
                job_write(J2);
                if ~isempty(J2.pid), run_lock('clear', J2.outpath, J2.pid); end
            end
            J2.progress = ''; J2.tail = {}; J = J2;
        end
    end
    jobs{k} = J;
end
jobs = jobs(~cellfun(@isempty, jobs));
if ~isempty(jobs)
    [~, i] = sort(cellfun(@(j) j.created, jobs, 'UniformOutput', false));
    jobs = jobs(flip(i));
end
% log tail and progress: running jobs and the latest finished ones
for k = 1:numel(jobs)
    if ~strcmp(jobs{k}.state, 'queued') && (k <= 15 || ismember(jobs{k}.state, {'starting' 'running'}))
        [jobs{k}.tail, jobs{k}.progress] = log_tail(jobs{k}.log, 12);
    end
end
end


function [tail, prog] = log_tail(fn, n)
tail = {}; prog = '';
if ~isfile(fn), return, end
try
    t = splitlines(string(fileread(fn)));
catch
    return
end
t = regexprep(t, '\[\x08|\]\x08', ''); % MATLAB's warning markers in a -batch diary
t = regexprep(t, '[\x00-\x08\x0B-\x1F]', '');
t = t(strlength(strtrim(t)) > 0);
t = t(~contains(t, ["statrobustfit" "robustfit" "opentoline" "In parallel_function" "Iteration limit"]));
tail = cellstr(t(max(1, end - n + 1):end));
m = regexp(t, '^(\w+): (\d+) of (\d+) files', 'tokens', 'once');
k = find(~cellfun(@isempty, m), 1, 'last');
if ~isempty(k), prog = sprintf('%s %s of %s', m{k}{1}, m{k}{2}, m{k}{3}); end
end


function J = job_read(fn)
J = jsondecode(fileread(fn));
if ~iscell(J.steps), J.steps = cellstr(string(J.steps)); end
if ~isfield(J, 'pid') || isempty(J.pid), J.pid = []; end
end


function job_write(J)
J = rmfield_if(J, {'progress' 'tail'});
tmp = [J.file '.tmp'];
fid = fopen(tmp, 'w');
if fid < 0, error('console_job:write', 'cannot write %s', tmp), end
fwrite(fid, jsonencode(J, 'PrettyPrint', true), 'char'); fclose(fid);
movefile(tmp, J.file, 'f');
end


function p = launcher_pid(J)
p = NaN;
f = strrep(J.file, '.json', '.launcher');
if isfile(f), p = str2double(fileread(f)); end
end


function a = pid_alive(pid)
[~, o] = system(sprintf('tasklist /FI "PID eq %d" /NH /FO CSV', round(pid)));
a = contains(o, sprintf('"%d"', round(pid)));
end


function s = rmfield_if(s, f)
s = rmfield(s, intersect(fieldnames(s), f));
end


function s = now_str()
s = char(datetime('now'), 'yyyy-MM-dd''T''HH:mm:ss');
end


function r = default_root()
r = "\\sosiknas1\Lab_data\Attune\cruise_data\";
end
