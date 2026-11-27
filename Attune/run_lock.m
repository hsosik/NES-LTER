function out = run_lock(action, outpath, varargin)
% function out = run_lock(action, outpath, ...)
% One process_cruise run at a time per output folder: a lock file
% <outpath>\process_cruise.lock (JSON: host, pid, user, started, steps, sample, job,
% token) held while a run is going.
%
%   L   = run_lock('check', outpath)         L.exists, L.alive, L.info, L.file, L.message
%   tok = run_lock('acquire', outpath, info) takes the lock (info: struct of steps,
%                                            sample, job); error process_cruise:locked
%                                            when a live run holds it. A lock left by
%                                            a process that is gone (same computer) is
%                                            taken over with a warning.
%         run_lock('release', outpath, tok)  removes the lock if it is still ours
%         run_lock('clear', outpath, pid)    removes the lock if process pid holds it
%                                            (console, after cancelling a run)
%
% A lock from another computer cannot be checked: it counts as live; delete the
% file by hand if that run is not going.
%
% Nov 2026, cruise console (phase 3: runs)

fn = fullfile(char(outpath), 'process_cruise.lock');
switch action
    case 'check'
        out = check(fn);
    case 'acquire'
        info = struct(); if ~isempty(varargin), info = varargin{1}; end
        L = check(fn);
        if L.exists && L.alive
            error('process_cruise:locked', '%s', L.message)
        elseif L.exists
            warning('process_cruise:stalelock', 'taking over a lock left by a run that is not going: %s', L.message)
        end
        tok = char(java.util.UUID.randomUUID());
        r = struct('host', host(), 'pid', feature('getpid'), 'user', getenv('USERNAME'), ...
            'started', char(datetime('now'), 'yyyy-MM-dd''T''HH:mm:ss'), 'token', tok);
        for f = fieldnames(info)', r.(f{1}) = info.(f{1}); end
        fid = fopen(fn, 'w');
        if fid < 0, error('process_cruise:lockwrite', 'cannot write %s', fn), end
        fwrite(fid, jsonencode(r), 'char'); fclose(fid);
        out = tok;
    case 'release'
        L = check(fn);
        if L.exists && isfield(L.info, 'token') && strcmp(L.info.token, varargin{1}), delete(fn); end
        out = [];
    case 'clear'
        L = check(fn);
        if L.exists && isfield(L.info, 'pid') && L.info.pid == varargin{1}, delete(fn); end
        out = [];
    otherwise
        error('run_lock:action', 'unknown action %s', action)
end
end


function L = check(fn)
L = struct('exists', false, 'alive', false, 'info', struct(), 'file', fn, 'message', '');
if ~isfile(fn), return, end
L.exists = true;
try
    L.info = jsondecode(fileread(fn));
catch
    L.alive = true; L.message = sprintf('unreadable lock %s; delete it if no run is going', fn);
    return
end
i = L.info;
same = isfield(i, 'host') && strcmpi(i.host, host());
if same && isfield(i, 'pid')
    L.alive = pid_alive(i.pid);
else
    L.alive = true; % another computer: cannot check
end
steps = ''; if isfield(i, 'steps'), steps = strjoin(cellstr(string(i.steps)), ' '); end
L.message = sprintf('a run (%s) started %s on %s (pid %d, %s) holds %s%s', steps, fo(i, 'started'), fo(i, 'host'), ...
    fo(i, 'pid', 0), fo(i, 'user'), fn, iif(L.alive, '; delete the lock only if that run is not going', ' but that process is gone'));
end


function a = pid_alive(pid)
[~, o] = system(sprintf('tasklist /FI "PID eq %d" /NH /FO CSV', round(pid)));
a = contains(o, sprintf('"%d"', round(pid)));
end


function h = host()
h = getenv('COMPUTERNAME');
if isempty(h), [~, h] = system('hostname'); h = strtrim(h); end
end


function v = fo(s, f, d)
if nargin < 3, d = ''; end
if isfield(s, f), v = s.(f); else, v = d; end
end


function v = iif(c, a, b)
if c, v = a; else, v = b; end
end
