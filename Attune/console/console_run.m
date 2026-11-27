function console_run(jobfile)
% function console_run(jobfile)
% Runs one cruise console job (console_job) in this MATLAB process: called by the
% process console_job('launch') starts (matlab -batch). Marks the job running,
% writes the output to the job's log (diary), runs process_cruise with the job's
% steps and options (process_cruise holds the cruise's run lock), marks the job done
% or failed, and refreshes that cruise's row of the status cache.
%
% Nov 2026, cruise console (phase 3: runs)

here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here))
J = console_job('read', jobfile);
J.state = 'running'; J.pid = feature('getpid'); J.host = getenv('COMPUTERNAME');
J.started = char(datetime('now'), 'yyyy-MM-dd''T''HH:mm:ss');
console_job('write', J);
lf = J.log;
if ~isfolder(fileparts(lf)), mkdir(fileparts(lf)); end
diary(lf);
fprintf('cruise console job %s: %s %s, steps %s\n', J.id, J.cruise, J.sample, strjoin(J.steps, ' '));
setenv('CRUISE_CONSOLE_JOB', J.id);
o = J.options;
args = {'Steps', string(J.steps), 'Sample', string(J.sample), 'NewOnly', logical(o.NewOnly), ...
    'Force', logical(o.Force), 'Review', logical(o.Review)};
if ~isempty(o.Parallel), args = [args {'Parallel', logical(o.Parallel)}]; end
if ~isempty(o.Files), args = [args {'Files', o.Files}]; end
if ~isempty(o.OutPath), args = [args {'OutPath', string(o.OutPath)}]; end
if ~isempty(o.Root), args = [args {'Root', string(o.Root)}]; end
try
    process_cruise(J.cruise, args{:});
    state = 'done'; msg = '';
catch ME
    state = 'failed'; msg = ME.message;
    fprintf(2, 'job failed: %s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
J = console_job('read', jobfile);
if ~strcmp(J.state, 'cancelled')
    J.state = state; J.message = msg;
    J.finished = char(datetime('now'), 'yyyy-MM-dd''T''HH:mm:ss');
    console_job('write', J);
end
fprintf('cruise console job %s: %s\n', J.id, J.state);
if isempty(o.OutPath) % a test run in another folder is not in the status cache
    try
        cruise_status_all(Cruises=string(J.cruise), Samples=string(J.sample), Root=string(o.Root), Show=false);
    catch ME
        fprintf('status not refreshed: %s\n', ME.message);
    end
end
diary off
end
