function txt = status_json(S, made, mode, busy, jobs)
% function txt = status_json(S, made, mode, busy, jobs)
% Cruise status (cruise_status results) as JSON text for the cruise console page:
% plain types only, datetimes as ISO strings, and per step the latest step record
% (process_state.mat: finished, ok, message).
%
% mode  'browser' (default: read-only page) or 'matlab' (page inside cruise_console:
%       refresh and open buttons go to MATLAB)
% busy  message shown at the top of the page (e.g. 'refreshing EN727 ...'); refresh
%       buttons are disabled while it is set
% jobs  cell array of console runs (console_job('list')), shown in the page's run list
%
% Nov 2026, cruise console

if nargin < 3, mode = 'browser'; end
if nargin < 4, busy = ''; end
if nargin < 5, jobs = {}; end
J = struct('made', dstr(made), 'mode', mode, 'busy', char(busy), 'cruises', [], 'jobs', {jobs});
C = cell(numel(S), 1);
for i = 1:numel(S)
    s = S(i);
    st = s.steps;
    for k = 1:numel(st)
        st(k).date = dstr(st(k).date);
        last = struct('finished', '', 'ok', [], 'message', '');
        if isstruct(s.last_steps) && isfield(s.last_steps, st(k).step)
            r = s.last_steps.(st(k).step);
            last.finished = dstr(r.finished); last.ok = logical(r.ok);
            if isstruct(r.info) && isfield(r.info, 'message'), last.message = char(string(r.info.message)); end
        end
        st(k).last = last;
    end
    br = s.bead_review;
    if isstruct(br) && isfield(br, 'saved'), br.saved = dstr(br.saved); end
    C{i} = struct('cruise', s.cruise, 'sample', s.sample, 'outpath', s.outpath, 'checked', dstr(s.checked), ...
        'needs_work', s.needs_work, 'archive', {cellstr(s.archive)}, 'links', s.links, 'counts', s.counts, ...
        'bead_review', br, 'steps', st);
end
J.cruises = C;
txt = jsonencode(J);
end


function v = dstr(t)
if isempty(t) || ~isdatetime(t) || isnat(t), v = ''; else, v = char(t, 'yyyy-MM-dd''T''HH:mm:ss'); end
end
