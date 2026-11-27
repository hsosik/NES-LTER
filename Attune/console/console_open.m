function msg = console_open(target, dryrun)
% function msg = console_open(target, dryrun)
% Cruise console 'open' buttons: opens an output folder in Explorer or a review page
% (index.html) in the default browser. msg: '' when opened, otherwise why not.
% dryrun true: only check (tests).
%
% Nov 2026, cruise console

if nargin < 2, dryrun = false; end
msg = "";
target = char(target);
if ~(isfolder(target) || isfile(target))
    msg = "not found: " + string(target);
    return
end
if ~(isfolder(target) || endsWith(lower(target), {'.html' '.htm'}))
    msg = "only folders and review pages are opened: " + string(target);
    return
end
if dryrun, return, end
try
    winopen(target);
catch ME
    msg = "could not open " + string(target) + ": " + ME.message;
end
end
