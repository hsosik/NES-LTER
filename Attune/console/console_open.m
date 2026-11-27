function msg = console_open(target, dryrun)
% function msg = console_open(target, dryrun)
% Cruise console 'open' buttons: opens an output folder in Explorer or a review page
% (index.html) in the default browser, or a run log in the default text viewer. msg: '' when opened, otherwise why not.
% dryrun true: only check (tests).
%
% Nov 2026, cruise console

if nargin < 2, dryrun = false; end
msg = "";
target = char(target);
if ~(isfolder(target) || isfile(target))
    msg = "not found: " + string(target);
    if ~isempty(regexp(target, '^[A-Za-z]:\\', 'once')) && ~isfolder(target(1:3))
        msg = msg + " (no drive " + string(target(1:2)) + " on this computer: map \\sosiknas1\Lab_data to it, or set the NAS drive letter at the top of the console)";
    end
    return
end
if ~(isfolder(target) || endsWith(lower(target), {'.html' '.htm' '.log'}))
    msg = "only folders, review pages and run logs are opened: " + string(target);
    return
end
if dryrun, return, end
try
    winopen(target);
catch ME
    msg = "could not open " + string(target) + ": " + ME.message;
end
end
