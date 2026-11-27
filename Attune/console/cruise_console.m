function fig = cruise_console(opt)
% function fig = cruise_console(Name=Value)
% Cruise console: processing state of every cruise (underway and preserved) in one
% window. The page (cruise_console_page.html) runs inside a MATLAB uihtml component;
% status comes from the cruise_status_all cache. In the window:
%   - overview grid: one row per cruise and sample type, a coloured state per step
%     (files, beads, bead review, classify, class review, volume, table, metadata);
%     hover a cell for the reason; filters for sample type and state
%   - details of the selected cruise: per step state, reason, output date and the
%     latest step record (process_state.mat), counts, archive
%   - buttons: refresh this cruise or all cruises (cruise_status_all, updates the
%     cache and the browser copy), open the output folder and the bead, class and
%     SYBR review pages
% The same page, read-only, is written next to the cache for browsers
% (<Root>\processing_status\index.html, console_publish).
% Launching process_cruise runs from the console: phase 3 (not built yet).
%
% Options
%   Refresh  true: scan all cruises before opening (default false: the cache; scans
%            when there is no cache)
%   Root     folder of cruise folders (default \\sosiknas1\Lab_data\Attune\cruise_data\)
%
% Nov 2026, cruise console (phase 2)

arguments
    opt.Refresh (1,1) logical = false
    opt.Root (1,1) string = ""
end
here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here))
root = opt.Root;
if root == "", root = "\\sosiknas1\Lab_data\Attune\cruise_data\"; end
cf = fullfile(root, 'processing_status', 'status.mat');
if opt.Refresh || ~isfile(cf)
    [~, S] = cruise_status_all(Root=opt.Root, Show=false);
end
o = load(cf, 'S', 'made');
st = struct('S', o.S, 'made', o.made, 'root', opt.Root);

fig = uifigure('Name', 'Attune cruise console', 'Position', [80 80 1500 860]);
g = uigridlayout(fig, [1 1], 'Padding', [0 0 0 0]);
h = uihtml(g, 'HTMLSource', fullfile(here, 'cruise_console_page.html'));
h.HTMLEventReceivedFcn = @(src, evt) console_event(src, evt);
h.UserData = st;
h.Data = status_json(st.S, st.made, 'matlab');
end


function console_event(h, evt)
% events from the page: 'refresh' {cruise, sample} ('' = all), 'open' path
switch evt.HTMLEventName
    case 'open'
        msg = console_open(char(evt.HTMLEventData));
        if msg ~= "", st = h.UserData; h.Data = status_json(st.S, st.made, 'matlab', msg); end
    case 'refresh'
        d = evt.HTMLEventData;
        st = h.UserData;
        cr = string(d.cruise); smp = string(d.sample);
        if cr == ""
            what = "all cruises (about 2 minutes)";
        else
            what = cr + " " + smp;
        end
        h.Data = status_json(st.S, st.made, 'matlab', "refreshing " + what + " ...");
        drawnow
        try
            st = console_refresh(st, cr, smp);
            h.UserData = st;
            h.Data = status_json(st.S, st.made, 'matlab');
        catch ME
            h.Data = status_json(st.S, st.made, 'matlab', "refresh failed: " + ME.message);
        end
end
end
