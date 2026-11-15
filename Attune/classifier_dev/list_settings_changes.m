function list_settings_changes(cruisedir)
% function list_settings_changes(cruisedir)
% Print each run of consecutive sample files with the same trigger settings,
% to line up step changes in class counts with instrument setting changes.
%
% Nov 2026, classifier redesign prototype

root = '\\sosiknas1\lab_data\Attune\cruise_data\';
load(fullfile(root, cruisedir, 'bead_calibrated', 'FCSfileinfo.mat'), 'FCSfileinfo')
F = sortrows(FCSfileinfo(~contains(FCSfileinfo.fcslist, {'bead' 'Bead'}), :), 'matdate_start');
key = strcat(F.trigger1_parameter, " ", string(F.trigger1_threshhold), " ", F.trigger2_parameter, " ", string(F.trigger2_threshhold));
starts = [1; find(key(2:end) ~= key(1:end-1)) + 1];
stops = [starts(2:end) - 1; height(F)];
for j = 1:numel(starts)
    fprintf('%-32s %s to %s (%d files)\n', key(starts(j)), datestr(F.matdate_start(starts(j)), 'dd-mmm HH:MM'), datestr(F.matdate_start(stops(j)), 'dd-mmm HH:MM'), stops(j) - starts(j) + 1);
end
end
