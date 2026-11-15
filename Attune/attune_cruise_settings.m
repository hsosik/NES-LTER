function s = attune_cruise_settings(cruise, root)
% function s = attune_cruise_settings(cruise, root)
% Per-cruise processing settings for process_cruise: folder names, underway
% source and extra exclusion patterns from attune_cruise_settings.csv, with
% the optical configuration (OD2 placement, FCB beads or not, filter layout,
% channel map) from classifier_dev/attune_config.
%
% cruise - cruise name, e.g. 'EN727' (case-insensitive)
% root   - optional, folder that holds the <yyyymmdd>_<CRUISE> cruise folders
%          (default \\sosiknas1\Lab_data\Attune\cruise_data\)
%
% s.cruise, s.folder, s.basepath
% s.underway   - underway source: blank in the csv means the NES-LTER API
% s.exclude    - cellstr; FCS files whose names contain any of these are not
%                classified (defaults plus the csv's exclude_extra, ';'-separated)
% s.notes      - csv notes (printed by process_cruise)
% s.config     - attune_config(cruise)
%
% To add a cruise: add a row to attune_cruise_settings.csv; check that
% attune_config handles its OD2 placement and FCB beads.
%
% Nov 2026, processing wrapper redesign

if nargin < 2 || isempty(root)
    root = [char([92 92]) 'sosiknas1\Lab_data\Attune\cruise_data\'];
end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'classifier_dev'))

opts = detectImportOptions(fullfile(here, 'attune_cruise_settings.csv'), 'TextType', 'string', 'Delimiter', ',');
opts = setvartype(opts, 'string');
T = readtable(fullfile(here, 'attune_cruise_settings.csv'), opts);
T = fillmissing(T, 'constant', "");

cruise = upper(string(cruise));
i = find(upper(T.cruise) == cruise, 1);
s.cruise = char(cruise);
if isempty(i)
    % not in the table yet: find the folder by name
    d = dir(fullfile(root, ['*_' char(cruise)]));
    d = d([d.isdir]);
    if numel(d) ~= 1
        error('attune_cruise_settings:folder', 'no unique cruise folder %s*_%s; add a row to attune_cruise_settings.csv', root, cruise)
    end
    warning('attune_cruise_settings:notInTable', '%s is not in attune_cruise_settings.csv; using defaults', cruise)
    s.folder = d.name;
    row = struct('underway_source', "", 'exclude_extra', "", 'notes', "");
else
    s.folder = char(T.folder(i));
    row = table2struct(T(i, :));
end
s.basepath = fullfile(root, s.folder);
if ~exist(s.basepath, 'dir')
    error('attune_cruise_settings:folder', 'cruise folder not found: %s', s.basepath)
end

if strlength(row.underway_source) > 0
    s.underway = char(row.underway_source);
else
    s.underway = ['https://nes-lter-api.whoi.edu/api/underway/' lower(char(cruise)) '.csv'];
end

s.exclude = {'fcb_bead' 'FCB_bead' 'bead' 'Bead' 'Cast' 'cast' '(lab test)' 'Dockwater' 'discrete' ...
    'Rinses' 'Filter config' 'Cultures' 'test'};
if strlength(row.exclude_extra) > 0
    s.exclude = [s.exclude cellstr(strtrim(split(row.exclude_extra, ';')))'];
end
s.notes = char(row.notes);
s.config = attune_config(s.cruise);
end
