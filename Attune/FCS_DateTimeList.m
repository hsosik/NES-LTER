function [ FCSfileinfo ] = FCS_DateTimeList(fcs_path, FCSfileinfo_name)
%function [ FCSfileinfo ] = FCS_DateTimeList(fcs_path, FCSfileinfo_name)
%input:
%   fcs_path - path to a directory of fcs files
%   FCSfileinfo_name - optional: existing FCSfileinfo.mat (full path). Files
%       already in it are kept as they are and only new files are read. Rows
%       without detector voltages (tables made before Nov 2026) get them from
%       the FCS headers, so pro_measured can be recomputed.
%output: table FCSfileinfo, one row per FCS file, in time order
%   fcslist                  file name
%   matdate_start, matdate_stop, date_processed, vol_analyzed
%   trigger1_parameter, trigger1_threshhold, trigger1_hv  (same for trigger2;
%       '' / NaN when there is no second trigger). Thresholds are in -H units.
%   hv_FSC ... hv_GL4        detector voltages (NaN if the detector is absent)
%   layout                   filter layout: 'v2019' before 19 Aug 2019, else 'current'
%                            (see classifier_dev\attune_config)
%   file_type                from the name: bead, cast, discrete, dockwater, rinse,
%                            culture, stained, test, or sample
%   QC_flag, QC_flowrates, QC_dataintegrity
%       QC_flowrates: median and std of SSC-A/SSC-H for events with both > 200;
%       QC_flag 1 (good) when median < 1.5, std < 2 and event times never decrease
%   pro_measured, trig_SSC_eq220, trig_BL3_eq340  see attune_pro_measured
%
%Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2019
% modified 4/26/2021 to make output into a table, and to include quality control data
% Nov 2026: append mode fixed (it used a column that does not exist); QC channels
%   found by name; detector voltages, layout and file_type columns; pro_measured
%   from trigger settings scaled to reference voltages (attune_pro_measured)
%   instead of "SSC trigger below 500", which marked the Euk/Syn alternating runs
%   of EN688, EN706 and AR78 as Pro and missed AR104's Pro runs.

if ~exist('fcs_path', 'var')
    fcs_path = uigetdir(pwd, 'Pick a Directory of FCS files');
elseif ~exist(fcs_path, 'dir')
    disp('WARNING: Directory not found. Check input path.')
    FCSfileinfo = [];
    return
end

fcslist = dir(fullfile(fcs_path, '*.fcs'));
fcslist = {fcslist.name}';

old = table;
if exist('FCSfileinfo_name', 'var') && ~isempty(FCSfileinfo_name) && exist(FCSfileinfo_name, 'file')
    old = load(FCSfileinfo_name, 'FCSfileinfo');
    old = old.FCSfileinfo;
    if isstruct(old), old = struct2table(old); end
    if ismember('filelist', old.Properties.VariableNames) && ~ismember('fcslist', old.Properties.VariableNames)
        old = renamevars(old, 'filelist', 'fcslist');
    end
    old = old(ismember(old.fcslist, fcslist), :); % files no longer in the folder are dropped
    fcslist = setdiff(fcslist, old.fcslist);
end

b = numel(fcslist);
R = cell(b, 1);
for ii = 1:b
    if rem(ii, 20) == 0
        disp([num2str(ii) ' of ' num2str(b)])
    end
    R{ii} = read_one(fcs_path, fcslist{ii});
end
new = vertcat(R{:});
if isempty(new), new = read_one('', ''); new(1, :) = []; end

if ~isempty(old)
    [old, new] = align_columns(old, new);
    FCSfileinfo = [old; new(:, old.Properties.VariableNames)];
    e = cellfun(@isempty, FCSfileinfo.file_type);
    FCSfileinfo.file_type(e) = cellfun(@file_type, FCSfileinfo.fcslist(e), 'UniformOutput', false);
    e = cellfun(@isempty, FCSfileinfo.layout) & isfinite(FCSfileinfo.matdate_start);
    FCSfileinfo.layout(e) = {'current'};
    FCSfileinfo.layout(e & FCSfileinfo.matdate_start < datenum(2019, 8, 19)) = {'v2019'};
else
    FCSfileinfo = new;
end

%detector voltages and trigger settings for rows from older tables
FCSfileinfo = attune_pro_measured(FCSfileinfo, fcs_path);

%put table in chronological order
FCSfileinfo = sortrows(FCSfileinfo, 'matdate_start');
end


function r = read_one(fcs_path, name)
% one row of FCSfileinfo; NaN times for an empty file (as before)
det = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
r = table({name}, NaN, NaN, NaT, NaN, {''}, NaN, NaN, {''}, NaN, NaN, ...
    'VariableNames', {'fcslist' 'matdate_start' 'matdate_stop' 'date_processed' 'vol_analyzed' ...
    'trigger1_parameter' 'trigger1_threshhold' 'trigger1_hv' 'trigger2_parameter' 'trigger2_threshhold' 'trigger2_hv'});
for i = 1:numel(det), r.(['hv_' det{i}]) = NaN; end
r.layout = {''};
r.file_type = {file_type(name)};
r.QC_flag = NaN; r.QC_flowrates = [NaN NaN]; r.QC_dataintegrity = NaN;
if isempty(name), return, end

[fcsdat, fcshdr] = fca_readfcs(fullfile(fcs_path, name));
if fcshdr.TotalEvents == 0, return, end
pn = {fcshdr.par.name};
r.matdate_start = datenum([fcshdr.date ', ' fcshdr.starttime]);
r.matdate_stop = datenum([fcshdr.date ', ' fcshdr.stoptime]);
r.date_processed = datetime(fcshdr.date);
r.vol_analyzed = fcshdr.VOL;
r.trigger1_parameter = {fcshdr.tr1_par};
r.trigger1_threshhold = fcshdr.tr1_level;
r.trigger1_hv = det_hv(fcshdr, extractAfter(fcshdr.tr1_par, '_'));
if isfield(fcshdr, 'tr2_par') && ~isempty(fcshdr.tr2_par)
    r.trigger2_parameter = {fcshdr.tr2_par};
    r.trigger2_threshhold = fcshdr.tr2_level;
    r.trigger2_hv = det_hv(fcshdr, extractAfter(fcshdr.tr2_par, '_'));
end
for i = 1:numel(det), r.(['hv_' det{i}]) = det_hv(fcshdr, det{i}); end
if r.matdate_start < datenum(2019, 8, 19), r.layout = {'v2019'}; else, r.layout = {'current'}; end

%quality control flags
ia = find(strcmp(pn, 'SSC-A'), 1); ih = find(strcmp(pn, 'SSC-H'), 1);
if ~isempty(ia) && ~isempty(ih)
    t = fcsdat(:, ih) > 200 & fcsdat(:, ia) > 200;
    q = fcsdat(t, ia) ./ fcsdat(t, ih);
    r.QC_flowrates = [median(q) std(q)];
end
r.QC_dataintegrity = min(diff(fcsdat(:, 1)));
r.QC_flag = double(r.QC_flowrates(2) < 2 && r.QC_flowrates(1) < 1.5 && r.QC_dataintegrity >= 0);
end


function hv = det_hv(fcshdr, det)
% voltage of detector det (e.g. 'SSC') from its -H parameter
hv = NaN;
i = find(strcmp({fcshdr.par.name}, [char(det) '-H']), 1);
if ~isempty(i) && ~isempty(fcshdr.par(i).hv), hv = fcshdr.par(i).hv; end
end


function t = file_type(name)
n = lower(name);
types = {'bead' {'bead'}; 'cast' {'cast'}; 'discrete' {'discrete'}; 'dockwater' {'dockwater'}; ...
    'rinse' {'rinse'}; 'culture' {'culture'}; 'stained' {'sybr' 'hbac'}; 'test' {'test' 'filter config'}};
t = 'sample';
for i = 1:size(types, 1)
    if contains(n, types{i, 2}), t = types{i, 1}; return, end
end
end


function [a, b] = align_columns(a, b)
% give each table the other's missing columns (NaN, '' or NaT)
[a, b] = add_missing(a, b);
[b, a] = add_missing(b, a);
end

function [a, b] = add_missing(a, b)
for v = setdiff(b.Properties.VariableNames, a.Properties.VariableNames)
    x = b.(v{1});
    if iscell(x), f = repmat({''}, height(a), size(x, 2));
    elseif isdatetime(x), f = NaT(height(a), size(x, 2));
    else, f = NaN(height(a), size(x, 2));
    end
    a.(v{1}) = f;
end
end
