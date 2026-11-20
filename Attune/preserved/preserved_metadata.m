function M = preserved_metadata(p, s)
% function M = preserved_metadata(p, s)
% Sampling metadata for each preserved sample of a cruise: process_cruise
% metadata step with Sample="preserved".
%
% CTD samples (cast > 0) are matched by cast and niskin to the bottle data
% (preserved_bottles: lat, lon, date sampled, depth, potential temperature,
% salinity, nearest station). Samples without a cast in the name (underway or
% bucket discretes, cast 0) are matched to the event log (s.elog_source):
%   - EXPORTS cruises (bottle_format 'exports'): by the event number in the
%     file name (SR1812 _nnnnn_UW, DY131 Exnnnnn); depth 6 m, cast/niskin -9999
%   - otherwise: the one event whose comment mentions FCM or bucket (EN657,
%     AR38); depth 5 m; temperature and salinity from the underway data
%     nearest in time
%
% M - table, one row per sample (cast, niskin, sample_label, sizefrac) with the
%     metadata; saved with the unmatched samples in <outpath>\BottleMetadata.mat.
%     The table step joins it to the per-Niskin table.
%
% Nov 2026, process_cruise preserved mode (from step 3 of
% Discrete_Samples\Process_Preserved_AllSteps_2026)

load([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo')
F = FCSfileinfo(~contains(FCSfileinfo.fcslist, p.filetype2exclude) & ~strcmp(FCSfileinfo.runtype, 'other'), :);
[~, i] = unique(F(:, {'cast' 'niskin' 'sample_label' 'sizefrac'}), 'rows', 'stable');
M = F(i, {'cast' 'niskin' 'sample_label' 'sizefrac' 'fcslist'});
M = renamevars(M, 'fcslist', 'example_file');

B = preserved_bottles(s);
M = outerjoin(M, B, 'Keys', {'cast' 'niskin'}, 'MergeKeys', true, 'Type', 'left');
M.source = repmat("bottle", height(M), 1);
M.source(isnat(M.date_sampled) & isnan(M.latitude)) = "";

uw = find(M.cast == 0);
if ~isempty(uw)
    try
        E = read_table(s.elog_source);
    catch ME
        warning('preserved_metadata:elog', 'no event log (%s): %s', s.elog_source, ME.message)
        E = table;
    end
    if ~isempty(E)
        if strcmpi(s.bottle_format, 'exports')
            M = match_exports(M, uw, E, s.cruise);
        else
            M = match_fcm_event(M, uw, E, s);
        end
    end
end

unmatched = M(M.source == "", {'cast' 'niskin' 'sample_label' 'example_file'});
save([p.outpath 'BottleMetadata.mat'], 'M', 'unmatched')
fprintf('metadata: %d samples, %d matched to bottles, %d to the event log, %d unmatched -> %sBottleMetadata.mat\n', ...
    height(M), sum(M.source == "bottle"), sum(M.source == "event log"), height(unmatched), p.outpath)
if ~isempty(unmatched)
    fprintf('  unmatched: %s\n', strjoin(string(unmatched.sample_label), ', '))
end
end


function M = match_exports(M, uw, E, cruise)
% event number in the file name -> RoeslerEvent in the event log
name = string(M.example_file(uw));
if startsWith(upper(cruise), 'SR')
    ev = "EX" + regexprep(string(regexp(name, '_\d{5}_UW', 'match', 'once')), '[_UW]', '');
else
    ev = upper(string(regexp(name, 'Ex\d{5}', 'match', 'once')));
end
[ok, j] = ismember(ev, upper(string(E.RoeslerEvent)));
k = uw(ok); j = j(ok);
M.date_sampled(k) = to_dt(E.dateTime8601(j));
M.latitude(k) = E.Latitude(j); M.longitude(k) = E.Longitude(j);
M.r2r_event(k) = string(E.R2R_Event(j));
M.depth_m(k) = 6; M.cast(k) = -9999; M.niskin(k) = -9999;
M.source(k) = "event log";
end


function M = match_fcm_event(M, uw, E, s)
% the single event logged as an FCM or bucket sample
if ~ismember('Comment', E.Properties.VariableNames), return, end
c = string(E.Comment); c(ismissing(c)) = "";
j = find(contains(c, 'FCM') | contains(c, 'bucket', 'IgnoreCase', true));
if numel(j) ~= 1
    warning('preserved_metadata:event', '%d events mention FCM or bucket; underway discretes not matched', numel(j))
    return
end
t = to_dt(E.dateTime8601(j));
M.date_sampled(uw) = t; M.latitude(uw) = E.Latitude(j); M.longitude(uw) = E.Longitude(j);
M.depth_m(uw) = 5; M.source(uw) = "event log";
try
    U = read_table(s.underway);
    ut = to_dt(U.date);
    [~, m] = min(abs(ut - t));
    tv = intersect({'tsg1_temperature' 'sbe48t'}, U.Properties.VariableNames, 'stable');
    sv = intersect({'tsg1_salinity' 'sbe45s'}, U.Properties.VariableNames, 'stable');
    if ~isempty(tv), M.potemp090c(uw) = U.(tv{1})(m); end
    if ~isempty(sv), M.salinity(uw) = U.(sv{1})(m); end
catch ME
    warning('preserved_metadata:underway', 'no underway T/S for the underway discretes: %s', ME.message)
end
end


function T = read_table(src)
if startsWith(src, 'http')
    T = webread(src, weboptions('Timeout', 60, 'ContentType', 'table'));
else
    T = readtable(src, 'Delimiter', ',');
end
end


function d = to_dt(x)
if isdatetime(x), d = x; if ~isempty(d.TimeZone), d.TimeZone = 'UTC'; d.TimeZone = ''; end, return, end
x = regexprep(strtrim(string(x)), '(\+00:?00|Z)$', '');
d = datetime(replace(x, 'T', ' '), 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
end
