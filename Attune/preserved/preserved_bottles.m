function B = preserved_bottles(s)
% function B = preserved_bottles(s)
% CTD bottle data for a cruise in one form, for matching preserved samples by
% cast and niskin. s = attune_cruise_settings(cruise); s.bottle_format picks
% the reader and s.bottle_source the file or URL:
%   'api'     NES-LTER API ctd/bottles (default), with nearest_station from
%             ctd/metadata
%   'exports' EXPORTS bottle files (SR1812, DY131, SG2105)
%   'tn368'   SPIROPA TN368 bottle table (.mat) plus niskin depths from Olga's
%             tn368_niskin_pressure_depth.txt
%   'hb1907'  OTZ HB1907 bottle csv (its temperature is in situ t090c, not
%             potential temperature)
%   'ar43'    OTZ AR43 bottle csv
%
% B - table: cast, niskin, date_sampled (datetime, UTC), latitude, longitude,
%     depth_m, potemp090c, salinity, nearest_station, r2r_event
%
% Nov 2026, process_cruise preserved mode (readers from step 3 of
% Discrete_Samples\Process_Preserved_AllSteps_2026)

src = s.bottle_source;
switch lower(s.bottle_format)
    case 'api'
        B = read_csv(src);
        try
            M = read_csv(regexprep(src, 'bottles', 'metadata'));
            if ismember('nearest_station', M.Properties.VariableNames)
                B = outerjoin(B, M(:, {'cast' 'nearest_station'}), 'Keys', 'cast', 'MergeKeys', true, 'Type', 'left');
            end
        catch ME
            warning('preserved_bottles:metadata', 'no cast metadata: %s', ME.message)
        end
        B = ren(B, {'date' 'depsm' 'sal00'}, {'date_sampled' 'depth_m' 'salinity'});
    case 'exports'
        B = read_csv(src);
        B = ren(B, {'Cast' 'BottleNo' 'Lat' 'Lon' 'depth' 'potTemp1' 'sal1' 'sdate'}, ...
            {'cast' 'niskin' 'latitude' 'longitude' 'depth_m' 'potemp090c' 'salinity' 'date_sampled'});
    case 'tn368'
        S = load(src);
        f = fieldnames(S); BTL = S.(f{1});
        if isfield(S, 'BTL'), BTL = S.BTL; end
        dfile = [char([92 92]) 'sosiknas1\Lab_data\SPIROPA\20190705_TN368\fromOlga\tn368_niskin_pressure_depth.txt'];
        t = importdata(dfile);
        bd = array2table(t.data, 'VariableNames', matlab.lang.makeValidName(t.colheaders));
        B = outerjoin(bd, BTL, 'LeftKeys', {'x_Cast' 'TargetDepth_m_'}, 'RightKeys', {'Cast' 'TargetDepth_m'}, ...
            'RightVariables', {'datetime' 'Latitude_decimalDeg' 'Longitude_decimalDeg' 'Potemp090C' 'Sal00'}, 'Type', 'left');
        B = unique(B, 'rows'); % BTL repeats rows for productivity entries
        B = ren(B, {'x_Cast' 'Bottle' 'Depth_m_' 'Potemp090C' 'Sal00' 'Latitude_decimalDeg' 'Longitude_decimalDeg' 'datetime'}, ...
            {'cast' 'niskin' 'depth_m' 'potemp090c' 'salinity' 'latitude' 'longitude' 'date_sampled'});
    case 'hb1907'
        B = read_csv(src);
        B = ren(B, {'depth' 'date' 't090c'}, {'depth_m' 'date_sampled' 'potemp090c'});
    case 'ar43'
        B = read_csv(src);
        B = ren(B, {'depth' 'date' 'sal00'}, {'depth_m' 'date_sampled' 'salinity'});
    otherwise
        error('preserved_bottles:format', 'unknown bottle_format "%s"', s.bottle_format)
end

% Seabird names used by several formats (HB1907 and AR43 too)
B = ren(B, {'depsm' 'sal00' 'date'}, {'depth_m' 'salinity' 'date_sampled'});

% one standard set of columns
num = {'cast' 'niskin' 'latitude' 'longitude' 'depth_m' 'potemp090c' 'salinity'};
for v = num
    if ~ismember(v{1}, B.Properties.VariableNames), B.(v{1}) = NaN(height(B), 1); end
    if ~isnumeric(B.(v{1})), B.(v{1}) = str2double(string(B.(v{1}))); end
end
for v = {'nearest_station' 'r2r_event'}
    if ~ismember(v{1}, B.Properties.VariableNames), B.(v{1}) = strings(height(B), 1); end
    B.(v{1}) = string(B.(v{1})); B.(v{1})(ismissing(B.(v{1}))) = "";
end
if ~ismember('date_sampled', B.Properties.VariableNames), B.date_sampled = NaT(height(B), 1); end
B.date_sampled = to_datetime(B.date_sampled);
B = B(:, [num(1:2) {'date_sampled'} num(3:end) {'nearest_station' 'r2r_event'}]);
B = B(isfinite(B.cast) & isfinite(B.niskin), :);
[~, i] = unique(B(:, {'cast' 'niskin'}), 'rows', 'stable');
B = B(i, :);
end


function T = read_csv(src)
if startsWith(src, 'http')
    T = webread(src, weboptions('Timeout', 60, 'ContentType', 'table'));
else
    T = readtable(src, 'Delimiter', ',', 'VariableNamingRule', 'modify');
end
end


function T = ren(T, from, to)
for i = 1:numel(from)
    if ismember(from{i}, T.Properties.VariableNames) && ~ismember(to{i}, T.Properties.VariableNames)
        T = renamevars(T, from{i}, to{i});
    end
end
end


function d = to_datetime(x)
% datetime (UTC, no time zone) from the date formats in the bottle files
if isdatetime(x)
    d = x; if ~isempty(d.TimeZone), d.TimeZone = 'UTC'; d.TimeZone = ''; end
    return
end
if isnumeric(x), d = datetime(x, 'ConvertFrom', 'datenum'); return, end
x = regexprep(strtrim(string(x)), '(\+00:?00|Z)$', '');
x = replace(x, 'T', ' ');
d = NaT(size(x));
for f = ["yyyy-MM-dd HH:mm:ss.SSS" "yyyy-MM-dd HH:mm:ss" "yyyy-MM-dd HH:mm" "dd/MM/yyyy HH:mm" "dd/MM/yyyy HH:mm:ss" "dd-MMM-yyyy HH:mm:ss"]
    todo = isnat(d) & strlength(x) > 0;
    if ~any(todo), break, end
    try
        d(todo) = datetime(x(todo), 'InputFormat', f);
    catch
    end
end
end
