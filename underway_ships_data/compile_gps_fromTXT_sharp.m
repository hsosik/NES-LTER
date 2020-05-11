%special case for NOAA Ship Sharp GPS data from ship

%S1_1802
bpath1 = '\\sosiknas1\Lab_data\LTER\NESLTER_broadscale\20181102_S1_1802\'; %cruise
bpath2  = 'ship_provided\HRS1809JP\SMS\'; %CNAV files under cruise path

inpath = [bpath1 bpath2];
f = dir([inpath '*.txt']);
outpath = [bpath1 'gps\'];
if ~exist(outpath, 'dir')
    mkdir(outpath)
end
vartypes = {'datetime' 'cell' 'double' 'double' 'double' };
varnames = {'date', 'time' 'matdate' 'lat' 'lon'};
gps = table('size', [1 5],'VariableTypes', vartypes, 'VariableNames', varnames);
for cc = 1:length(f)
    disp(cc)
    %T = readtable([inpath f(cc).name], 'FileType', 'text', 'HeaderLines',3, 'Format', '%{MM/dd/yyyy}D%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
    T = readtable([inpath f(cc).name], 'FileType', 'text');
    temp = char(T{:,2});
    mdate = datenum(T{:,2}) + datenum(T{:,3});
    if iscell(T{:,3})
        mdate = datenum(T{:,2}) + datenum(datetime(T{:,3})) - floor(datenum(datetime));
    end
    lat = T{:,9};
    if iscell(lat)
        lat = str2num(cell2mat(lat));
    end
    latsign = ones(size(mdate)); lonsign = latsign;
    latsign(strmatch('S', T{:,8})) = -1;
    lat = latsign.*lat;
    lon = T{:,12};
    if iscell(lon)
        lon = str2num(cell2mat(lon));
    end
    lonsign(strmatch('W', T{:,11})) = -1;
    lon = lonsign.*lon;
    gps = [gps; table(T{:,2}, cellstr(char(T{:,3})), mdate, lat, lon, 'VariableNames', varnames)];
end
gps(1,:) = [];
ii = find(gps.lat == -99 & gps.lon == -99);
gps(ii,:) = [];

save([outpath 'gps_all'], 'gps')
disp('results saved:')
disp([outpath 'gps_all'])
