function gps = compile_gps_fromRAW(cruise)
% function gps = compile_gps_fromRAW(cruise)
% function to compile GPS info for NOAA ships Bigelow, Okeanos Explorer, and
% Gunter
% See notes at bottom about source of raw files with GPS for each ship
%
% Heidi M. Sosik, Woods Hole Oceanographic Institution, May 2020
% used to compile metadata for underway IFCB sampling
%%
basepath = '\\sosiknas1\Lab_data\LTER\NESLTER_broadscale\';
temp = dir([basepath '*' cruise]);
bpath1 = [basepath temp.name filesep];

ship = cruise(1:2);
switch ship
    case 'HB'
        bpath2  = 'nodc_noaa_download\POSMV\'; %CNAV files under cruise path
        filestr = 'POSMV-INGGA*.Raw';
    case 'EX'
        bpath2  = 'NCEI_nav\'; %CNAV files under cruise path
        filestr = 'CNAV*.Raw';
    case 'GU'
        bpath2 = 'nodc_noaa_download\GPS_Data\'; 
        filestr = 'GP150-NMEA-GPGGA*.Raw';
end

inpath = [bpath1 bpath2];
f = dir([inpath filestr]);
outpath = [bpath1 'gps\'];
if ~exist(outpath, 'dir')
    mkdir(outpath)
end
vartypes = {'datetime' 'cell' 'double' 'double' 'double' };
varnames = {'date', 'time' 'matdate' 'lat' 'lon'};
gps = table('size', [1 5],'VariableTypes', vartypes, 'VariableNames', varnames);
for cc = 1:length(f)
    switch ship
        case 'EX'
            T = readtable([inpath f(cc).name], 'FileType', 'text', 'HeaderLines',3, 'Format', '%{MM/dd/yyyy}D%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
        otherwise
            T = readtable([inpath f(cc).name],'FileType', 'text', 'Format', '%{MM/dd/yyyy}D%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');    
    end                  
    temp = char(T{:,2});
    mdate = datenum(T{:,1}) + datenum(0,0,0,str2num(temp(:,1:2)),str2num(temp(:,4:5)),str2num(temp(:,7:end)));
    lat = char(regexprep(T{:,5}, char(0),'')); %POSMV files from HB have occasional odd spaces (char(0)) in lat/lon records
    latsign = ones(size(mdate)); lonsign = latsign;
    latsign(strmatch('S', regexprep(T{:,6}, char(0),''))) = -1;
    lat = latsign.*(str2num(lat(:,1:2)) + str2num(lat(:,3:end))/60);
    lon = char(regexprep(T{:,7}, char(0), '')); %POSMV files from HB have occasional odd spaces (char(0)) in lat/lon records
    lonsign(strmatch('W', regexprep(T{:,8}, char(0),''))) = -1;
    lon = lonsign.*(str2num(lon(:,1:3)) + str2num(lon(:,4:end))/60);
    gps = [gps; table(T{:,1}, T{:,2}, mdate, lat, lon, 'VariableNames', varnames)];
end
gps(1,:) = [];

save([outpath 'gps_all'], 'gps')
disp('results saved:')
disp([outpath 'gps_all'])

end


% if needed copy *RAW files from here (NCEI), example for HB1701 (you can paste this path in windows explorer):
%ftp://ftp.nodc.noaa.gov/nodc/archive/arc0110/0164796/1.1/data/0-data/HB-0_2017-07-20-143703/Posmv/
%start here and search for cruise abbreviation: https://www.nodc.noaa.gov/archivesearch/catalog/search/search.page
%if that doesn't work, search ship name (e.g., bigelow) and dates
% EX cruises
% 1. get CNAV* files from here (NCEI), example for EX1305 (you can paste this path in windows explorer):
% ftp://ftp.nodc.noaa.gov/nodc/archive/arc0062/0113335/1.1/data/0-data/vessel/Ship_Navigation_Data/NAV/
% copy to sosiknas1 (e.g., \\sosiknas1\Lab_data\LTER\NESLTER_broadscale\20130824_EX1305\NCEI_nav\)
% 
% HB cruises
% 1. Get POSMV-INGGA-RAW*.RAW files from NODC ftp site (https://www.nodc.noaa.gov/archivesearch/catalog/search/search.page)
% 	search for text = bigelow and dates 
% copy to Sosiknas1 (e.g., \\sosiknas1\Lab_data\LTER\NESLTER_broadscale\20170211_HB1701\nodc_noaa_download\POSMV\)
% 
% GU cruises
% 1. get GP150-NMEA-GPGGA*.Rar files from NODC ftp site
% copy to sosiknas1 (e.g., \\sosiknas1\Lab_data\LTER\NESLTER_broadscale\20151013_GU1506\nodc_noaa_download\GPS_Data\)
