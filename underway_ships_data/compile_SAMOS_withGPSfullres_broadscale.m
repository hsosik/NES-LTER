function uw = compile_SAMOS_withGPSfullres_broadscale(cruise)
% function uw = compile_SAMOS_withGPSfullres_broadscale(cruise)
% function to compile NOAA ship's underway data from SAMOS netcdf files
% including ship specific collation of high resolution GPS data, if needed
% 
% compilation needed to prepare IFCB metadata records
%
%Heidi M. Sosik, Woods Hole Oceanographic Institution, May 2020

basepath = '\\sosiknas1\Lab_data\LTER\NESLTER_broadscale\';
temp = dir([basepath '*' cruise]);
bpath1 = [basepath temp.name filesep]; clear temp
ncpath = [bpath1 '\samos_netcdf\netcdf\'];
outfile = [bpath1 '\compiled_underway\'];
if ~exist(outfile, 'dir')
    mkdir(outfile)
end
outfile = [outfile cruise 'uw_compiled'];

if ~strncmp('PC',cruise,2)
    %load the gps table compiled from ship's data with compile_gps_from*.m 
    load([bpath1 'gps\gps_all.mat']); 
end
%%
ncfiles = dir([ncpath '*.nc']);
ncfiles = {ncfiles.name}';
%keep only the latest version of each day
temp = char(ncfiles);
v = str2num(temp(:,15:19));
temp = cellstr(temp(:,6:13));
temp2 = unique(temp);
keep = ones(size(temp));
for cc = 1:length(temp2)
    ii = strmatch(temp2(cc), temp);
    keep(ii(v(ii)<max(v(ii)))) = 0;
end
info = ncinfo([ncpath ncfiles{1}]);
varname = {info.Variables.Name};
varname = setdiff(varname,{'flag', 'history'});
uw = table;
for count = 1:length(ncfiles)
    T = table;
    for vcount = 1:length(varname)
        T.(varname{vcount})  = ncread([ncpath ncfiles{count}],varname{vcount});
        T.Properties.VariableUnits{varname{vcount}} = ncreadatt([ncpath ncfiles{count}],varname{vcount}, 'units');
        T.Properties.VariableDescriptions{varname{vcount}} = ncreadatt([ncpath ncfiles{count}],varname{vcount}, 'long_name');   
    end
    uw = [uw; T];
end
uw.matdate = datenum(datenum('1-1-1980') + double(uw.time)/60/24);

if exist('gps', 'var')
    ind = NaN(size(uw.matdate));
    mdate_gps = datenum(gps.matdate);
    for count = 1:length(uw.matdate)
        [dd,tt] = min(abs(uw.matdate(count)-mdate_gps));
        if dd < 2/60/24 %2 minutes as days
            ind(count) = tt;
        end
    end

    temp = NaN(size(ind));
    uw.latitude_fullres = temp;
    uw.longitude_fullres = temp;
    uw.mdate_fullres = temp;
    nind = find(~isnan(ind));
    uw.latitude_fullres(nind) = gps.lat(ind(nind));
    uw.longitude_fullres(nind) = gps.lon(ind(nind));
    uw.mdate_fullres(nind) = mdate_gps(ind(nind));
    uw = movevars(uw,{'latitude_fullres', 'longitude_fullres'},'Before',1);
    uw.Properties.VariableNames{'lat'} = 'lat_SAMOS';
    uw.Properties.VariableNames{'lon'} = 'lon_SAMOS';
else
    uw = movevars(uw,{'lat', 'lon'},'Before',1); 
    uw.Properties.VariableNames{'lat'} = 'latitude_fullres';
    uw.Properties.VariableNames{'lon'} = 'longitude_fullres';
    uw.Properties.VariableNames{'matdate'} = 'mdate_fullres';
end
notes = {'Heidi Sosik, WHOI, produced with compile_SAMOS_withGPSfullres_broadscale.m from downloaded SAMOS netcdf files and appended higher resolution lat, lon from raw data provided by ship or NCEI download'};
save(outfile, 'uw', 'notes')
disp('results saved: ') 
disp(outfile)

end