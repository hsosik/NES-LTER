function [IFCB_match_btl_results] =  IFCB_match_btl(IFCB_files, cast, niskin, bottle_data)
% function IFCB_match = IFCB_match_uw(IFCB_files, uw);
% input a list of IFCB filenames and a table of underway ship's info
% (mdate, lat, lon)
% output a table of closest match lat/lon for the IFCB files (interpolate
% if larger than 5 minute gap in the underway data)
% Heidi M. Sosik, Woods Hole Oceanographic Instition, Decemeber 2020
%
warning off

IFCB_mdate = IFCB_file2date(cellstr(IFCB_files));
iso8601format = 'yyyy-mm-dd hh:MM:ss';
btl_mdate = datenum(bottle_data.date, iso8601format);
t = bottle_data.Properties.VariableNames;
ilat = find(contains(t, 'latitude'));
ilon = find(contains(t, 'longitude'));
idepth = find(contains(t, 'dep'));
IFCB_match_btl_results = table('Size', [length(IFCB_files) size(bottle_data,2)],...
    'VariableTypes', varfun(@class,bottle_data,'OutputFormat','cell'), 'VariableNames', bottle_data.Properties.VariableNames);
IFCB_match_btl_results.lat(1) = NaN; 
IFCB_match_btl_results.lon(1) = NaN;
IFCB_match_btl_results.datetime(1) = {'null'}; 
numericVars = varfun(@isnumeric,IFCB_match_btl_results,'output','uniform');
numericVars = find(numericVars);
nanvec = num2cell(NaN(1,numel(numericVars)));

for count = 1:length(IFCB_files)
    ind = find(bottle_data.cast == cast(count) & bottle_data.niskin == niskin(count));
    if isempty(ind)
        %disp('No match up with bottle file info:')
        %disp([IFCB_files(count) cast(count) niskin(count)])
        %keyboard
        IFCB_match_btl_results(count,numericVars) = nanvec;
        IFCB_match_btl_results.cruise(count) = bottle_data.cruise(1);
        IFCB_match_btl_results.cast(count) = cast(count);
        IFCB_match_btl_results.niskin(count) = niskin(count);
        %IFCB_match_btl_results.datetime(count) = NaN;
    else
        IFCB_match_btl_results(count,1:end-3) = bottle_data(ind,:);
        IFCB_match_btl_results.lat(count) = bottle_data.(t{ilat})(ind);
        IFCB_match_btl_results.lon(count) = bottle_data.(t{ilon})(ind);
        IFCB_match_btl_results.datetime(count) = bottle_data.date(ind);
    end    
end

IFCB_match_btl_results = addvars(IFCB_match_btl_results, IFCB_files, 'before', 1, 'NewVariableNames',{'pid'});
IFCB_match_btl_results = addvars(IFCB_match_btl_results, IFCB_match_btl_results.(t{idepth}), 'NewVariableNames',{'depth'});

end
