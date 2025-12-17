function [IFCB_match] =  IFCB_match_uw(IFCB_files, IFCB_mdate, uw)
% function IFCB_match = IFCB_match_uw(IFCB_files, uw);
% input a list of IFCB filenames and a table of underway ship's info
% (mdate, lat, lon)
% output a table of closest match lat/lon for the IFCB files (interpolate
% if larger than 5 minute gap in the underway data)
% Heidi M. Sosik, Woods Hole Oceanographic Instition, Decemeber 2020
%
warning off

%IFCB_mdate = IFCB_file2date(cellstr(IFCB_files));
iso8601format = 'yyyy-mm-dd hh:MM:ss';
if ismember('mdate_fullres', uw.Properties.VariableNames)
    uw_mdate = uw.mdate_fullres; %case for NESLTER_broadscale with added full resolution gps
elseif ismember('matdate', uw.Properties.VariableNames) %case for NESLTER_broadscale only SAMOS
    uw_mdate = uw.matdate; %case for NESLTER_broadscale
else
    %uw_mdate = datenum(uw.date, iso8601format);
    uw_mdate = datenum(strcat(char(uw.date_gmt), char(uw.time_gmt)),'yyyy/mm/ddHH:MM:ss.FFF'); %api2
end
t = uw.Properties.VariableNames;
%ilat = find(contains(t, 'latitude'));
%ilon = find(contains(t, 'longitude'));
ilat = find(contains(t, 'lat'));
ilon = find(contains(t, 'lon'));
uw_lat = uw.(t{ilat(1)});
uw_lon = uw.(t{ilon(1)});
IFCB_match(1,:) = uw(1,:);
IFCB_match.lat(1:length(IFCB_mdate)) = NaN(size(IFCB_mdate)); 
IFCB_match.lon(1:length(IFCB_mdate)) = NaN(size(IFCB_mdate));
nnind = find(~isnan(uw_mdate));

for count = 1:length(IFCB_mdate)
    [m,ia] = min(abs(IFCB_mdate(count)-uw_mdate));
    if m < 5/60/24 %5 minutes as days
            IFCB_match(count,1:end-2) = uw(ia,:);
        IFCB_match.lat(count) = uw_lat(ia);
        IFCB_match.lon(count) = uw_lon(ia);        
    else
        if ia < length(uw_mdate) && ia > 1 && m<2/24 %otherwise no match (IFCB file after end of uw data or before beginning) OR gap longer than 2 h
            if IFCB_mdate(count) > uw_mdate(ia) %closest to end of gap
                it = ia;
            else %closest to start of gap
                it = ia-1;
            end
            step = floor((uw_mdate(it+1)- uw_mdate(it))*24*60); %one minute interpolation
            [lat,lon] = track2(uw_lat(it), uw_lon(it),uw_lat(it+1), uw_lon(it+1),[],[], step);
            it2 = round((IFCB_mdate(count)-uw_mdate(it))/(uw_mdate(it+1)-uw_mdate(it))*step); %index of closest interpolated minute
            IFCB_match.lat(count) = lat(it2);
            IFCB_match.lon(count) = lon(it2);
            IFCB_match(count,2:end-2) = array2table(interp1(uw_mdate(nnind), double(uw{nnind,2:end}), IFCB_mdate(count)));
            if ismember('latitude_fullres', IFCB_match.Properties.VariableNames)
                IFCB_match.latitude_fullres(count) = IFCB_match.lat(count);
                IFCB_match.longitude_fullres(count) = IFCB_match.lon(count);
            end
            disp(['CHECK interpolation file: ' char(IFCB_files(count))])
%             keyboard
        else
            %IFCB_match(count,1:size(uw,2)) = array2table(NaN(size(uw(1,:))));
            %IFCB_match(count,1:size(uw,2)-1) = array2table(NaN(size(uw(1,1:end-1)))); %skip date cell on end??
            IFCB_match{count,strcmp(IFCB_match.Properties.VariableTypes, 'double')} = NaN;
        end
    end
end
IFCB_match = addvars(IFCB_match, IFCB_files, 'before', 1, 'NewVariableNames',{'pid'});
IFCB_match = addvars(IFCB_match, IFCB_mdate, 'NewVariableNames',{'mdate'});

end
