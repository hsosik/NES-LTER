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
uw_mdate = datenum(uw.date, iso8601format);
t = uw.Properties.VariableNames;
ilat = find(contains(t, 'latitude'));
ilon = find(contains(t, 'longitude'));
uw_lat = uw.(t{ilat(1)});
uw_lon = uw.(t{ilon(1)});
IFCB_match = uw(1,:);
IFCB_match.lat(1) = NaN; 
IFCB_match.lon(1) = NaN;

for count = 1:length(IFCB_mdate)
    [m,ia] = min(abs(IFCB_mdate(count)-uw_mdate));
    if m < 5/60/24 %5 minutes as days
        IFCB_match(count,1:end-2) = uw(ia,:);
        IFCB_match.lat(count) = uw_lat(ia);
        IFCB_match.lon(count) = uw_lon(ia);        
    else
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
        IFCB_match(count,1:end-2) = interp1(uw_mdate, uw, IFCB_mdate(count));
        disp('CHECK interpolation')
        keyboard
    end
end
IFCB_match = addvars(IFCB_match, IFCB_files, 'before', 1, 'NewVariableNames',{'pid'});
IFCB_match = addvars(IFCB_match, IFCB_mdate, 'NewVariableNames',{'mdate'});

end
