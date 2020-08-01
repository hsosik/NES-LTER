cruises = {'AR22' 'AR24A' 'AR24B' 'AR24C' 'EN608' 'AR28A' 'AR28B' 'EN617'...
    'AR31A' 'AR31B' 'AR31C' 'AR32' 'EN627' 'AR34A' 'AR34B' 'EN644'};  %'AR16' 'AR34B' 'AR38' 'AR39A' 'AR39A'

%cruises = {'EN608' 'AR28B' 'EN617' 'AR31B'};  %'AR16' 'AR34B' 'AR38' 'AR39'
cruises = {'EN608' 'EN617' 'EN627' 'EN644'};
cruises = {'EN617' 'EN644' 'EN649'};

ubase = '\\sosiknas1\IFCB_data\NESLTER_transect\match_up\';

match_uw = table;
nut = table;
match_cast = table;
slist = {'pid' 'cruise' 'lat' 'lon' 'temperature' 'salinity' 'mdate'};
slist_cast = {'pid' 'cruise' 'lat' 'lon' 't090c' 'sal00' 'depth' 'mdate'};

 for count1 = 1:length(cruises)
     disp(cruises{count1})
     if strmatch('EN627', cruises{count1})
         orig_var = {'tsg2_temperature' 'tsg2_salinity'};
     elseif strmatch('EN', cruises{count1}(1:2))
         orig_var = {'tsg1_temperature' 'tsg1_salinity'};
     else
         orig_var = {'sbe48t' 'sbe45s'};
     end
     u = load([ubase 'NESLTER_transect_' cruises{count1} '_uw_match.mat']);
     u.IFCB_match_uw_results.Properties.VariableNames(orig_var) = {'temperature' 'salinity'};
     [~,ia] = ismember(slist, u.IFCB_match_uw_results.Properties.VariableNames);
     match_uw = [match_uw; u.IFCB_match_uw_results(:,ia)];
     if exist([ubase 'NESLTER_transect_' cruises{count1} '_cast_match.mat'], 'file')
        u = load([ubase 'NESLTER_transect_' cruises{count1} '_cast_match.mat']);
        u.IFCB_match_btl_results.mdate = datenum(u.IFCB_match_btl_results.datetime, 'yyyy-mm-dd hh:MM:SS+00:00');
        [~,ia] = ismember(slist_cast, u.IFCB_match_btl_results.Properties.VariableNames);
        match_cast = [match_cast; u.IFCB_match_btl_results(:,ia)];
     end
     n = webread(['https://nes-lter-data.whoi.edu/api/nut/' cruises{count1} '.csv']);
     n.alternate_sample_id = []; %move this column since the type doesn't match between all cruises
     nut = [nut; n];
 end
nut.mdate = datenum(nut.date, 'yyyy-mm-dd hh:MM:ss+00:00');
% FUDGE for bad point for now
ind = find(strcmp(match_uw.cruise, 'EN627') & match_uw.salinity > 34);
match_uw.temperature(ind) = NaN;
match_uw.salinity(ind) = NaN;

clear u ubase ia ind count1 n slist* orig_var 
