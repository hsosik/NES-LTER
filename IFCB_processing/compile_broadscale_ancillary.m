ubase = '\\sosiknas1\IFCB_data\NESLTER_broadscale\match_up\';
filelist = dir([ubase 'NESLTER*uw_match.mat']);
%ubase = '\\sosiknas1\IFCB_data\URI_ecomon\match_up\';
%filelist = dir([ubase 'URI*uw_match.mat']);
filelist = {filelist.name}';
cruises = split(filelist,'_');
cruises = cruises(:,3);
%cruises = cruises(3);
clear temp
match_uw = table;
slist = {'pid' 'cruise' 'lat' 'lon' 'temperature' 'salinity' 'mdate'};

for count1 = 1:length(cruises)
     disp(cruises{count1})
      if strmatch('GU1905', cruises{count1})
          orig_var = {'TSG-Ext-SBE38-Temp' 'TSG-Salinity'};
      elseif strncmp('AR', cruises{count1},2)
          orig_var = {'sbe48t' 'sbe45s'};
      else
          orig_var = {'T' 'SSPS'};
      end
%      elseif strmatch('EN', cruises{count1}(1:2))
%          orig_var = {'tsg1_temperature' 'tsg1_salinity'};
%      else
%          orig_var = {'sbe48t' 'sbe45s'};
%      end
     u = load([ubase filelist{count1}]);
     if ~ismember(orig_var(2), u.IFCB_match_uw_results.Properties.VariableNames) %no salinity data
         u.IFCB_match_uw_results.(orig_var{2})(:) = NaN;
     end
     u.IFCB_match_uw_results.cruise(:) = cruises(count1);
     u.IFCB_match_uw_results.Properties.VariableNames(orig_var) = {'temperature' 'salinity'};
     [~,ia] = ismember(slist, u.IFCB_match_uw_results.Properties.VariableNames);
     match_uw = [match_uw; u.IFCB_match_uw_results(:,ia)];
end

clear u ia ind count1 n slist* orig_var 
save([ubase 'compiledTS_tables'], 'match*', 'cruises')
disp('results saved:')
disp([ubase 'compiledTS_tables'])