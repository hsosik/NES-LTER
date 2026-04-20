% function to read discrete bottled data from the API
% although it is unnecessary to save the nutrient data file to the local
% computer, I did it still... just in case there might be occassions where internet is unavailable 

% Bofu Zheng
% zhengbofuzju@gmail.com
% USE EN715 as an example

%% 
clear all
cruise_name = 'EN715';
cruise_name_low = 'en715';


%% load botttled nitrate data
options = weboptions('ContentType', 'table', 'Timeout', 30);
nutrient = webread(['https://nes-lter-data.whoi.edu/api/nut/',cruise_name,'.csv'], options);
% nut_varb_name = nutrient.Properties.VariableNames;
% sta_ind = find(strcmp('nearest_station', nut_varb_name));
% sta_name = nutrient(:,sta_ind);
cruisename = table2array(nutrient(1,1));
cruisename = cruisename{1};
CTD_bottle.cruise_name = cruisename;
CTD.cast = table2array(nutrient(:,2)); 
CTD_bottle.niskin = table2array(nutrient(:,3)); 
CTD_bottle.time = datenum(table2array(nutrient(:,4)),'yyyy-mm-dd HH:MM:SS');
CTD_bottle.lat  = table2array(nutrient(:,5)); 
CTD_bottle.lon  = table2array(nutrient(:,6));
CTD_bottle.d    = table2array(nutrient(:,7)); 
CTD_bottle.n    = table2array(nutrient(:,10)); 

save(['/Users/warrbob/Desktop/WHOI/research/sunaQC/',cruise_name,'/SUNA/',cruise_name,'_bottle.mat'],'CTD_bottle')



