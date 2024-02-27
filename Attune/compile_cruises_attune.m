basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\';
filebase = '\bead_calibrated\AttuneTable_uw_match.mat';

yr = 2018;
for yr = 2018:2020
cruise = dir([basepath num2str(yr) '*']);

yrTable = table;
for ii = 1:length(cruise)
  tempTable = table;
  f = [basepath cruise(ii).name filebase];
  if exist(f, 'file')
    temp = load(f);
  if ~isempty(strfind(cruise(ii).name, 'EN'))
    tempTable.lat = temp.Attune_uw_match.gps_furuno_latitude;
    tempTable.lon = temp.Attune_uw_match.gps_furuno_longitude;
    tempTable.temperature = temp.Attune_uw_match.tsg1_temperature;
    tempTable.salinity = temp.Attune_uw_match.tsg1_salinity;
  elseif ~isempty(strfind(cruise(ii).name, 'AR'))
    tempTable.lat = temp.Attune_uw_match.dec_lat;
    tempTable.lon = temp.Attune_uw_match.dec_lon;
    tempTable.temperature = temp.Attune_uw_match.sbe48t;
    tempTable.salinity =  temp.Attune_uw_match.sbe45s;
  elseif ~isempty(strfind(cruise(ii).name, 'TN'))
    tempTable.lat = temp.Attune_uw_match.lat_tsg;
    tempTable.lon = temp.Attune_uw_match.lon_tsg;  
    tempTable.temperature = temp.Attune_uw_match.t1;
    tempTable.salinity = temp.Attune_uw_match.s;
  end
  tempTable.mdate = temp.Attune_uw_match.StartDate;
  if strmatch('Syn_count', temp.Attune_uw_match.Properties.VariableNames)
    tempTable.Syn_count = temp.Attune_uw_match.Syn_count;
    tempTable.Syn_carbon = temp.Attune_uw_match.Syn_carbon;  
  else
    tempTable.Syn_count = temp.Attune_uw_match.(" Syn_count ");
    tempTable.Syn_carbon = temp.Attune_uw_match.(" Syn_carbon ");
  end
  tempTable.VolAnalyzed_ml = temp.Attune_uw_match.VolAnalyzed_ml;
  yrTable = [yrTable; tempTable(temp.good,:)];
  else
      disp([cruise(ii).name ' missing'])
  end
end
eval(['yrTable' num2str(yr) ' = yrTable;'])
end

%%

yrTable = yrTable2018;

lat_smooth = round(yrTable.lat,1);
ilat = unique(lat_smooth);
ilat = ilat(~isnan(ilat));
dv = datevec(yrTable.mdate);
Z2 = yrTable.Syn_count./yrTable.VolAnalyzed_ml;
tt = find(ilat<39.5 | ilat>41.5);
ilat(tt) = [];
Z2(tt) = [];
lat_smooth(tt) = [];
yrTable(tt,:) = [];
ind = find(yrTable.lon < -70.883+.24 & yrTable.lon > -70.883-.24); 



%%
%For 2021 site review, Q3 presentation by Rachel
tstr = {'Jan-Mar' 'Apr-Jun' 'Jul-Sep' 'Oct-Dec'};
ymax = [2e4 1e4 1.2e5 1.2e5];
%tstr = {'Feb-Mar' 'Apr-May' 'Jun-Aug' 'Sep-Nov'};
%for yr = 2018 %:2020
figure 
set(gcf,'position', [360 60 500 600])
tl = tiledlayout(4,1, 'TileSpacing', 'compact');
ccmonth = {[1:3] [4:6] [7:9] [10:12]};
%ccmonth = {[2:3] [4:5] [6:8] [9:11]};
ax1 = nexttile;
for cc = 1:length(ccmonth)-1 %2:2:9
    %subplot(5,1,cc/2)
    ii = find((dv(ind,2)>=ccmonth{cc}(1) & dv(ind,2)<=ccmonth{cc}(end)) & dv(ind,1) == yr);
    
    boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    set(gca, 'xdir', 'rev', 'xticklabel', [])
    a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
    yl = ylim; ylim([0 yl(2)]); %auto y-scale
    ylim([0 ymax(cc)]), yl = ylim;
    text(5, yl(2)*.8, tstr{cc}, 'fontsize', 14)
    %text(5, yl(2)*.7, [a(1,:) '-' a(2,:)], 'fontsize', 14)
    %ylim([0 1e5]) %all same y scale
    %text(5, 12e5, [a(1,:) '-' a(2,:)], 'fontsize', 14)
    nexttile
end

cc = cc+1;
%cc = 10;
%ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
ii = find((dv(ind,2)>=ccmonth{cc}(1) & dv(ind,2)<=ccmonth{cc}(end)) & dv(ind,1) == yr);    
boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
set(gca, 'xdir', 'rev')
yl = ylim; ylim([0 yl(2)]);
ylim([0 ymax(cc)]); yl = ylim;
a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
text(5, yl(2)*.8, tstr{cc}, 'fontsize', 14)
set(gca, 'XTickLabelRotation',90)
set(gcf, 'paperposition', [.25 .25 6 10.5])

ylabel(tl, '\itSynechococcus\rm concentration (ml^{-1})')
xlabel(tl, 'Latitude')
title(ax1, yr, 'fontsize', 14)
%end

%%
dv2018 = datevec(yrTable2018.mdate);
dv2019 = datevec(yrTable2019.mdate);
dv2020 = datevec(yrTable2020.mdate);
figure
X = [yrTable2018.temperature; yrTable2019.temperature; yrTable2020.temperature];
Y = [yrTable2018.Syn_count./yrTable2018.VolAnalyzed_ml; yrTable2019.Syn_count./yrTable2019.VolAnalyzed_ml; yrTable2020.Syn_count./yrTable2020.VolAnalyzed_ml];
C = [dv2018(:,2); dv2019(:,2); dv2020(:,2)];
scatter(X,Y,10,C, 'filled')
% scatter(yrTable2018.temperature, yrTable2018.Syn_count./yrTable2018.VolAnalyzed_ml,10,dv2018(:,2), 'filled')
% hold on
% scatter(yrTable2019.temperature, yrTable2019.Syn_count./yrTable2019.VolAnalyzed_ml,10,dv2019(:,2), 'filled')
% scatter(yrTable2020.temperature, yrTable2020.Syn_count./yrTable2020.VolAnalyzed_ml,10,dv2020(:,2), 'filled')
cbh = colorbar
set(cbh, 'xdir', 'rev')
caxis([1 12])
set(cbh, 'ticklabels', {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'}')
title('Cross-shelf observations 2018-2020')
axis([-2 28 20 1e6])
set(gca, 'yscale', 'log', 'box', 'on')
ylabel('\itSynechococcus\rm concentration (ml^{-1})', 'fontsize', 14)
xlabel('Temperature (\circC)', 'fontsize', 14)
