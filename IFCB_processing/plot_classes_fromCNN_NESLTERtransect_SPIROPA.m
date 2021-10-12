cruises = {'AR22' 'AR24A' 'AR24B' 'AR24C' 'EN608' 'AR28A' 'AR28B' 'EN617'...
    'AR31A' 'AR31B' 'AR31C' 'AR32' 'AR44' 'EN627' 'AR34A' 'AR34B' 'AR38' 'AR39a'...
    'AR39B' 'EN644' 'EN649' 'EN655' 'EN657' 'EN661'};  %'AR16' 'AR48A' 'AR48B'  
cruises2 = {'AR29' 'RB1904' 'TN368' };

ubase = '\\sosiknas1\IFCB_data\NESLTER_transect\match_up\';
ubase2 = '\\sosiknas1\IFCB_data\SPIROPA\match_up\';
%ubase = 'c:\work\IFCB_products\NESLTER_transect\match_up\';

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
%  %   n = webread(['https://nes-lter-data.whoi.edu/api/nut/' cruises{count1} '.csv']);
%  %   n.alternate_sample_id = []; %move this column since the type doesn't match between all cruises
%  %   nut = [nut; n];
 end
% %nut.mdate = datenum(nut.date, 'yyyy-mm-dd hh:MM:ss+00:00');
ind = find(strcmp(match_uw.cruise, 'EN627') & match_uw.salinity > 34);
match_uw.temperature(ind) = NaN;
match_uw.salinity(ind) = NaN;
%%
for count1 = 1:length(cruises2)
     disp(cruises2{count1})
     if strmatch('AR29', cruises2{count1})
         orig_var = {'SBE48T' 'SBE45S'};
     else
         orig_var = {'T' 'SSPS'};
     end
     u = load([ubase2 'SPIROPA_' cruises2{count1} '__newdashboard_USEMEuw_match.mat']);
     u.IFCB_match_uw_results.Properties.VariableNames(orig_var) = {'temperature' 'salinity'};
     [~,ia] = ismember(slist, u.IFCB_match_uw_results.Properties.VariableNames);
     match_uw = [match_uw; u.IFCB_match_uw_results(:,ia)];
     if exist([ubase2 'SPIROPA_' cruises2{count1} '__newdashboard_USEMEcast_match.mat'], 'file')
        u = load([ubase2 'SPIROPA_' cruises2{count1} '__newdashboard_USEMEcast_match.mat']);  
        u.IFCB_match_btl_results.Properties.VariableNames = lower(u.IFCB_match_btl_results.Properties.VariableNames);
        if strmatch('AR29', cruises2{count1})
            uind = find(u.IFCB_match_btl_results.cast == 106);
            if isempty(u.IFCB_match_btl_results.datetime{uind(1)})
                u.IFCB_match_btl_results.Depth_m(uind) = 2;
                u.IFCB_match_btl_results.datetime(uind) = {'2018-04-25 04:13:51+00:00'};
            end
        end    
        u.IFCB_match_btl_results.mdate = datenum(u.IFCB_match_btl_results.datetime, 'yyyy-mm-dd hh:MM:SS+00:00');
        [~,ia] = ismember(slist_cast, u.IFCB_match_btl_results.Properties.VariableNames);
        match_cast = [match_cast; u.IFCB_match_btl_results(:,ia)];
     end
 end
%%

s2017 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017.mat');
s2018 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018.mat');
s2019 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019.mat');
s2020 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2020.mat');
s2021 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2021.mat');
s2018b = load('\\sosiknas1\IFCB_products\SPIROPA\summary\summary_biovol_allHDF_min20_2018.mat');
s2019b = load('\\sosiknas1\IFCB_products\SPIROPA\summary\summary_biovol_allHDF_min20_2019.mat');
tag5 = repmat(cellstr(''),size(s2018b.meta_data,1),1);
s2018b.meta_data = addvars(s2018b.meta_data, tag5, 'After', 'tag4', 'NewVariableNames', 'tag5' );
s2018b.meta_data.cast = cellstr(num2str(s2018b.meta_data.cast));
s2018b.meta_data.tag4 = cellstr(num2str(s2018b.meta_data.tag4));
tag5 = repmat(cellstr(''),size(s2019b.meta_data,1),1);
s2019b.meta_data = addvars(s2019b.meta_data, tag5, 'After', 'tag4', 'NewVariableNames', 'tag5' );
s2019b.meta_data.cast = cellstr(num2str(s2019b.meta_data.cast));
s2019b.meta_data.tag4 = cellstr(num2str(s2019b.meta_data.tag4));
%s2017 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017.mat');
%s2018 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018.mat');
%s2019 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019.mat');

IFCBsum = table;
slist = {'filelist' 'classcount' 'meta_data' 'classbiovol' 'mdate'};
for count = 1:length(slist)
    s = slist{count}; IFCBsum.(s) = [s2017.(s); s2018.(s); s2019.(s); s2020.(s); s2021.(s); s2018b.(s); s2019b.(s)];
end

class2use = s2017.class2use;

match_uw = sortrows(match_uw, 'mdate');
[a,b] = ismember(match_uw.pid, IFCBsum.filelist);
%
match_cast = sortrows(match_cast, 'mdate');
[a_cast,b_cast] = ismember(match_cast.pid, IFCBsum.filelist);


temp = find(a==0);
if ~isempty(temp)
    disp('Skipping unmatched files: ')
    disp(match_uw.pid(temp))
    disp('Skipped above unmatched files')
    match_uw(temp,:) = [];
    b(a==0) = [];
end

temp = find(a_cast==0);
if ~isempty(temp)
    disp('Skipping unmatched files: ')
    disp(match_cast.pid(temp))
    disp('Skipped above unmatched files')
    match_cast(temp,:) = [];
    b_cast(a_cast==0) = [];
end

group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
%%7 April 2021, including the line below was a mistake on plots for NES-LTER annual meeting 2021, P-n accidentally excluded from diatoms...
%%group_table.CNN_classlist(strmatch('Pseudo-nitzschia', group_table.CNN_classlist)) = {'Pseudo_nitzschia'};
[~,ia,ib] = intersect(group_table.CNN_classlist, class2use);
diatom_ind = ib(find(group_table.Diatom(ia)));
dino_ind = ib(find(group_table.Dinoflagellate(ia)));

dv = datevec(match_uw.mdate);
yd_vec = match_uw.mdate-datenum(dv(:,1),1,0);

Z2 = (sum(IFCBsum.classbiovol(b,diatom_ind),2)./IFCBsum.meta_data.ml_analyzed(b));

Z2all = IFCBsum.classbiovol(b,diatom_ind)./IFCBsum.meta_data.ml_analyzed(b);
ind = find(match_uw.lon < -70.883+.24 & match_uw.lon > -70.883-.24); 

for ii = 1:12, numyrs(ii) = length(unique(dv(ind(find(dv(ind,2) == ii)),1))); end
lat_smooth = round(match_uw.lat,1);
ilat = unique(lat_smooth);
ilat = ilat(~isnan(ilat));

%%
for ii = 1:12, numyrs(ii) = length(unique(dv(ind(find(dv(ind,2) == ii)),1))); end
yy = 2017:2021;
bins = linspace(3,6.5,100);
figure, tl = tiledlayout(11,1, 'TileSpacing', 'compact')
for iii = 1:11
    nexttile
    hh = NaN(length(yy), length(bins)-1);
    for cc = 1:length(yy)
        tt = find(dv(ind,2) == iii & dv(ind,1) == yy(cc));
        hh(cc,:) = histcounts(log10(Z2(ind(tt))),bins);
    end
    bar(bins(1:end-1),hh, 'stacked', 'edgecolor', 'none'), set(gca, 'xticklabel', [])
    yl = ylim;
    text(3.1, yl(2)*.8, datestr(datenum(2020,iii,1), 'mmm'))
end
ylabel(tl, 'Frequency of samples (transect 2017-2020)')
xlabel(tl, {'Diatom biovolume concentration';  '(log10 \mum^3 ml^{-1})'})
set(gca, 'xticklabelmode', 'auto')
nexttile(6)
legend(num2str(yy'), 'location', 'southeast')
set(gcf, 'paperposition', [.25 .25 4 10.5])
set(gcf, 'position', [200 50 400 600])

print('c:\work\diatom_bv_month', '-dpng')

%%
figure
boxplot(Z2(ind),dv(ind,2), 'whisker', 3, 'notch', 'on', 'datalim', [0 1.5e6], 'extrememode', 'compress', 'outliersize', 4)
xlim([.5 12.5])
set(gca, 'xtick', 1:12)
set(gca, 'xticklabel', num2str(get(gca, 'xtick')'))
set(gca, 'position', [.13 .11 .775 .7])
tt = regexprep(cellstr(num2str(histcounts(dv(:,2),1:13)')), ' ', '');
text(1:12,2e6*ones(12,1),tt, 'fontsize', 10, 'horizontalalignment', 'center', 'color', 'b')
for ii = 1:12, numyrs(ii) = length(unique(dv(ind(find(dv(ind,2) == ii)),1))); end
text(1:12,2.15e6*ones(12,1),num2str(numyrs'), 'fontsize', 10, 'horizontalalignment', 'center', 'color', 'b')
set(gca, 'xticklabel', ['JFMAMJJASOND']')
ylabel('Diatom biovolume concentration (\mum^3 ml^{-1})')

print('c:\work\diatom_bv_monthly_box', '-dpng')

%%
%[ mdate_mat, y_mat, yearlist, yd ] = timeseries2ydmat( match_uw.mdate, Z2 );
%for ii = 1:12, y_month(ii) = nanmean(Z2(dv(:,2)==ii)); end
%ilat = 39.2:.05:41.5;
%for ii = 1:length(ilat)-1, iii = find(match_uw.lat(ind) >ilat(ii) & match_uw.lat(ind)<=ilat(ii+1)); Z2latmd(ii) = nanmedian(Z2(ind(iii)));end

%figure
%plot(match_uw.lat(ind), Z2(ind), '.')
%hold on
%plot(ilat(1:end-1)+.05, Z2latmd,'-', 'linewidth', 3)
%set(gca, 'xdir', 'rev')

%figure
%plot(1:366, nanmean(y_mat,2), 'linewidth', 3)
%ilat = (39.2:.05:41.5)';

% for ii = 1:12, numyrs(ii) = length(unique(dv(ind(find(dv(ind,2) == ii)),1))); end
% lat_smooth = round(match_uw.lat,1);
% ilat = unique(lat_smooth);
% ilat = ilat(~isnan(ilat));
figure
t = boxplot(Z2(ind),lat_smooth(ind), 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'inline');
set(gca, 'xdir', 'rev')
yl = ylim; ylim([0 yl(2)])
ylabel('Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel('Latitude')
print('c:\work\diatom_bv_lat_box', '-dpng')

LTER = [41 11.8 70 53; 41 1.8 70 53; 40 51.8 70 53; 40 41.8 70 53; 40 30.8 70 53; 40 21.8 70 53; 40 13.6 70 53; ...
    40 08.2 70 46.5; 40 5.9 70 53; 39 56.4 70 53; 39 46.4 70 53];
LTER_labels = {'L1'; 'L2'; 'L3'; 'L4'; 'L5'; 'L6'; 'L7'; 'L8'; 'L9'; 'L10'; 'L11'};
LTERdeg(:,1) = LTER(:,1)+LTER(:,2)/60;
LTERdeg(:,2) = LTER(:,3)+LTER(:,4)/60;

figure
hh = histcounts(lat_smooth(ind),ilat(~isnan(ilat)));
bar(ilat(2:end), hh), set(gca, 'xdir', 'rev', 'xtick', ilat, 'XTickLabelRotation', 90)
set(gca, 'xdir', 'rev', 'xgrid', 'on')
xlim([39.1 41.55])
ylabel({'Number of'; 'samples'})
set(gcf, 'position', [360 530 560 200])
set(gca, 'position', [.13 .1386 .775 .4])
hold on
text(LTERdeg(:,1), 1500*ones(size(LTER_labels)), 'v', 'fontsize', 8, 'verticalalignment', 'bottom', 'fontname', 'arial')
text(LTERdeg(:,1), 1700*ones(size(LTER_labels)), LTER_labels, 'fontsize', 8, 'verticalalignment', 'bottom', 'horizontalalignment', 'left')
print('c:\work\samplenum_by_lat', '-dpng')

%%

figure 
%set(gcf,'paperposition', [.5 .5 5 10])
tl = tiledlayout(5,1, 'TileSpacing', 'compact')
for cc = 2:2:9
    %subplot(5,1,cc/2)
    nexttile
    ii = find(dv(ind,2)==cc |dv(ind,2)==cc+1); 
    %boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on');
    boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    set(gca, 'xdir', 'rev', 'xticklabel', [])
    a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
    yl = ylim; ylim([0 yl(2)]); %auto y-scale
    text(5, yl(2)*.7, [a(1,:) '-' a(2,:)], 'fontsize', 14)
    %ylim([0 1.83e6]) %all same y scale
    %text(5, 12e5, [a(1,:) '-' a(2,:)], 'fontsize', 14)
end
nexttile
cc = 10;
ii = find(dv(ind,2)==cc |dv(ind,2)==cc+1); 
boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
set(gca, 'xdir', 'rev')
%yl = ylim; ylim([0 yl(2)]);
ylim([0 1.83e6]);
a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
text(5, 12e5, [a(1,:) '-' a(2,:)], 'fontsize', 14)
set(gca, 'XTickLabelRotation',90)
set(gcf, 'paperposition', [.25 .25 6 10.5])
print('c:\work\diatom_bv_lat_by_2mon_box', '-dpng')

ylabel(tl, 'Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel(tl, 'Latitude')

%%
tstr = {'Feb-Mar' 'Apr-May' 'Jun-Jul' 'Aug-Sep' 'Oct-Nov'};
for yy = 2017:2021
figure 
set(gcf,'position', [360 60 500 600])
tl = tiledlayout(5,1, 'TileSpacing', 'compact')
for cc = 2:2:9
    %subplot(5,1,cc/2)
    nexttile
    ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
    boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    set(gca, 'xdir', 'rev', 'xticklabel', [])
    a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
    yl = ylim; ylim([0 yl(2)]); %auto y-scale
    text(5, yl(2)*.7, tstr{cc/2}, 'fontsize', 14)
    %text(5, yl(2)*.7, [a(1,:) '-' a(2,:)], 'fontsize', 14)
    %ylim([0 1.83e6]) %all same y scale
    %text(5, 12e5, [a(1,:) '-' a(2,:)], 'fontsize', 14)
end
nexttile
cc = 10;
ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
set(gca, 'xdir', 'rev')
%yl = ylim; ylim([0 yl(2)]);
ylim([0 1.83e6]);
a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
text(5, 12e5, tstr{cc/2}, 'fontsize', 14)
set(gca, 'XTickLabelRotation',90)
set(gcf, 'paperposition', [.25 .25 6 10.5])

ylabel(tl, 'Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel(tl, 'Latitude')
title(tl, yy)
end

%%
%%
tstr = {'Feb-Mar' 'Apr-May' 'Jun-Jul' 'Aug-Sep' 'Oct-Nov'};
for yy = 2017:2020
figure 
set(gcf,'position', [360 20 500 600])
tl = tiledlayout(5,1, 'TileSpacing', 'compact');
%tZ2 = Z2all(:,strmatch('Hemiaulus', class2use(diatom_ind)));
c2plot = 'Guinardia_flaccida';
c2plot = 'Guinardia_delicatula';
tZ2 = sum(Z2all(:,strmatch(c2plot, class2use(diatom_ind))),2);
for cc = 2:2:9
    %subplot(5,1,cc/2)
    nexttile
    ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
    boxplot([tZ2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    set(gca, 'xdir', 'rev', 'xticklabel', [])
    a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
    yl = ylim; ylim([0 yl(2)]); %auto y-scale
    %text(5, yl(2)*.7, tstr{cc/2}, 'fontsize', 14)
    %text(5, yl(2)*.7, [a(1,:) '-' a(2,:)], 'fontsize', 14)
    ylim([0 1e6]) %all same y scale
    text(7, .8e6, tstr{cc/2}, 'fontsize', 14)
end
nexttile
cc = 10;
ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
boxplot([tZ2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
set(gca, 'xdir', 'rev')
%yl = ylim; ylim([0 yl(2)]);
ylim([0 1.83e6]);
a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
text(5, 12e5, tstr{cc/2}, 'fontsize', 14)
set(gca, 'XTickLabelRotation',90)
set(gcf, 'paperposition', [.25 .25 6 10.5])

ylabel(tl, 'Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel(tl, 'Latitude')
title(tl, [c2plot ' ' num2str(yy)], 'interpreter', 'none')
end

%%
ii2018 = find(dv(ind,2)==2 & dv(ind,1) == 2018);
[~, is2018] = sort(sum(Z2all(ind(ii2018),:)), 'descend');
class2use(diatom_ind(is2018(1:10)))

ii2019 = find(dv(ind,2)==2 & dv(ind,1) == 2019);
[~, is2019] = sort(sum(Z2all(ind(ii2019),:)), 'descend');

ii2020 = find(dv(ind,2)==2 & dv(ind,1) == 2020);
[~, is2020] = sort(sum(Z2all(ind(ii2020),:)), 'descend');

top10ind = unique([is2018(1:10) is2019(1:10) is2020(1:10)]);
class2use(diatom_ind(unique([is2018(1:10) is2019(1:10) is2020(1:10)])))

cind2 = strmatch('Thalassiosira', class2use(diatom_ind));
cind = strmatch('Guinardia', class2use(diatom_ind));

figure
tl = tiledlayout(3,1, 'TileSpacing', 'compact');
ind2plot = {is2019(1) cind cind2};
for cc = 1:2
    nexttile
    boxplot([sum(Z2all(ind(ii2019),ind2plot{cc}),2); NaN(size(ilat))],[lat_smooth(ind(ii2019)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    set(gca, 'xdir', 'rev', 'xticklabel', [])
    ylim([0 .5e6]) %all same y scale
end
nexttile
cc = 3;
boxplot([sum(Z2all(ind(ii2019),ind2plot{cc}),2); NaN(size(ilat))],[lat_smooth(ind(ii2019)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
set(gca, 'xdir', 'rev')
ylim([0 .5e6]) %all same y scale
set(gca, 'XTickLabelRotation',90)
%set(gcf, 'paperposition', [.25 .25 6 10.5])

ylabel(tl, 'Biovolume concentration (\mum^3 ml^{-1})')
xlabel(tl, 'Latitude')
title(tl, 'Feb 2019, top 3 diatom genera')

%%
ii2019Aug = find(dv(ind,2)==8 & dv(ind,1) == 2019);
figure
tl = tiledlayout(3,1, 'TileSpacing', 'compact');
ind2plot = strmatch('Hemiaulus', class2use(diatom_ind))
    nexttile
    boxplot([sum(Z2all(ind(ii2019Aug),ind2plot),2); NaN(size(ilat))],[lat_smooth(ind(ii2019Aug)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    set(gca, 'xdir', 'rev')
    ylim([0 2e6]) %all same y scale

set(gca, 'XTickLabelRotation',90)
%set(gcf, 'paperposition', [.25 .25 6 10.5])

ylabel({'Biovolume concentration'; '(\mum^3 ml^{-1})'})
xlabel('Latitude')
title(tl, 'Aug 2019, \itHemiaulus')
%for 2021 NES LTER site review report
figure
set(gcf, 'position', [350 300 550 180])
boxplot([sum(Z2all(ind(ii2019Aug),ind2plot),2); NaN(size(ilat))],[lat_smooth(ind(ii2019Aug)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
set(gca, 'xdir', 'rev')
ylim([0 2e6])
set(gca, 'XTickLabelRotation',90)
ylabel({'Biovolume concentration'; '(\mum^3 ml^{-1})'})
xlabel('Latitude')
xlim([8 26])
text(25.5, 1.1e6, {'\itHemiaulus\rm'; 'August 2019'})

%%
T_smooth = round(match_uw.temperature);
iT = unique(T_smooth);
iT = iT(~isnan(iT));

figure
t = boxplot(Z2(ind),T_smooth(ind), 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on');
yl = ylim; ylim([0 yl(2)])
ylabel('Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel('Temperature (\circC)')
set(gca, 'yscale', 'log')
ylim([1e3 2e6])
%print('c:\work\diatom_bv_T_box', '-dpng')

%%
figure
tt = find(dv(ind,1) == 2018 & dv(ind,2) < 3);
semilogy(match_uw.temperature(ind(tt)), Z2(ind(tt)), '.'), hold on
ylim([1e3 2e6])
xlim([0 30])
ylabel('Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel('Temperature (\circC)')
legend('2018 winter')
print('c:\work\lter\diatom_bv_T1', '-dpng')
%pause
tt = find(dv(ind,1) == 2018 & (dv(ind,2) == 7 | dv(ind,2) ==8));
semilogy(match_uw.temperature(ind(tt)), Z2(ind(tt)), '.')
legend('2018 winter', '2018 summer', 'location', 'southwest')
print('c:\work\lter\diatom_bv_T2', '-dpng')
%pause
tt = find(dv(ind,1) == 2019 & dv(ind,2) < 3);
semilogy(match_uw.temperature(ind(tt)), Z2(ind(tt)), '.')
legend('2018 winter', '2018 summer', '2019 winter', 'location', 'southwest')
print('c:\work\lter\diatom_bv_T3', '-dpng')
tt = find(dv(ind,1) == 2019 & (dv(ind,2) == 7 | dv(ind,2) ==8));
semilogy(match_uw.temperature(ind(tt)), Z2(ind(tt)), 'y.')
legend('2018 winter', '2018 summer', '2019 winter', '2019 summer', 'location', 'southwest')
print('c:\work\lter\diatom_bv_T4', '-dpng')
tt = find(dv(ind,1) == 2020 & dv(ind,2) < 3);
semilogy(match_uw.temperature(ind(tt)), Z2(ind(tt)), '.')
legend('2018 winter', '2018 summer', '2019 winter', '2019 summer', '2029 winter', 'location', 'southwest')
print('c:\work\lter\diatom_bv_T5', '-dpng')
tt = find(dv(ind,1) == 2020 & (dv(ind,2) == 7 | dv(ind,2) ==8));
semilogy(match_uw.temperature(ind(tt)), Z2(ind(tt)), '.')
legend('2018 winter', '2018 summer', '2019 winter', '2019 summer', '2020 winter', '2020 summer', 'location', 'southwest')
print('c:\work\lter\diatom_bv_T6', '-dpng')
tt = find(dv(ind,1) == 2021 & dv(ind,2) < 3);
semilogy(match_uw.temperature(ind(tt)), Z2(ind(tt)), '.')
legend('2018 winter', '2018 summer', '2019 winter', '2019 summer', '2020 winter', '2020 summer', '2021 winter', 'location', 'southwest')
print('c:\work\lter\diatom_bv_T7', '-dpng')

figure
scatter(match_uw.temperature(ind), Z2(ind),10, dv(ind,2))
set(gca, 'yscale', 'log')
ylim([1e3 2e6])
xlim([0 30])
ylabel('Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel('Temperature (\circC)')
cbh = colorbar;
colormap hsv
cmap = colormap;
cmap = cmap(1:22:end,:);
colormap(cmap)
caxis([1 13])
set(cbh, 'Ticks', 1.5:12.5, 'TickLabels', ['JFMAMJJASOND']', 'Direction', 'rev')
set(gca, 'box', 'on')

print('c:\work\lter\diatom_bv_T_all', '-dpng')

%%
%ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy);
ii = find((dv(ind,2)==11 | dv(ind,2)==11));%& dv(ind,1) == 2019); 
ii = find(dv(ind,2));
figure
scatter(match_uw.salinity(ind(ii)), Z2(ind(ii)),10, dv(ind(ii),2), 'filled')
%scatter(match_uw.salinity(ind(ii)), sum(Z2(ind(ii),:),2),10, dv(ind(ii),2), 'filled')
set(gca, 'yscale', 'log')
ylim([1e3 2e6])
xlim([30 36])
ylabel('Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel('Salinity')
cbh = colorbar;
colormap hsv
cmap = colormap;
cmap = cmap(1:22:end,:);
colormap(cmap)
caxis([1 13])
set(cbh, 'Ticks', 1.5:12.5, 'TickLabels', ['JFMAMJJASOND']', 'Direction', 'rev')
set(gca, 'box', 'on')


%%
%%
%ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
ii = find(dv(ind,2));
figure
scatter(match_uw.lat(ind(ii)), match_uw.salinity(ind(ii)), 10, dv(ind(ii),2), 'filled')
set(gca, 'yscale', 'log', 'xdir', 'rev')
%ylim([1e3 2e6])
ylim([30 37])
ylabel('Salinity')
xlabel('Latitude')
cbh = colorbar;
colormap hsv
cmap = colormap;
cmap = cmap(1:22:end,:);
colormap(cmap)
caxis([1 13])
set(cbh, 'Ticks', 1.5:12.5, 'TickLabels', ['JFMAMJJASOND']', 'Direction', 'rev')
set(gca, 'box', 'on')



%%
tstr = {'Feb-Mar' 'Apr-May' 'Jun-Jul' 'Aug-Sep' 'Oct-Nov'};
for yy = 2017:2021
figure 
set(gcf,'position', [360 60 500 600])
tl = tiledlayout(5,1, 'TileSpacing', 'compact')
for cc = 2:2:9
    %subplot(5,1,cc/2)
    nexttile
    ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
    %boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    plot(match_uw.lat(ind(ii)), match_uw.salinity(ind(ii)), '.')
    set(gca, 'xdir', 'rev', 'xticklabel', [])
    a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
    ylim([31.5 36.5]), xlim([39.6 41.6])
    text(41.5, 36, tstr{cc/2}, 'fontsize', 14)
    line(xlim, [34.5 34.5],'color', 'r')
    line(xlim, [35 35],'color', 'r')
    grid
end
nexttile
cc = 10;
ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
plot(match_uw.lat(ind(ii)), match_uw.salinity(ind(ii)), '.')
set(gca, 'xdir', 'rev')
a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
ylim([31.5 36.5])
xlim([39.6 41.6])
text(41.5, 36, tstr{cc/2}, 'fontsize', 14)
set(gcf, 'paperposition', [.25 .25 6 10.5])
grid
ylabel(tl, 'Salinity')
xlabel(tl, 'Latitude')
title(tl, yy)

line(xlim, [34.5 34.5],'color', 'r')
line(xlim, [35 35],'color', 'r')
end

%%
%For 2021 site review, Q3 presentation by Rachel
tt = find(ilat<39.5 | ilat>41.5);
ilat(tt) = [];
Z2(tt) = [];
lat_smooth(tt) = [];
match_uw(tt,:) = [];
dv = datevec(match_uw.mdate);
ind = find(match_uw.lon < -70.883+.24 & match_uw.lon > -70.883-.24); 

tstr = {'Jan-Mar' 'Apr-Jun' 'Jul-Sep' 'Oct-Dec'};
%tstr = {'Feb-Mar' 'Apr-May' 'Jun-Aug' 'Sep-Nov'};
for yy = 2018:2020
figure 
set(gcf,'position', [360 60 500 600])
tl = tiledlayout(4,1, 'TileSpacing', 'compact')
ccmonth = {[1:3] [4:6] [7:9] [10:12]};
%ccmonth = {[2:3] [4:5] [6:8] [9:11]};
ax1 = nexttile
for cc = 1:length(ccmonth)-1 %2:2:9
    ii = find((dv(ind,2)>=ccmonth{cc}(1) & dv(ind,2)<=ccmonth{cc}(end)) & dv(ind,1) == yy);
    
    boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    set(gca, 'xdir', 'rev', 'xticklabel', [])
    a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
    yl = ylim; ylim([0 yl(2)]); %auto y-scale
    ylim([0 1.83e6]); yl = ylim;
    text(6, yl(2)*.7, tstr{cc}, 'fontsize', 14)
    %xlim([5 27.5])
    %text(5, yl(2)*.7, [a(1,:) '-' a(2,:)], 'fontsize', 14)
    %ylim([0 1.83e6]) %all same y scale
    %text(5, 12e5, [a(1,:) '-' a(2,:)], 'fontsize', 14)
    nexttile
end
cc = cc+1;
%cc = 10;
%ii = find((dv(ind,2)==cc | dv(ind,2)==cc+1) & dv(ind,1) == yy); 
ii = find((dv(ind,2)>=ccmonth{cc}(1) & dv(ind,2)<=ccmonth{cc}(end)) & dv(ind,1) == yy);    
boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
set(gca, 'xdir', 'rev')
%yl = ylim; ylim([0 yl(2)]);
ylim([0 1.83e6]);
a = datestr(datenum(2020, unique(dv(ind(ii),2)),1), 'mmm');
text(6, 12e5, tstr{cc}, 'fontsize', 14)
set(gca, 'XTickLabelRotation',90)
set(gcf, 'paperposition', [.25 .25 6 10.5])
%xlim([5 27.5])

ylabel(tl, 'Diatom biovolume concentration (\mum^3 ml^{-1})')
xlabel(tl, 'Latitude')
title(ax1, yy, 'fontsize', 14)
end
