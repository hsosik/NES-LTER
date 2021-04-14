cruises = {'AR22' 'AR24A' 'AR24B' 'AR24C' 'EN608' 'AR28A' 'AR28B' 'EN617'...
    'AR31A' 'AR31B' 'AR31C' 'AR32' 'EN627' 'AR34A' 'AR34B' 'AR38' 'AR39a'...
    'AR39B' 'EN644'};  %'AR16' 'AR44' 'EN649' 'EN655' 'EN657' 'EN661'

%cruises = {'EN608' 'AR28B' 'EN617' 'AR31B'};  %'AR16' 'AR34B' 'AR38' 'AR39'
%cruises = {'EN608' 'EN617' 'EN627' 'EN644'};
%cruises = {'EN617' 'EN644'};

ubase = '\\sosiknas1\IFCB_data\NESLTER_transect\match_up\';
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

s2017 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017.mat');
s2018 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018.mat');
s2019 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019.mat');
%s2017 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017.mat');
%s2018 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018.mat');
%s2019 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019.mat');

IFCBsum = table;
slist = {'filelist' 'classcount' 'meta_data' 'classbiovol' 'mdate'};
for count = 1:length(slist)
    s = slist{count}; IFCBsum.(s) = [s2017.(s); s2018.(s); s2019.(s)];
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
%%
%%
%using the m_map toolbox: https://www.eoas.ubc.ca/~rich/map.html
%cc = strmatch('Guinardia_delicatula', class2use, 'exact');
%cc = strmatch('Chrysochromulina', class2use, 'exact');
cc = strmatch('Hemiaulus', class2use, 'exact');
m_proj('Mollweide','long',[-72.5 -69.5],'lat',[39.5 43]);
[ASIT(1), ASIT(2)] = m_ll2xy(-70.5667,41.325);
%for cc1 = 1:length(ciliate_ind)
%cc = ciliate_ind(cc1);
for cc = cc %1:length(class2use)
    for count = 1:length(cruises)
        ind = strmatch(cruises{count}, match_uw.cruise);
        figure
        m_usercoast('gumby3','patch',[.8 .8 .8],'edgecolor','k');
        m_grid('box','fancy','tickdir','out','fontsize', 10, 'fontname','Times New Roman');
        %m_elev('contour',[-200 -200],'edgecolor','k');
        m_tbase('contour',[-200 -200],'edgecolor','k');
        [X,Y] = m_ll2xy(match_uw.lon(ind), match_uw.lat(ind)); %convert positions of IFCB manual points to map coordinates
        Z = log10(sum(IFCBsum.classcount(b(ind),cc),2)./IFCBsum.meta_data.ml_analyzed(b(ind)));
        ii = find(isinf(Z));
        plot(X(ii),Y(ii), 'ok', 'markersize',4)
        hold on
        Z(ii) = NaN;
        scatter(X,Y,20,Z,'filled') %microg per liter
        cbh = colorbar;
        caxis([log10(.5) log10(100)])
        set(cbh, 'ytick', log10([1 10 100]))
        set(cbh, 'yticklabel', 10.^(get(cbh, 'ytick')))
        set(cbh, 'position', [.76 .12 .03 .6])
        title(cbh, {'chains' ; 'ml^{-1}'}, 'fontsize', 14)
        plot(ASIT(1), ASIT(2), 'rpentagram', 'markersize', 20, 'markerfacecolor', 'r')
        text(ASIT(1), ASIT(2)*1.003, '  MVCO', 'fontsize', 14, 'color', 'r')
        title([cruises{count} ' ' datestr(min(floor(match_uw.mdate(ind))),'dd mmm') '-' datestr(max(floor(match_uw.mdate(ind))),'dd mmm yyyy')])
    end
end

%%
cc = strmatch('Hemiaulus', class2use);
%cc = strmatch('Guinardia_delicatula', class2use, 'exact');
%cc = strmatch('Guinardia_delicatula', class2use);
figure
for count = 1:length(cruises)
    ind = strmatch(cruises{count}, match_uw.cruise);
    subplot(4,4,count)
    plot(match_uw.lat(ind), IFCBsum.classcount(b(ind),cc)./IFCBsum.meta_data.ml_analyzed(b(ind)), '*', 'markersize', 6)
    plot(match_uw.lat(ind), sum(IFCBsum.classcount(b(ind),cc)./IFCBsum.meta_data.ml_analyzed(b(ind)),2), '*', 'markersize', 6)
    title(cruises(count))
    set(gca, 'xdir', 'reverse', 'xlim', [39.5 41.5])
    xlabel('Latitude', 'fontsize', 14)
    ylabel('Concentration (ml^{-1})', 'fontsize', 14)
    ylim([0 120])
end
%pause
%%
figure
cc = strmatch('Hemiaulus', class2use);
%cc = strmatch('Chaetoceros', class2use,'exact');
%cc = strmatch('Thalassiosira', class2use,'exact');
ind = strmatch('EN644', match_cast.cruise);
s = scatter(match_cast.lat(ind), match_cast.depth(ind), 20, log10(IFCBsum.classcount(b_cast(ind),cc)./IFCBsum.meta_data.ml_analyzed(b_cast(ind))))
set(gca, 'xdir', 'rev', 'ydir', 'rev')
set(s,'MarkerFaceColor', get(s,'markeredgecolor'))
colorbar

%figure
%ind = strmatch('EN644', match_uw.cruise);
%[count_sort,a] = sort(sum(IFCBsum.classcount(b(ind),:)), 'descend')
%[class2use(a(1:20))' cellstr(num2str(count_sort(1:20)'))]
%pause

%%
for cc = 1:length(class2use)
    figure(2), clf
    clim = ([-1 max([.1 max(log10(IFCBsum.classcount(b,cc)./IFCBsum.meta_data.ml_analyzed(b)))])]);
    tsdiagram([30 37], [-2 28], 10)
    hold on    
    scatter(match_uw.salinity, match_uw.temperature, 10, log10(IFCBsum.classcount(b,cc)./IFCBsum.meta_data.ml_analyzed(b)+.001))
    caxis(clim)
    cbh = colorbar;
    set(cbh, 'yticklabel', round(100*10.^(get(cbh, 'ytick')))/100)
    title([class2use{cc} '  cells or chains per ml'], 'interpreter','none')
    %print(['c:\work\lter\IFCB_classes\TS\' class2use{cc}], '-dpng')
    pause
end

return

%%
figure
count = strmatch('EN644', cruises)
ind = strmatch(cruises{count}, match_uw.cruise);
subplot(1,3,1)
scatter(match_uw.lon(ind), match_uw.lat(ind), 10, match_uw.temperature(ind))
colorbar
subplot(1,3,2)
scatter(match_uw.lon(ind), match_uw.lat(ind), 10, match_uw.salinity(ind))
caxis([26 inf])
colorbar
subplot(1,3,3)
scatter(match_uw.temperature(ind), match_uw.salinity(ind), 10, day(datetime(match_uw.mdate(ind), 'ConvertFrom', 'datenum'), 'dayofyear'))
ylim([26 36])
colorbar

%%
cstr = 'EN617';
figure
ind = strmatch(cstr, match_uw.cruise);
yyaxis left, plot(match_uw.lat(ind), match_uw.temperature(ind), '.-'), ylim([19 28])
yyaxis right, plot(match_uw.lat(ind), match_uw.salinity(ind), '.-'), ylim([30 35])
set(gca, 'xdir', 'reverse')
title(cstr)

%%
figure
nind = find(ismember(nut.cruise, cstr) & nut.depth <5);
yyaxis left, plot(nut.latitude(nind), nut.nitrate_nitrite(nind), '.-')
hold on
yyaxis left, plot(nut.latitude(nind), nut.ammonium(nind), '*-')
yyaxis left, plot(nut.latitude(nind), nut.phosphate(nind), '+-')
ylim([0 .5])
yyaxis right, plot(nut.latitude(nind), nut.silicate(nind), '.-')
ylim([0 3])
set(gca, 'xdir', 'reverse')
title(cstr)

%%
for count = 1:length(class2use)
    ind = find(match_uw.mdate < datenum('2018-3-1')); %EN608
    plot(match_uw.lat(ind), IFCBsum.classcount(b(ind),count)./IFCBsum.ml_analyzed(b(ind)), '*b', 'markersize', 6)
    hold on
    ind = find(match_uw.mdate > datenum('2018-7-1') & match_uw.mdate < datenum('2018-9-1')); %EN617
    plot(match_uw.lat(ind), IFCBsum.classcount(b(ind),count)./IFCBsum.ml_analyzed(b(ind)), 'or','markerfacecolor', 'r', 'markersize', 4)
    ind = find(match_uw.mdate > datenum('2019-1-1') & match_uw.mdate < datenum('2019-3-1')); %EN627
    plot(match_uw.lat(ind), IFCBsum.classcount(b(ind),count)./IFCBsum.ml_analyzed(b(ind)), 'sb','markerfacecolor', 'b', 'markersize', 4)
    ind = find(match_uw.mdate > datenum('2019-7-1')); %EN644
    plot(match_uw.lat(ind), IFCBsum.classcount(b(ind),count)./IFCBsum.ml_analyzed(b(ind)), '^r', 'markersize', 6)
    title(class2use(count), 'interpreter', 'none', 'fontsize', 14)
    set(gca, 'xdir', 'reverse', 'xlim', [39.5 41.5])
    xlabel('Latitude', 'fontsize', 14)
    ylabel('Concentration (ml^{-1})', 'fontsize', 14)
    legend('EN608', 'EN617', 'EN627', 'EN644')
    axis square
    pause
    clf
end
%%
%plot for NES-LTER SC lightning talk May 2020

group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
group_table.CNN_classlist(strmatch('Pseudo-nitzschia', group_table.CNN_classlist)) = {'Pseudo_nitzschia'};
[~,ia,ib] = intersect(group_table.CNN_classlist, class2use);
diatom_ind = ib(find(group_table.Diatom(ia)));
%[~,exclude_ind] = intersect(class2use, {'Bacillariophyceae' 'Licmophora' 'Nanoneis' 'Thalassiosira' 'Thalassiosira_TAG_external_detritus'});
%diatom_ind = setdiff(diatom_ind, exclude_ind);
%dino_ind = ib(find(group_table.Dinoflagellate(ia)));
%%
figure, hold on
Z2 = (sum(IFCBsum.classbiovol(b,diatom_ind),2)./IFCBsum.meta_data.ml_analyzed(b));
ind = find(strcmp(match_uw.cruise, 'EN608') & match_uw.lat < 41.2);
plot(match_uw.temperature(ind), Z2(ind),'ob', 'markerfacecolor', 'b', 'markersize', 4)
ind = find(strcmp(match_uw.cruise, 'EN617') & match_uw.lat < 41.2);
plot(match_uw.temperature(ind), Z2(ind), '^r', 'markerfacecolor', 'r', 'markersize', 4)
ind = find(strcmp(match_uw.cruise, 'EN627') & match_uw.lat < 41.2);
plot(match_uw.temperature(ind), Z2(ind), '^c', 'markerfacecolor', 'c', 'markersize', 4)
xlabel('Temperature (\circC)', 'fontsize', 18)
ylabel('Biovolume (\mum^3 ml^{-1})', 'fontsize', 18)
text(10, 1e6, 'Diatoms', 'fontsize', 18)
ylim([0 15e5])
set(gca, 'yscale', 'log')
legend('EN608', 'EN617', 'EN627', 'location', 'southwest')
legend('Winter 2018', 'Summer 2018', 'Winter 2019', 'location', 'south')
xlim([0 30])

ind = find(strcmp(match_uw.cruise, 'EN644') & match_uw.lat < 41.2);
plot(match_uw.temperature(ind), Z2(ind),'^', 'color', [.93 .69 .13], 'markerfacecolor', [.93 .69 .13], 'markersize', 4)

legend('Winter 2018', 'Summer 2018', 'Winter 2019', 'Summer 2019', 'location', 'south')
%legend('EN608', 'EN617', 'EN627', 'EN644', 'location', 'southwest')
%% 
% NES-LTER 2021 annual meeting
Z2 = (sum(IFCBsum.classbiovol(b,diatom_ind),2)./IFCBsum.meta_data.ml_analyzed(b));
ind = find(match_uw.lat>41.25 & match_uw.lon <-71.1); %points towards Narragansett Bay
Z2(ind) = NaN;
dv = datevec(match_uw.mdate);
yd_vec = match_uw.mdate-datenum(dv(:,1),1,0);
figure
scatter(match_uw.temperature, Z2, 20,yd_vec, '.')
text(10, 1e6, 'Diatoms', 'fontsize', 18)
ylim([0 15e5])
set(gca, 'yscale', 'log')
xlim([0 30])
xlabel('Temperature (\circC)', 'fontsize', 18)
ylabel('Biovolume (\mum^3 ml^{-1})', 'fontsize', 18)
cb = colorbar; colormap jet
set(cb, 'xtick', datenum(0,1:12,1))
datetick(cb,'x', 'mmm', 'keepticks')

%%
figure
ind = find(dv(:,1) == 2017);
scatter(match_uw.temperature(ind), Z2(ind), 20,yd_vec(ind), '.')
title(dv(ind(1),1))
text(2, 1500, 'Diatoms', 'fontsize', 18)
set(gca, 'yscale', 'log')
ylim([1e3 6e6])
xlim([0 30])
xlabel('Temperature (\circC)', 'fontsize', 18)
ylabel('Biovolume (\mum^3 ml^{-1})', 'fontsize', 18)
cb = colorbar; colormap hsv
set(cb, 'xtick', datenum(0,1:12,1))
datetick(cb,'x', 'mmm', 'keepticks')
%%
figure
plot(match_uw.lat, Z2,'.')
set(gca, 'xdir', 'rev')
%%

group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
[~,ia,ib] = intersect(group_table.CNN_classlist, class2use);
diatom_ind = ib(find(group_table.Diatom(ia)));
notalive_ind = [ib(find(group_table.OtherNotAlive(ia))); ib(find(group_table.IFCBArtifact(ia)))];
alive_ind = 1:length(class2use); alive_ind(notalive_ind) = [];
alive_ind(strmatch( 'Phaeocystis', class2use(alive_ind))) = []; %seems to be detritus for en644
alive_ind(strmatch( 'unclassified', class2use(alive_ind))) = [];
%%
ncp617 = load('C:\work\LTER\Stanley_data\ncplterEn617.mat');
ncp644 = load('C:\work\LTER\Stanley_data\ncplterEn644.mat');

%%
%%
%cc = strmatch('Guinardia_delicatula', class2use, 'exact');
cc = strmatch('Hemiaulus', class2use, 'exact');
for cc = cc %1:length(class2use)
    figure(1), clf
    clim = ([-1 max([.1 max(log10(IFCBsum.classcount(b,cc)./IFCBsum.meta_data.ml_analyzed(b)))])]);
    for count = 1:length(cruises)
        ind = strmatch(cruises{count}, match_uw.cruise);
        subplot(4,4,count)
        scatter(match_uw.lon(ind), match_uw.lat(ind), 20, log10(IFCBsum.classcount(b(ind),cc)./IFCBsum.meta_data.ml_analyzed(b(ind))+.001))
        title([cruises{count} ' ' datestr(min(floor(match_uw.mdate(ind))),'dd mmm') '-' datestr(max(floor(match_uw.mdate(ind))),'dd mmm yyyy')])
        set(gca, 'ylim', [39.5 41.5], 'xlim', [-71.5 -70.5])
        caxis(clim)
        %xlabel('Latitude', 'fontsize', 14)
        %ylabel('Concentration (ml^{-1})', 'fontsize', 14)
    end
    %caxis([-1 max([1 max(log10(IFCBsum.classcount(b,cc)./IFCBsum.ml_analyzed(b)))])])
    cbh = colorbar;
    set(cbh, 'position', [.92 .11 .01 .3]) %.0074 .2158
    set(cbh, 'yticklabel', round(100*10.^(get(cbh, 'ytick')))/100)
    text(-75.5,50.6, [class2use{cc} '  cells or chains per ml'],'interpreter', 'none', 'fontsize', 20)
    set(gcf, 'Position', get(0, 'Screensize'));
   % print(['c:\work\lter\IFCB_classes\' class2use{cc}], '-dpng')
    pause
end

%%
%simple matlab mapping
%cc = strmatch('Guinardia_delicatula', class2use);
cc = strmatch('Hemiaulus', class2use, 'exact');
for cc = cc %1:length(class2use)
    %figure(1), clf
    %clim = ([-1 max([.1 max(log10(IFCBsum.clasascount(b,cc)./IFCBsum.ml_analyzed(b)))])]);
    ASIT = [41.325 -70.5667];
    for count = 1:length(cruises)
        ind = strmatch(cruises{count}, match_uw.cruise);
     %   subplot(4,4,count)
        figure
        latlim = [39.5 43]; %[39.5 41.5];
        lonlim = [-72 -69]; %[-71.5 -70.5];
        ax = usamap(latlim, lonlim)
        geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
        %scatterm(meta.latitude(bm), meta.longitude(bm),10,log10(sum(classC(am,b),2)/1000./ml_analyzed(am)), 'filled') %microg per liter
        scatterm(match_uw.lat(ind), match_uw.lon(ind), 20, log10(sum(IFCBsum.classcount(b(ind),cc),2)./IFCBsum.meta_data.ml_analyzed(b(ind))+.001), 'filled')        %caxis([log10(.3) log10(250)])
        cbh = colorbar;
        %caxis([log10(.1) log10(50)])
        caxis([log10(.1) log10(100)])
        set(cbh, 'ytick', log10([ 1 10 100]))
        %set(cbh, 'ytick', log10([.5 5 50]))
        set(cbh, 'yticklabel', 10.^(get(cbh, 'ytick')))
        plotm(ASIT(1), ASIT(2), 'rpentagram', 'markersize', 20, 'markerfacecolor', 'r')
        title([cruises{count} ' ' datestr(min(floor(match_uw.mdate(ind))),'dd mmm') '-' datestr(max(floor(match_uw.mdate(ind))),'dd mmm yyyy')])
        end
    %caxis([-1 max([1 max(log10(IFCBsum.classcount(b,cc)./IFCBsum.ml_analyzed(b)))])])
    %cbh = colorbar;
    %set(cbh, 'position', [.92 .11 .01 .3]) %.0074 .2158
    %set(cbh, 'yticklabel', round(100*10.^(get(cbh, 'ytick')))/100)
    %text(-75.5,50.6, [class2use{cc} '  cells or chains per ml'],'interpreter', 'none', 'fontsize', 20)
    %set(gcf, 'Position', get(0, 'Screensize'));
    %print(['c:\work\lter\IFCB_classes\' class2use{cc}], '-dpng')
%    pause
end
