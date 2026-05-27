%run save_groupsTable_byCruise_IFCBattune_merged.m to get sumT (table)
%load("c:\work\IFCB_products\NESLTER_transect\summary\Attune_IFCB_merge_class_summary_uw_allcruises")
p = '\\sosiknas1\Lab_data\Attune\cruise_data\IFCB_Attune_merge\summary_files\';
load([p 'Attune_IFCB_merge_class_summary_uw_allcruises.mat'])
sumT = phytoC_all;
%just along the main transect
ii2use = (~isnan(sumT.phytoC_total) & sumT.longitude < -70.883+.24 & sumT.longitude > -70.883-.24 & sumT.latitude > 39.5 & sumT.latitude < 41.35); %omit too far south of L11, omit too far north into NB

cStr_list = unique(sumT.cruise);
%synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';

%%
lat_step = 0.2;
lat_min = 39.8; lat_max = 41.5;
ix = (lat_min:lat_step:lat_max)';
[u,v] = discretize(sumT.latitude, ix-lat_step/2);
sumT.lat_bin = NaN(size(sumT.latitude));
sumT.lat_bin((~isnan(u))) = ix(u(~isnan(u)));

%NOTE subtract 30 days to offset by 1 mon, i.e., winter=Dec-Jan-Feb...
season_offset = 30; %days
sumT.season = quarter(sumT.dt+season_offset);
sumT.doy = day(sumT.dt, 'dayofyear');
sumT.year = year(sumT.dt);

%f_vr = @(x)variance_ratio(x);
f_vr = @(x)loreau_vr(x);
sumTvr = sumT(ii2use,ismember(sumT.Properties.VariableNames, {'cruise' 'dt' 'season' 'doy' 'year' 'lat_bin' 'temp' 'sal'}));
sumTvr.phytoC_comp = [sumT(ii2use,:).phytoC_0_2 sumT(ii2use,:).phytoC_2_20 sumT(ii2use,:).phytoC_20_100 ];
sumTvr.phytoC_comp = sumTvr.phytoC_comp; 
sumTvr.phytoC_comp_cube_root = sumTvr.phytoC_comp.^(1/3); %cube root
sumTvr = sortrows(sumTvr, "dt");

tbl_cruise = groupsummary(sumTvr, ["cruise" "lat_bin"], "mean", ["dt" "season" "phytoC_comp"]);
tbl_cruise.mean_phytoC_comp_cube_root = tbl_cruise.mean_phytoC_comp.^(1/3);
tbl_cruiseVR = groupsummary(sumTvr, ["cruise" "lat_bin" ], f_vr,"phytoC_comp_cube_root");
tbl_cruiseVR.mean_season = tbl_cruise.mean_season;
tbl_cruise_seasonVR = groupsummary(tbl_cruise, ["mean_season" "lat_bin"], f_vr, "mean_phytoC_comp_cube_root");
tbl_cruise_allVR = groupsummary(tbl_cruise, ["lat_bin"], f_vr, "mean_phytoC_comp_cube_root");
tbl_cruiseVR_season = groupsummary(tbl_cruiseVR, ["mean_season" "lat_bin"], "mean", "fun1_phytoC_comp_cube_root");
%temp = groupsummary(sumTvr, "lat_bin", f_vr, "phytoC_comp");
%tbl_cruiseVR_season = groupsummary(tbl_cruiseVR, ["mean_season" "lat_bin"], "mean", "fun1_phytoC_comp");

tbl_season_year = groupsummary(sumTvr, ["year" "season" "lat_bin"], "mean", ["dt" "season" "phytoC_comp"], 'IncludeEmptyGroups',true);
tbl_season_year.mean_phytoC_comp_cube_root = tbl_season_year.mean_phytoC_comp.^(1/3);
%tbl_season_year t = groupsummary(sumTvr, ["year" "season" "lat_bin"], @(x)nanmean(x,1), ["dt" "phytoC_comp"], 'IncludeEmptyGroups',true);
%tbl_season_year.mean_phytoC_comp = tbl_season_year.fun1_phytoC_comp;
tbl_season_yearVR = groupsummary(sumTvr, ["year" "season" "lat_bin"], f_vr,"phytoC_comp_cube_root", 'IncludeEmptyGroups',true);
tbl_seasonVR = groupsummary(tbl_season_year, ["season" "lat_bin"], f_vr,"mean_phytoC_comp_cube_root");
tbl_season_year.mean_phytoC_comp_cube_root12 = tbl_season_year.mean_phytoC_comp_cube_root(:,1:2);
tbl_season_year.mean_phytoC_comp_cube_root23 = tbl_season_year.mean_phytoC_comp_cube_root(:,2:3);
tbl_season_year.mean_phytoC_comp_cube_root13 = tbl_season_year.mean_phytoC_comp_cube_root(:,[1 3]);

tempi = ~isnan(tbl_season_year.mean_season);
tbl_season_allVR = groupsummary(tbl_season_year(tempi,:), ["lat_bin"], f_vr, ["mean_phytoC_comp_cube_root" "mean_phytoC_comp_cube_root12" "mean_phytoC_comp_cube_root23" "mean_phytoC_comp_cube_root13"]);

tblStats_all = grpstats(sumT(ii2use,:), "lat_bin" , ["nanmean" "nanstd" "nanmedian"],"DataVars",["doy" "temp" "sal" sumT.Properties.VariableNames(strncmp('phyto', sumT.Properties.VariableNames,5))]);
tblStats_cruise = grpstats(sumT(ii2use,:), ["cruise" "lat_bin" ], ["nanmean" "nanstd" "nanmedian"],"DataVars",["doy" "temp" "sal" sumT.Properties.VariableNames(strncmp('phyto', sumT.Properties.VariableNames,5))]);
tblStats_season = grpstats(sumT(ii2use,:), ["season" "lat_bin" ], ["nanmean" "nanstd" "nanmedian"],"DataVars",["doy" "temp" "sal" sumT.Properties.VariableNames(strncmp('phyto', sumT.Properties.VariableNames,5))]);
tblStats_cruise.season = quarter(tblStats_cruise.nanmean_doy+season_offset);
tblStats_cruise1_season2 = grpstats(tblStats_cruise, ["season" "lat_bin" ], ["mean" "std" "median"],"DataVars",["nanmean_doy" "nanmean_temp" "nanmean_sal" tblStats_cruise.Properties.VariableNames(contains(tblStats_cruise.Properties.VariableNames,'phyto'))]);

notes = {'saved from phytoC_transect_stats.m'}
%save([synchrony_path 'Attune_IFCB_pico_nano_micro_seasonal_transect_means_latbin'], 'tbl_season_year', 'tbl_season_allVR', 'f_vr', 'notes');
save([p 'Attune_IFCB_pico_nano_micro_seasonal_transect_means_latbin'], 'tbl_season_year', 'tbl_season_allVR', 'f_vr', 'notes');

cruiseStat = groupsummary(tbl_cruise,["cruise"], "mean", {'mean_dt'}); %save cruise start date
%save([synchrony_path 'Attune_IFCB_pico_nano_micro_seasonal_transect_means_latbin_pt1'], 'tbl_season_year', 'tbl_season_allVR', 'f_vr', 'notes', 'tblStats_cruise', 'cruiseStat');
save([p 'Attune_IFCB_pico_nano_micro_seasonal_transect_means_latbin_pt1'], 'tbl_season_year', 'tbl_season_allVR', 'f_vr', 'notes', 'tblStats_cruise', 'cruiseStat');

return
%%
%seasonal mean of cruise means
var2plot = strcat('mean_nanmean_', {'phytoC_total' 'phytoC_0_2' 'phytoC_2_20' 'phytoC_20_100'});
varStr = {'Total phyto C' 'Pico C' 'Nano C' 'Micro C'};
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
for v = 1:length(var2plot)
    nexttile
    gs = gscatter(tblStats_cruise1_season2.lat_bin, tblStats_cruise1_season2.(var2plot{v}),tblStats_cruise1_season2.season);
    set(gs, 'linestyle', '-', LineWidth=2)
    set(gca, 'xdir', 'rev')
    ylabel([varStr{v} ' \mug l^{-1}'])
    if v>1
        set(legend, 'visible', 'off')
    else
        set(legend, 'string', sStr)
    end
    yl(v,:) = ylim;
end
title(tl, 'Seasonal means of cruise means')

var2plot = regexprep(var2plot, 'mean_nanmean', 'std_nanmean');
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
for v = 1:length(var2plot)
    nexttile
    gs = gscatter(tblStats_cruise1_season2.lat_bin, tblStats_cruise1_season2.(var2plot{v}),tblStats_cruise1_season2.season);
    set(gs, 'linestyle', '-', LineWidth=2)
    set(gca, 'xdir', 'rev')
    ylabel([varStr{v} ' \mug l^{-1}'])
    if v>1
        set(legend, 'visible', 'off')
    else
        set(legend, 'string', sStr)
    end
    ylim(yl(v,:))
end
title(tl, 'Seasonal standard deviation of cruise means')
%%
if 0
var2plot = {'median_phytoC_total' 'median_phytoC_0_2' 'median_phytoC_2_20' 'median_phytoC_20_100'}
figure
tiledlayout(2,4)
for ii = 1:4
    varStr = var2plot{ii};
    nexttile
    gs = gscatter(tblStats_season.lat_bin, tblStats_season.(varStr),tblStats_season.season);
    set(gs, 'linestyle', '-', LineWidth=2)
    set(gca, 'xdir', 'rev')
    %set(legend, 'string', {'winter' 'spring' 'summer' 'fall'})
end
for ii = 1:4
    varStr = regexprep(var2plot{ii}, 'median', 'std');
    nexttile
    gs = gscatter(tblStats_season.lat_bin, tblStats_season.(varStr),tblStats_season.season);
    set(gs, 'linestyle', '-', LineWidth=2)
    set(gca, 'xdir', 'rev')
    %set(legend, 'string', {'winter' 'spring' 'summer' 'fall'})
end
end
%%
%cruise means by season
var2plot = strcat('nanmean_', {'phytoC_total' 'phytoC_0_2' 'phytoC_2_20' 'phytoC_20_100'});
%var2plot = strcat('nanmean_', {'phytoC_total' 'phytoC_0_3' 'phytoC_3_10' 'phytoC_10_100'});
varStr = {'Total phyto C' 'Pico C' 'Nano C' 'Micro C'};
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
for s = 1:4
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
for v = 1:length(var2plot)
    nexttile
    ii = (tblStats_cruise.season == s); gs = gscatter(tblStats_cruise.lat_bin(ii), tblStats_cruise.(var2plot{v})(ii),tblStats_cruise.cruise(ii));
    set(gs, 'linestyle', '-', LineWidth=2)
    set(gca, 'xdir', 'rev')
    ylabel([varStr{v} ' \mug l^{-1}'])
    if v>1
        set(legend, 'visible', 'off')
    end
end
title(tl, sStr{s})
print(['C:\work\IFCB_products\NESLTER_transect\summary\phyto_size_composition_underway_by_cruise_' sStr{s}], '-dpng')
end

%%
% cruise means vs T or S in latitude bin 
if 0
var2plot = strcat('nanmean_', {'phytoC_total' 'phytoC_0_2' 'phytoC_2_20' 'phytoC_20_100'});
varStr = {'Total phyto C' 'Pico C' 'Nano C' 'Micro C'};
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
for s = 1:4
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
for v = 1:length(var2plot)
    nexttile
    ii = (tblStats_cruise.season == s); gs = gscatter(tblStats_cruise.nanmean_sal(ii), tblStats_cruise.(var2plot{v})(ii),tblStats_cruise.cruise(ii));
   % set(gs, 'linestyle', '-', LineWidth=2)
   % set(gca, 'xdir', 'rev')
    ylabel([varStr{v} ' \mug l^{-1}'])
    if v>1
        set(legend, 'visible', 'off')
    end
end
title(tl, sStr{s})
end
end
%%
%cruise VR by season
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
for s = 1:4
    nexttile
    ii = (tbl_cruise. mean_season == s & tbl_cruiseVR.GroupCount>2); 
    gs = gscatter(tbl_cruiseVR.lat_bin(ii), tbl_cruiseVR.fun1_phytoC_comp_cube_root(ii), tbl_cruiseVR.cruise(ii));
    set(gs, 'linestyle', '-', LineWidth=2)
    set(gca, 'xdir', 'rev')
    title(sStr{s})
    ylim([0 1])
    lh = line(xlim, [1 1]./size(sumTvr.phytoC_comp,2), 'linestyle', '--', 'color', 'k');
    set( get( get( lh, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set(legend, 'fontsize', 6, 'Location', 'best')
end
%ylabel(tl, 'Variance ratio')

%%
%season VR by year
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
for s = 1:4
    nexttile
    ii = (tbl_season_year.season == s & (tbl_season_yearVR.GroupCount>2 | tbl_season_yearVR.GroupCount==0)); 
    gs = gscatter(tbl_season_yearVR.lat_bin(ii), tbl_season_yearVR.fun1_phytoC_comp_cube_root(ii), tbl_season_yearVR.year(ii));
    set(gs, 'linestyle', '-', LineWidth=2)
    set(gca, 'xdir', 'rev')
    title(sStr{s})
    ylim([0 1])
    lh = line(xlim, [1 1]./size(sumTvr.phytoC_comp,2), 'linestyle', '--', 'color', 'k');
    set( get( get( lh, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set(legend, 'fontsize', 6, 'Location', 'best')
end
ylabel(tl, 'Variance ratio')

%%
%seasonal VR of cruise means
varStr = {'VR of cruise means'};
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
%figure
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
nexttile
gs = gscatter(tbl_cruise_seasonVR.lat_bin, tbl_cruise_seasonVR.fun1_mean_phytoC_comp_cube_root,tbl_cruise_seasonVR.mean_season);
set(gs, 'linestyle', '-', LineWidth=2)
set(gca, 'xdir', 'rev')
ylabel(varStr)
set(legend, 'string', sStr)
%yl(v,:) = ylim;
ylim([0 1])
lh = line(xlim, [1 1]./size(sumTvr.phytoC_comp,2), 'linestyle', '--', 'color', 'k');
set( get( get( lh, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
title('Interannual VR by season')

%%
%mean of VR by season
varStr = {'Mean of VR'};
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
%figure
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
nexttile
gs = gscatter(tbl_cruiseVR_season.lat_bin, tbl_cruiseVR_season.mean_fun1_phytoC_comp_cube_root,tbl_cruiseVR_season.mean_season);
set(gs, 'linestyle', '-', LineWidth=2)
set(gca, 'xdir', 'rev')
ylabel(varStr)
set(legend, 'string', sStr)
ylim([0 1])
lh = line(xlim, [1 1]./size(sumTvr.phytoC_comp,2), 'linestyle', '--', 'color', 'k');
set( get( get( lh, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
title('Mean VR by season')

%%
%season means by year
var2plot = 'mean_phytoC_comp';
varStr = {'Total phyto C' 'Pico C' 'Nano C' 'Micro C'};
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
for s = 1:4
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
for v = 1:4 %length(var2plot)
    nexttile
    if v == 1
        X = sum(tbl_season_year.mean_phytoC_comp,2);
    else
        X = tbl_season_year.mean_phytoC_comp(:,v-1);
    end
    ii = (tbl_season_year.season == s); gs = gscatter(tbl_season_year.lat_bin(ii), X(ii),tbl_season_year.year(ii));
    set(gs, 'linestyle', '-', LineWidth=2)
    set(gca, 'xdir', 'rev')
    ylabel([varStr{v} ' \mug l^{-1}'])
    if v>1
        set(legend, 'visible', 'off')
    end
end
title(tl, sStr{s})
end

%%
%seasonal VR of season means
varStr = {'VR of season means'};
sStr = {'Winter' 'Spring' 'Summer' 'Fall'};
%figure
figure('Position',[50 250 1230 290])
tl = tiledlayout(1,4);
nexttile
gs = gscatter(tbl_seasonVR.lat_bin, tbl_seasonVR.fun1_mean_phytoC_comp_cube_root,tbl_seasonVR.season);
set(gs, 'linestyle', '-', LineWidth=2)
set(gca, 'xdir', 'rev')
ylabel(varStr)
set(legend, 'string', sStr)
%yl(v,:) = ylim;
ylim([0 1])
lh = line(xlim, [1 1]/3, 'linestyle', '--', 'color', 'k');
set( get( get( lh, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
title('Interannual VR by season')

%%
% VR of seasonal means
figure
plot(tbl_season_allVR.lat_bin, tbl_season_allVR.fun1_mean_phytoC_comp_cube_root, '.-', 'LineWidth',2, 'MarkerSize',20)
set(gca, 'xdir', 'rev')
ylabel('VR, Loreau')
ylim([0 1])
xlim([39.5 41.5])
lh = line(xlim, [1 1]./3, 'linestyle', '--', 'color', 'k');
set( get( get( lh, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
hold on
plot(tbl_season_allVR.lat_bin, tbl_season_allVR.fun1_mean_phytoC_comp_cube_root13, '.-', 'LineWidth',2, 'MarkerSize',20)
plot(tbl_season_allVR.lat_bin, tbl_season_allVR.fun1_mean_phytoC_comp_cube_root12, '.-', 'LineWidth',2, 'MarkerSize',20)
plot(tbl_season_allVR.lat_bin, tbl_season_allVR.fun1_mean_phytoC_comp_cube_root23, '.-', 'LineWidth',2, 'MarkerSize',20)
legend({'pico-nano-micro' 'pico-micro' 'pico-nano' 'nano-micro'})
title('VR of seasonal means by latitude bin')
%%
% VR of seasonal means
figure
plot(tbl_season_allVR.lat_bin, tbl_season_allVR.fun1_mean_phytoC_comp_cube_root, '.-', 'LineWidth',2, 'MarkerSize',20)
set(gca, 'xdir', 'rev')
ylabel('VR, Loreau')
ylim([0 1])
xlim([39.5 41.5])
lh = line(xlim, 0.3*[1 1], 'linestyle', '--', 'color', 'k');
set( get( get( lh, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
title('VR of pico-nano-micro seasonal means by latitude bin')
%%
figure('position', [350 40 560 600])
tl  = tiledlayout(4,1);
nexttile
gs = gscatter(tbl_season_year.year+tbl_season_year.season/4-.25,tbl_season_year.mean_phytoC_comp_cube_root(:,1), tbl_season_year.lat_bin, winter(length(ix)-1));
set(gs, 'linestyle', '-', LineWidth=2)
ylabel('PicoC (\mug l^{-1})')
lh = legend('Location', 'eastoutside'); delete(lh)
nexttile
gs = gscatter(tbl_season_year.year+tbl_season_year.season/4-.25,tbl_season_year.mean_phytoC_comp_cube_root(:,2), tbl_season_year.lat_bin, winter(length(ix)-1));
set(gs, 'linestyle', '-', LineWidth=2)
ylabel('NanoC (\mug l^{-1})')
lh = legend('Location', 'eastoutside'); delete(lh)
nexttile
gs = gscatter(tbl_season_year.year+tbl_season_year.season/4-.25,tbl_season_year.mean_phytoC_comp_cube_root(:,3), tbl_season_year.lat_bin, winter(length(ix)-1));
set(gs, 'linestyle', '-', LineWidth=2)
ylabel('MicroC (\mug l^{-1})')
lh = legend('Location', 'eastoutside'); delete(lh)
nexttile
gs = gscatter(tbl_season_year.year+tbl_season_year.season/4-.25,sum(tbl_season_year.mean_phytoC_comp,2), tbl_season_year.lat_bin, winter(length(ix)-1));
set(gs, 'linestyle', '-', LineWidth=2)
ylabel('Total C (\mug l^{-1})')
legend('Location', 'eastoutside')
title(tl,'Seasonal means by latitude bin, cube root')

%%
latleg = cellstr(num2str(unique(tbl_season_year.lat_bin(~isnan(tbl_season_year.lat_bin)))))'
nlat = length(ix)-1;
nlat = length(ix);
pico = reshape(tbl_season_year.mean_phytoC_comp(:,1),nlat,size(tbl_season_year,1)/nlat)';
pico(:,end) = []; %this is the nan lat_bin column
%vr_pico_shelf = variance_ratio(pico(~isnan(sum(pico,2)),:));
vr_pico_shelf = loreau_vr(pico(~isnan(sum(pico,2)),:));
nano = reshape(tbl_season_year.mean_phytoC_comp(:,2),nlat,size(tbl_season_year,1)/nlat)';
nano(:,end) = []; %this is the nan lat_bin column
%vr_nano_shelf = variance_ratio(nano(~isnan(sum(nano,2)),:));
vr_nano_shelf = loreau_vr(nano(~isnan(sum(nano,2)),:));
micro = reshape(tbl_season_year.mean_phytoC_comp(:,3),nlat,size(tbl_season_year,1)/nlat)';
micro(:,end) = []; %this is the nan lat_bin column
%vr_micro_shelf = variance_ratio(micro(~isnan(sum(micro,2)),:));
vr_micro_shelf = loreau_vr(micro(~isnan(sum(micro,2)),:));

%total = pico+nano+micro;
total = (pico.^3+nano.^3+micro.^3).^(1/3); %case for cube roots values

figure
tl = tiledlayout(3,1);
nexttile
t = pico;
colororder(winter(nlat))
plot(2018:.25:2023-.25, t,'.-', 'linewidth', 1)
hold on
%plot(2018:.25:2023-.25, sum(t,2),'k.-', 'linewidth',1, 'markersize', 10)
plot(2018:.25:2023-.25, sum(t.^3,2).^(1/3),'k.-', 'linewidth',1, 'markersize', 10)
ylabel('PicoC (\mug l^{-1})')
lh = legend('Location', 'eastoutside'); delete(lh)
nexttile
t = nano;
colororder(winter(nlat))
plot(2018:.25:2023-.25, t,'.-', 'linewidth', 1)
hold on
plot(2018:.25:2023-.25, sum(t.^3,2).^(1/3),'k.-', 'linewidth',1, 'markersize', 10)
ylabel('NanoC (\mug l^{-1})')
lh = legend([latleg 'Total'],'Location', 'eastoutside');
nexttile
t = micro;
colororder(winter(nlat))
plot(2018:.25:2023-.25, t,'.-', 'linewidth', 1)
hold on
plot(2018:.25:2023-.25, sum(t.^3,2).^(1/3),'k.-', 'linewidth',1, 'markersize', 10)
ylabel('MicroC (\mug l^{-1})')
lh = legend('Location', 'eastoutside'); delete(lh)
