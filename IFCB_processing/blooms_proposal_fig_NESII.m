%synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';
p = '\\sosiknas1\Lab_data\Attune\cruise_data\IFCB_Attune_merge\summary_files\';
%save this from phytoC_transect_stats with latbin step = 0.1;
load([p 'Attune_IFCB_pico_nano_micro_seasonal_transect_means_latbin_pt1'])

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
%print(['C:\work\IFCB_products\NESLTER_transect\summary\phyto_size_composition_underway_by_cruise_' sStr{s}], '-dpng')
end

%%
% figure for NES II proposal
cruiseStat.newStr=cellstr(datestr(cruiseStat.mean_mean_dt, 'yyyy-mmm'));
varStr2 = regexprep(varStr, ' C', 'plankton C')
figure('Position',[50 250 900 230])
tl = tiledlayout(1,3);
nexttile
s = 2; %spring
v = 3; %nano
ii = (tblStats_cruise.season == s); gs = gscatter(tblStats_cruise.lat_bin(ii), tblStats_cruise.(var2plot{v})(ii),tblStats_cruise.cruise(ii));
set(gs, 'linestyle', '-', LineWidth=2)
set(gca, 'xdir', 'rev')
ylabel([varStr2{v} ' \mug l^{-1}'])
Str = get(legend, 'string');
[~, ia] = ismember(Str,cruiseStat.cruise);
set(legend, 'string', cruiseStat.newStr(ia),'fontsize', 7);
nexttile
s = 3; %summer
v = 4; %micro
ii = (tblStats_cruise.season == s); gs = gscatter(tblStats_cruise.lat_bin(ii), tblStats_cruise.(var2plot{v})(ii),tblStats_cruise.cruise(ii));
set(gs, 'linestyle', '-', LineWidth=2)
set(gca, 'xdir', 'rev')
ylabel([varStr2{v} ' \mug l^{-1}'])
ylim([0 100])
Str = get(legend, 'string');
[~, ia] = ismember(Str,cruiseStat.cruise);
set(legend, 'string', cruiseStat.newStr(ia), 'fontsize', 7);
nexttile
s = 2; %spring
v = 4; %nano
ii = (tblStats_cruise.season == s); gs = gscatter(tblStats_cruise.lat_bin(ii), tblStats_cruise.(var2plot{v})(ii),tblStats_cruise.cruise(ii));
set(gs, 'linestyle', '-', LineWidth=2)
set(gca, 'xdir', 'rev')
ylabel([varStr2{v} ' \mug l^{-1}'])
ylim([0 100])
Str = get(legend, 'string');
[~, ia] = ismember(Str,cruiseStat.cruise);
set(legend, 'string', cruiseStat.newStr(ia),'fontsize', 7);
%print(['C:\work\IFCB_products\NESLTER_transect\summary\phyto_size_composition_underway_by_cruise_NESII_proposal'], '-dpng', '-r600')
