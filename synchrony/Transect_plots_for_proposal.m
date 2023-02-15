synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';
load([synchrony_path 'Attune_pico_nano_micro_seasonal_transect_means'])

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
lat_bins = unique(tbl_season_year.lat_bin);
lat_bins(isnan(lat_bins)) = [];
bin2use = 8; %PICK a bin: 1 = most offshore, 8 = most inshore
ilat = find(tbl_season_year.lat_bin==lat_bins(bin2use));
figure
plot(tbl_season_year.year(ilat)+(tbl_season_year.season(ilat)-1)/4, tbl_season_year.mean_phytoC_comp_cube_root(ilat,:), '.-', 'linewidth', 2, 'markersize', 20)
ylabel('Carbon \mug l^{-1}, cube root')
legend('Pico','Nano', 'Micro')
title([num2str(lat_bins(bin2use)) '\circN'])