% Generalized variance ratio calculations and plotting for various MVCO FCB and IFCB data input options
%
% Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2023

synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';
%
%case2plot = 'Syn_Euk2_Euk10';
%case2plot = 'top20diatom';
case2plot = 'Syn_Gdel';
load([synchrony_path case2plot '_filled'])
plot_partials = false;
titleStr = class2use;
switch case2plot
    case 'Syn_Euk2_Euk10'
    case 'top20diatom'
        titleStr = 'Top 20 diatom types at MVCO';
    case 'Syn_Gel'
end

%%

if plot_partials % this is slow to plot so skip if desired
    for ii = 1:length(class2use)
        figure, subplot(2,1,1),plotPartialDependence(m.(class2use{ii}), 'doy')
        subplot(2,1,2), plotPartialDependence(m.(class2use{ii}), 'year')
        print([synchrony_path filesep 'MVCO_phytoplankton_plots' filesep 'GAM_partial_effects' filesep class2use{ii} '_GAM_partial_effects'], '-dpng')
    end
end
%%
for ii = 1:length(class2use)
    figure('position', [50 200 1200 400])
    plot(Tt.Time, yfit_nogaps.(class2use{ii}), '-','linewidth', 1.5)
    hold on
    plot(Tt.Time, Tt.(class2use{ii}), '.')
    ylabel([(class2use{ii}) ' carbon conc, cube root'], 'interpreter', 'none')
    set(gca, 'xtick', datetime(2003:2023,1,1), 'xgrid', 'on')
    print([synchrony_path filesep 'MVCO_phytoplankton_plots' filesep 'GAM_fits2timeseries' filesep class2use{ii} '_GAM_fit2timeseries'], '-dpng')   
end

%%
%Now get various VR values
ii = ~isnan(sum(Tt{:,class2use},2)); %find the rows with no NaNs
VR_total_gaps = variance_ratio(Tt{ii,class2use});
VR_total_filled = variance_ratio(Tt_nogaps{:,class2use})
VR_doy_effect = variance_ratio(pd_doy{:,class2use})
VR_year_effect = variance_ratio(pd_year{:,class2use})
figure
bar([VR_total_gaps VR_total_filled VR_doy_effect VR_year_effect])
line(xlim, [1 1], 'color', 'k', 'linestyle', '--')
ylim([0 length(class2use)])
set(gca, 'xTickLabel', {'Total (w/gaps)' 'Total (filled)' 'Year day effect' 'Year effect'})
ylabel('VR, Schluter')
title(titleStr, 'interpreter', 'none')
print([synchrony_path filesep 'MVCO_phytoplankton_plots' filesep char(join(titleStr, '--')) '_VR_Schluter'], '-dpng')   

%Now get various VR values, Loreau
ii = ~isnan(sum(Tt{:,class2use},2)); %find the rows with no NaNs
VR_L_total_gaps = loreau_vr(Tt{ii,class2use});
VR_L_total_filled = loreau_vr(Tt_nogaps{:,class2use})
VR_L_doy_effect = loreau_vr(pd_doy{:,class2use})
VR_L_year_effect = loreau_vr(pd_year{:,class2use})

figure
bar([VR_L_total_gaps VR_L_total_filled VR_L_doy_effect VR_L_year_effect])
ylim([0 1])
line(xlim, [1 1]/length(class2use), 'color', 'k', 'linestyle', '--')
set(gca, 'xTickLabel', {'Total (w/gaps)' 'Total (filled)' 'Year day effect' 'Year effect'})
ylabel('VR, Loreau')
title(titleStr, 'interpreter', 'none')
print([synchrony_path filesep 'MVCO_phytoplankton_plots' filesep char(join(titleStr, '--')) '_VR_Loreau'], '-dpng')   

%%
%full time-scale specific VR
[phi_sigma,freq] = tsvr_all_sigma(Tt_nogaps{:,class2use});
figure
semilogx(1./freq, phi_sigma)
line(xlim, [1 1], 'color', 'k', 'linestyle', '--')
line([365 365], ylim, 'color', 'r', 'linewidth', 1)
xlabel('Days')
ylabel('VR')
title(titleStr, 'interpreter', 'none')
print([synchrony_path filesep 'MVCO_phytoplankton_plots' filesep char(join(titleStr, '--')) '_tsvr'], '-dpng')   

%%
%now do some time-scale specific ranges
thresholds = [7 30 90 180 365]; %day edges: 1-7, 7-30..., >365
[vr_scales] = tsvr_scale(Tt_nogaps{:,class2use},thresholds);
figure
bar(vr_scales)
line(xlim, [1 1], 'color', 'k', 'linestyle', '--')
ylabel('VR')
xlabel('Day range')
set(gca, 'xTickLabel', {'<7' '7-30' '30-90' '90-180' '180-365' '>365'})
title(titleStr, 'interpreter', 'none')
print([synchrony_path filesep 'MVCO_phytoplankton_plots' filesep char(join(titleStr, '--')) '_tsvr_ranges'], '-dpng')   
