synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';
load([synchrony_path 'top20diatom_filled'])
%%
if true % this is slow to plot so skip if desired
    for ii = 1:length(class2use)
        figure, subplot(2,1,1),plotPartialDependence(m.(class2use{ii}), 'doy')
        subplot(2,1,2), plotPartialDependence(m.(class2use{ii}), 'year')
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
set(gca, 'xTickLabel', {'Total (w/gaps)' 'Total (filled)' 'Year day effect' 'Year effect'})
ylabel('VR, Schluter')
%title(class2use, 'interpreter', 'none')
title('Top 20 diatom types at MVCO')

%Now get various VR values
ii = ~isnan(sum(Tt{:,class2use},2)); %find the rows with no NaNs
VR_total_gaps = loreau_vr(Tt{ii,class2use});
VR_total_filled = loreau_vr(Tt_nogaps{:,class2use})
VR_doy_effect = loreau_vr(pd_doy{:,class2use})
VR_year_effect = loreau_vr(pd_year{:,class2use})

figure
bar([VR_total_gaps VR_total_filled VR_doy_effect VR_year_effect])
ylim([0 1])
line(xlim, [1 1]*.5, 'color', 'k', 'linestyle', '--')
set(gca, 'xTickLabel', {'Total (w/gaps)' 'Total (filled)' 'Year day effect' 'Year effect'})
ylabel('VR, Loreau')
title('Top 20 diatom types at MVCO')


%%
%pairwise VR matrices
VR_mat = NaN(length(class2use)+1);
VR_mat_doy = VR_mat;
VR_mat_year = VR_mat;
for ia = 1:length(class2use)
    for ib = 1:length(class2use)
        VR_mat(ia,ib) = variance_ratio(Tt_nogaps{:,class2use([ia ib])});
        VR_mat_doy(ia,ib) = variance_ratio(pd_doy{:,class2use([ia ib])}); 
        VR_mat_year(ia,ib) = variance_ratio(pd_year{:,class2use([ia ib])});
    end
end
figure
pcolor(VR_mat)
set(gca, 'YTick',1.5:20.5, 'YTickLabel' ,regexprep(class2use, '_', ' '))
set(gca, 'xTick',1:20, 'xTickLabel' ,regexprep(class2use, '_', ' '), 'XTickLabelRotationMode', 'auto')
axis square
colorbar
caxis([0 2])
title('Total VR (no gaps)')
figure
pcolor(VR_mat_doy)
set(gca, 'YTick',1.5:20.5, 'YTickLabel' ,regexprep(class2use, '_', ' '))
set(gca, 'xTick',1:20, 'xTickLabel' ,regexprep(class2use, '_', ' '), 'XTickLabelRotationMode', 'auto')
axis square
colorbar
caxis([0 2])
title('VR year day effect')
figurei
pcolor(VR_mat_year)
set(gca, 'YTick',1.5:20.5, 'YTickLabel' ,regexprep(class2use, '_', ' '))
set(gca, 'xTick',1:20, 'xTickLabel' ,regexprep(class2use, '_', ' '), 'XTickLabelRotationMode', 'auto')
axis square
colorbar
caxis([0 2])
title('VR year effect')
%%
%%
%full time-scale specific VR
[phi_sigma,freq] = tsvr_all_sigma(Tt_nogaps{:,class2use});
figure
semilogx(1./freq, phi_sigma)
line(xlim, [1 1], 'color', 'k', 'linestyle', '--')
line([365 365], ylim, 'color', 'r', 'linewidth', 1)
xlabel('Days')
ylabel('VR')
%title(class2use, 'interpreter', 'none')
title('Top 20 diatom types at MVCO')

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
%title(class2use, 'interpreter', 'none')
title('Top 20 diatom types at MVCO')
return

%below here is how to load time series data, fit GAM and save gap-filled time series
%%
% from MVCO_compile_hourly_daily_from_class_tables.m
load(['C:\work\IFCB_products\MVCO\summary_v4\MVCO_carbon_class_group_summary'])

Cday{:,2:end} = Cday{:,2:end}./Cday.ml/1000;
Cday_opt{:,2:end} = Cday_opt{:,2:end}./Cday_opt.ml/1000;
Cday_adhoc{:,2:end} = Cday_adhoc{:,2:end}./Cday_adhoc.ml/1000;
Cday_group_opt{:,2:end} = Cday_group_opt{:,2:end}./Cday_group_opt.ml/1000;

group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
[~,ia,ib] = intersect(group_table.CNN_classlist, Cday.Properties.VariableNames);
diatom_ind = ib(find(group_table.Diatom_noDetritus(ia)));

t = nanmean(Cday_opt{:,diatom_ind});
[~,si] = sort(t, 'descend');
Ttop = table;
Ttop.index = (1:length(diatom_ind))';
Ttop.class = Cday_opt.Properties.VariableNames(diatom_ind(si))';
Ttop.meanC_opt = t(si)';
t = nanmean(Cday_adhoc{:,diatom_ind});
Ttop.meanC_adhoc = t(si)';

top20diatom = Ttop.class(1:20);
Cconc_day_top20diatom = Cday(:,top20diatom);
top20diatom = regexprep(top20diatom, '-', '');
class2use = top20diatom';
Cconc_day_top20diatom = renamevars(Cconc_day_top20diatom, 'Pseudo-nitzschia', 'Pseudonitzschia');

synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';

%make a temporary table for the VR calculations
Tt = Cconc_day_top20diatom;
Tt{:,:} = (Tt{:,:}).^(1/3); %cube root transform to be closer to normal
Tt.doy = day(Tt.Time,"dayofyear"); %predictor variable 1
Tt.year = year(Tt.Time); %predictor variable 2
Tt = retime(Tt, 'daily'); %no gaps now
Tt = timetable2table(Tt); %ftirgam does not take timetables as input

%fill the gaps with a GAM fit
for ii = 1:length(class2use)
    disp(class2use(ii))
    m.(class2use{ii}) = fitrgam(Tt,class2use{ii}, 'PredictorNames', {'doy' 'year'},'CategoricalPredictors', [2],'Interactions', 'all', 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions',struct('UseParallel',true), 'Verbose',false, 'MaxPValue', .05, 'HyperparameterOptimizationOptions',struct('ShowPlots',false));
    yfit_nogaps.(class2use{ii}) = predict(m.(class2use{ii}), Tt(:,["doy" "year"]));
end
%here's the whole predicted time series
yfit_nogaps = struct2table(yfit_nogaps);

%now make a gap-filled time series for VR input
Tt_nogaps = Tt(:,['Time' class2use]);
for ii = 1:length(class2use)
    iii = isnan(Tt.(class2use{ii}));
    Tt_nogaps.(class2use{ii})(iii) = yfit_nogaps.(class2use{ii})(iii);
end

%save the GAM-derived partial effects for day of year and year
%these could be VR inputs for different scales
pd_doy = table;
pd_year = table;
for ii = 1:length(class2use)
    disp(class2use(ii))
    [pd, x_doy] = partialDependence(m.(class2use{ii}), ["doy"]);
    pd_doy.(class2use{ii}) = pd';
    [pd,x_year] = partialDependence(m.(class2use{ii}), ["year"]);
    pd_year.(class2use{ii}) = pd';
end
pd_doy.doy = x_doy;
pd_year.year = x_year;

%
save([synchrony_path 'top20diatom_filled'], 'Cconc_day_top20diatom', 'm', 'pd_doy', 'pd_year', 'yfit_nogaps', 'Tt', 'Tt_nogaps', 'class2use')
