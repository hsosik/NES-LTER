synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';
load([synchrony_path 'Syn_Gdel_filled'])
class2use = {'Guinardia_delicatula' 'Syn'};

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
%Now get various VR values, Schluter
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
title(class2use, 'interpreter', 'none')

%Now get various VR values, Loreau
ii = ~isnan(sum(Tt{:,class2use},2)); %find the rows with no NaNs
VR_L_total_gaps = loreau_vr(Tt{ii,class2use});
VR_L_total_filled = loreau_vr(Tt_nogaps{:,class2use})
VR_L_doy_effect = loreau_vr(pd_doy{:,class2use})
VR_L_year_effect = loreau_vr(pd_year{:,class2use})

figure
bar([VR_L_total_gaps VR_L_total_filled VR_L_doy_effect VR_L_year_effect])
ylim([0 1])
line(xlim, [1 1]*.5, 'color', 'k', 'linestyle', '--')
set(gca, 'xTickLabel', {'Total (w/gaps)' 'Total (filled)' 'Year day effect' 'Year effect'})
ylabel('VR, Loreau')
title(class2use, 'interpreter', 'none')

%%
%full time-scale specific VR
[phi_sigma,freq] = tsvr_all_sigma(Tt_nogaps{:,class2use});
figure
semilogx(1./freq, phi_sigma)
line(xlim, [1 1], 'color', 'k', 'linestyle', '--')
line([365 365], ylim, 'color', 'r', 'linewidth', 1)
xlabel('Days')
ylabel('VR')
title(class2use, 'interpreter', 'none')

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
title(class2use, 'interpreter', 'none')

return
%% 
%below here is how to load time series data, fit GAM and save gap-filled time series

%load(['C:\work\IFCB_products\MVCO\summary_v4\MVCO_carbon_class_group_summary'])
synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';

load([synchrony_path 'MVCO_carbon_class_group_summary'])
clear Cday Cday_adhoc Cday_group_opt
%fudge for now due to bad bins on this date
%ii = (Cday_opt.Time==datetime('27-feb-2014'));
tt = datetime({'08-jan-2008' '26-dec-2009' '27-dec-2009' '7-jul-2012' '27-feb-2014' '14-jul-2014' '15-jul-2014' '16-jul-2014' '9-Jun-2020' '10-Jun-2020' '11-Jun-2020' '29-Apr-2021' '30-Apr-2021' '1-May-2021' '2-May-2021' '3-May-2021' '4-May-2021' '14-oct-2021' '15-oct-2021' '16-oct-2021'})
%ii = (Cday.Time==datetime('27-feb-2014'));
[~,ii] = intersect(Cday_opt.Time, tt);
Cday_opt(ii,:) = [];

class2use1 = {'Guinardia_delicatula'};
Tt_Cconc1 = timetable(Cday_opt.Time,Cday_opt.(class2use1{1})./Cday_opt.ml/1000, 'VariableNames',class2use1);

synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';
load([synchrony_path 'FCB_compiledC_tables_daily'])
class2use2 = {'Syn'};
Tt_Cconc2 = timetable(TTsumC_day_conc.Time, TTsumC_day_conc.(class2use2{1}), 'VariableNames',class2use2);

class2use = [class2use1 class2use2];
%Tt = retime(Tt_Cconc1, Tt_Cconc2);
Tt = retime(Tt_Cconc1, Tt_Cconc2, Tt_Cconc1.Time); %start at first Gdel time to not extrapolate back for IFCB
Tt = retime(Tt, 'daily'); %no gaps now

Tt{:,:} = (Tt{:,:}).^(1/3); %cube root transform to be closer to normal
Tt.doy = day(Tt.Time,"dayofyear"); %predictor variable 1
Tt.year = year(Tt.Time); %predictor variable 2
Tt = timetable2table(Tt); %ftirgam does not take timetables as input

%fill the gaps with a GAM fit
yfit_nogaps = table; %this will be the whole GAM-predicted time series
for ii = 1:length(class2use)
    disp(class2use(ii))
    m.(class2use{ii}) = fitrgam(Tt,class2use{ii}, 'PredictorNames', {'doy' 'year'}, 'CategoricalPredictors', [2], 'Interactions', 'all', 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions',struct('UseParallel',true), 'Verbose',false, 'MaxPValue', .05, 'HyperparameterOptimizationOptions',struct('ShowPlots',false));
    yfit_nogaps.(class2use{ii}) = predict(m.(class2use{ii}), Tt(:,["doy" "year"]));
end

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

save([synchrony_path 'Syn_Gdel_filled'], 'm', 'pd_doy', 'pd_year', 'yfit_nogaps', 'Tt', 'Tt_nogaps', 'class2use')
disp('saved results:')
disp([synchrony_path 'Syn_Gdel_filled'])
