% Let's try an example of variance ratio and time-scale dependent variance
% ratios for three components of the phytoplankton in the MVCO time series
% Syn, Euk<=2microns (pico), Euk2-10microns (small nano)
%
% Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2023

synchrony_path = 'C:\Users\heidi\Woods Hole Oceanographic Institution\Michael G Neubert - NES_Synchrony_Data\';
load([synchrony_path 'FCB_compiledC_tables_daily'])

%subtract to get the 2-10 micron Euk size fraction
TTsumC_day_conc.Euk2_10microns = TTsumC_day_conc.Eukleq10microns-TTsumC_day_conc.Eukleq2microns;

class2use = {'Syn' 'Eukleq2microns' 'Euk2_10microns'};
%make a temporary table for the VR calculations
Tt = TTsumC_day_conc;
Tt{:,:} = (Tt{:,:}).^(1/3); %cube root transform to be closer to normal
Tt.doy = day(Tt.Time,"dayofyear"); %predictor variable 1
Tt.year = year(Tt.Time); %predictor variable 2
Tt = timetable2table(Tt); %ftirgam does not take timetables as input

%fill the gaps with a GAM fit
for ii = 1:length(class2use)
    disp(class2use(ii))
    m.(class2use{ii}) = fitrgam(Tt,class2use{ii}, 'PredictorNames', {'doy' 'year'}, 'Interactions', 'all', 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions',struct('UseParallel',true), 'Verbose',false, 'MaxPValue', .05);
    yfit_nogaps.(class2use{ii}) = predict(m.(class2use{ii}), Tt(:,["doy" "year"]));
end
%here's the whole gap-filled time series for VR input
yfit_nogaps = struct2table(yfit_nogaps);
%%
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
VR_total = variance_ratio(yfit_nogaps{:,:})
VR_doy_effect = variance_ratio(pd_doy{:,class2use})
VR_year_effect = variance_ratio(pd_year{:,class2use})

figure
bar([VR_total VR_doy_effect VR_year_effect])
line(xlim, [1 1], 'color', 'k', 'linestyle', '--')
set(gca, 'xTickLabel', {'Total' 'Year day effect' 'Year effect'})
ylabel('VR')
title(class2use, 'interpreter', 'none')

%%
%full time-scale specific VR
[phi_sigma,freq] = tsvr_all_sigma(yfit_nogaps{:,:}')
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
[vr_scales] = tsvr_scale(yfit_nogaps{:,:}',thresholds)
figure
bar(vr_scales)
line(xlim, [1 1], 'color', 'k', 'linestyle', '--')
ylabel('VR')
xlabel('Day range')
set(gca, 'xTickLabel', {'<7' '7-30' '30-90' '90-180' '180-365' '>365'})
title(class2use, 'interpreter', 'none')