load('\\sosiknas1\Lab_data\MVCO\FCB\summary\FCB_compiledC_tables.mat')
%remove rows with NaN for Syn
TTsumC(isnan(TTsumC.Syn),:) = [];
TTsumC_day = retime(TTsumC, 'daily', 'sum');

TTsumC_day_conc = TTsumC_day;
TTsumC_day_conc{:,:} = TTsumC_day{:,:}./TTsumC_day.ml/1000; %carbon conc microgram per ml
TTsumC_day_conc = removevars(TTsumC_day_conc,'ml');

%climatology by yearday
TTsumC_day_conc.doy = day(TTsumC_day_conc.Time, 'dayofyear');
TTsumC_conc_yearday = groupsummary(TTsumC_day_conc,'doy',"mean");
figure
plot(TTsumC_conc_yearday.doy, TTsumC_conc_yearday.mean_Syn, '.-')

%class2use = {"Syn" "("Euk<=2microns"")" "Euk<=3microns" "Euk<=10microns"};
%%
class2use = {'Syn' '("Euk<=2microns")' '("Euk2_10microns")'};
Tt = TTsumC_day_conc;
Tt{:,:} = (Tt{:,:}).^(1/3);
Tt.doy = day(Tt.Time,"dayofyear");
Tt.year = year(Tt.Time);
Tt = timetable2table(Tt);
clear ypred ypred_noInt m

for ii = 2:length(class2use)
    disp(class2use(ii))
    m.(class2use{ii}) = fitrgam(Tt,class2use{ii}, 'PredictorNames', {'doy' 'year'}, 'Interactions', 'all', 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions',struct('UseParallel',true), 'Verbose',false, 'MaxPValue', .05);
%    m_noInt.(class2use{ii}) = fitrgam(Tt,class2use{ii}, 'PredictorNames', {'doy' 'year'}, 'Interactions', 0, 'OptimizeHyperparameters', 'all-univariate', 'HyperparameterOptimizationOptions',struct('UseParallel',true), 'Verbose',false, 'MaxPValue', .05);
    ypred.(class2use{ii}) = predict(m.(class2use{ii}), Tt(:,["doy" "year"]));
    ypred_noInt.(class2use{ii}) = predict(m.(class2use{ii}), Tt(:,["doy" "year"]),'IncludeInteractions',false);
end
ypred = struct2table(ypred);
ypred_noInt = struct2table(ypred_noInt);

%%
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

figure, subplot(2,1,1),plotPartialDependence(m.(class2use{ii}), 'doy'), subplot(2,1,2), plotPartialDependence(m.(class2use{ii}), 'year')


%%
figure
plot(Tt.Time, Tt.Syn, '.')
hold on
plot(Tt.Time, ypred.Syn, '-')
plot(Tt.Time, ypred_noInt.Syn, '-')
