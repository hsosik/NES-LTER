
clear all; close all;
%Load FCB data
load('\\sosiknas1\Lab_data\MVCO\FCB\summary\FCB_compiledC_tables.mat')

%Load Attune data
load('\\sosiknas1\Lab_data\Attune\cruise_data\MVCO\preserved\outputs\SummaryTable.mat')

%load EPICS data?
%load('\\sosiknas1\Lab_data\MVCO\EPICS_MVCO_FCM_analysis\wksp_1May2014.mat')
%make EPICS workspace date usable
%dateplot = datetime(datevec(matdate), 'Format', 'uuuu-MM-dd HH:mm:ss');
load('\\sosiknas1\Lab_data\MVCO\EPICS_MVCO_FCM_analysis\classify_FCMs.mat')
tempEPICSdate = cell2mat(FCMs_classified(:,3));
dateplot = datetime(datevec(tempEPICSdate), 'Format', 'uuuu-MM-dd HH:mm:ss');
EPICSsynconc=cell2mat(FCMs_classified(:,6))./(cell2mat(FCMs_classified(:,9))*0.1);
EPICSeukconc=cell2mat(FCMs_classified(:,7))./(cell2mat(FCMs_classified(:,9))*0.1);

figure
plot(TTsumnum.Time, TTsumnum.Syn./TTsumnum.ml,'.', 'Color',[.5 .5 .5], 'MarkerSize', 2)
hold on
plot(CNTable.date_sampled, CNTable.syn_per_ml, '.', 'Color',[1 0 0], 'MarkerSize', 9)
plot(dateplot, EPICSsynconc,'.', 'Color', [.6 .2 .2], 'MarkerSize', 9)
legend('FCB', 'Attune', 'EPICS2014');
title('Syn');
ylabel('Syn per ml')
set(gcf, 'Position', [87 600 1200 400])
grid
datetick('keeplimits')

figure
plot(TTsumnum.Time, TTsumnum.Euk_total./TTsumnum.ml, '.', 'Color',[.5 .5 .5], 'MarkerSize', 2)
hold on
plot(CNTable.date_sampled, CNTable.euk_per_ml, '.', 'Color',[0 1 0], 'MarkerSize', 9)
plot(dateplot, EPICSeukconc, '.', 'Color',[.2 .6 .4], 'MarkerSize', 9)
legend('FCB', 'Attune', 'EPICS2014');
title('Euks');
ylabel('Euks per ml')
set(gcf, 'Position', [87 100 1200 400])
grid
datetick('keeplimits')

