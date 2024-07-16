
clear all; close all;
%Load FCB data
load('\\sosiknas1\Lab_data\MVCO\FCB\summary\FCB_compiledC_tables.mat')

%Load Attune data
load('\\sosiknas1\Lab_data\Attune\cruise_data\MVCO\preserved\outputs\SummaryTable.mat')

%load EPICS data?
load('\\sosiknas1\Lab_data\MVCO\EPICS_MVCO_FCM_analysis\wksp_1May2014.mat')
%make EPICS workspace date usable
dateplot = datetime(datevec(matdate), 'Format', 'uuuu-MM-dd HH:mm:ss');


figure
plot(TTsumnum.Time, TTsumnum.Syn./TTsumnum.ml, 'k.', 'MarkerSize', 3)
hold on
plot(CNTable.date_sampled, CNTable.syn_per_ml, 'r.', 'MarkerSize', 8)
plot(dateplot, synconc, 'g.', 'MarkerSize', 8)
legend('FCB', 'Attune', 'EPICS2014');
title('Syn');
ylabel('Syn per ml')
set(gcf, 'Position', [87 600 1200 400])
grid
datetick('keeplimits')

figure
plot(TTsumnum.Time, TTsumnum.Euk_total./TTsumnum.ml, 'k.', 'MarkerSize', 3)
hold on
plot(CNTable.date_sampled, CNTable.euk_per_ml, 'r.', 'MarkerSize', 8)
%plot(dateplot, synconc, 'g.', 'MarkerSize', 8)
legend('FCB', 'Attune');
title('Euks');
ylabel('Euks per ml')
set(gcf, 'Position', [87 100 1200 400])
grid
datetick('keeplimits')

