
clear all; close all;
%Load FCB data
load('\\sosiknas1\Lab_data\MVCO\FCB\summary\FCB_compiledC_tables.mat')

%Load Attune data
load('\\sosiknas1\Lab_data\Attune\cruise_data\MVCO\preserved\outputs\SummaryTable.mat')

%or 
AttuneT = readtable('\\sosiknas1\Lab_data\Attune\EDI_data_packages\MVCO_FCMdiscrete\nes-lter-fcm-discrete-mvco.csv');

%load EPICS data?
%load('\\sosiknas1\Lab_data\MVCO\EPICS_MVCO_FCM_analysis\wksp_1May2014.mat')
%make EPICS workspace date usable
%dateplot = datetime(datevec(matdate), 'Format', 'uuuu-MM-dd HH:mm:ss');
load('\\sosiknas1\Lab_data\MVCO\EPICS_MVCO_FCM_analysis\classify_FCMs.mat')
tempEPICSdate = cell2mat(FCMs_classified(:,3));
dateplot = datetime(datevec(tempEPICSdate), 'Format', 'uuuu-MM-dd HH:mm:ss');
EPICSsynconc=cell2mat(FCMs_classified(:,6))./(cell2mat(FCMs_classified(:,9))*0.1);
EPICSeukconc=cell2mat(FCMs_classified(:,7))./(cell2mat(FCMs_classified(:,9))*0.1);

% temp_month = month(CNTable.date_sampled);
% summer = find(temp_month == 7 | temp_month == 8 | temp_month == 9);
% winter = find(temp_month == 1 | temp_month == 2 | temp_month == 3);
% spring = find(temp_month == 4 | temp_month == 5 | temp_month == 6);
% fall = find(temp_month == 10 | temp_month == 11 | temp_month == 12);

figure
plot(TTsumnum.Time, TTsumnum.Syn./TTsumnum.ml,'.', 'Color',[.5 .5 .5], 'MarkerSize', 2)
hold on
plot(AttuneT.date_sampled, AttuneT.syn_cells_per_ml, '.', 'Color',[1 0 0], 'MarkerSize', 9)
plot(dateplot, EPICSsynconc,'.', 'Color', [.6 .2 .2], 'MarkerSize', 9)
plot(miraRun.datetime, miraRun.syn_conc, '.', 'Color', [1 .5 .5], 'MarkerSize', 9)
legend('FCB', 'Attune by Emily', 'EPICS2014', 'Attune by Mira');
title('Syn');
ylabel('Syn per ml')
set(gcf, 'Position', [87 600 1200 400])
grid
datetick('keeplimits')

figure
plot(TTsumnum.Time, TTsumnum.Euk_total./TTsumnum.ml, '.', 'Color',[.5 .5 .5], 'MarkerSize', 2)
hold on
plot(AttuneT.date_sampled, AttuneT.redeuk_leq_20um_cells_per_ml, '.', 'Color',[0 1 0], 'MarkerSize', 9)
plot(dateplot, EPICSeukconc, '.', 'Color',[.2 .6 .4], 'MarkerSize', 9)
legend('FCB', 'Attune', 'EPICS2014');
title('Euks');
ylabel('Euks per ml')
set(gcf, 'Position', [87 100 1200 400])
grid
datetick('keeplimits')

figure

hold on
plot(miraRun.datetime, miraRun.hbac_conc, 'Marker','.', 'Color', [.2 .5 1], 'MarkerSize', 9)
plot(CNTable.date_sampled, CNTable.bac_per_ml,'Marker', '.', 'Color',[0 0 1], 'MarkerSize', 9)
legend('Attube by Mira', 'Attune by Emily');
title('Bacteria');
ylabel('Bacteria per ml')
set(gcf, 'Position', [87 100 1200 400])
grid
datetick('keeplimits')

figure
subplot(3,1,3); plot(CNTable.date_sampled, CNTable.bac_per_ml,'Marker', '.', 'Color',[0 0 1], 'MarkerSize', 9)
hold on
title('Bacteria')
ylabel('Heterotrophic bactera per ml')
set(gcf, 'Position', [87 100 1200 400])
grid
datetick('keeplimits')
xlim([AttuneT.date_sampled(13) AttuneT.date_sampled(end)])
hold on
subplot(3,1,1); 
plot(TTsumnum.Time, TTsumnum.Euk_total./TTsumnum.ml, '.', 'Color',[.4 .4 .4], 'MarkerSize', 3)
hold on
plot(AttuneT.date_sampled, AttuneT.redeuk_leq_20um_cells_per_ml, 'Marker', '.', 'Color',[0 1 0], 'MarkerSize', 9)
title('Euks');
ylabel('Euks per ml')
grid
datetick('keeplimits')
xlim([AttuneT.date_sampled(13) AttuneT.date_sampled(end)])
subplot(3,1,2); 
hold on
plot(TTsumnum.Time, TTsumnum.Syn./TTsumnum.ml,'.', 'Color',[.4 .4 .4], 'MarkerSize', 3)
plot(AttuneT.date_sampled, AttuneT.syn_cells_per_ml, 'Color',[1 0 0], 'Marker', '.', 'MarkerSize', 9)
title('Syn');
ylabel('Syn per ml')
grid
datetick('keeplimits')
xlim([AttuneT.date_sampled(13) AttuneT.date_sampled(end)])

%plots without MVCO data
% figure
% plot(TTsumnum.Time, TTsumnum.Syn./TTsumnum.ml,'.', 'Color',[.5 .5 .5], 'MarkerSize', 2)
% hold on
% plot(CNTable.date_sampled, CNTable.syn_per_ml, '.', 'Color',[1 0 0], 'MarkerSize', 9)
% plot(dateplot, EPICSsynconc,'.', 'Color', [.6 .2 .2], 'MarkerSize', 9)
% legend('FCB', 'Attune', 'EPICS2014');
% title('Syn');
% ylabel('Syn per ml')
% set(gcf, 'Position', [87 600 1200 400])
% grid
% datetick('keeplimits')
% 
% figure
% plot(TTsumnum.Time, TTsumnum.Euk_total./TTsumnum.ml, '.', 'Color',[.5 .5 .5], 'MarkerSize', 2)
% hold on
% plot(CNTable.date_sampled, CNTable.euk_per_ml, '.', 'Color',[0 1 0], 'MarkerSize', 9)
% plot(dateplot, EPICSeukconc, '.', 'Color',[.2 .6 .4], 'MarkerSize', 9)
% legend('FCB', 'Attune', 'EPICS2014');
% title('Euks');
% ylabel('Euks per ml')
% set(gcf, 'Position', [87 100 1200 400])
% grid
% datetick('keeplimits')
% 
% figure
% subplot(3,1,3); plot(CNTable.date_sampled, CNTable.bac_per_ml,'Marker', '.', 'Color',[0 0 1], 'MarkerSize', 9)
% title('Bacteria')
% ylabel('Heterotrophic bactera per ml')
% set(gcf, 'Position', [87 100 1200 400])
% grid
% datetick('keeplimits')
% xlim([CNTable.date_sampled(11) CNTable.date_sampled(end)])
% hold on
% subplot(3,1,1); plot(CNTable.date_sampled, CNTable.euk_per_ml, 'Marker', '.', 'Color',[0 1 0], 'MarkerSize', 9)
% title('Euks');
% ylabel('Euks per ml')
% grid
% datetick('keeplimits')
% xlim([CNTable.date_sampled(11) CNTable.date_sampled(end)])
% subplot(3,1,2); plot(CNTable.date_sampled, CNTable.syn_per_ml, 'Color',[1 0 0], 'Marker', '.', 'MarkerSize', 9)
% title('Syn');
% ylabel('Syn per ml')
% grid
% datetick('keeplimits')
% xlim([CNTable.date_sampled(11) CNTable.date_sampled(end)])
% 



