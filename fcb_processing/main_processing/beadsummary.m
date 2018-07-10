%beadpath = 'c:\work\syn_lab\data\processed\beads\';

if exist([beadpath  'beadresults.mat']),
    eval(['delete ' beadpath  'beadresults.mat'])
end;

filelist = dir([savepath '*beads*.*']);

date = datenum(cat(1,filelist.date));
[temp, fileorder] = sort(date);
clear temp date

filelist = filelist(fileorder);

beadsall = [];
for filenum = 1:length(filelist),
    load([beadpath filelist(filenum).name])
    beadsall = [beadsall; beadresults];
    clear beadresults
end;

t = find(beadsall(:,1));
beadsall = beadsall(t,:);
beadresults = beadsall;
clear file* datapath ans t

nind = find(isnan(beadresults(:,20)));  %these are old format files without syringe volume
%beadresults(nind,1) = beadresults(nind,1)/3600/24;  %seconds to days

figure(99), clf
subplot(321)
plot(beadresults(:,1), beadresults(:,5), '.')
%tax = axis;
%text(tax(1)*1.1, tax(3)*1.2, [num2str((max(beadresults(:,5)) - min(beadresults(:,5)))./max(beadresults(:,5))*100) '% drop'])
ylabel('PE')
title('Bead summary')
datetick('x', 6, 'keepticks', 'keeplimits')
subplot(322)
plot(beadresults(:,1), beadresults(:,6), 'r.')
ylabel('FLS')
datetick('x', 6, 'keepticks', 'keeplimits')
subplot(323)
plot(beadresults(:,1), beadresults(:,7), 'g.')
ylabel('CHL')
xlabel('Yearday')
datetick('x', 6, 'keepticks', 'keeplimits')
subplot(324)
plot(beadresults(:,1), beadresults(:,8), 'm.')
ylabel('SSC')
xlabel('Yearday')
datetick('x', 6, 'keepticks', 'keeplimits')
subplot(325)
plot(beadresults(:,1), beadresults(:,4)./(beadresults(:,3)*0.05/60), 'b.')  %time based
hold on
plot(beadresults(:,1), beadresults(:,4)./beadresults(:,20), 'ro') %vol based
ylabel('Beads ml^{-1}')
xlabel('Yearday')
datetick('x', 6, 'keepticks', 'keeplimits')
if classplotflag
    pause
end
eval(['save ' beadpath  'beadresults beadresults beadtitles'])