function [datmerged, fit] = fcbmergeproc(dat, fittouse)
%renamed to fcbmergeproc (from cytosubproc3), 3/03 Heidi
%modified to accept input fit for case of beads without good overlap
PEcut = 100;
FLScut = 100;
CHLcut = 50;
SSCcut = 50;
if isnan(fittouse)
    t = find(dat(:,1) > 0 & dat(:,1) < PEcut & dat(:,2) > 0);
    %fitPE  = [polyfit(dat(t,1), dat(t,2), 1) length(t)]; %1 = slope, 2 = int, 3 = count
    if length(t) > 3, temp = corrcoef(dat(t,1), dat(t,2)); fitPE  = [polyfit(dat(t,1), dat(t,2), 1) length(t) temp(2)]; end; %1 = slope, 2 = int, 3 = count, 4 = r
    t = find(dat(:,3) > 0 & dat(:,3) < FLScut & dat(:,4) > 0);
    if length(t) > 3, temp = corrcoef(dat(t,3), dat(t,4)); fitFLS = [polyfit(dat(t,3), dat(t,4), 1) length(t) temp(2)];  end;
    %fitFLS = [polyfit(dat(t,3), dat(t,4), 1) length(t)]; 
    t = find(dat(:,5) > 0 & dat(:,5) < CHLcut & dat(:,6) > 0);
    if length(t) > 3, temp = corrcoef(dat(t,5), dat(t,6)); fitCHL = [polyfit(dat(t,5), dat(t,6), 1) length(t) temp(2)];  end;
    %fitCHL = [polyfit(dat(t,5), dat(t,6), 1) length(t)]; 
    t = find(dat(:,7) > 0 & dat(:,7) < SSCcut & dat(:,8) > 0);
    if length(t) > 3, temp = corrcoef(dat(t,7), dat(t,8)); fitSSC = [polyfit(dat(t,7), dat(t,8), 1) length(t) temp(2)];   end;
    %fitSSC = [polyfit(dat(t,7), dat(t,8), 1) length(t)]; 
else
    fitPE = [fittouse(1,:) NaN NaN];
    fitFLS = [fittouse(2,:) NaN NaN];
    fitCHL = [fittouse(3,:) NaN NaN];
    fitSSC = [fittouse(4,:) NaN NaN];
end;

%use half of cut value as merge point
datmerged(:,1) = dat(:,1)*fitPE(1)+fitPE(2);
t = find(dat(:,1) < PEcut/2 & dat(:,2) < PEcut/2*fitPE(1)+fitPE(2));  %unambiguously offscale low on low gain...
datmerged(t,1) = dat(t,2);  %reassign stuff on baseline to be high gain
t = find((dat(:,1) <= PEcut/2 & dat(:,2) >= PEcut/2*fitPE(1)+fitPE(2)) | (dat(:,1) >= PEcut/2 & dat(:,2) <= PEcut/2*fitPE(1)+fitPE(2)));  %points with one signal above and one signal below merge point
indtemp = 1:2:length(t);  %odd indices
datmerged(t(indtemp),1) = dat(t(indtemp),2); %"dither" by reassigning every other ambiguous point

datmerged(:,2) = dat(:,3)*fitFLS(1)+fitFLS(2);
t = find(dat(:,3) < FLScut/2 & dat(:,4) < FLScut/2*fitFLS(1)+fitFLS(2));  %unambiguously offscale low on low gain...
datmerged(t,2) = dat(t,4);  %reassign stuff on baseline to be high gain
t = find((dat(:,3) <= FLScut/2 & dat(:,4) >= FLScut/2*fitFLS(1)+fitFLS(2)) | (dat(:,3) >= FLScut/2 & dat(:,4) <= FLScut/2*fitFLS(1)+fitFLS(2)));  %points with one signal above and one signal below merge point
indtemp = 1:2:length(t);  %odd indices
datmerged(t(indtemp),2) = dat(t(indtemp),4); %"dither" by reassigning every other ambiguous point

datmerged(:,3) = dat(:,5)*fitCHL(1)+fitCHL(2);
t = find(dat(:,5) < CHLcut/2 & dat(:,6) < CHLcut/2*fitCHL(1)+fitCHL(2));  %unambiguously offscale low on low gain...
datmerged(t,3) = dat(t,6);  %reassign stuff on baseline to be high gain
t = find((dat(:,5) <= CHLcut/2 & dat(:,6) >= CHLcut/2*fitCHL(1)+fitCHL(2)) | (dat(:,5) >= CHLcut/2 & dat(:,6) <= CHLcut/2*fitCHL(1)+fitCHL(2)));  %points with one signal above and one signal below merge point
indtemp = 1:2:length(t);  %odd indices
datmerged(t(indtemp),3) = dat(t(indtemp),6); %"dither" by reassigning every other ambiguous point

datmerged(:,4) = dat(:,7)*fitSSC(1)+fitSSC(2);
t = find(dat(:,7) < SSCcut/2 & dat(:,8) < SSCcut/2*fitSSC(1)+fitSSC(2));  %unambiguously offscale low on low gain...
datmerged(t,4) = dat(t,8);  %reassign stuff on baseline to be high gain
t = find((dat(:,7) <= SSCcut/2 & dat(:,8) >= SSCcut/2*fitSSC(1)+fitSSC(2)) | (dat(:,7) >= SSCcut/2 & dat(:,8) <= SSCcut/2*fitSSC(1)+fitSSC(2)));  %points with one signal above and one signal below merge point
indtemp = 1:2:length(t);  %odd indices
datmerged(t(indtemp),4) = dat(t(indtemp),8); %"dither" by reassigning every other ambiguous point

fit = [fitPE; fitFLS; fitCHL; fitSSC];

if 0,   %make plots or not
figure(1)
clf
subplot(221)
plot(dat(:,1), dat(:,2), '.')
hold on
eval(['fplot(''x * ' num2str(fitPE(1)) ' + ' num2str(fitPE(2)) ''', [0, 100], ''color'', ''r'')'])
axis([0 200 0 20000])
ylabel('PE high gain')
xlabel('PE low gain')
subplot(222)
plot(dat(:,3), dat(:,4), '.')
hold on
eval(['fplot(''x * ' num2str(fitFLS(1)) ' + ' num2str(fitFLS(2)) ''', [0, 100], ''color'', ''r'')'])
axis([0 200 0 20000])
ylabel('FLS high gain')
xlabel('FLS low gain')
subplot(223)
plot(dat(:,5), dat(:,6), '.')
hold on
eval(['fplot(''x * ' num2str(fitCHL(1)) ' + ' num2str(fitCHL(2)) ''', [0, 50], ''color'', ''r'')'])
axis([0 200 0 20000])
ylabel('CHL high gain')
xlabel('CHL low gain')
subplot(224)
plot(dat(:,7), dat(:,8), '.')
hold on
eval(['fplot(''x * ' num2str(fitSSC(1)) ' + ' num2str(fitSSC(2)) ''', [0, 50], ''color'', ''r'')'])
axis([0 200 0 20000])
ylabel('SSC high gain')
xlabel('SSC low gain')

figure(2)
clf
fittemp = fitFLS;
dattoplot = dat(:,7:8); %1:2 PE, 3:4 FLS, 5:6 CHL, 7:8 SSC
datmergedtoplot = datmerged(:,4);  %1 = PE, 2 = FLS, 3 = CHL, 4 = SSC
maxchannel = 16383*fittemp(1)+fittemp(2)+10;
bins = 10.^(0:log10(maxchannel+1)/256:log10(maxchannel));  %make 256 log spaced bins
lowhist = histc(dattoplot(:,1)*fittemp(1)+fittemp(2),bins);
highhist = histc(dattoplot(:,2),bins);
mergedhist = histc(datmergedtoplot, bins);
%subplot(211)
bar(bins, [highhist,lowhist, mergedhist], 'histc')
set(gca, 'xscale', 'log')
axis([100 1e5 0 5000])
legend('bins', 'low', 'hi', 'merg')
title('SSC')
figure(1)
pause

end; %if 1 for plots