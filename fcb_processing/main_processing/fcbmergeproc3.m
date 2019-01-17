function [datmerged, fit] = fcbmergeproc2(dat, fittouse)
%renamed to fcbmergeproc (from cytosubproc3), 3/03 Heidi
%modified to accept input fit for case of beads without good overlap
% April 12 2005 - fcbmergeproc2 modified from fcbmergeproc to deal with
% addition of SSC2 (new lower HV PMT, also with 2 gains)
% January 2015, add check to see if no good PE fit (in case use PE high gn only)
% December 2018, update from fcbmergproc2 adding >0 lower limits for high
%   gain channels before polyfit, necessary to deal with high noise levels
%   apparent in 2017 FCB1 results, Heidi

global special_flag

PEcut = 300; %200;
FLScut = 300; %200;
CHLcut = 200; %100;
SSC1cut = 200; %150;
SSC2cut = 200; %100;
SSC1_2acut = 4000; %3000;
SSC1_2bcut = 250; %150;

LowGain_min = 200;

fitPE = [NaN NaN NaN NaN]; fitFLS = fitPE; fitCHL = fitPE; fitSSC1 = fitPE; fitSSC1_2a = fitPE; fitSSC1_2b = fitPE;
datmerged = NaN.*ones(size(dat,1),6);
if isnan(fittouse)
    t = find(dat(:,1) > 0 & dat(:,1) < PEcut & dat(:,2) > PEcut);
    if length(t) > 3, temp = corrcoef(dat(t,1), dat(t,2)); fitPE  = [polyfit(dat(t,1), dat(t,2), 1) length(t) temp(2)]; end; %1 = slope, 2 = int, 3 = count, 4 = r
    t = find(dat(:,3) > 0 & dat(:,3) < FLScut & dat(:,4) > FLScut);
    if length(t) > 3, temp = corrcoef(dat(t,3), dat(t,4)); fitFLS = [polyfit(dat(t,3), dat(t,4), 1) length(t) temp(2)];  end;
    t = find(dat(:,5) > 0 & dat(:,5) < CHLcut & dat(:,6) > CHLcut*4);
    if length(t) > 3, temp = corrcoef(dat(t,5), dat(t,6)); fitCHL = [polyfit(dat(t,5), dat(t,6), 1) length(t) temp(2)];  end;
    t = find(dat(:,7) > 0 & dat(:,7) < SSC1cut & dat(:,8) > SSC1cut);
    if length(t) > 3, temp = corrcoef(dat(t,7), dat(t,8)); fitSSC1 = [polyfit(dat(t,7), dat(t,8), 1) length(t) temp(2)];   end;
else
    fitPE = [fittouse(1,:) NaN NaN];
    fitFLS = [fittouse(2,:) NaN NaN];
    fitCHL = [fittouse(3,:) NaN NaN];
    fitSSC1 = [fittouse(4,:) NaN NaN];
    fitSSC1_2a = [fittouse(5,:) NaN NaN];
    fitSSC1_2b = [fittouse(6,:) NaN NaN];
end;

if ~isnan(fitPE(1))
%use half of cut value as merge point
datmerged(:,1) = dat(:,1)*fitPE(1)+fitPE(2);
t = find(dat(:,1) < PEcut/2 & dat(:,2) < PEcut/2*fitPE(1)+fitPE(2));  %unambiguously offscale low on low gain...
datmerged(t,1) = dat(t,2);  %reassign stuff on baseline to be high gain
t = find((dat(:,1) <= PEcut/2 & dat(:,2) >= PEcut/2*fitPE(1)+fitPE(2)) | (dat(:,1) >= PEcut/2 & dat(:,2) <= PEcut/2*fitPE(1)+fitPE(2)));  %points with one signal above and one signal below merge point
indtemp = 1:2:length(t);  %odd indices
datmerged(t(indtemp),1) = dat(t(indtemp),2); %"dither" by reassigning every other ambiguous point
t = find(dat(:,1)< LowGain_min);
datmerged(t,1) = dat(t,2);  %Make sure low low gain signals stay as high gain (not dithered) Dec 2018
else
    disp('using PE high gain only')
    datmerged(:,1) = dat(:,2); 
end

datmerged(:,2) = dat(:,3)*fitFLS(1)+fitFLS(2);
t = find(dat(:,3) < FLScut/2 & dat(:,4) < FLScut/2*fitFLS(1)+fitFLS(2));  %unambiguously offscale low on low gain...
datmerged(t,2) = dat(t,4);  %reassign stuff on baseline to be high gain
t = find((dat(:,3) <= FLScut/2 & dat(:,4) >= FLScut/2*fitFLS(1)+fitFLS(2)) | (dat(:,3) >= FLScut/2 & dat(:,4) <= FLScut/2*fitFLS(1)+fitFLS(2)));  %points with one signal above and one signal below merge point
indtemp = 1:2:length(t);  %odd indices
datmerged(t(indtemp),2) = dat(t(indtemp),4); %"dither" by reassigning every other ambiguous point
t = find(dat(:,3)< LowGain_min);
datmerged(t,2) = dat(t,4);  %Make sure low low gain signals stay as high gain (not dithered) Dec 2018

datmerged(:,3) = dat(:,5)*fitCHL(1)+fitCHL(2);
t = find(dat(:,5) < CHLcut/2 & dat(:,6) < CHLcut/2*fitCHL(1)+fitCHL(2));  %unambiguously offscale low on low gain...
datmerged(t,3) = dat(t,6);  %reassign stuff on baseline to be high gain
t = find((dat(:,5) <= CHLcut/2 & dat(:,6) >= CHLcut/2*fitCHL(1)+fitCHL(2)) | (dat(:,5) >= CHLcut/2 & dat(:,6) <= CHLcut/2*fitCHL(1)+fitCHL(2)));  %points with one signal above and one signal below merge point
indtemp = 1:2:length(t);  %odd indices
datmerged(t(indtemp),3) = dat(t(indtemp),6); %"dither" by reassigning every other ambiguous point
t = find(dat(:,5)< LowGain_min);
datmerged(t,3) = dat(t,6);  %Make sure low low gain signals stay as high gain (not dithered) Dec 2018

datmerged(:,5) = dat(:,7)*fitSSC1(1)+fitSSC1(2);
t = find(dat(:,7) < SSC1cut/2 & dat(:,8) < SSC1cut/2*fitSSC1(1)+fitSSC1(2));  %unambiguously offscale low on low gain...
datmerged(t,5) = dat(t,8);  %reassign stuff on baseline to be high gain
t = find((dat(:,7) <= SSC1cut/2 & dat(:,8) >= SSC1cut/2*fitSSC1(1)+fitSSC1(2)) | (dat(:,7) >= SSC1cut/2 & dat(:,8) <= SSC1cut/2*fitSSC1(1)+fitSSC1(2)));  %points with one signal above and one signal below merge point
indtemp = 1:2:length(t);  %odd indices
datmerged(t(indtemp),5) = dat(t(indtemp),8); %"dither" by reassigning every other ambiguous point
t = find(dat(:,7)< LowGain_min);
datmerged(t,5) = dat(t,8);  %Make sure low low gain signals stay as high gain (not dithered) Dec 2018

%NOW merge over merged SSC1 with high gain SSC2
fitSSC1_2a = [NaN NaN NaN NaN];
if isnan(fittouse)
    t = find(dat(:,10) > 1000 & dat(:,10) < SSC1_2acut & datmerged(:,5) > 1e5 & datmerged(:,5) < 8e5); %5000*10);
    if length(t) > 100, 
        temp = corrcoef(dat(t,10), datmerged(t,5));
        fitSSC1_2a = [polyfit(dat(t,10), datmerged(t,5), 1) length(t) temp(2)]; 
    end;
else
    fitSSC1_2a = [fittouse(5,:) NaN NaN];
end;

if special_flag == 1, %case for 2009 FCB2 bad SSC2 data
    if fitSSC1_2a(1) > 220 | fitSSC1_2a(4) < 0.95,
        fitSSC1_2a = [NaN NaN NaN NaN];
    end;
end;

if ~isnan(fitSSC1_2a(1)),
    datmerged(:,6) = dat(:,10)*fitSSC1_2a(1)+fitSSC1_2a(2);
    t = find(dat(:,10) < SSC1_2acut/2 & datmerged(:,5) < SSC1_2acut/2*fitSSC1_2a(1)+fitSSC1_2a(2));  %unambiguously offscale low on low gain...
    datmerged(t,6) = datmerged(t,5);  %reassign stuff on baseline to be high gain
    t = find((dat(:,10) <= SSC1_2acut/2 & datmerged(:,5) >= SSC1_2acut/2*fitSSC1_2a(1)+fitSSC1_2a(2)) | (dat(:,10) >= SSC1_2acut/2 & datmerged(:,5) <= SSC1_2acut/2*fitSSC1_2a(1)+fitSSC1_2a(2)));  %points with one signal above and one signal below merge point
    indtemp = 1:2:length(t);  %odd indices
    datmerged(t(indtemp),6) = datmerged(t(indtemp),5); %"dither" by reassigning every other ambiguous point
else  %case where no good SSC2 results, use merged SSC1 for final
    datmerged(:,4) = datmerged(:,5); 
    disp('using only SSC1')
end;

%NOW merge over merged (SSC1 & high gain SSC2) and low gain SSC
fitSSC1_2b = [NaN NaN NaN NaN];
if isnan(fittouse)
    t = find(dat(:,9) > 0 & dat(:,9) < SSC1_2bcut & datmerged(:,6) > 1e4*10);
    if length(t) > 100 & ~isnan(fitSSC1_2a(1)), %only do this merge if enough points and SSC1_2a was good
        temp = corrcoef(dat(t,9), datmerged(t,6));
        fitSSC1_2b = [polyfit(dat(t,9), datmerged(t,6), 1) length(t) temp(2)];  
    end;
else
    fitSSC1_2b = [fittouse(6,:) NaN NaN];
end;
if ~isnan(fitSSC1_2b(1)), 
    datmerged(:,4) = dat(:,9)*fitSSC1_2b(1)+fitSSC1_2b(2);
    t = find(dat(:,9) < SSC1_2bcut/2 & datmerged(:,6) < SSC1_2bcut/2*fitSSC1_2b(1)+fitSSC1_2b(2));  %unambiguously offscale low on low gain...
    datmerged(t,4) = datmerged(t,6);  %reassign stuff on baseline to be high gain
    t = find((dat(:,9) <= SSC1_2bcut/2 & datmerged(:,6) >= SSC1_2bcut/2*fitSSC1_2b(1)+fitSSC1_2b(2)) | (dat(:,9) >= SSC1_2bcut/2 & datmerged(:,6) <= SSC1_2bcut/2*fitSSC1_2b(1)+fitSSC1_2b(2)));  %points with one signal above and one signal below merge point
    indtemp = 1:2:length(t);  %odd indices
    datmerged(t(indtemp),4) = datmerged(t(indtemp),6); %"dither" by reassigning every other ambiguous point
elseif ~isnan(fitSSC1_2a(1))  %case where no SSC2 low results
    datmerged(:,4) = datmerged(:,6); %use merged (SSC1 & SSC2 hi) as final\
    disp('using only SSC1 and SSC2 hi')
end;

fit = [fitPE; fitFLS; fitCHL; fitSSC1; fitSSC1_2a; fitSSC1_2b];

if 0  %make plots or not
figure(1)
clf
subplot(231)
plot(dat(:,1), dat(:,2), '.')
hold on
eval(['fplot(''x * ' num2str(fitPE(1)) ' + ' num2str(fitPE(2)) ''', [0, ' num2str(PEcut) '], ''color'', ''r'')'])
axis([0 600 0 20000])
ylabel('PE high gain')
xlabel('PE low gain')
subplot(232)
plot(dat(:,3), dat(:,4), '.')
hold on
eval(['fplot(''x * ' num2str(fitFLS(1)) ' + ' num2str(fitFLS(2)) ''', [0 ' num2str(FLScut) '], ''color'', ''r'')'])
axis([0 600 0 20000])
ylabel('FLS high gain')
xlabel('FLS low gain')
subplot(233)
plot(dat(:,5), dat(:,6), '.')
hold on
eval(['fplot(''x * ' num2str(fitCHL(1)) ' + ' num2str(fitCHL(2)) ''', [0 ' num2str(CHLcut) '], ''color'', ''r'')'])
axis([0 600 0 20000])
ylabel('CHL high gain')
xlabel('CHL low gain')
subplot(234)
plot(dat(:,7), dat(:,8), '.')
hold on
eval(['fplot(''x * ' num2str(fitSSC1(1)) ' + ' num2str(fitSSC1(2)) ''', [0 ' num2str(SSC1cut) '], ''color'', ''r'')'])
axis([0 600 0 20000])
ylabel('SSC1 high gain')
xlabel('SSC1 low gain')

subplot(235)
plot(dat(:,10), datmerged(:,5), '.')
hold on
eval(['fplot(''x * ' num2str(fitSSC1_2a(1)) ' + ' num2str(fitSSC1_2a(2)) ''', [0 ' num2str(SSC1_2acut) '], ''color'', ''r'')'])
axis([0 10000 0 1e6])
xlabel('SSC2 high gain')
ylabel('SSC1, merged')

subplot(236)
plot(dat(:,9), datmerged(:,6), '.')
hold on
eval(['fplot(''x * ' num2str(fitSSC1_2b(1)) ' + ' num2str(fitSSC1_2b(2)) ''', [0 ' num2str(SSC1_2bcut) '], ''color'', ''r'')'])
axis([0 500 0 4e6])
ylabel('SSC1&SSC2hi merged')
xlabel('SSC2 low')

figure(2)
clf
fittemp = fitSSC1;
if length(fittouse) > 1,
    fittemp = fittouse(5,:);
end;
dattoplot = dat(:,7:8); %1:2 PE, 3:4 FLS, 5:6 CHL, 7:8 SSC1, 9:10 SSC2
datmergedtoplot = datmerged(:,5);  %1 = PE, 2 = FLS, 3 = CHL, 4 = SSC1_2, 5 = SSC1, 6 = SSC2
maxchannel = 16383*fittemp(1)+fittemp(2)+10;
bins = 10.^(0:log10(maxchannel+1)/256:log10(maxchannel));  %make 256 log spaced bins
lowhist = histc(dattoplot(:,1)*fittemp(1)+fittemp(2),bins);
highhist = histc(dattoplot(:,2),bins);
mergedhist = histc(datmergedtoplot, bins);
%subplot(211)
if ~isnan(bins)
    bar(bins, [highhist,lowhist, mergedhist], 'histc')
end
set(gca, 'xscale', 'log')
axis([100 1e5 0 5000])
legend('low', 'hi', 'merg', 'bins')
title('SSC1')
figure(1)
figure(3)
loglog(datmerged(:,4), datmerged(:,3), '.'), axis([10 1e7 1e2 1e6])
pause %(.1)
%keyboard
end; %if 1 for plots