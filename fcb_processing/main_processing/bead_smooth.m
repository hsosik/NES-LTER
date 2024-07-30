function bead_smooth(year2do)
%year2do = 2021;
monstr = 'Jan';
switch year2do
    case {2003, 2006}
        monstr = 'May';
    case {2004, 2005}
        monstr = 'Apr';
    case 2007
        monstr = 'Mar';
end
beadpath = ['\\sosiknas1\Lab_data\MVCO\FCB\MVCO_' monstr num2str(year2do) '\data\processed\beads\'];
% beadpath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Apr2004\data\processed\beads\';
% beadpath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Apr2005\data\processed\beads\';
% beadpath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_May2006\data\processed\beads\';
% beadpath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Mar2007\data\processed\beads\';
% beadpath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2010\data\processed\beads\';

gpath = regexprep(beadpath, 'beads', 'grouped');
g = load([gpath 'groupsum']);

plotflag = 1;
load([beadpath 'beadresults'])

[ss is]=sort(beadresults(:,1)); %happens for 2006 where out of time sync?
beadresults=beadresults(is,:);

%smooth mean of bead SSC with 3-point moving average:
% but first, check for outliers and other bead anomalies:
if plotflag ==1
    figure(6)
    subplot(2,1,1,'replace')
    plot(beadresults(:,1),beadresults(:,13),'.-')
%    keyboard
end

%known bead outliers:
switch year2do
    case 2003
       ind=find(beadresults(:,13) < 0.9e4); %outlier
       beadresults(ind,13)=NaN;
    case 2009
        ind=find(beadresults(:,13) > 4.05e4); %outlier
        beadresults(ind,13)=NaN;
    case 2010
        ind=find(beadresults(:,1) < datenum('1-aug-2010') & beadresults(:,13) < 2.5e4); %outlier
        beadresults(ind,13)=NaN;
    case 2011
        ind=find(beadresults(:,13) > 7e4); %outlier
        beadresults(ind,13)=NaN; 
    case 2013
        ind=find(beadresults(:,13) > 10e4); %outlier
        beadresults(ind,13)=NaN;
    case 2016
        ind=find(beadresults(:,13) > 20e4); %outlier
        beadresults(ind,13)=NaN;
end

sm_bead_avgSSC2=mvco_running_average(beadresults(:,1),beadresults(:,13:14),3,2); %running average smoothing function that takes into account breaks in FCB deployments
hold on
plot(beadresults(:,1), sm_bead_avgSSC2, '.-')

%insert constants across gaps
tt = find(diff(beadresults(:,1))>2);
for count = length(tt):-1:1
    t = tt(count);
    thalf = mean(beadresults(t:t+1,1));
    beadresults = [beadresults(1:t,:); [thalf-1e-6 beadresults(t,2:end)]; [thalf+1e-6 beadresults(t+1,2:end)]; beadresults(t+1:end,:)];
    sm_bead_avgSSC2 = [sm_bead_avgSSC2(1:t,:); sm_bead_avgSSC2(t); sm_bead_avgSSC2(t+1); sm_bead_avgSSC2(t+1:end)];
end


beadmatchSSCsmooth = interp1(beadresults(:,1), sm_bead_avgSSC2, g.cellresultsall(:,1));

save([gpath 'groupsum'], 'beadmatchSSCsmooth', '-append')

figure
plot(beadresults(:,1),beadresults(:,13),'.-k')
hold on
plot(beadresults(:,1),sm_bead_avgSSC2,'.-g')
plot(g.cellresultsall(:,1), beadmatchSSCsmooth, '.-b')

return
figure
plot(beadresults(:,1),beadresults(:,13)/15,'.-k')
hold on
plot(g.cellresultsall(:,1), g.cellSSCmodeall(:,1), 'c.-')
plot(g.cellresultsall(:,1), g.cellSSCmodeall(:,1)./g.beadmatchall(:,5).*nanmean(g.beadmatchall(:,5)), 'r.-')
plot(g.cellresultsall(:,1), g.beadmatchall(:,5)/15, 'b.-')

