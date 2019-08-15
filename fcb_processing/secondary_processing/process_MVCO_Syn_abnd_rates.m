%Summarize and process MVCO Synechococcus abundance, model division rate estimates, net growth and loss rate

% MAKE SURE CONNECTED TO SOSIKNAS FIRST!!!

addpath /Users/kristenhunter-cevera/mvco_tools
addpath ~/NES-LTER/fcb_processing/secondary_processing/
addpath ~/NES-LTER/fcb_processing/miscellaneous/
%/Volumes/Lab_data/MVCO/FCB/Syn_divrate_model/

%% Synechococcus cell abundance, fluorescence and size:

beadplotflag=1; %for QC of beads, SSC, PE and CHL fluorescence
abnd_plotflag=1; %for QC abundance

allmatdate=[];
allsynconc=[];
allbeads=[];
allsynSSC=[];
allsynSSCmode=[];
allsynPE=[];
allsynPEmode=[];
allsynCHL=[];
allsynCHLmode=[];
allsynvol=[];
allsynvolmode=[];
FCBnum=[];
%%
for year2do=2012
    
    switch year2do
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    
    %use the most up-to-date processing:
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/data/processed/grouped/groupsum.mat'])
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/data/processed/beads/beadresults.mat'])
    [ss is]=sort(beadresults(:,1)); %sometimes data points are out of order...
    beadresults=beadresults(is,:);

    to_use=exclude_data(cellresultsall,year2do); %ALSO REMOVES NANS!
    
    %known bead outliers:
    switch year2do
        case 2003
            ind=find(beadresults(:,13) < 0.9e4); %SSC outlier
            beadresults(ind,:)=NaN;
        case 2009
            ind=find(beadresults(:,13) > 4.05e4); %SSC outlier
            beadresults(ind,:)=NaN;
        case 2010
            ind=find(beadresults(:,10) > 2e4); %PE outlier
            beadresults(ind,:)=NaN;
        case 2011
            ind=find(beadresults(:,13) > 7e4); %SSC outlier
            beadresults(ind,:)=NaN;
        case 2013
            ind=find(beadresults(:,13) > 10e4); %SSC outlier
            beadresults(ind,:)=NaN;
            ind=find(beadresults(:,10) > 1.8e4); %PE outlier
            beadresults(ind,:)=NaN;
        case 2016
            ind=find(beadresults(:,13) > 8e5); %SSC outlier
            beadresults(ind,:)=NaN;
            ind=find(beadresults(:,10) > 4e4); %PE outlier
            beadresults(ind,:)=NaN;
        case 2017
            ind=find(beadresults(:,10) > 1.5e4); %PE/CHL outlier
            beadresults(ind,:)=NaN;
    end
    
    % smooth bead mean ?
    sm_bead_avgSSC=mvco_running_average(beadresults(:,1),beadresults(:,13),3,1); %running average smoothing function that takes into account breaks in FCB deployments
    sm_bead_avgPE=mvco_running_average(beadresults(:,1),beadresults(:,10),3,1); %running average smoothing function that takes into account breaks in FCB deployments
    sm_bead_avgCHL=mvco_running_average(beadresults(:,1),beadresults(:,12),3,1); %running average smoothing function that takes into account breaks in FCB deployments
    
    %make a 'beadmatch' equivalent: we will use the bead SSC mean...and
    %make sure that for a given day, the same value is used...
    
    % INTERPOLATED BEAD VALUES:
    beadmatch=mvco_interpolation(beadresults(:,1),sm_bead_avgSSC,cellresultsall(to_use,1),1); %one day as a gap
    beadPEmatch=mvco_interpolation(beadresults(:,1),sm_bead_avgPE,cellresultsall(to_use,1),1);
    beadCHLmatch=mvco_interpolation(beadresults(:,1),sm_bead_avgCHL,cellresultsall(to_use,1),1);
    
    %and plot to check!
    if beadplotflag==1
        figure(13), clf, set(gcf,'position',[164    55   965   930])
        subplot(3,1,1,'replace'), hold on
        h1=plot(beadresults(:,1),beadresults(:,13),'.--');
        h2=plot(beadresults(:,1),sm_bead_avgSSC,'o');
        h3=plot(cellresultsall(to_use),beadmatch,'.');
        %check the syn!
        sc=max(beadresults(:,13))./max(cellSSCall(to_use,1)./beadmatch);
        h5=plot(cellresultsall(setxor(to_use,1:length(cellresultsall)),1),quantile(beadresults(:,13),0.9),'.','color',[0.5 0.5 0.5]); %shows where data is missing
        h4=plot(cellresultsall(to_use,1),sc*cellSSCall(to_use,1)./beadmatch,'.-'); %scaled SSC
        plot(cellresultsall(to_use,1),sc*cellSSCmodeall(to_use,1)./beadmatch,'.','color',[0 0 0.5]) %scaled SSC
        datetick('x','mm/dd')
        ylabel('SSC')
        legend([h1(1);h2(1);h3(1);h4(1);h5(1)],'bead data','smoothed bead data','matched data','scaled Syn data','unused data','location','Eastoutside')
        title(num2str(year2do))
    
        subplot(3,1,2,'replace'), hold on
        plot(beadresults(:,1),beadresults(:,10),'.--')
        plot(beadresults(:,1),sm_bead_avgPE,'o')
        plot(cellresultsall(to_use),beadPEmatch,'.')
        plot(cellresultsall(setxor(to_use,1:length(cellresultsall)),1),quantile(beadresults(:,10),0.9),'.','color',[0.5 0.5 0.5]) %shows where data is missing
        sc=max(beadresults(:,10))./max(cellPEall(to_use,1)./beadPEmatch);
        plot(cellresultsall(to_use,1),sc*cellPEall(to_use,1)./beadPEmatch,'.-') %scaled PE
        plot(cellresultsall(to_use,1),sc*cellPEmodeall(to_use,1)./beadPEmatch,'.','color',[0 0 0.5])
        datetick('x','mm/dd')
        ylabel('PE')
    
        subplot(3,1,3,'replace'), hold on
        plot(beadresults(:,1),beadresults(:,12),'.--')
        plot(beadresults(:,1),sm_bead_avgCHL,'o')
        plot(cellresultsall(to_use),beadCHLmatch,'.')
        plot(cellresultsall(setxor(to_use,1:length(cellresultsall)),1),quantile(beadresults(:,12),0.9),'.','color',[0.5 0.5 0.5]) %shows where data is missing
        sc=max(beadresults(:,12))./max(cellCHLall(to_use,1)./beadCHLmatch);
        plot(cellresultsall(to_use,1),sc*cellCHLall(to_use,1)./beadCHLmatch,'.-') %scaled CHL
        plot(cellresultsall(to_use,1),sc*cellCHLmodeall(to_use,1)./beadCHLmatch,'.','color',[0 0 0.5])
        datetick('x','mm/dd')
        ylabel('CHL')
%     
        figure(14) %raw values...
        subplot(2,2,1,'replace')
        plot(cellPEmodeall(:,1),cellPEall(:,1),'.','color',[1 0.5 0])
        title('PE'), xlabel('mode'), ylabel('mean')
        subplot(2,2,2,'replace')
        plot(cellCHLmodeall(:,1),cellCHLall(:,1),'.','color',[1 0 0])
        title('CHL'), xlabel('mode'), ylabel('mean')
        subplot(2,2,3,'replace')
        plot(cellSSCmodeall(:,1),cellSSCall(:,1),'.','color',[0 0.8 0])
        title('SSC'), xlabel('mode'), ylabel('mean')
        subplot(2,2,4,'replace')
        plot(cellPEall(:,1),cellCHLall(:,1),'.')
        xlabel('PE')
        ylabel('CHL')
        
        keyboard
    end
    
    if abnd_plotflag==1
        
        figure(6)
        plot(cellresultsall(:,1),cellNUMall(:,1)./cellresultsall(:,3),'.-')
        datetick('x','mm/dd')
        
        keyboard
    end
    %normalize if needed and record:
    allmatdate=[allmatdate; cellresultsall(to_use,1)];
    allsynconc=[allsynconc; cellNUMall(to_use,1)./cellresultsall(to_use,3)]; %syn cell counts
    
    allsynSSC=[allsynSSC; cellSSCall(to_use,1)./beadmatch]; %SSC
    allsynSSCmode=[allsynSSCmode; cellSSCmodeall(to_use,1)./beadmatch];
    
    allsynvol = [allsynvol; cytosub_SSC2vol(cellSSCall(to_use,1)./beadmatch)]; %volume calc
    allsynvolmode = [allsynvolmode; cytosub_SSC2vol(cellSSCmodeall(to_use,1)./beadmatch)];
    
    allsynPE=[allsynPE; cellPEall(to_use,1)./beadPEmatch]; %PE
    allsynPEmode=[allsynPEmode; cellPEmodeall(to_use,1)./beadPEmatch];
    
    allsynCHL=[allsynCHL; cellCHLall(to_use,1)./beadCHLmatch]; %CHL
    allsynCHLmode=[allsynCHLmode; cellCHLmodeall(to_use,1)./beadCHLmatch];
    
    FCBnum=[FCBnum; FCBnumberall(to_use)];
end


%DO NOT SHIFT! BETTER TO WORK WITH DATA AS UTC!!!

%smooth the abundance data over 48 hours, with data separated by 1 day treated as separate chunks:
synrunavg=mvco_running_average(allmatdate, allsynconc,48,1);

clearvars -except synrunavg allmatdate allsynconc allsynSSC* allsynPE* allsynvol* allsynCHL* FCBnum

%% bin syn cell abundance:

%from average values:
[time_syn_ns_dy, daily_syn_ns] = timeseries2ydmat(allmatdate, allsynconc); %raw, unsmoothed abundance
[time_syn, daily_syn, synyears, ydmu] = timeseries2ydmat(allmatdate, synrunavg); %smoothed abundance

%and weekly:
[weekly_syn, time_syn_wk, yd_wk] = ydmat2weeklymat(daily_syn, synyears);
[syn_avg_wk, syn_std_wk] = dy2wkmn_climatology(daily_syn, synyears);

%daily avg:
syn_avg=nanmean(daily_syn,2);
syn_avg(end)=NaN; %this is because only one year with a leap year could be included, so the average isn't really good....
syn_std=nanstd(daily_syn,0,2);

%median:
syn_med=nanmedian(daily_syn,2);
[syn_med_wk] = dy2wkmn_medclimatology(daily_syn, synyears);

%% Division rates from the model-----------------------------------------------------------------------------------------------------------------------------------------------

allgrowthrates=[];
modelallmatdate=[];
allMR=[];

for year2do=2003:2018
    
    switch year2do
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    
    rootpath='/Volumes/Lab_data/MVCO/FCB/';   
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/model/output_June2019/mvco_14par_dmn_' num2str(year2do) '.mat'])       
 
    [days2redo, days2exclude]=exclude_modeldata(year2do);
    
    if ~isempty(days2exclude)
        days2exclude=str2num(cell2mat(days2exclude(:,1)));
        to_use=find(ismember(modelresults(:,1),days2exclude)==0);
    else
        to_use=find(~isnan(modelresults(:,1)) & modelresults(:,1)~=0); %just in case ;)
    end
    
    %     clf,
    %     plot(modelresults(:,1),modelresults(:,17),'r.','markersize',20)
    %     hold on
    %     plot(modelresults(to_use,1),modelresults(to_use,17),'.','color',[0.2081    0.1663    0.5292],'markersize',20)
    %     title(num2str(year2do))
    %     pause
    
    allgrowthrates=[allgrowthrates; modelresults(to_use,17)];
    modelallmatdate=[modelallmatdate; modelresults(to_use,1)];
    allMR=[allMR; modelresults];
    
end

%
% remove days that don't have a rate (model run didn't work...)
% ind=find(modelallmatdate==0);
% ii=setdiff(1:length(modelallmatdate),ind);
% modelallmatdate=modelallmatdate(ii);
% allgrowthrates=allgrowthrates(ii);
% allMR=allMR(ii,:);

%% and those above 2 per day?
ind=find(allgrowthrates > 2);
allgrowthrates(ind)=NaN;

% % exclude spurious mu's based on temperature data:
% [mdate_mu, daily_mu, yearlist, ydmu ] = timeseries2ydmat(modelallmatdate, allgrowthrates);
% daily_mu=[nan(366,1) daily_mu]; %just if don't have 2003 data yet...
% jj=find((daily_mu > 0.5 & Tday < 8) | (daily_mu > 0.4 & Tday < 4)) % | (daily_mu > 0.3 & Tday < 6));
% daily_mu(jj)=nan;

clearvars 'modelresults*' 'allmodelruns*' ii ind MR filelist filename %remove individual years from workspace
%clearvars -except synrunavg allmatdate allsynconc allgrowthrates allMR modelallmatdate

% bin division rates
[time_mu, daily_mu, muyears] = timeseries2ydmat(modelallmatdate, allgrowthrates);

% %exclude suspicious mu's baesd on temperature?
[smu, new_mu_est, no_new_est]=replace_suspicious_mus(time_mu,daily_mu,0); %1 or 0 for plotflag!

%%


if exist('/Volumes/Lab_data/MVCO/','dir')
    rootpath='/Volumes/Lab_data/MVCO/FCB/';
    load(fullfile(rootpath,'/Syn_and_MVCO_packaged_data/mvco_envdata_current.mat'),'Tbeam_corr')
else
    rootpath='\\sosiknas\Lab_data\MVCO\FCB\';
   load(fullfile(rootpath,'\Syn_and_MVCO_packaged_data/mvco_envdata_current.mat'),'Tbeam_corr');
end
%load('~/MVCO_light_at_depth/syn_data_analysis/mvco_envdata_29Jan2018.mat', 'Tbeam_corr') %mvco_envdata_16Aug2016.mat

figure, plot(Tbeam_corr,daily_mu,'.','color',[0.4 0.4 0.4])
hold on, plot(Tbeam_corr(new_mu_est(:,1)),daily_mu(new_mu_est(:,1)),'rp')
plot(Tbeam_corr(no_new_est(:,1)),daily_mu(no_new_est(:,1)),'bp')
% plot(Tbeam_corr(new_mu_est(:,1)),new_mu_est(:,2),'kp')
plot(Tbeam_corr(no_new_est(:,1)),no_new_est(:,3),'cp')

%% Here, have the choice to replace suspicious mus with close estimates, but for now, let's just set those to nan:
daily_mu(new_mu_est(:,1))=nan; %new_mu_est(:,2);
daily_mu(no_new_est(:,1))=nan;
%%
% adjust size for mu matrix that may not have all the years yet:
[yearsmissing, iy] = setdiff(synyears,muyears);
if ~isempty(yearsmissing)
    daily_mu(:,iy)=nan(size(daily_mu,1),length(iy));
end

mu_avg=nanmean(daily_mu,2);
mu_std=nanstd(daily_mu,0,2);

% The script saprse_weeklybin.m  bins the data according to calendar weeks. If however, a
%week only has 1 observation, then that obs is tacked onto the
%previous or following week, depending on which has closest data.


%[weekly_mu, time_mu_wk] = ydmat2weeklymat(daily_mu, synyears);
[time_mu_wk, weekly_mu]=sparse_weeklybin(time_mu,daily_mu,synyears);
[mu_avg_wk, mu_std_wk] = dy2wkmn_climatology(daily_mu, synyears);
[mu_med_wk] = dy2wkmn_medclimatology(daily_mu, synyears);
% mu_avg_wk=nanmean(weekly_mu,2);
% mu_std_wk=nanstd(weekly_mu,0,2);


%% Compute net growth rate and loss rates:
%-----------------------------------------------------------------------------------------------------------------------------------------------
%find the syn abundance for each day that matches around dawn and use this to calculate net growth rate

%Find dawn hour from solar.mat files for each year:

dawnhours=[];
for year2do=2003:2018;
    
    switch year2do
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    
    if ismember(year2do,[2005:2007 2010:2013])
        eval(['load ' rootpath 'MVCO_' filelabel num2str(year2do) '/model/solar' num2str(year2do) '_w_buoy.mat'])
    else
    eval(['load ' rootpath 'MVCO_' filelabel num2str(year2do) '/model/solar' num2str(year2do) '.mat'])
    end
    
    %eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' year2do num2str(year) '/model/solar' num2str(year) '.mat;']) %Change path here!
    dawnhours=[dawnhours; dawn];
    
end

[~, ia]=unique(dawnhours(:,1));
dawnhours=dawnhours(ia,:);
[~, daily_dawn, ~, ~ ] = timeseries2ydmat(dawnhours(:,1), dawnhours(:,2));
dawn_avg=nanmedian(daily_dawn,2);

%Can also find it from the dielstarthr for each day in the setupdays files... see original climatology_data.m around June 2016

%% ----------------------------------------------------------------------------------------------------------
%For each day in the syn abundance data, see if there is a corresponding model day,
%record index for matching day and find the time period that best matches the
%diel period to use to calculate net growth rate:

daylist=unique(floor(allmatdate));
dawnavg=zeros(length(daylist),5);
dawn_ind=cell(length(daylist),1);
netind=[];

for j=1:length(daylist)
    
    day=daylist(j);
    
    %average over just a few hours around dawn if possible:
    dw=find(dawnhours(:,1)==day);
    dielst=dawnhours(dw,2);
    
    if isempty(dw) | isnan(dielst) %have not found a dawn value or has a nan value
        %use closest day value:
        dii=find(dawnhours(:,1) >= day-2 & dawnhours(:,1) <= day+2);
        if ~isempty(dii)
            dielst=nanmean(dawnhours(dii,2));
        else %if still empty, use climatology average:
            yd=find_yearday(day);
            dielst=dawn_avg(yd);
        end
        
    end
    
    %Now that we have an estimate of dawn for each day - Are there measurements near dawn?
    ind=find(allmatdate(:,1) >= day+((dielst(1)-2)/24) & allmatdate(:,1) <= day+((dielst(1)+2)/24));
    
    if ~isempty(ind)
        dawnavg(j,:)=[day nanmean(synrunavg(ind)) length(ind) 24*(allmatdate(ind(end))-allmatdate(ind(1))) dielst(1)];
        dawn_ind(j)={ind};
    else %if not, then need to skip...
        dawnavg(j,:)=[day nan(1,3) dielst];
    end
    
    
end

%% ---------------------------------------------------------------------------------------------------------
% Calculate the net growth rate
%there are a few ways to do this, temporarily have settled on taking the
%average cell conc around dawn (+/- 3 hours) and using this as the syn abnd
%for each day to do the net growth rate calc:

netmu_avg=[dawnavg(1:end-1,1) log(dawnavg(2:end,2)./dawnavg(1:end-1,2)) diff(dawnavg(:,1)) diff(dawnavg(:,5))]; %[day    net mu    diff in time     diff in dawnhr]

%exclude any gap longer than a day:
jj=find(netmu_avg(:,3) > 1);
netmu_avg(jj,2)=nan;

%Find the days that having matching model runs:
[tt, id, im]=intersect(netmu_avg(:,1),allMR(:,1));

%THERE ARE SOME MODEL DAYS THAT DO NOT HAVE A NET GROWTH RATE -
%LIKEY, THESE ARE DAYS WHERE DAWN HOURS CAN BE MISSING!!!
%BUT NEED TO CHECK THIS!

%INDEED, all days come back as either missing 1-5 of first dawn hours, or
%abundance was screened out later due to other problems...

% Pause here for sanity figures!
% lightblue=[191 	239 	255]./[255 255 255];
%
% figure
% plot(allmatdate,synrunavg,'k.-')
% hold on
% for j=1:length(dawnavg)
%     plot(allmatdate(dawn_ind{j}),synrunavg(dawn_ind{j}),'r.--','markersize',10)
% end
%
% temp=netmu_avg(id,:);
% qq=find(isnan(temp(:,2)));
% missing_days=temp(qq,1);
%
% for q=1:length(missing_days);
%
%     day=missing_days(q);
%     w1=find(dawnavg(:,1)==day);
%     dielst=dawnavg(w1,5);
%
%     w2=find(netmu_avg(:,1)==day);
%
%     %How many hours of data are past dawn?
%     %  ind=find(allmatdate(:,1) > day+((dielst-2)/24) & allmatdate(:,1) < day+((dielst+3)/24)+1);
%     %  fprintf('Length of ind: %f\n',length(ind))
%
%     f1=fill([day; day+1; day+1; day]+dielst/24,[0 0 5e5 5e5],lightblue);hold on;
%     set(f1,'linestyle', 'none'), xlim([day-3 day+3]),  uistack(f1,'bottom')
%     ylim([0.5*dawnavg(w1+1,2) 1.2*dawnavg(w1+1,2)])
%
%     datetick('x','mm/dd/yy','keeplimits')
%     disp([num2str(q) ' out of ' num2str(length(missing_days))])
%     pause
%
% end

%% Bin net growth rate data:

[time_net, daily_net, netyears] = timeseries2ydmat(netmu_avg(:,1), netmu_avg(:,2));
net_avg=nanmean(daily_net,2);
net_std=nanstd(daily_net,0,2);

[weekly_net, time_net_wk] = ydmat2weeklymat(daily_net, netyears );
[net_avg_wk, net_std_wk] = dy2wkmn_climatology(daily_net, netyears);
[net_med_wk] = dy2wkmn_medclimatology(daily_net, synyears);
%% Back calculate loss rates:
%---------------------------------------------------------------------------------------------------------

daily_loss=daily_mu-daily_net;
loss_avg=nanmean(daily_loss,2);
loss_std=nanstd(daily_loss,0,2);

[time_loss_wk, weekly_loss]=sparse_weeklybin(time_net,daily_loss,synyears);
[loss_avg_wk, loss_std_wk] = dy2wkmn_climatology(daily_loss, synyears);
[loss_med_wk] = dy2wkmn_medclimatology(daily_loss, synyears);

%% we need volume and fluorescence to match up to light, which was calc on local time:

%from hourly mean values (already bead normalized):
% [time_PE, daily_PE] = timeseries2ydmat(allmatdate-4/24, allsynPE); %Syn PE fluorescence
% [time_SSC, daily_SSC] = timeseries2ydmat(allmatdate-4/24, allsynSSC); %SSC, bead normalized
% [time_vol, daily_vol] = timeseries2ydmat(allmatdate-4/24, allsynvol); %Cell volume from SSC-bead normalized
% [time_CHL, daily_CHL] = timeseries2ydmat(allmatdate-4/24, allsynCHL); %Syn CHL fluorescence
% 
% PE_avg=nanmean(daily_PE,2);
% PE_med=nanmedian(daily_PE,2);
% SSC_avg=nanmean(daily_SSC,2);
% SSC_med=nanmedian(daily_SSC,2);
% vol_avg=nanmean(daily_vol,2);
% vol_med=nanmedian(daily_vol,2);
% CHL_avg=nanmean(daily_CHL,2);
% CHL_med=nanmedian(daily_CHL,2);
% 
% [PE_avg_wk, PE_std_wk] = dy2wkmn_climatology(daily_PE, synyears);
% [CHL_avg_wk, CHL_std_wk] = dy2wkmn_climatology(daily_CHL, synyears);
% [SSC_avg_wk, SSC_std_wk] = dy2wkmn_climatology(daily_SSC, synyears);
% [vol_avg_wk, vol_std_wk] = dy2wkmn_climatology(daily_vol, synyears);
% 
% %from hourly mode values:
% [time_PEmode, daily_PEmode] = timeseries2ydmat(allmatdate-4/24, allsynPEmode); %Syn PE fluorescence
% [time_SSCmode, daily_SSCmode] = timeseries2ydmat(allmatdate-4/24, allsynSSCmode); %smoothed abundance
% [time_volmode, daily_volmode] = timeseries2ydmat(allmatdate-4/24, allsynvolmode); %smoothed abundance
% [time_CHLmode, daily_CHLmode] = timeseries2ydmat(allmatdate-4/24, allsynCHLmode); %Syn CHL fluorescence
% 
% PE_mode_avg=nanmean(daily_PEmode,2);
% SSC_mode_avg=nanmean(daily_SSCmode,2);
% vol_mode_avg=nanmean(daily_volmode,2);
% CHL_mode_avg=nanmean(daily_CHLmode,2);
% 
% [PE_mode_avg_wk, PE_mode_std_wk] = dy2wkmn_climatology(daily_PEmode, synyears);
% [SSC_mode_avg_wk, SSC_mode_std_wk] = dy2wkmn_climatology(daily_SSCmode, synyears);
% [vol_mode_avg_wk, vol_mode_std_wk] = dy2wkmn_climatology(daily_volmode, synyears);
% [CHL_mode_avg_wk, CHL_mode_std_wk] = dy2wkmn_climatology(daily_CHLmode, synyears);

%medians and max/min from the mode values:
% [time_vol_Q, daily_vol_Q10] = timeseries2ydmat_quantile(allmatdate, allsynvolmode, .10);
% [time_vol_Q, daily_vol_Q90] = timeseries2ydmat_quantile(allmatdate, allsynvolmode, .90);
% [time_vol_Q, daily_vol_Q50] = timeseries2ydmat_quantile(allmatdate, allsynvolmode, .50);
%[time_vol_Q, daily_vol_max] = timeseries2ydmat_quantile(allmatdate-4/24, allsynvolmode, 1);
%% So...we also want the minimum cell volume and then the corresponding PE fluorscence that goes with it as metrics:

unqdays=unique(floor(allmatdate));
minvol=nan(length(unqdays),3); maxvol=nan(length(unqdays),3);  minPE=nan(length(unqdays),3); maxPE=nan(length(unqdays),3);
for q=1:length(unqdays)
    day=unqdays(q);
    ww=find(dawnhours(:,1)==day);
    if isempty(ww) %use median value of dawn instead
        jj=find_yearday(day);
        dawn=dawn_avg(jj);
        disp('Using average value of dawn')
    else
        dawn=dawnhours(ww,2);
    end
    
    qq=find(allmatdate >= day+dawn/24 & allmatdate < day+1+dawn/24); %should be dawn to dawn window...
    if length(qq) > 22  %we also really only want mins from complete days...
        
        [mm, im]=min(allsynvolmode(qq));
        minvol(q,1)=mm; %min volu
        minvol(q,2)=qq(im);  %index back for date and other matrices
        minvol(q,3)=(allmatdate(qq(im))-day-dawn/24)*24; %hours after or before dawn
        
        [mm, im]=max(allsynvolmode(qq));
        maxvol(q,1)=mm; %min volu
        maxvol(q,2)=qq(im);  %index back for date and other matrices
        maxvol(q,3)=(allmatdate(qq(im))-day-dawn/24)*24; %hours after or before dawn
        
        [mm, im]=min(allsynPEmode(qq));
        minPE(q,1)=mm; %min volu
        minPE(q,2)=qq(im);  %index back for date and other matrices
        minPE(q,3)=(allmatdate(qq(im))-day-dawn/24)*24; %hours after or before dawn
        
        [mm, im]=max(allsynPEmode(qq));
        maxPE(q,1)=mm; %min volu
        maxPE(q,2)=qq(im);  %index back for date and other matrices
        maxPE(q,3)=(allmatdate(qq(im))-day-dawn/24)*24; %hours after or before dawn
        
    end
end


%% IF WANT MEASUREMENT AROUND DAWN FOR EACH DAY:
%One could also imagine just taking the metrics of the before dawn hour of each day...
%find predawn hour of each day:
unqdays=unique(floor(allmatdate));
dawnPE=nan(length(unqdays),4);
for q=1:length(unqdays)
    day=unqdays(q);
    ww=find(dawnhours(:,1)==day);
    if isempty(ww) %use median value of dawn instead
        jj=find_yearday(day);
        dawn=dawn_avg(jj);
        disp('Using average value of dawn')
    else
        dawn=dawnhours(ww,2);
    end
    
    qq=find(allmatdate >= day+dawn/24-1.5/24 & allmatdate <= day+dawn/24+1.5/24); %3 hours around dawn
    
    if ~isempty(qq)
        [mm, im]=min(allsynPE(qq)); %min(allsynvolmode(qq));
        dawnPE(q,1)=mm; %min mean PE value around dawn
        dawnPE(q,2)=qq(im);  %index back for date and other matrices
        dawnPE(q,3)=(allmatdate(qq(im))-day)*24-dawn; %hours after or before 'dawn'
        dawnPE(q,4)=mean(allsynPE(qq));
    end
end



%%
minvol=minvol(~isnan(minvol(:,1)),:); %remove nan's
maxvol=maxvol(~isnan(maxvol(:,1)),:);
minPE=minPE(~isnan(minPE(:,1)),:); %min PE over day
dawnPE=dawnPE(~isnan(dawnPE(:,1)),:); %min PE around dawn

%keyboard
%% into matrices:
[time_volmin,daily_vol_min]=timeseries2ydmat_quantile(allmatdate(minvol(:,2)),minvol(:,1),0); %at this point, just using the script to bin into matrix for timeseries
[time_volmax,daily_vol_max]=timeseries2ydmat_quantile(allmatdate(maxvol(:,2)),maxvol(:,1),0); %at this point, just using the script to bin into matrix for timeseries

%normalize PE by volume:
[time_PEmin,daily_PE_min]=timeseries2ydmat(allmatdate(minPE(:,2)),minPE(:,1)); %at this point, just using the script to bin into matrix for timeseries
[time_PEmin,daily_PEminvol_ratio]=timeseries2ydmat(allmatdate(minPE(:,2)),allsynPEmode(minPE(:,2))./allsynvolmode(minPE(:,2))); %at this point, just using the script to bin into matrix for timeseries

%PE around dawn:

%MEAN of MEAN PE's
%RATIO THOUGH IS MIN PE AROUND DAWN divided by that volue...
[time_PE_dawn,daily_PE_dawn]=timeseries2ydmat(allmatdate(dawnPE(:,2)),dawnPE(:,4)); 
[time_PEmin_dawn,daily_PEvol_ratio_dawn]=timeseries2ydmat(allmatdate(dawnPE(:,2)),dawnPE(:,1)./allsynvolmode(dawnPE(:,2))); 

%%
%[time_vol_Q, daily_vol_min] = timeseries2ydmat_quantile(allmatdate-4/24, allsynvolmode, 0);
%[time_vol_Q, daily_vol_max] = timeseries2ydmat_quantile(allmatdate-4/24, allsynvolmode, 1);

% [time_PE, daily_PE_min] = timeseries2ydmat_quantile(allmatdate-4/24, allsynPEmode,0); %Syn PE fluorescence
% [time_PE, daily_PE_max] = timeseries2ydmat_quantile(allmatdate-4/24, allsynPEmode,1); %Syn PE fluorescence

vol_avgmax=nanmean(daily_vol_max,2);
vol_avgmin=nanmean(daily_vol_min,2);

%PE_avgmax=nanmean(daily_PE_max,2);
PE_avgmin=nanmean(daily_PE_min,2);
PE_avgdawn=nanmean(daily_PE_dawn,2);
PEminvol_ratio_avg=nanmean(daily_PEminvol_ratio,2);
PEdawnvol_ratio_avg=nanmean(daily_PEvol_ratio_dawn,2);

[PEmin_wkmed] = dy2wkmn_medclimatology(daily_PE_min, synyears);
[PEdawn_wkmed] = dy2wkmn_medclimatology(daily_PE_dawn, synyears);
%[PEmax_wkmed] = dy2wkmn_medclimatology(daily_PE_max, synyears);
[volmin_wkmed] = dy2wkmn_medclimatology(daily_vol_min, synyears);

%%

% eval(['save /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/syndata_TEMP_' dd(1:2) dd(4:6) dd(8:end) '.mat *syn* *mu* *net* *loss* *PE* *SSC* *vol* *CHL* allgrowthrates allMR modelallmatdate allmatdate'])

if exist('/Volumes/Lab_data/MVCO/','dir')
    savepath=fullfile('/Volumes/Lab_data/MVCO/FCB/Syn_and_MVCO_packaged_data/');
else
    savepath=fullfile('\\sosiknas\Lab_data\MVCO\FCB\Syn_and_MVCO_packaged_data\');
end

dd=date;
eval(['save ' savepath 'syndata_current.mat *syn* *mu* *net* *loss* *PE* *SSC* *vol* *CHL* allgrowthrates allMR modelallmatdate allmatdate FCBnum'])
%save a backup copy, if you'd like to return to this later:
eval(['save ' fullfile(savepath,'older_products/') 'syndata_' dd(1:2) dd(4:6) dd(8:end) '.mat *syn* *mu* *net* *loss* *PE* *SSC* *vol* *CHL* allgrowthrates allMR modelallmatdate allmatdate FCBnum'])

