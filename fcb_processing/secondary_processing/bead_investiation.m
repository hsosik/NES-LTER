%bead investigation:

%looks at whether or not syn volumes drop if bead values also drop:
clear all
close all

beadSSC=[];
beadres=[];
allmatdate=[];
allmatdateSSC=[];
allmatdatePE=[];


allsynSSC=[];
allsynSSCmode=[];
allsynPE=[];
allsynPEmode=[];
allsynCHL=[];
allsynCHLmode=[];
allsynvol=[];
allsynvolmode=[];
FCBnum=[];

%% load in bead and other data:

addpath ~/NES-LTER/fcb_processing/secondary_processing/
addpath ~/NES-LTER/fcb_processing/miscellaneous/
addpath /Users/kristenhunter-cevera/mvco_tools

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
    
    %use the most up-to-date processing:
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/data/processed/grouped/groupsum.mat'])
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/data/processed/beads/beadresults.mat'])
    [ss is]=sort(beadresults(:,1)); %sometimes data points are out of order...
    beadresults=beadresults(is,:);
    
    [to_use,to_use_SSC,to_use_PE]=exclude_data(cellresultsall,year2do); %ALSO REMOVES NANS!
    
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
            ind=find(floor(beadresults(:,1))==datenum('5-9-13') | floor(beadresults(:,1))==datenum('5-8-13')); %bad bead values for this day
            beadresults(ind,:)=NaN;
            ind=find(beadresults(:,13) > 10e4); %SSC outlier
            beadresults(ind,:)=NaN;
            ind=find(beadresults(:,10) > 1.8e4); %PE outlier
            beadresults(ind,:)=NaN;
        case 2015
            ind1=find(floor(beadresults(:,1))>=datenum('5-19-15') &floor(beadresults(:,1))<=datenum('5-21-15')); %bad bead values for these days
            ind2=find(beadresults(:,1)>= (datenum('5-24-15')+23/24) & beadresults(:,1)<=datenum('5-30-15')); %bad bead values for these days
            ind3=find(beadresults(:,1)>= (datenum('June-6-15')+12/24) & beadresults(:,1)<=datenum('June-8-15')); %just weird bead values for these days
            beadresults([ind1;ind2;ind3],:)=NaN;
        case 2016
            ind=find(floor(beadresults(:,1))==datenum('9-23-16') | floor(beadresults(:,1))==datenum('9-24-16')); %bad bead values for these days
            beadresults(ind,:)=NaN;
            ind=find(floor(beadresults(:,1))==datenum('10-11-16')); %bad bead values for these days
            beadresults(ind,:)=NaN;
            ind=find(beadresults(:,13) > 8e5); %SSC outlier
            beadresults(ind,:)=NaN;
            ind=find(beadresults(:,10) > 4e4); %PE outlier
            beadresults(ind,:)=NaN;
        case 2017
            ind1=find(beadresults(:,1)>=datenum('8-28-17') & beadresults(:,1)<=datenum('8-29-17')+12/24); %did not find beads!
            ind2=find(floor(beadresults(:,1))==datenum('8-30-17')); %did not find beads!
            ind3=find(floor(beadresults(:,1))==datenum('9-1-17')); %did not find beads!
            ind4=find(floor(beadresults(:,1))==datenum('9-5-17')); %did not find beads!
            beadresults([ind1;ind2;ind3;ind4],:)=NaN;
            ind=find(beadresults(:,10) > 1.5e4); %PE/CHL outlier
            beadresults(ind,:)=NaN;
        case 2018
            ind=find(floor(beadresults(:,1))==datenum('7-30-18')); %did not find beads!
            beadresults(ind,:)=NaN;
    end
    
    % smooth bead mean ?
    sm_bead_avgSSC=mvco_running_average(beadresults(:,1),beadresults(:,13),3,1); %running average smoothing function that takes into account breaks in FCB deployments
    sm_bead_avgPE=mvco_running_average(beadresults(:,1),beadresults(:,10),3,1); %running average smoothing function that takes into account breaks in FCB deployments
    sm_bead_avgCHL=mvco_running_average(beadresults(:,1),beadresults(:,12),3,1); %running average smoothing function that takes into account breaks in FCB deployments
    
    %make a 'beadmatch' equivalent: we will use the bead SSC mean...and
    %make sure that for a given day, the same value is used...
    
    % INTERPOLATED BEAD VALUES:
    beadmatch=mvco_interpolation(beadresults(:,1),sm_bead_avgSSC,cellresultsall(to_use_SSC,1),1); %one day as a gap
    beadPEmatch=mvco_interpolation(beadresults(:,1),sm_bead_avgPE,cellresultsall(to_use_PE,1),1);
    beadCHLmatch=mvco_interpolation(beadresults(:,1),sm_bead_avgCHL,cellresultsall(to_use,1),1);
    
 
    
    %normalize if needed and record:
    allmatdate=[allmatdate; cellresultsall(to_use,1)];
    allmatdateSSC=[allmatdateSSC; cellresultsall(to_use_SSC,1)];
    allmatdatePE=[allmatdatePE; cellresultsall(to_use_PE,1)];
    
    
    beadSSC=[beadSSC; beadmatch];
    beadres=[beadres; beadresults];
    
    allsynSSC=[allsynSSC; cellSSCall(to_use_SSC,1)]; %SSC
    allsynSSCmode=[allsynSSCmode; cellSSCmodeall(to_use_SSC,1)];
    
    allsynvol = [allsynvol; cytosub_SSC2vol(cellSSCall(to_use_SSC,1)./beadmatch)]; %volume calc
    allsynvolmode = [allsynvolmode; cytosub_SSC2vol(cellSSCmodeall(to_use_SSC,1)./beadmatch)];
    
    allsynPE=[allsynPE; cellPEall(to_use_PE,1)./beadPEmatch]; %PE
    allsynPEmode=[allsynPEmode; cellPEmodeall(to_use_PE,1)./beadPEmatch];
    
    allsynCHL=[allsynCHL; cellCHLall(to_use,1)./beadCHLmatch]; %CHL
    allsynCHLmode=[allsynCHLmode; cellCHLmodeall(to_use,1)./beadCHLmatch];
    
    FCBnum=[FCBnum; FCBnumberall(to_use)];
end

%%

figure, 
plot(beadres(:,1),beadres(:,13),'.')
datetick('x','yyyy')
set(gca,'xgrid','on')
%%
figure, 
plot(allmatdateSSC,beadSSC,'.')
datetick('x','yyyy')
set(gca,'xgrid','on')
%%

figure
plot(allmatdateSSC,allsynSSCmode./beadSSC,'.')
datetick('x','yyyy')
set(gca,'xgrid','on')
%%
figure
[years,~,~,~,~,~]=datevec(allmatdateSSC);
scatter(beadSSC,allsynSSCmode./beadSSC,20,years,'filled')
jj=find(beadSSC > 2e4 & beadSSC < 3e4 & allsynSSCmode./beadSSC < 0.035);
hold on
plot(beadSSC(jj),allsynSSCmode(jj)./beadSSC(jj),'k*')
%%
figure,
scatter(beadSSC,allsynSSCmode,20,years,'filled')
xlabel('Bead SSC')
ylabel('Syn SSC')