%qualityControl_attune
basepath = '\\sosiknas1\Lab_data\Attune\EN608';

cpath = [basepath '\ExportedFCS\'];
outpath = [basepath '\Summary'];
load([outpath '\Attune'])

% Spike Check
%finds indices of points with excessive spikes
 synBad = dataqc_spiketest(Attune.Count.SynTotal,0);
 synSpike = find(synBad == 1);
 
 %finds indices of points with excessive spikes
eukBad = dataqc_spiketest(Attune.Count.EukTotal,0);
eukSpike = find(eukBad == 1);

% Unexpected CV/ Ratio
%set threshold
t_CV = 80; %identified visually 80% spread is bad in the y direction
t_Count = 500;

synRatio = find(Attune.Count.SynYCV >= t_CV & Attune.Count.SynTotal >= t_Count);
eukRatio = find(Attune.Count.SynYCV >= t_CV & Attune.Count.EukTotal >= t_Count);

%% Summary Plot of Tests

figure('units','normalized','outerposition',[0 0 1 1])
[~,ii] = sort(Attune.FCSfileinfo.matdate_start)
SynConc = Attune.Count.SynTotal(ii)./Attune.vol_analyzed(ii)
EukConc = Attune.Count.EukTotal(ii)./Attune.vol_analyzed(ii)

subplot(4,1,1)
plot(SynConc*1000, 'b.-','LineWidth',1)
hold on
plot(EukConc*1000, 'g.-','LineWidth',1)
xlim([0 length(SynConc)])
lh = legend('\itSynechococcus', 'Small eukaryotes','location', 'northwest');
title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
title('Raw Data')
set(gca, 'fontsize', 12)


% Results of Ratio Test
subplot(4,1,2)
plot(SynConc*1000, 'b.-','LineWidth',1)
hold on
plot(EukConc*1000, 'g.-','LineWidth',1)
plot(synRatio,SynConc(synRatio)*1000,'r*','MarkerSize',10,'LineWidth',1)
plot(eukRatio,EukConc(eukRatio)*1000,'r*','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
xlim([0 length(SynConc)])
lh = legend('\itSynechococcus','Small Euks','Unexpected CV Ratio', 'location', 'northwest');
title('Unexpected CV Ratio Check')
set(gca, 'fontsize', 12)


subplot(4,1,3)
plot(SynConc*1000 ,'b.-','LineWidth',1)
hold on
plot(synRatio, SynConc(synRatio)*1000, 'rx','MarkerSize',10,'LineWidth',1)
xlim([0 length(SynConc)])
lh = legend('\itSynechecoccus','Excessive Spike', 'location', 'northwest');
title('Excessive Syn Spike Check')
set(gca, 'fontsize', 12)

subplot(4,1,4)
plot(EukConc*1000, 'g.-','LineWidth',1)
hold on
plot(eukRatio,EukConc(eukRatio)*1000, 'rx','MarkerSize',10,'LineWidth',1)
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
lh = legend('Small eukaryotes','Excessive Spike', 'location', 'northwest');
title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
title('Excessive Euk Spike Check')
set(gca, 'fontsize', 12)
suptitle('Concentrations and Quality Control Test Results')

%% Expert Review

filename_reviewlist = vertcat(Attune.FCSfileinfo.filelist(synSpike)...
    ,Attune.FCSfileinfo.filelist(eukSpike)...
    ,Attune.FCSfileinfo.filelist(synRatio)...
    ,Attune.FCSfileinfo.filelist(eukRatio));

%fcs_path = '\\sosiknas1\Lab_data\Attune\EN608\ExportedFCS\'

reviewlist = cellstr(strcat(cpath,filename_reviewlist));
reviewlist2 = erase(reviewlist," ");
Attune.qc.review = zeros(length(Attune.FCSfileinfo.filelist),1);

for ii = 1:length(reviewlist2)
    [~,fcshdr,fcsdatscaled] =fca_readfcs(char(reviewlist2(ii)));
    fig = uifigure('Name','PE Signal','Position',[10 10 1000 800]);
    ax = uiaxes(fig,'Position',[30 30 790 770]);
    loglog(ax,fcsdatscaled(:,11),fcsdatscaled(:,19),'.')
    ax.XLim = [10^2  10^6];
    ax.YLim = [10^2  10^6];
    ax.XLabel.String = 'FSC-H';
    ax.YLabel.String  = 'GL1-H';
    Attune.qc.index = strmatch(char(filename_reviewlist(ii)),char(Attune.FCSfileinfo.filelist));
    bg = uibuttongroup(fig,'Position',[850 400 100 100],'Title','Options','SelectionChangedFcn',@(bg,event) bselection(bg,event,Attune));
        r1 = uiradiobutton(bg,'Text','Bad','Position',[10 60 150 15]);
        r2 = uiradiobutton(bg,'Text','Concerning','Position',[10 38 150 15]);
        r3 = uiradiobutton(bg,'Text','Good','Position',[10 15 150 15]);
              bg.Visible = 'on';
        btn = uibutton(fig,'push','Text','Next',...
               'Position',[850, 218, 100, 22],...
               'ButtonPushedFcn', @(btn,event) closefig(fig));
           btn.Visible = 'on';
        pause
end
% Create the function for the ButtonPushedFcn callback
         
function [Attune] = bselection(bg,event,Attune)
display(['Previous: ', event.OldValue.Text]);
display(['Current: ', event.NewValue.Text]);
switch event.NewValue.Text
    case 'Bad'
        Attune.qc.review(Attune.qc.index) = 3;
        display('index changed to 3');
        display('------------------');
    case 'Concerning'
        Attune.qc.review(Attune.qc.index) = 2;
        display('index changed to 2');
        display('------------------');
    case 'Good'
        Attune.qc.review(Attune.qc.index) = 0 ;
        display('index changed to 0');
        display('------------------');
end
end

function closefig(fig)
    delete(fig)
end