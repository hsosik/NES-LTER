cd C:\Users\mps48\Documents\GitHub\NES-LTER\Attune

% There are two levels of quality control:
% 1. Primary Level (asseses quality of the data)
%     Good
%     Bad
%     Questionable/Suspicious  
% 2.Secondary Level (provides justification for quality assesement with QC tests)
%     Unexpected CV and count Ratio for Synchecoccus
%     Excessive Spike Check on Synecheccous Data
%     Excessive Spike Check on Small Eukaryotes Data
%     Expert Review
% input SynConc EukConc SynYcv EukYcv SynCount EukCount from compile_attune
%
% output  spike_syn spike_exp   (same length as synConc, EukConc and CV value)

%% Loading in statistical data created by compile_attune
if ~exist('\\sosiknas1\Lab_data\Attune\EN608\Summary\compiled_stats.mat','file')
    open compile_attune
else
load \\sosiknas1\Lab_data\Attune\EN608\Summary\compiled_stats.mat;
end

%write something to check for compiled_stats and to run compile_attune if
%it is not found
%% "Excessive Spike Check"
%this returns a list of 0 for failing the spike test (spike present) and a
%1 if the point passes the spike check (no spike)

%checking for spikes in the Syn Concentration Data
%the second argument is the specificity of the test

spikeSyn = dataqc_spiketest(SynConc,0);

%checking for spikes in the Euk Concentration Data
%the second argument is the specificity of the test
spikeEuk = dataqc_spiketest(EukConc,0);

index_spikeSyn = find(spikeSyn == 1);
index_spikeEuk= find(spikeEuk == 1);


%% Unexpected CV and Count Ratios for Syn Data
% this creates a vector of ones the same length as the SynConc vector
unexpcv = zeros(length(SynConc),1);

%setting the threshold for expected Cv and Count values
t_CV = 80;
t_Count = 600;

%this indexes any data points that are above the threshold (bad points)
index_threshold = find(SynYcv >= t_CV & SynCount >= t_Count);

%this gives those points (indexed as bad) a value of zero, indicating they
%have failed the test
unexpcv(index_threshold) = 1

%% Summary Plot of Tests

figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,1,1)
plot(SynConc*1000, 'b.-')
hold on
plot(EukConc*1000, 'g.-')
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
lh = legend('\itSynechococcus', 'Small eukaryotes','location', 'northwest');
title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
title('Raw Data')
set(lh, 'fontsize', 10)

subplot(4,1,2)
plot(SynConc*1000, 'b.-')
hold on
plot(index_threshold,SynConc(index_threshold)*1000,'r*','MarkerSize',10)
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
lh = legend('\itSynechococcus','Unexpected CV Ratio', 'location', 'northwest');
title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
title('Unexpected CV Ratio Check')
set(lh, 'fontsize', 10)
%plot(index_threshold,EukConc(index_threshold)*1000,'go','MarkerSize',5,'MarkerFaceColor',[.49 1 .63])

subplot(4,1,3)
plot(SynConc*1000, 'b.-')
hold on
plot(index_spikeSyn, SynConc(index_spikeSyn)*1000, 'rv','MarkerSize',10)
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
lh = legend('\itSynechococcus','Syn Spike', 'location', 'northwest');
title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
title('Excessive Syn Spike Check')
set(lh, 'fontsize', 10)

subplot(4,1,4)
plot(EukConc*1000, 'g.-')
hold on
plot(index_spikeEuk,EukConc(index_spikeEuk)*1000, 'rv','MarkerSize',10)
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
lh = legend('Small eukaryotes','Euk Spikes', 'location', 'northwest');
title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
title('Excessive Euk Spike Check')
set(lh, 'fontsize', 10)
suptitle('Concentrations and Quality Control Test Results')

%% For Further Expert Review

filename_reviewlist = vertcat(fcsfile_syn(index_spikeSyn), fcsfile_syn(index_spikeEuk),fcsfile_syn(index_threshold));
fcs_path = '\\sosiknas1\Lab_data\Attune\EN608\ExportedFCS\'
exp = zeros(length(SynConc),1);

for ii = 1:length(filename_reviewlist)
    [~,fcshdr,fcsdatscaled] =fca_readfcs(char(fullfile(fcs_path,filename_reviewlist(ii))));
    fig = uifigure('Name','PE Signal','Position',[10 10 1000 800]);
    ax = uiaxes(fig,'Position',[30 30 790 770]);
    loglog(ax,fcsdatscaled(:,11),fcsdatscaled(:,19),'.')
    ax.XLim = [10^2  10^6]
    ax.YLim = [10^2  10^6]
    ax.XLabel.String = 'FSC-H'
    ax.YLabel.String  = 'GL1-H'
    index = strmatch(char(filename_reviewlist(ii)),char(fcsfile_syn))
    F = struct('exp',exp, 'index',index);
    bg = uibuttongroup(fig,'Position',[850 400 100 100],'Title','Options','SelectionChangedFcn',@(bg,event) bselection(bg,event,F));
        r1 = uiradiobutton(bg,'Text','Bad','Position',[10 60 150 15]);
        r2 = uiradiobutton(bg,'Text','Concerning','Position',[10 38 150 15]);
        r3 = uiradiobutton(bg,'Text','Good','Position',[10 15 150 15]);
              bg.Visible = 'on';
              pause
end
  
%% "Visual ID"
%Visually I was able to identify these suspicious data points
index_visual = [512, 668, 1020, 1510, 1924, 2035, 2124, 2126, 2206, 2320, 2476, 2562];
index_visual= index_visual';

expreview = zeros(length(SynConc),1);
expreview(index_visual)= 1;

%%
function [F] = bselection(bg,event,F)
display(['Previous: ', event.OldValue.Text]);
display(['Current: ', event.NewValue.Text]);
switch event.NewValue.Text
    case 'Bad'
        exp(F.index) = 3;
        display('index changed to 3');
        display('------------------');
    case 'Concerning'
        exp(F.index) = 2;
        display('index changed to 2');
        display('------------------');
    case 'Good'
        exp(F.index) = 0 ;
        display('index changed to 0');
        display('------------------');
end
end