%this sets the current directory to the file with the matlab scripts for
%processing the attune data.
basepath = '\\sosiknas1\Lab_data\Attune\EN608\';
%Input: Summary from code that process attune ran

% if exist(Attune) == 1
%     load [basepath 'Summary\Attune']
% else 
%     disp('process_attune needs to run first')
%     process_attune(basepath)
% end
%run process attune if it hasn't been run yet

%% "Excessive Spike Check"
%this returns a list of 0 for failing the spike test (spike present) and a
%1 if the point passes the spike check (no spike)

%checking for spikes in the Syn Concentration Data
%the second argument is the specificity of the test

spikeSyn = dataqc_spiketest(x,0);

%checking for spikes in the Euk Concentration Data
%the second argument is the specificity of the test
index_x = find(x == 1);


%% Unexpected CV and Count Ratios for Syn Data
% this creates a vector of ones the same length as the SynConc vector

YCV = std(x)./mean(x)

% figure
% plot(SynYcv, SynCount, 'r.')

unexpcv = zeros(length(x),1);

%setting the threshold for expected Cv and Count values
t_CV = 80;
t_Count = 600;

%this indexes any data points that are above the threshold (bad points)
index_threshold = find(YCV >= t_CV & SynCount >= t_Count);

%this gives those points (indexed as bad) a value of zero, indicating they
%have failed the test
unexpcv(index_threshold) = 1

%% Summary Plot of Tests

figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,1,1)
plot(x*1000, 'b.-','LineWidth',3)
% hold on
% plot(EukConc*1000, 'g.-','LineWidth',3)
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
lh = legend('\itSynechococcus', 'Small eukaryotes','location', 'northwest');
% title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
% title('Raw Data')
set(gca, 'fontsize', 25)

subplot(4,1,2)
plot(x*1000, 'b.-','LineWidth',3)
hold on
plot(index_threshold,x(index_threshold)*1000,'r*','MarkerSize',20,'LineWidth',2)
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
lh = legend('\itSynechococcus','Unexpected CV Ratio', 'location', 'northwest');
% title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
% title('Unexpected CV Ratio Check')
set(gca, 'fontsize', 25)
%plot(index_threshold,EukConc(index_threshold)*1000,'go','MarkerSize',5,'MarkerFaceColor',[.49 1 .63])

subplot(4,1,3)
figure
plot(SynConc*1000 ,'b.-','LineWidth',3)
hold on
plot(index_spikeSyn, SynConc(index_spikeSyn)*1000, 'r*','MarkerSize',20,'LineWidth',2)
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
lh = legend('\itSmall Eukaryotes','Excessive Spike', 'location', 'northwest');
% title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
% title('Excessive Syn Spike Check')
set(lh, 'fontsize', 24)

subplot(4,1,4)
figure
plot(EukConc*1000, 'g.-','LineWidth',3)
hold on
plot(index_spikeSyn,EukConc(index_spikeSyn)*1000, 'r*','MarkerSize',20,'LineWidth',3)
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
lh = legend('Small eukaryotes','Euk Spikes', 'location', 'northwest');
% title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
% title('Excessive Euk Spike Check')
set(gca, 'fontsize', 25)
suptitle('Concentrations and Quality Control Test Results')

%% For Further Expert Review

filename_reviewlist = vertcat(fcsfile_syn(index_spikeSyn), fcsfile_syn(index_spikeEuk),fcsfile_syn(index_threshold));
%fcs_path = '\\sosiknas1\Lab_data\Attune\EN608\ExportedFCS\'
fcs_path = '\\sosiknas1\Backup\SPIROPA\20180414_AR29\Attune\FCSexport'
exp = zeros(length(SynConc),1);x

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
SynConc_plot = SynConc;

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