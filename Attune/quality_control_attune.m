% There are two levels of quality control:
% 1. Primary Level (asseses quality of the data)
%     Good
%     Bad
%     Questionable/Suspicious  
QC1 = [];
% 2.Secondary Level (provides justification for quality assesement with QC tests)
%     Unexpected CV and count Ratio
%     Excessive Spike Check
%     Expert Review
QC2 = blanks(length(SynConc));
% input SynConc EukConc SynYcv EukYcv SynCount EukCount
% output QC1 QC2 spike exprev unexp   (same length as synConc, EukConc and CV value)

%% Loading in statistical data created by compile_attune

load \\sosiknas1\Lab_data\Attune\EN608\Summary\compiled_stats.mat;

%write something to check for compiled_stats and to run compile_attune if
%it is not found

%% "Excessive Spike Check"
%this returns a list of 0 for failing the spike test (spike present) and a
%1 if the point passes the spike check (no spike)

%checking for spikes in the Syn Concentration Data
%the second argument is the specificity of the test
spikeSyn = dataqc_spiketest(SynConc,.0001)

%checking for spikes in the Euk Concentration Data
%the second argument is the specificity of the test
spikeEuk = dataqc_spiketest(EukConc,.0001)

index_spikeSyn = find(spikeSyn == 0)
index_spikeEuk= find(spikeEuk == 0)

%% Unexpected CV and Count Ratios for Syn Data
% this creates a vector of ones the same length as the SynConc vector
unexpcv = ones(length(SynConc),1);

%setting the threshold for expected Cv and Count values
t_CV = 80;
t_Count = 600;

%this indexes any data points that are above the threshold (bad points)
index_threshold = find(SynYcv >= t_CV & SynCount >= t_Count);

%this gives those points (indexed as bad) a value of zero, indicating they
%have failed the test
unexpcv(index_threshold) = 0


%% "Expert Review"
%Visually I was able to identify these suspicious data points
index_visual = [512, 668, 1020, 1510, 1924, 2035, 2124, 2126, 2206, 2320, 2476, 2562];
index_visual= index_visual'

expreview = ones(length(SynConc),1);
expreview(index_visual)= 0

%% Summary Plot and Final Output

QC1 = unexpcv + spikeSyn + spikeEuk + expreview

%QC2

%% Summary Plot of Tests vs. Expert Review

figure
plot(SynConc*1000, '.-')
hold on
plot(EukConc*1000, '.-')
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
lh = legend('\itSynechococcus', 'Small eukaryotes', 'location', 'northwest');
title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
set(lh, 'fontsize', 14)
plot(index_visual,SynConc(index_visual)*1000,'md','MarkerSize',10)
plot(index_visual, EukConc(index_visual)*1000, 'md','MarkerSize',10)
plot(index_threshold,SynConc(index_threshold)*1000,'go','MarkerSize',5,'MarkerFaceColor',[.49 1 .63])
plot(index_threshold,EukConc(index_threshold)*1000,'go','MarkerSize',5,'MarkerFaceColor',[.49 1 .63])
% plot(index_spikeSyn, SynConc(index_spikeSyn)*1000, 'r*','MarkerSize',10)
plot(index_spikeEuk,EukConc(index_spikeEuk)*1000, 'r*','MarkerSize',10)
title('Concentrations and Quality Control Test Results')
lh = legend('\itSynechococcus', 'Small eukaryotes','Expert Review','Expert Review','Unexpected Ratio Syn','Unexpected Ratio  Euk','Synechococcus Spike Check', 'location', 'northwest');


%% Plotting CV vs. Counts for Unexpected CV and Count Ratio 
%plotting the Syn counts vs. %CV values
figure
plot(SynYcv, SynCount,'r.')
xlim([0 max(SynYcv)])
xlabel('%CV')
ylabel('Count')
title('CV vs. Count for Synchecoccus Data')


%% Plotting Spikes and Spike Test Results
figure
plot(SynConc*1000, '.-')
hold on
plot(EukConc*1000, '.-')
xlim([0 length(SynConc)])
ylabel('Cell concentration (ml^{-1})')
xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
lh = legend('\itSynechococcus', 'Small eukaryotes', 'location', 'northwest');
title('Concentrations and Spike Test Results')
set(lh, 'fontsize', 14)
plot(spikeSyn*25000,'r:')
plot(spikeEuk*25000, 'b:')
