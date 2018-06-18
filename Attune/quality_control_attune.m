% There are two levels of quality control:
% 1. Primary Level (asseses quality of the data)
%     Good
%     Bad
%     Questionable/Suspicious    
% 2.Secondary Level (provides justification for quality assesement with QC tests)
%     Unexpected CV and count Ratio
%     Excessive Spike Check
%     Expert Review
%
% input synConc eukConc 
% output spike exprev unexp(same length as synConc, EukConc and CV value)

%% Loading in statistical data created by compile_attune

fpath = '\\sosiknas1\Lab_data\Attune\EN608\Summary\';
outpath = '\\sosiknas1\Lab_data\Attune\EN608\Summary\';

load \\sosiknas1\Lab_data\Attune\EN608\Summary\compiled_stats.mat;

%write something to check for compiled_stats and to run compile_attune if
%it is not found

%% Plotting Concentration over time
% figure
% plot(SynConc*1000, '.-')
% hold on
% plot(EukConc*1000, '.-')
% xlim([0 2563])
% ylabel('Cell concentration (ml^{-1})')
% xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
% lh = legend('\itSynechococcus', 'Small eukaryotes', 'location', 'northwest');
% title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
% set(lh, 'fontsize', 14)
% 
% %% Plotting Count over time
% figure
% plot(SynCount, '.-')
% hold on
% plot(EukCount, '.-')
% xlim([0 length(SynCount)])
% ylabel('Cell Count')
% xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
% lh = legend('\itSynechococcus', 'Small eukaryotes', 'location', 'northwest');
% title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
% set(lh, 'fontsize', 14)
% 
% %% "Expert Review"
% %Visually I was able to identify these suspicious data points
% index_visual = [512, 668, 1020, 1510, 1924, 2035, 2124, 2126, 2206, 2320, 2476, 2562];
% index_visual= index_visual'
% 
% expreview = ones(length(SynConc),1);
% expreview(index_visual)= 0
% 
% figure
% plot(SynConc*1000, '.-')
% hold on
% plot(EukConc*1000, '.-')
% xlim([0 2563])
% ylabel('Cell concentration (ml^{-1})')
% xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
% lh = legend('\itSynechococcus', 'Small eukaryotes', 'location', 'northwest');
% title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
% set(lh, 'fontsize', 14)
% plot(index_visual,SynConc(index_visual)*1000,'md','MarkerSize',10)
% plot(index_visual, EukConc(index_visual)*1000, 'md','MarkerSize',10)
% plot(index_threshold,SynConc(index_threshold)*1000,'go','MarkerSize',5,'MarkerFaceColor',[.49 1 .63])
% plot(index_threshold,EukConc(index_threshold)*1000,'go','MarkerSize',5,'MarkerFaceColor',[.49 1 .63])

%% "Excessive Spike Check"

open dataqc_spiketest.m

spikeSyn = dataqc_spiketest(SynConc,.001)
spikeEuk = dataqc_spiketest(EukConc,.0001)

%% Plotting CV vs. Counts for Unexpected CV and Count Ratio 
%plotting the Syn counts vs. %CV values
figure
plot(SynYcv, SynCount,'r.')
xlim([0 max(SynYcv)])
xlabel('%CV')
ylabel('Count')
title('Syn')

%% Setting Threshold for CV and Count
unexpcv = ones(length(SynConc),1);
t_CV = 80;
t_Count = 600;
index_threshold = find(SynYcv >= t_CV & SynCount >= t_Count);
unexpcv(index_threshold) = 0

%%
figure
plot(SynConc*1000, '.-')
hold on
plot(EukConc*1000, '.-')
xlim([0 2563])
ylabel('Cell concentration (ml^{-1})')
xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018')
lh = legend('\itSynechococcus', 'Small eukaryotes', 'location', 'northwest');
title('onshore              \leftarrow                     offshore            \rightarrow                    onshore')
set(lh, 'fontsize', 14)

plot(spikeSyn*25000,'r:')

%%
figure 
plot(EukConc*1000, '.-')
hold on
plot(spikeEuk*25000, 'b:')

QC1 = unexpcv + spikeSyn + spikeEuk + expreview

