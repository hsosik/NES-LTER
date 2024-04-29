close all;
clear all;

load('\\sosiknas1\Lab_data\Attune\cruise_data\20220806_EN688\discrete_live\outputs\SummaryTable.mat');
live = CNTable;
load('\\sosiknas1\Lab_data\Attune\cruise_data\20220806_EN688\preserved\outputs\SummaryTable.mat');
preserved = CNTable;
load('\\sosiknas1\Lab_data\Attune\cruise_data\20220806_EN688\preserved\original_outputs\SummaryTable.mat');
original = CNTable;

clear CNTable;

 

live_casts = unique(live.Cast);
preserved_casts = unique(preserved.Cast);
compare = intersect(preserved_casts, live_casts);

 figure
for i = 1:length(compare)
    single_cast_live = find(live.Cast == compare(i));
    subplot(2,4,i), plot(live.syn_per_ml(single_cast_live), live.depth_m(single_cast_live), 'r*', 'LineStyle','-')
    hold on
    single_cast_preserved = find(preserved.Cast == compare(i));
    subplot(2,4,i), plot(preserved.syn_per_ml(single_cast_preserved), preserved.depth_m(single_cast_preserved), 'k+', 'LineStyle','--')
    set(gca, 'YDir', 'reverse')
    ylim([0 100])
    xlim([0 200000])
    ylabel('Depth(m)')
    xlabel('Syn per ml')
    title({'cast' compare(i)})
    legend('live', 'Pres;CHL runs')
end


  figure
for i = 1:length(compare)
    single_cast_live = find(live.Cast == compare(i));
    subplot(2,4,i), plot(live.pro_per_ml(single_cast_live), live.depth_m(single_cast_live), 'm*', 'LineStyle','-')
    hold on
    single_cast_preserved = find(preserved.Cast == compare(i));
    subplot(2,4,i), plot(preserved.pro_per_ml(single_cast_preserved), preserved.depth_m(single_cast_preserved), 'k+', 'LineStyle','--')
    set(gca, 'YDir', 'reverse')
    ylim([0 100])
    xlim([0 500000])
    ylabel('Depth(m)')
    xlabel('Pro per ml')
    title({'cast' compare(i)})
    legend('live', 'preserved')
end

   figure
for i = 1:length(compare)
    single_cast_live = find(live.Cast == compare(i));
    subplot(2,4,i), plot(live.euk_per_ml(single_cast_live), live.depth_m(single_cast_live), 'g*', 'LineStyle','-')
    hold on
    single_cast_preserved = find(preserved.Cast == compare(i));
    subplot(2,4,i), plot(preserved.euk_per_ml(single_cast_preserved), preserved.depth_m(single_cast_preserved), 'k+', 'LineStyle','--')
    set(gca, 'YDir', 'reverse')
    ylim([0 100])
    xlim([0 20000])
    ylabel('Depth(m)')
    xlabel('Euk per ml')
    title({'cast' compare(i)})
    legend('live', 'preserved')
end

 figure
for i = 1:length(compare)
    single_cast_original = find(original.Cast == compare(i));
    subplot(2,4,i), plot(original.syn_per_ml(single_cast_original), original.depth_m(single_cast_original), 'r*', 'LineStyle','-')
    hold on
    single_cast_preserved = find(preserved.Cast == compare(i));
    subplot(2,4,i), plot(preserved.syn_per_ml(single_cast_preserved), preserved.depth_m(single_cast_preserved), 'k+', 'LineStyle','--')
    single_cast_live = find(live.Cast == compare(i));
    subplot(2,4,i), plot(live.syn_per_ml(single_cast_live), live.depth_m(single_cast_live), 'b*', 'LineStyle', ':')
    set(gca, 'YDir', 'reverse')
    ylim([0 100])
    xlim([0 200000])
    ylabel('Depth(m)')
    xlabel('Syn per ml')
    title({'cast' compare(i)})
    legend('Pres;PE runs', 'Pres; CHL runs', 'Live')
end

% figure
% for i = 1:length(compare)
%     single_cast_live = find(live.Cast == compare(i));
%     single_cast_preserved = find(preserved.Cast == compare(i));
%     plot(live.pro_per_ml(single_cast_live), live.depth_m(single_cast_live), 'g*', 'LineStyle','-')
%     hold on
%     
%     plot(preserved.pro_per_ml(single_cast_preserved), preserved.depth_m(single_cast_preserved), 'r+', 'LineStyle','--')
%     set(gca, 'YDir', 'reverse')
%     ylim([0 100])
%     ylabel('Depth(m)')
%     xlabel('Pro per ml')
%     title({'cast' compare(i)})
% end