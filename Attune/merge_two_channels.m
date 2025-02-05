%This function will convert bead values to volumes separately on each
%scattering channel GL1 and SSC,
%make plots to see which SSC channel is more suitable for a particle, and
%chose size cut off.
%Output a new version of classfiles

%%Inputs: 
%   outpath is outpath from wrapper script. usually 'basepath/bead_calibrated'
%   classpath is where class files are 
%   DIM is dimension of signal being used to estimate volume, either 'A' or 'H'
%   for ssc-a or ssc-h (area or height of signal) 
%   OD2setting is either 'GL1', 'SSC', or 'None', depending no where filter
%   was set during the cruise. 

function merge_two_channels(outpath, classpath, DIM, OD2setting)
%DIM should be either 'A' or 'H' for ssc-a or ssc-h
%ssc_ch_num = 3; %3 ssc-a, 12 ssc-h

figpath = [classpath 'merge']; %within class files we save calibraiton informatino and some figures
classlist = dir([classpath, '*.mat']);

%figpath = [saverpath '/two_channel'];
if ~exist(figpath, 'dir')
    mkdir(figpath)
end

load([outpath, 'beadstat_2021.mat'])
fpath = regexprep(outpath, 'bead_calibrated', 'FCS');

%DIM == 'A';

for counti = 1:length(classlist);
   
    %load corresponding fcs file
    filename = [fpath, regexprep(classlist(counti).name, '.mat', '.fcs')];
    [fcsdat, fcshdr] = fca_readfcs(filename);

    ssc_ch_num = strmatch(['SSC-' DIM], {fcshdr.par.name});
    gl1_ch_num = strmatch(['GL1-' DIM], {fcshdr.par.name});
    gl1_vals = fcsdat(:, gl1_ch_num);
    ssc_value = fcsdat(:,ssc_ch_num);

    if DIM == 'A'
        negA_ind = ssc_value<=0; %keep track of particles for which area is neagtive. 
        ssch_ch_num = strmatch(['SSC-H'], {fcshdr.par.name});
        ssc_value(negA_ind) = fcsdat(negA_ind, ssch_ch_num); %replace negative A values with H value as proxy
    end



file_hv = fcshdr.par(ssc_ch_num).hv;
bead_value = [beadstat_2021.OD2_hv beadstat_2021.OD2centers(:,2)];
bead_value_to_convert = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value on GL1
bead_value = [beadstat_2021.OD2_hv beadstat_2021.NoOD2centers(:,2)];
bead_value = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value on SSC

overlap = [.4 .7];
overlap2 = [.6 1];
overlap3 = [.8 1.2];
temp = gl1_vals./bead_value_to_convert;
diff = (gl1_vals./bead_value_to_convert - ssc_value./bead_value);
above = find((((temp - ssc_value./bead_value) > 0)) & (temp < overlap(2)) & (temp > overlap(1)));
below = find((((temp - ssc_value./bead_value) < 0)) & (temp < overlap(2)) & (temp > overlap(1)));

above2 = find((((temp - ssc_value./bead_value) > 0)) & (temp < overlap(2)) & (temp > overlap(1)));
below2 = find((((temp- ssc_value./bead_value) < 0)) & (temp < overlap(2)) & (temp > overlap(1)));

above3 = find((((temp - ssc_value./bead_value) > 0)) & (temp < overlap3(2)) & (temp > overlap3(1)));
below3 = find((((temp- ssc_value./bead_value) < 0)) & (temp < overlap3(2)) & (temp > overlap3(1)));






     

figure(98)
set(gcf, 'Position', [87 100 1000 800])
        subplot(3,2,1) 
        plot(gl1_vals./bead_value_to_convert, ssc_value./bead_value, '.', 'MarkerSize', 2)
        hold on
        %plot(gl1_vals(above)./bead_value_to_convert, ssc_value(above)./bead_value, 'r.', 'MarkerSize', 2)
        %plot(gl1_vals(below)./bead_value_to_convert, ssc_value(below)./bead_value, 'r.', 'MarkerSize', 2)
        xlabel('GL1 - low sensitivity; bead normalized')
        ylabel('SSC - high sensitivity; bead normalized')
        hold on 
        %axis([-.5 10 -.5 10])
        axis([-.1 3 -.1 3])
        line(xlim,xlim, 'color', 'r', 'linewidth', 2)
        line([overlap(1) overlap(1)], ylim, 'linewidth', 2)
        line([overlap(2) overlap(2)], ylim, 'linewidth', 2)
        line([overlap2(1) overlap2(1)], ylim,  'linewidth', 1, 'lineStyle', '--', 'Color', 'r')
        line([overlap2(2) overlap2(2)], ylim, 'linewidth', 1, 'lineStyle', '--', 'Color', 'r')
        line([overlap3(1) overlap3(1)], ylim, 'linewidth', 2, 'lineStyle', ':', 'Color', 'g')
        line([overlap3(2) overlap3(2)], ylim, 'linewidth', 2, 'lineStyle', ':', 'Color', 'g')
        grid

%         subplot(2,2,2); 
%         plot(gl1_vals./bead_value_to_convert, ssc_value./bead_value, '.', 'MarkerSize', 4)
%         xlabel('GL1 - low sensitivity; bead normalized')
%         ylabel('SSC - high sensitivity; bead normalized')
%         hold on 
%         axis([.4 1.4 .4 1.4])
%         line(xlim,xlim, 'color', 'r', 'linewidth', 2)
%         line([.6 .6], ylim, 'linewidth', 2)
%         line([1 1], ylim, 'linewidth', 2)
%         grid

        subplot(3,2,3);
        plot(ssc_value./bead_value, (gl1_vals./bead_value_to_convert - ssc_value./bead_value), '.', 'MarkerSize', 2)
        xlabel('SSC - high sensitivity; bead normalized')
        ylabel('bead normalized difference')
       
          line([overlap(1) overlap(1)], ylim, 'linewidth', 2)
        line([overlap(2) overlap(2)], ylim, 'linewidth', 2)
        xlim([.1 2])
        ylim([-1 1])
        set(gca, 'xscale', 'log')
        line([.1 2], [0 0], 'color', 'r', 'linewidth', 2)
        line([overlap2(1) overlap2(1)], ylim,  'linewidth', 1, 'lineStyle', '--', 'Color', 'r')
        line([overlap2(2) overlap2(2)], ylim, 'linewidth', 1, 'lineStyle', '--', 'Color', 'r')
        line([overlap3(1) overlap3(1)], ylim, 'linewidth', 2, 'lineStyle', ':', 'Color', 'g')
        line([overlap3(2) overlap3(2)], ylim, 'linewidth', 2, 'lineStyle', ':', 'Color', 'g')
        grid

        subplot(3,2,5);
        plot(gl1_vals./bead_value_to_convert, (gl1_vals./bead_value_to_convert - ssc_value./bead_value), '.', 'MarkerSize', 2)
        xlabel('GL1 - low sensitivity; bead normalized')
        ylabel('bead normalized difference')
        
         line([overlap(1) overlap(1)], ylim, 'linewidth', 2)
        line([overlap(2) overlap(2)], ylim, 'linewidth', 2)
        xlim([.1 2])
        ylim([-1 1])
        set(gca, 'xscale', 'log')
        line([.01 2], [0 0], 'color', 'r', 'linewidth', 2)
        line([overlap2(1) overlap2(1)], ylim,  'linewidth', 1, 'lineStyle', '--', 'Color', 'r')
        line([overlap2(2) overlap2(2)], ylim, 'linewidth', 1, 'lineStyle', '--', 'Color', 'r')
        line([overlap3(1) overlap3(1)], ylim, 'linewidth', 2, 'lineStyle', ':', 'Color', 'g')
        line([overlap3(2) overlap3(2)], ylim, 'linewidth', 2, 'lineStyle', ':', 'Color', 'g')
        grid

       
subplot(3,2,2)

    hist((diff([below; above])))
    %hist((diff(above)))
    line([0 0], ylim, 'linewidth', 2, 'Color', 'r')
    title(overlap)
xlim([-.6 .6])
title(['bead normalized difference' string(overlap)])

    subplot(3,2,4)
     hist((diff([below2; above2])))
    %hold on
    %hist((diff(above2)))
    line([0 0], ylim, 'linewidth', 2, 'Color', 'r')
    title(overlap2)
    xlim([-.6 .6])
    title(['bead normalized difference' string(overlap2)])

    subplot(3,2,6)
    hist((diff([below3; above3])))
    %hold on
    %hist((diff(above3)))
    line([0 0], ylim, 'linewidth', 2, 'Color', 'r')
     title(overlap3)
     xlim([-.6 .6])
     title(['bead normalized difference' string(overlap3)])



        print(figure(98), fullfile(figpath, regexprep(classlist(counti).name, '.mat', '.png')), '-dpng')
        clf(98)

%         volume = 10.^(1.24*log10(scatter_value./bead_value) + 1.064); %based on linear fit to scaled ssch on OD2 filter March 2019
%         volumestring = '10.^(1.24*log10(scatter_value./bead_value) + 1.064';
%      
%         %treat H values differently from A values
%         volume(negA_ind) = 10.^(1.4225*log10(scatter_value(negA_ind)./bead_value) + 1.1432); 
% 
% 
%  %now convert modified scatter_values to volumes & save results
%   volume = 10.^(1.24*log10(scatter_value./bead_value) + 1.064); %based on linear fit to scaled ssch on OD2 filter March 2019
%   volumestring = '10.^(1.24*log10(scatter_value./bead_value) + 1.064';
%      
%   %treat H values differently from A values
%   volume(negA_ind) = 10.^(1.4225*log10(scatter_value(negA_ind)./bead_value) + 1.1432);


end