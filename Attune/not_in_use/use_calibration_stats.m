% This function takes the output of get_calibration_stats and uses them to 
%   a) convert GL1 values to SSC values for all triggers in the FCS files 
%   b) save a handful of figures to the "onefit" subfolder of the calibration 
%       directory for quality control, to show how average linear regression
%       applies to data it wasn't fit to. 
%   c) convert SSC to Volume estimates for each trigger and append volume to
%       class files. 

%Inputs: 
%   outpath is outpath from wrapper script. usually 'basepath/bead_calibrated'
%   classpath is where class files are 
%   DIM is dimension of signal being used to estimate volume, either 'A' or 'H'
%   for ssc-a or ssc-h (area or height of signal) 


function use_calibration_stats(outpath, classpath, DIM)
%DIM should be either 'A' or 'H' for ssc-a or ssc-h
%ssc_ch_num = 3; %3 ssc-a, 12 ssc-h

fpath = regexprep(outpath, 'bead_calibrated', 'FCS');

saverpath = [classpath 'calibration']; %within class files we save calibraiton informatino and some figures
classlist = dir([classpath, '*.mat']);
load([outpath, 'beadstat.mat'])
load([classpath '/calibration/table.mat'])

figpath = [saverpath '/onefit'];
if ~exist(figpath, 'dir')
    mkdir(figpath)
end

%get average linear model statistics from table, only for files where
%quality check flag wasn't flagged
joint_table.intercept(joint_table.intercept == -Inf | joint_table.intercept == Inf ) = NaN;
joint_table.slope(joint_table.slope == -Inf | joint_table.slope == Inf) = NaN;

intercept = nanmean(joint_table.intercept(joint_table.qc == 0));
slope = nanmean(joint_table.slope(joint_table.qc == 0));
R_bound = nanmean(joint_table.rightbound(joint_table.qc == 0));


for counti = 1:length(classlist)
    
    %load corresponding fcs file
    filename = [fpath, regexprep(classlist(counti).name, '.mat', '.fcs')];
    [fcsdat,fcshdr] = fca_readfcs(filename);
    
    
    %find channels we want to use
    ssc_ch_num = strmatch(['SSC-' DIM], {fcshdr.par.name});
    gl1_ch_num = strmatch(['GL1-' DIM], {fcshdr.par.name});
    gl1_vals = fcsdat(:, gl1_ch_num);
    ssc_value = fcsdat(:,ssc_ch_num);
    
    file_hv = fcshdr.par(ssc_ch_num).hv; %heidi
    
    %get bead value of interest for this file
    if DIM == 'A'
        bead_value = bead_med_ssca(bead_med_ssca(:,1)==file_hv,2); %bead_value
    elseif DIM == 'H'
        bead_value = bead_med_ssch(bead_med_ssch(:,1)==file_hv,2); %bead_value
    end
    
    
    %replace negative values with zero
    ssc_value(ssc_value<0) =0;
    l_bound = 5e2;
    r_bound = R_bound;
    
    %replace high values with estimates from low sensitivity channel
    new_ssc_vals = ssc_value;
    new_ssc_vals(gl1_vals>r_bound) = 10.^[intercept + slope.*(log10(gl1_vals(gl1_vals>r_bound)))];
    
    %dither according to heidi's example
    t = find((gl1_vals <= r_bound & ssc_value >= 10.^(intercept + slope.*(log10(r_bound)))) | (gl1_vals >= r_bound & ssc_value <= 10.^(intercept + slope.*(log10(r_bound)))));
    new_ssc_vals(t(1:2:end)) = 10.^[intercept + slope.*(log10(gl1_vals(t(1:2:end))))];
    
    scatter_value = new_ssc_vals;
    
    
    if ~rem(counti,50) %save some images occasionally
        disp([num2str(counti) ' of ' num2str(length(classlist))])
        
        figure(98)
        loglog(gl1_vals, ssc_value, '.')
        xlabel('GL1 - low sensitivity')
        ylabel('SSC - high sensitivity')
        hold on
        plot([1 max(gl1_vals)], 10.^[intercept (intercept + slope*log10(max(gl1_vals)))])
        loglog(gl1_vals(t), ssc_value(t), '.b')
        loglog(gl1_vals, new_ssc_vals, 'g.')
        
        print(figure(98), fullfile(figpath, regexprep(classlist(counti).name, '.mat', '.png')), '-dpng')
        clf(98)
    end
    
    %now convert modified scatter_values to volumes & save results
    if DIM == 'A'
        volume = 10.^(1.357*log10(scatter_value./bead_value) + 1.315); %based on linear fit to scaled ssca on OD2 filter March 2019
        volumestring = '10.^(1.357*log10(scatter_value./bead_value) + 1.315)';
    elseif DIM == 'H'
        volume = 10.^(1.4225*log10(scatter_value./bead_value) + 1.1432);
        volumestring = '10.^(1.4225*log10(scatter_value./bead_value) + 1.1432)';
        %or possibly 
        %volume = 10.^(1.2232*log10(ssca./fcb_mean) + 1.0868);
    end
    vol_notes = {strcat('calibrated: ', string(datetime())); strcat('using SSC-', DIM, ' and GL1 Linear Fit: right bound ', num2str(R_bound), 'intercept ', num2str(intercept), 'slope ', num2str(slope));
        volumestring};
    save([classpath classlist(counti).name], 'volume', 'ssc_value', 'vol_notes', '-append')
        

end

end





