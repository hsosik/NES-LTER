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
%   OD2setting is either 'GL1', 'SSC', or 'None', depending no where filter
%   was set during the cruise. 

function use_calibration_stats_linear(outpath, classpath, DIM, OD2setting)
%DIM should be either 'A' or 'H' for ssc-a or ssc-h
%ssc_ch_num = 3; %3 ssc-a, 12 ssc-h


fpath = regexprep(outpath, 'bead_calibrated', 'FCS');

saverpath = [classpath 'calibration_l']; %within class files we save calibraiton informatino and some figures
classlist = dir([classpath, '*.mat']);

figpath = [saverpath '/onefit'];
if ~exist(figpath, 'dir')
    mkdir(figpath)
end

calibrate = 0; 
if strcmp(OD2setting, 'GL1')
if exist([classpath '/calibration_l/table.mat'])
    calibrate = 1;
    load([classpath '/calibration_l/table.mat'])
    %get average linear model statistics from table, only for files where
    %quality check flag wasn't flagged
    joint_table.intercept(joint_table.intercept == -Inf | joint_table.intercept == Inf ) = NaN;
    joint_table.slope(joint_table.slope == -Inf | joint_table.slope == Inf) = NaN;

    intercept = nanmean(joint_table.intercept(joint_table.qc == 0));
    slope = nanmean(joint_table.slope(joint_table.qc == 0));
    R_bound = nanmean(joint_table.rightbound(joint_table.qc == 0));
end
end

if exist([outpath, 'beadstat_2021.mat'])
    load([outpath, 'beadstat_2021.mat'])
else 
    beadtype = 'PT only'; 
end

for counti = 1:length(classlist)
   
    %load corresponding fcs file
    filename = [fpath, regexprep(classlist(counti).name, '.mat', '.fcs')];
    [fcsdat, fcshdr] = fca_readfcs(filename);

    
    %find channels we want to use
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
    
    %get median 1micron bead center for this cruise over relevant hv
    %settings for this file 
    if strcmp(beadtype, 'FCB') 
        if strcmp(OD2setting, 'SSC') %on SSC OD2 cruises we can't convert to high sensitivity channel 
            %Only Od2 measurements are available, use those 
            bead_value = [beadstat_2021.OD2_hv beadstat_2021.OD2centers(:,2)];
            bead_value = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value
            
        elseif strcmp(OD2setting, 'GL1')
            %use OD2 measurements to project to NoOD2 values
            bead_value = [beadstat_2021.OD2_hv beadstat_2021.OD2centers(:,2)];
            bead_value_to_convert = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value on GL1
            bead_value = [beadstat_2021.OD2_hv beadstat_2021.NoOD2centers(:,2)];
            bead_value = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value on SSC
            
        else %OD2setting is 'none'
            %use noOD2 value
            bead_value = [beadstat_2021.NoOD2_hv beadstat_2021.NoOD2centers(:,2)];
            bead_value = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value
            
        end
    elseif strcmp(beadtype, 'PT') % no FCB measurements, but we have PT beads measured on seawater settings
        %get 2nd cluster values, 3.2 um beads. 
         if strcmp(OD2setting, 'SSC') %on SSC OD2 cruises we can't convert to high sensitivity channel 
            %Only Od2 measurements are available, use those 
            bead_value = [beadstat_2021.OD2_hv beadstat_2021.OD2centers(:,2)];
            bead_value = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value
            
        elseif strcmp(OD2setting, 'GL1')
            %use OD2 measurements to project to NoOD2 values
            bead_value = [beadstat_2021.OD2_hv beadstat_2021.OD2centers(:,2)];
            bead_value_to_convert = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value on GL1
            bead_value = [beadstat_2021.OD2_hv beadstat_2021.NoOD2centers(:,2)];
            bead_value = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value on SSC
            
        else %OD2setting is 'none'
            %use noOD2 value
            bead_value = [beadstat_2021.NoOD2_hv beadstat_2021.NoOD2centers(:,2)];
            bead_value = nanmedian(bead_value(bead_value(:,1)==file_hv,2)); %bead_value for 3.2 micron beads 
            
         end
        
        
         %now convert beadvalue of 3.2 um bead to what it would be for 1 um
         %bead? 
         bead_value_1 = 10.^(log10(bead_value) - (log10(3.2/1)./1.336)); %slope from bead experiment March 2019 fit to unscaled attune, neutral settings
         bead_value = bead_value_1; 
    end
    
    
    %if using SSC to GL1 calibration
    if calibrate == 1
        l_bound = 5e2;
        r_bound = R_bound;
        
        %replace high values with estimates from low sensitivity channel
        new_ssc_vals = ssc_value;
        new_ssc_vals(gl1_vals>r_bound) = [intercept + slope.*(gl1_vals(gl1_vals>r_bound))];
        
        %dither according to heidi's example
        t = find((gl1_vals <= r_bound & ssc_value >= intercept + slope.*(r_bound)) | (gl1_vals >= r_bound & ssc_value <= (intercept + slope.*(r_bound))));
        new_ssc_vals(t(1:2:end)) = [intercept + slope.*(gl1_vals(t(1:2:end)))];
        
        scatter_value = new_ssc_vals;
        %also convert bead value if it is large
        if bead_value_to_convert>r_bound
            bead_value = [intercept + slope.*(bead_value_to_convert)];
        end %if not, bead_value is SSC no filter measurement
          
    if ~rem(counti,50) %save some images occasionally
        disp([num2str(counti) ' of ' num2str(length(classlist))])
        
        figure(98)
        plot(gl1_vals, ssc_value, '.')
        xlabel('GL1 - low sensitivity')
        ylabel('SSC - high sensitivity')
        hold on
        plot([1 max(gl1_vals)], [intercept (intercept + slope*max(gl1_vals))])
        plot(gl1_vals(t), ssc_value(t), '.b')
        plot(gl1_vals, new_ssc_vals, 'g.')
        axis([0 1e5 0 6e6])
        
        print(figure(98), fullfile(figpath, regexprep(classlist(counti).name, '.mat', '.png')), '-dpng')
        clf(98)
    end
    
    else %not calibration for OD2setting == GL1

        scatter_value = ssc_value;  

        if strcmp(OD2setting, 'SSC')
        % do we want to make the adjustment as learned from hv experiments
        if file_hv < 230
        %yes we do! 
         %"forward" projection
            P1  =  0.0004172 + file_hv *(-3.483e-06) + (7.284e-09) * file_hv^2 ;
            scatter_value = (scatter_value.^2).*P1 + scatter_value.*0.05; 
            %do the same thing to our bead value
            bead_value = (bead_value.^2).*P1 + bead_value.*0.05; 
        end
        
        end
        
    end
    
    
    %now convert modified scatter_values to volumes & save results
    if DIM == 'A'
        volume = 10.^(1.24*log10(scatter_value./bead_value) + 1.064); %based on linear fit to scaled ssch on OD2 filter March 2019
        volumestring = '10.^(1.24*log10(scatter_value./bead_value) + 1.064';
    elseif DIM == 'H'
        volume = 10.^(1.4225*log10(scatter_value./bead_value) + 1.1432);
        volumestring = '10.^(1.4225*log10(scatter_value./bead_value) + 1.1432)'; %based on linear fit to scaled ssca on OD2 filter March 2019
        %below is a conversion we were using at one point? 
        %volume = 10.^(1.2232*log10(ssca./fcb_mean) + 1.0868);
        
        %treat H values differently from A values
        volume(negA_ind) = 10.^(1.24*log10(scatter_value(negA_ind)./bead_value) + 1.064); 
    end
    
    if calibrate == 1
    vol_notes = {strcat('calibrated: ', string(datetime())); strcat('using SSC-', DIM, ' and GL1 Linear Scale Fit: right bound ', num2str(R_bound), 'intercept ', num2str(intercept), 'slope ', num2str(slope));
        volumestring};
    else 
        vol_notes = {strcat('calibrated: ', string(datetime())); strcat('using SSC-', DIM, 'No Linear Fit to GL1   ', volumestring)};
    end
    save([classpath classlist(counti).name], 'volume', 'ssc_value', 'vol_notes', 'bead_value', 'negA_ind', 'file_hv', '-append' )
        

    
end



end





