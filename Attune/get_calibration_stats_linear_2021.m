% this function first looks through a folder of FCS files 
% for relationship between SSC low senstivity side scattering channel 
% and GL1 high sensitivity side scattering channel by fitting a linear
% regression to every Nth file. It saves plots to demonstrate the fit of
% each model and also saves the coefficients for each model to a table 
% which is saved to class/calibration/table.mat 

% Modified 11/2 so that data is plotted on linear, rather than log scale. 


%Inputs: 
% outpath is the path being saved to by the wrapper function. usually 
% basepath/bead_calibrated
% classpath is the location of the class files
% stepsize is the number of files to skip. Stepsize = 1 will fit a new
% linear regresion to every fsc file for which there is a class file. 
%  DIM is dimension of signal being used to estimate volume, either 'A' or 'H'
%   for ssc-a or ssc-h (area or height of signal) 




function get_calibration_stats_linear_2021(outpath, classpath, stepsize, DIM)
%DIM should be either 'A' or 'H' for ssc-a or ssc-h
         %ssc_ch_num = 3; %3 ssc-a, 12 ssc-h 
fpath = regexprep(outpath, 'bead_calibrated', 'FCS');

saverpath = [classpath 'calibration']; 
if ~exist([classpath 'calibration'], 'dir')
    mkdir([classpath 'calibration'])
end

 classlist = dir([classpath, '*.mat']); 
    %load([outpath, 'beadstat_2021.mat'])
    joint_table = []; 
    
     for counti = 1:stepsize:length(classlist)
         qc_warning = 0; 

        %load corresponding fcs file
         filename = [fpath, regexprep(classlist(counti).name, '.mat', '.fcs')];
         [fcsdat,fcshdr] = fca_readfcs(filename);
         
                  
        ssc_ch_num = strmatch(['SSC-' DIM], {fcshdr.par.name});
        gl1_ch_num = strmatch(['GL1-' DIM], {fcshdr.par.name});

         file_hv = fcshdr.par(ssc_ch_num).hv; %heidi    
    
    
        gl1_vals = fcsdat(:, gl1_ch_num);
        ssc_value = fcsdat(:,ssc_ch_num); 
   
        if DIM == 'A'
           ssch_ch_num = strmatch(['SSC-H'], {fcshdr.par.name});
           ssc_value(ssc_value<=0) = fcsdat(ssc_value<=0, ssch_ch_num); %replace negative A values with H value as proxy
        end
                    
        l_bound = 1e3;
        r_bound = 9e3; %fixed
        if r_bound < l_bound; 
            qc_warning =1; 
        end
    
        ind_to_fit = find(gl1_vals>l_bound & gl1_vals<r_bound); %not a super robust way to choose the range to fit
   
        LM = fitlm(gl1_vals(ind_to_fit), ssc_value(ind_to_fit));     
        intercept = LM.Coefficients.Estimate(1); 
        slope = LM.Coefficients.Estimate(2); 
    
        if isnan(intercept) %if data is bad, linear model doesn't work 
            qc_warning = 1; 
            scatter_value = fcsdat(:,ssc_ch_num); 
        elseif LM.Rsquared.Adjusted < .8
            qc_warning = 1; 
        end
    
    new_ssc_vals = ssc_value;
    new_ssc_vals(gl1_vals>r_bound) = 10.^[intercept + slope.*(gl1_vals(gl1_vals>r_bound))]; 
    
    %dither according to heidi's example
    t = find((gl1_vals <= r_bound & ssc_value >= 10.^(intercept + slope.*(r_bound))) | (gl1_vals >= r_bound & ssc_value <= 10.^(intercept + slope.*(r_bound)))); 
    new_ssc_vals(t(1:2:end)) = 10.^[intercept + slope.*(gl1_vals(t(1:2:end)))];
    
    calibrate_info = table; 
    calibrate_info.ssc_ch_num = ssc_ch_num; 
    calibrate_info.qc = qc_warning; 
    calibrate_info.intercept = intercept; 
    calibrate_info.slope = slope; 
    calibrate_info.rsquared = LM.Rsquared.Adjusted; 
    calibrate_info.numpoints = length(ind_to_fit); 
    calibrate_info.rightbound = r_bound; 
    
    joint_table = [joint_table; calibrate_info];

        figure(98)
        plot(gl1_vals, ssc_value, '.')
        xlabel('GL1 - low sensitivity')
        ylabel('SSC - high sensitivity')
        hold on 
        plot([1 max(gl1_vals)], [intercept (intercept + slope*(max(gl1_vals)))])
        plot(gl1_vals(t), ssc_value(t), '.b')
        plot(gl1_vals, new_ssc_vals, 'g.')
        plot(gl1_vals(ind_to_fit), ssc_value(ind_to_fit), 'r.')
        title(num2str(qc_warning))
        axis([0 1e5 0 1e6])
        print(figure(98), fullfile(saverpath, regexprep(classlist(counti).name, '.mat', '.png')), '-dpng')
        clf(98)
                 
        
     end
     
     save([saverpath '/table.mat'], 'joint_table')
     
end

     