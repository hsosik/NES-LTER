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

function calibrate_two_channel(outpath, classpath, DIM, OD2setting)
%DIM should be either 'A' or 'H' for ssc-a or ssc-h
%ssc_ch_num = 3; %3 ssc-a, 12 ssc-h

saverpath = [classpath 'calibration']; %within class files we save calibraiton informatino and some figures
classlist = dir([classpath, '*.mat']);

figpath = [saverpath '/two_channel'];
if ~exist(figpath, 'dir')
    mkdir(figpath)
end

load([outpath, 'beadstat_2021.mat'])


for counti = 1:length(classlist)
   
    %load corresponding fcs file
    filename = [fpath, regexprep(classlist(counti).name, '.mat', '.fcs')];
    [fcsdat, fcshdr] = fca_readfcs(filename);

    sc_ch_num = strmatch(['SSC-' DIM], {fcshdr.par.name});
    gl1_ch_num = strmatch(['GL1-' DIM], {fcshdr.par.name});
    gl1_vals = fcsdat(:, gl1_ch_num);
    ssc_value = fcsdat(:,ssc_ch_num);

    if DIM == 'A'
        negA_ind = ssc_value<=0; %keep track of particles for which area is neagtive. 
        ssch_ch_num = strmatch(['SSC-H'], {fcshdr.par.name});
        ssc_value(negA_ind) = fcsdat(negA_ind, ssch_ch_num); %replace negative A values with H value as proxy
    end
