%% Make a movie of attune data given FCS files and class files

%inputs : 
% fpath = string path of folder where FCS files are, 
    % (path should end with a \)
% classpath = string path of folder where class files are, end in \
% framemaker = string of function name we want to use to plot data. 
    % may change from cruise to cruise, but we will try to make new 
    % framemaker scripts when needed rather than editing old ones. 

 %Summary : 
 % This script will look through class file folder, sort files in 
 % chronological order, and  make a movie with the data from all the files 
 % (This assumes that files we wish to exclude (e.g. grazer experiments) 
 % have already been excluded from the class assingment step)
 % This script will open fcs & class files, and pass data to framemaker
 % function. Then save the frames to the ouput movie Attune_cyto_vid.avi 
 % which will be saved into the classpath. 
 
 %This script assumes FCSfileinfo.mat exists one dir up from class files
 
 
 function [] = attune_lter_moviemaker(fpath, classpath, OD2setting, framemaker, channels, stepsize)
%  For example:
%attune_lter_moviemaker('\\sosiknas1\Lab_data\Attune\cruise_data\20180414_AR29\FCS\', '\\sosiknas1\Lab_data\Attune\cruise_data\20180414_AR29\bead_calibrated\class\', 'GL1', 'make_movieframe_EN649')

if ~exist('stepsize','var')
     % third parameter does not exist, so default it to something
      stepsize = 1;
end
 
%% First, some file management 
% get list of all class files
temp = dir([classpath, '*.mat']); 
classfilelist = {temp.name}; 
clear temp 

%also make corresponding list of FCS files of interest
filelist = regexprep(classfilelist,'mat', 'fcs'); 

% initialize movie output
v = VideoWriter([classpath 'Attune_cyto_vid_new.avi']); 
v.FrameRate = 10; 
open(v)

%let's make sure the files are run in chronological order
%assume FCSfileinfo.mat exists one dir up from class files
temp_path = regexprep(classpath, 'class\', '');
temp = load([temp_path 'FCSfileinfo']);
[~,a,b] = intersect(filelist, temp.FCSfileinfo.fcslist);
matdate = NaN(size(filelist));
QC_flags = matdate; 
matdate(a) = temp.FCSfileinfo.matdate_start(b);
QC_flags(a) = temp.FCSfileinfo.QC_flag(b); 
[matdate,sort_ind] = sort(matdate);
filelist = filelist(sort_ind);
classfilelist = classfilelist(sort_ind);
QC_flags = QC_flags(sort_ind); 


%% now go through list and make frames to save to video
for count = 1:stepsize:length(filelist)
    if ~rem(count,10)
        disp([num2str(count) ' of ' num2str(length(filelist))])
    end
    
    filename = [filelist{count}];
    disp(filename)
    
    %Load in both raw fcs file and class assignments
    [fcsdat,fcshdr] = fca_readfcs([fpath, filename]);
    load([classpath, classfilelist{count}], 'class')
    
    %% one step of data cleaning
    % assign scattering channels
    if strcmp(OD2setting, 'GL1')
        ssc_name = 'GL1'; 
    else
        ssc_name = 'SSC'; 
    end
     ssc_ch = strmatch([ssc_name '-A'], {fcshdr.par.name});
     ssch = strmatch([ssc_name '-H'], {fcshdr.par.name});

    % correct for negative SSC-A values
    cf = fitlm(fcsdat(fcsdat(:,ssch)<1000,ssch), fcsdat(fcsdat(:,ssch)<1000, ssc_ch), 'Intercept', false);
    cf = cf.Coefficients.Estimate;
    fcsdat(fcsdat(:,ssc_ch)<0, ssc_ch) = cf*fcsdat(fcsdat(:,ssc_ch)<0, ssch);
    
    % call the function to make the plots and getframe     
    eval(['Frame = ', framemaker, '(fcsdat, fcshdr, class, channels, QC_flags(count))']);
            
    % add to movie
    writeVideo(v, Frame);

    
    clear class fcsdat fcshdr
  
end

close(v) 
disp(['Result file saved:' classpath 'Attune_cyto_vid.avi'])

end
