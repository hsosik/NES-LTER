%This script is designed to be able to process a cruise of attune data 
%from start to finish or with each step optional. 
% This is pretty inefficient, because we have to open each of the
% fcs files in a loop multiple times. But it is totally modular. 

%Steps: 
%1. Generate FCSfileinfo 
        % this is a table of file names and dates, with volume and a
        % quality flag. For this table, QC_flag = 1 is good. 
%2. make new class files
        % this creates a new folder of files named after the .fcs files 
        % but filled with .mat files of class assingments for each particle. 
        % There is an option of making the class cytogram movies at this
        % step. 
%3. Assign beads to make beadstats and beadplots
        % generate beadplots amd beadstats.mat within bead_calibrated directory 
        % beadstats includes all relevant bead statistics and setup info
%4-5. Apply calibration to add volume to class files 
        % Get conversion for GL1 to SSC and then SSC to volume, save volume
        % values to Class files. 
%6. Generate attune table
        % Attune Table has final collection of filenames, volumes sampled, 
        % counts of different classes, and -- if volumes have been estimated- biovolumes. 
%7. Make a movie
        % option to make a movie from an existing folder of class files 
%8. Match attune table to underway data
        % use rest api or a local spreadsheet to get environmnetal data 
%9. Make standardized volume tables 
        %discretize syn and euk data into volume bins useful for division
        %rate estimation. 
        % Also, this step standardizes environmental variable names so
        % tables can later be merged across cruises. 

%Inputs: 
    %basepath is the full path to the cruise directory with FCS files
    %within it 
% e.g. basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20180720_EN617'
%
%basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20190705_TN368\'
%basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20191005_AR39\'
%basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20201013_EN657\'
%basepath = '\\vdm\PublicData\EN668_Sosik\Sosik-provided_data\Attune\'; 

%basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20180131_EN608\';
%basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20220216_AT46\';
%basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20190512_RB1904\';

function process_wrapper_2021(basepath)

%% choose which steps to do (1 yes, 0 skip) and adjust inputs as necessary
%it assigns these variables whether or not the step is running. 

step1 = 0; %Generate FCSfileinfo

step2 = 0; %make new class files
    dont_overwrite_volumes = 0; %change classes without changing volume estimates
    assign_class_function = 'assign_class_EN706'; 
    filetype2exclude = {'fcb_bead'; 'FCB_bead'; 'bead';  'Cast'; '(lab test)'; 'Dockwater'; 'discrete'; 'Rinses'; "Filter config"; "Cultures"; "cast"; "test"; "08Aug2023"}; % "Dilution";'test'; needed for Step2
    OD2setting = 'GL1'; %where was the OD2 filter on this cruise? 'SSC', 'GL1', or 'None' 
    
    appendonly = 0; %set to 1 if we don't want to change any existing class files.
    
    makemovieasyougo = 0; %option to make things more efficient. 
    framemaker = 'make_movieframe_density';
    stepsize = 1; %controls resolution of movie
    moviechannels = 'late'; %{'BL3-H', 'GL2-H', 'GL1-H', 'GL2-H'}; %parameter numbers for euk X euk Y synX and SynY polygons if framemaker is general
            %typically this is GL1-H for older cruises and GL2-H for new
    

step3 = 0; %Assign beads to make beadstats table and bead plots
    beadfiles2include = {'FCB_bead'};
    beadtype = 'FCB';   %'PT';%'PT';%
    %check OD2setting above in step 2 settings
    
step4 = 0; %set up calibration, only if OD2setting is 'GL1'
    SSCDIM = 'A'; %needed for Step 4 & 5, SSCDIM = 'A' or 'H'
    
step5 = 0; %apply calibration to add volume to class files 
    %Check SSCDIM above anpd OD2setting
    
step6 = 0 ; %Generate attune table

step7 = 0; %Make a movie out of class files after the fact. 
    %Check moviechannels, framemaker and stepsize above. 

step8 = 1; %match underway
    uw_fullname = 'https://nes-lter-data.whoi.edu/api/underway/en706.csv'; %path to find underway environmental data 

step9 = 0; %make standardized volume table and make quality control plot
    

%% Nothing below this section should change between cruises!
% Everything should be adjustable by making a different assign class or
% framemaker etc. 

%% some file structure setup  
fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'bead_calibrated' filesep];
beadfigpath = [outpath filesep 'bead_plots_2021'];
classpath = [outpath 'class' filesep];

if ~exist(outpath, 'dir')
    mkdir(outpath)
end
if ~exist(beadfigpath, 'dir')
    mkdir(beadfigpath)
end

if ~exist(classpath, 'dir')
    mkdir(classpath)
end

%% Save variables for steps being used

save([outpath '\Processing_variables.mat'], 'basepath')
if step2
    step2vars = {dont_overwrite_volumes, assign_class_function, filetype2exclude, OD2setting, appendonly, makemovieasyougo};
    save([outpath '\Processing_variables.mat'], 'step2vars', '-append')
end
if step3
    step3vars = {beadfiles2include, beadtype, OD2setting};
    save([outpath '\Processing_variables.mat'], 'step3vars', '-append')
end
if step4
    step4vars = {SSCDIM}
    save([outpath '\Processing_variables.mat'], 'step4vars', '-append')
end
if step5
    step5vars = {SSCDIM, OD2setting}
    save([outpath '\Processing_variables.mat'], 'step5vars', '-append')

end
if step7 | (step2 & makemovieasyougo) 
    step7vars = {makemovieasyougo, framemaker, moviechannels, stepsize};
    save([outpath '\Processing_variables.mat'], 'step7vars', '-append')
end
if step8 
    step8vars = {uw_fullname};
    save([outpath '\Processing_variables.mat'], 'step8vars', '-append')
end

%% STEP 1
if step1
    [FCSfileinfo] = FCS_DateTimeList(fpath); 
    save([outpath 'FCSfileinfo.mat'], 'FCSfileinfo'); 
else 
    load([outpath filesep 'FCSfileinfo.mat'])
end

%% STEP 2
if step2
    %assign class, save class files, with option to make movies
    step2function(basepath, assign_class_function, filetype2exclude, FCSfileinfo, makemovieasyougo, framemaker, OD2setting, appendonly, moviechannels, dont_overwrite_volumes)
end

%% STEP 3
if step3
   clear FCSfileinfo %need to get back to version that hasn't been cut down for class files in case step 2 was run
   load([outpath filesep 'FCSfileinfo.mat'])

   if strcmp(OD2setting, 'GL1')
       bead_ch_names = {'GL1-A', 'BL3-H', 'GL3-H'}; 
   else
       bead_ch_names = {'SSC-A', 'BL3-H', 'GL2-H'}; 
   end
    %First should be scattering channel, second chlorophyl 

   %process_beads_only(outpath, bead_ch_names, FCSfileinfo, beadfiles2include)
   process_beads_PT_adjust_2(basepath, FCSfileinfo, beadfiles2include, beadtype, OD2setting) %this line works for FCB bead
    
end
%% STEPS 4 - 5
%size calibration, can be redone without reassigning classes if bead processing is adjusted
if step4
    if strcmp(OD2setting, 'GL1')
        get_calibration_stats_linear_2021(outpath, classpath, 50, SSCDIM) %A means ssch_ch_num is ssc-a
    end
end
if step5 
    use_calibration_stats_linear(outpath, classpath, SSCDIM, OD2setting)  
end



%% STEP 6
if step6
    generate_attune_table(classpath, [outpath 'FCSfileinfo.mat'])
    %this function will generate attune table for files with class files
    %only, generally beads are removed
end

%% STEP 7 
% make a movie 
if step7
    attune_lter_moviemaker(fpath, classpath, OD2setting, framemaker, moviechannels, stepsize)
end

%% STEP 8
% match underway data 
if step8
    Attune_uw_match = match_Attune_underway_LTER([outpath 'AttuneTable.mat'],uw_fullname); 
end


%% STEP 9
% make standardized volume tables for division rate estimates and quality
% control %makes products that Bethany used for all of her Thesis work. 
% Also standardizes the names of the environmental data from step 8. 
if step9
    get_cruise_voldists_fromEDItable2(basepath)
    Plot_Voldists_function(basepath, outpath)
end

end


function step2function(basepath, assign_class_function, filetype2exclude, FCSfileinfo, makemovieasyougo, framemaker, OD2setting, appendonly, moviechannels, dont_overwrite_volumes) 

fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'bead_calibrated' filesep];
classpath = [outpath 'class' filesep];

if strcmp(OD2setting, 'GL1')
    ssc_name = 'GL1'; 
else
    ssc_name = 'SSC'; 
end

 %remove filetype2exclude from FCSfileinfo and generate filelist 
        %this part has been updated  for Table structure FCS filinfo 
   for iii = 1:length(filetype2exclude)    
         t = contains(FCSfileinfo.fcslist, filetype2exclude{iii});
         if ~isempty(t)
              FCSfileinfo(t, :) = []; 
         end
    end
    clear t iii 
   
    if appendonly %in this case remove elements from list that alredy have class files 
        fcslist = regexprep(FCSfileinfo.fcslist, '.fcs', '.mat'); %fcslist 
        classlist = dir(classpath); %existing class files 
        classlist = {classlist(:).name};
        FCSfileinfo(ismember(fcslist, classlist), :) = []; %remove them
    end
     filelist = FCSfileinfo.fcslist;
     
    %make sure its in chronological order
    matdate = FCSfileinfo.matdate_start;
    QC_flags = FCSfileinfo.QC_flag; 
    [matdate,sort_ind] = sort(matdate);
    filelist = filelist(sort_ind);
    QC_flags = QC_flags(sort_ind); 

    
    if makemovieasyougo
        v = VideoWriter([classpath 'Attune_cyto_vid.avi']); 
        v.FrameRate = 10; 
        open(v)
    end
    
    %keyboard
    %now go through files of interest, assign classes, and save results
    for count = 1:length(filelist)
         if ~rem(count,10)
            disp([num2str(count) ' of ' num2str(length(filelist))])
         end
        filename = [fpath filelist{count}];
        disp(filename)
        [fcsdat,fcshdr] = fca_readfcs(filename);
        [~,fname] = fileparts(filename);
        class = eval([assign_class_function '( fcsdat, fcshdr, 0, fname, FCSfileinfo.QC_flag(count), FCSfileinfo.matdate_start(count) );']); 
        clear fname
        notes = ['Class 1= Euk, Class 2 = Syn, Class 3 = lowPEeuks, Class 4 = hiPEeuks, Class 5 = Syn_euk_coincident1, Class 0 = junk; Cell volume in cubic microns;',  assign_class_function, string(datetime)];
        
        if dont_overwrite_volumes 
            save([classpath regexprep(filelist{count}, '.fcs', '')], 'class', 'notes', '-append') %I think class files wont have volume yet if we don't have bead statistics to calibrate 
        else
            save([classpath regexprep(filelist{count}, '.fcs', '')], 'class', 'notes') %I think class files wont have volume yet if we don't have bead statistics to calibrate 
        end
        
        if makemovieasyougo
            % assign scattering channels
            ssc_ch = strmatch([ssc_name '-A'], {fcshdr.par.name});
            ssch = strmatch([ssc_name '-H'], {fcshdr.par.name});
            
            % correct for negative SSC-A values
            cf = fitlm(fcsdat(fcsdat(:,ssch)<1000,ssch), fcsdat(fcsdat(:,ssch)<1000, ssc_ch), 'Intercept', false);
            cf = cf.Coefficients.Estimate;
            fcsdat(fcsdat(:,ssc_ch)<0, ssc_ch) = cf*fcsdat(fcsdat(:,ssc_ch)<0, ssch);
            
            % call the function to make the plots and getframe
            eval(['Frame = ', framemaker, '(fcsdat, fcshdr, class, moviechannels, QC_flags(count))']);
            
            % add to movie
            writeVideo(v, Frame);
            
        end
        
    end
    
    if makemovieasyougo
        close (v)
        disp(['Result video saved:' classpath 'Attune_cyto_vid.avi'])
    end
end

