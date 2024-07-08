function process_wrapper_2024(cruise, steps2do)
% function process_wrapper_2024(cruise, steps2do)
% e.g., 
%   process_wrapper_2024('EN661', [1 3]) %for just steps 1 and 3
% or 
%   process_wrapper_2024('EN661', [1, 3:10]) %for all steps at once
%
% February 2024 (Heidi M. Sosik, WHOI)
% update from process_wrapper_2021 (from Bethany); 
% handle defaults for all steps, with special cases coded for
% early non-standard cruises; change step to a vector for ease of input
%1. Generate FCSfileinfo
%2. Test class assignments
%3. make new class files
%4. Assign beads to make beadstats and beadplots
%5-6. Apply calibration to add volume to class files 
%7. Generate attune table
%8. Make a movie
%9. Match attune table to underway data
%10. Make standardized volume tables 

%This function is designed to be able to process a cruise of attune data 
%from start to finish or with each step optional. 
% This is pretty inefficient, because we have to open each of the
% fcs files in a loop multiple times. But it is totally modular. 

%Steps: 
%1. Generate FCSfileinfo 
        % this is a table of file names and dates, with volume and a
        % quality flag. For this table, QC_flag = 1 is good. 
%2. Test class assignments....
%3. make new class files
        % this creates a new folder of files named after the .fcs files 
        % but filled with .mat files of class assingments for each particle. 
        % There is an option of making the class cytogram movies at this
        % step. 
%4. Assign beads to make beadstats and beadplots
        % generate beadplots amd beadstats.mat within bead_calibrated directory 
        % beadstats includes all relevant bead statistics and setup info
%5-6. Apply calibration to add volume to class files 
        % Get conversion for GL1 to SSC and then SSC to volume, save volume
        % values to Class files. 
%7. Generate attune table
        % Attune Table has final collection of filenames, volumes sampled, 
        % counts of different classes, and -- if volumes have been estimated- biovolumes. 
%8. Make a movie
        % option to make a movie from an existing folder of class files 
%9. Match attune table to underway data
        % use rest api or a local spreadsheet to get environmnetal data 
%10. Make standardized volume tables 
        %discretize syn and euk data into volume bins useful for division
        %rate estimation. 
        % Also, this step standardizes environmental variable names so
        % tables can later be merged across cruises. 


%make step vector
step = zeros(1,10);
step(steps2do) = 1;
basepath_temp =  '\\sosiknas1\Lab_data\Attune\cruise_data\';
temp = dir([basepath_temp '*' cruise]);
if ~isempty(temp)
    p.basepath = [basepath_temp temp.name filesep];
else
    disp('No directory for cruise:')
    disp([basepath_temp '*' cruise])
    keyboard
end
%step(1) = 0; %Generate FCSfileinfo

%step(2) = 0; %make new class files
    p.dont_overwrite_volumes = 0; %change classes without changing volume estimates
    p.assign_class_function = ['assign_class_' cruise]; %'assign_class_AR43'; 
    p.filetype2exclude = {'fcb_bead'; 'FCB_bead'; 'bead';  'Cast'; '(lab test)'; 'Dockwater'; 'discrete'; 'Rinses'; "Filter config"; "Cultures"; "cast"; "test"; "08Aug2023"}; % "Dilution";'test'; needed for Step2
%    OD2setting = 'GL1'; %where was the OD2 filter on this cruise? 'SSC', 'GL1', or 'None'    
    p.appendonly = 0; %set to 1 if we don't want to change any existing class files.
    p.makemovieasyougo = 0; %option to make things more efficient. 
    p.framemaker = 'make_movieframe_density';
    p.stepsize = 1; %controls resolution of movie
    %moviechannels = 'late'; %{'BL3-H', 'GL2-H', 'GL1-H', 'GL2-H'}; %parameter numbers for euk X euk Y synX and SynY polygons if framemaker is general
            %typically this is GL1-H for older cruises and GL2-H for new

%step(3) = 0; %Assign beads to make beadstats table and bead plots
    p.beadfiles2include = {'FCB_bead'};
%    beadtype = 'FCB';   %'PT';%'PT';%
    %check OD2setting above in step 2 settings
    
%step(4) = 0; %set up calibration, only if OD2setting is 'GL1'
    p.SSCDIM = 'A'; %needed for Step 4 & 5, SSCDIM = 'A' or 'H'
    
%step(5) = 0; %apply calibration to add volume to class files 
    %Check SSCDIM above anpd OD2setting
    
%step(6) = 0 ; %Generate attune table

%step(7) = 0; %Make a movie out of class files after the fact. 
    %Check moviechannels, framemaker and stepsize above. 

%step(8) = 0; %match underway
%    uw_fullname = 'https://nes-lter-data.whoi.edu/api/underway/en688.csv'; %path to find underway environmental data 

%step(9) = 1; %make standardized volume table and make quality control plot

%SPECIAL CASES FOR EARLY CRUISES
%where was the OD2 filter on this cruise? 'SSC', 'GL1', or 'None'
if ismember(cruise, {'EN608' 'AR28B' 'AR29' 'EN617' 'AR31A'})
    p.OD2setting = 'SSC';
elseif ismember(cruise, {'EN627' 'AR34B' 'RB1904'})
    p.OD2setting = 'None';
else %DEFAULT
    p.OD2setting = 'GL1';     
end

if ismember(cruise, {'EN608' 'AR28B' 'AR29' 'EN617' 'AR31A' 'EN627' 'AR34B' 'RB1904' 'TN368'})
    p.moviechannels = 'early';
else
    p.moviechannels = 'late'; 
end

if ismember(cruise, {'EN617' 'AR31A' 'EN627'})
    p.beadtype = 'PT';
else
    p.beadtype = 'FCB';
end

switch cruise
    case 'AR29'
        p.uw_fullname = '\\sosiknas1\Lab_data\SPIROPA\20180414_AR29\underway\procAR29_underway_all.mat';
    case 'RB1904' 
        p.uw_fullname = '\\sosiknas1\Lab_data\SPIROPA\20180503_RB1904\compiled_underway\rb1904_uw_compiled.mat';
    case 'TN368' 
        p.uw_fullname = '\\sosiknas1\Lab_data\SPIROPA\20190705_TN368\compiled_underway\tn368_uw_compiled.mat';
    case 'AR43' 
        p.uw_fullname = '\\sosiknas1\Lab_data\OTZ\20200311_AR43\underway\proc\ar43_underway.csv';  
    otherwise %DEFAULT, NES LTER api
        p.uw_fullname = ['https://nes-lter-data.whoi.edu/api/underway/' lower(cruise) '.csv'];
end

%% Nothing below this section should change between cruises!
% Everything should be adjustable by making a different assign class or
% framemaker etc. 

%% some file structure setup  
fpath = [p.basepath filesep 'FCS' filesep];
outpath = [p.basepath filesep 'bead_calibrated' filesep];
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

save([outpath '\Processing_variables.mat'], 'p')
% if step(3)
%     step2vars = {dont_overwrite_volumes, assign_class_function, filetype2exclude, OD2setting, appendonly, makemovieasyougo};
%     save([outpath '\Processing_variables.mat'], 'step2vars', '-append')
% end
% if step(4)
%     step3vars = {beadfiles2include, beadtype, OD2setting};
%     save([outpath '\Processing_variables.mat'], 'step3vars', '-append')
% end
% if step(5)
%     step4vars = {SSCDIM}
%     save([outpath '\Processing_variables.mat'], 'step4vars', '-append')
% end
% if step(6)
%     step5vars = {SSCDIM, OD2setting}
%     save([outpath '\Processing_variables.mat'], 'step5vars', '-append')
% 
% end
% if step(8) | (step(3) & makemovieasyougo) 
%     step7vars = {makemovieasyougo, framemaker, moviechannels, stepsize};
%     save([outpath '\Processing_variables.mat'], 'step7vars', '-append')
% end
% if step(9) 
%     step8vars = {uw_fullname};
%     save([outpath '\Processing_variables.mat'], 'step8vars', '-append')
% end

%% STEP 1
if step(1)
    [FCSfileinfo] = FCS_DateTimeList(fpath); 
    save([outpath 'FCSfileinfo.mat'], 'FCSfileinfo'); 
else 
    load([outpath filesep 'FCSfileinfo.mat'])
end

%% STEP 2
if step(2)
    test_class_assignments(p)
end

%% STEP 3
if step(3)
    %assign class, save class files, with option to make movies
    makeClassFiles(p, FCSfileinfo)
end

%% STEP 4
if step(4)
   clear FCSfileinfo %need to get back to version that hasn't been cut down for class files in case step 2 was run
   load([outpath filesep 'FCSfileinfo.mat'])

   if strcmp(p.OD2setting, 'GL1')
       bead_ch_names = {'GL1-A', 'BL3-H', 'GL3-H'}; 
   else
       bead_ch_names = {'SSC-A', 'BL3-H', 'GL2-H'}; 
   end
    %First should be scattering channel, second chlorophyl 

   %process_beads_only(outpath, bead_ch_names, FCSfileinfo, beadfiles2include)
   process_beads_PT_adjust_2(p.basepath, FCSfileinfo, p.beadfiles2include, p.beadtype, p.OD2setting) %this line works for FCB bead
    
end
%% STEPS 5 - 6
%size calibration, can be redone without reassigning classes if bead processing is adjusted
if step(5)
    if strcmp(p.OD2setting, 'GL1')
        get_calibration_stats_linear_2021(outpath, classpath, 50, p.SSCDIM) %A means ssch_ch_num is ssc-a
    end
end
if step(6) 
    use_calibration_stats_linear(outpath, classpath, p.SSCDIM, p.OD2setting)  
end



%% STEP 7
if step(7)
    generate_attune_table(classpath, [outpath 'FCSfileinfo.mat'])
    %this function will generate attune table for files with class files
    %only, generally beads are removed
end

%% STEP 8 
% make a movie 
if step(8)
    attune_lter_moviemaker(fpath, classpath, p.OD2setting, p.framemaker, p.moviechannels, p.stepsize)
end

%% STEP 9
% match underway data 
if step(9)
    Attune_uw_match = match_Attune_underway_LTER([outpath 'AttuneTable.mat'],p.uw_fullname); 
end


%% STEP 10
% make standardized volume tables for division rate estimates and quality
% control %makes products that Bethany used for all of her Thesis work. 
% Also standardizes the names of the environmental data from step 8. 
if step(10)
    get_cruise_voldists_fromEDItable2(p.basepath)
    Plot_Voldists_function(p.basepath, outpath)
end

end


function makeClassFiles(p, FCSfileinfo) 

fpath = [p.basepath filesep 'FCS' filesep];
outpath = [p.basepath filesep 'bead_calibrated' filesep];
classpath = [outpath 'class' filesep];

if strcmp(p.OD2setting, 'GL1')
    ssc_name = 'GL1'; 
else
    ssc_name = 'SSC'; 
end

 %remove filetype2exclude from FCSfileinfo and generate filelist 
        %this part has been updated  for Table structure FCS filinfo 
   for iii = 1:length(p.filetype2exclude)    
         t = contains(FCSfileinfo.fcslist, p.filetype2exclude{iii});
         if ~isempty(t)
              FCSfileinfo(t, :) = []; 
         end
    end
    clear t iii 
   
    if p.appendonly %in this case remove elements from list that alredy have class files 
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

    
    if p.makemovieasyougo
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
        class = eval([p.assign_class_function '( fcsdat, fcshdr, 0, fname, FCSfileinfo.QC_flag(count), FCSfileinfo.matdate_start(count) );']); 
        clear fname
        notes = ['Class 1= Euk, Class 2 = Syn, Class 3 = lowPEeuks, Class 4 = hiPEeuks, Class 5 = Syn_euk_coincident1, Class 0 = junk; Cell volume in cubic microns;',  p.assign_class_function, string(datetime)];
        
        if p.dont_overwrite_volumes 
            save([classpath regexprep(filelist{count}, '.fcs', '')], 'class', 'notes', '-append') %I think class files wont have volume yet if we don't have bead statistics to calibrate 
        else
            save([classpath regexprep(filelist{count}, '.fcs', '')], 'class', 'notes') %I think class files wont have volume yet if we don't have bead statistics to calibrate 
        end
        
        if p.makemovieasyougo
            % assign scattering channels
            ssc_ch = strmatch([ssc_name '-A'], {fcshdr.par.name});
            ssch = strmatch([ssc_name '-H'], {fcshdr.par.name});
            
            % correct for negative SSC-A values
            cf = fitlm(fcsdat(fcsdat(:,ssch)<1000,ssch), fcsdat(fcsdat(:,ssch)<1000, ssc_ch), 'Intercept', false);
            cf = cf.Coefficients.Estimate;
            fcsdat(fcsdat(:,ssc_ch)<0, ssc_ch) = cf*fcsdat(fcsdat(:,ssc_ch)<0, ssch);
            
            % call the function to make the plots and getframe
            eval(['Frame = ', p.framemaker, '(fcsdat, fcshdr, class, moviechannels, QC_flags(count))']);
            
            % add to movie
            writeVideo(v, Frame);
            
        end
        
    end
    
    if p.makemovieasyougo
        close (v)
        disp(['Result video saved:' classpath 'Attune_cyto_vid.avi'])
    end
end

