clear
close all

for year2do = 2006:2018

    disp(num2str(year2do))
    % addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/ %has cytosub_SSC2vol.m 
    
    if ~isempty(strfind(computer,'WIN'))
        envpath='\\sosiknas1\lab_data\mvco\EnvironmentalData\';
        rootpath='\\sosiknas1\lab_data\MVCO\FCB\';
        model_path=fullfile('~/NES-LTER/phyto-division-rate-model/');
        addpath(fullfile('~/NES-LTER/fcb_processing/miscellaneous/')); %contains helpful scripts :)
    elseif ~isempty(strfind(computer,'GLNXA64'))
        envpath='/mnt/Lab_data/MVCO/EnvironmentalData/';
        rootpath='/mnt/Lab_data/MVCO/FCB/';
        model_path=fullfile('~/Documents/NES-LTER/phyto-division-rate-model/');
        addpath(fullfile('~/Documents/NES-LTER/fcb_processing/miscellaneous/')); 
    else
        envpath='/Volumes/Lab_data/MVCO/EnvironmentalData/';
        rootpath='/Volumes/Lab_data/MVCO/FCB/';
        model_path=fullfile('~/NES-LTER/phyto-division-rate-model/');
        addpath(fullfile('~/NES-LTER/fcb_processing/miscellaneous/'));
    end
    
    do_Solar = 0;
    buoy_flag=0;
    solarplotflag=0;
    
    do_setupdays = 0;
    folder_tag='June2019';
    do_setupdays_movie = 0;
    
    do_model = 1;
    components=1; %1 or 2...
    do_modelfit_movie = 0;
    redo_model=0;
    
    
    switch year2do
        case 2003
            datapath = fullfile(rootpath,'MVCO_May2003/');
        case 2004
            datapath = fullfile(rootpath,'MVCO_Apr2004/');
        case 2005
            datapath = fullfile(rootpath,'MVCO_Apr2005/');
        case 2006
            datapath = fullfile(rootpath,'MVCO_May2006/');
        case 2007
            datapath = fullfile(rootpath,'MVCO_Mar2007/');
        otherwise
            datapath =fullfile(rootpath,['MVCO_Jan' num2str(year2do) '/']);
    end

    if do_Solar
        solarsavepath=fullfile(datapath, '/model/');
        if ~exist(solarsavepath, 'dir'), mkdir(solarsavepath), end
        getSolar     
    end

    if do_setupdays
        mergedpath0 = fullfile([datapath,'data/processed/grouped/merged/']); %the '0' designation is important - mergedpath is a save variable name in the files that will be downloaded....
        groupedpath = fullfile(datapath,'data/processed/grouped/');
        beadpath=fullfile(datapath, 'data/processed/beads/');
        modelinputpath = fullfile(datapath,[ 'model/input_beadmean_' folder_tag '/']);
        if ~exist(modelinputpath, 'dir'), mkdir(modelinputpath), end %where daily input will go....
        plotflag=1;
        setup_days_all
    end

    % make movies of input days:
    if do_setupdays_movie
        modelinputpath = fullfile(datapath,[ 'model/input_beadmean_' folder_tag '/']);
        MVCO_setup_days_QC
    end

    %run the model:
    if do_model
        
        disp(num2str(year2do))
        pathname=fullfile(datapath, ['model/input_beadmean_' folder_tag '/']);
        
        switch components
            case 2
                savepath=fullfile(datapath,['model/output_' folder_tag '/']); %Change path to sosiknas here!
                if ~exist(savepath,'dir'), mkdir(savepath), end % make directory if doesn't exist yet
                
                addpath(model_path)
                addpath(fullfile(model_path,'model_src/two_subpopulation'))
                
                call_to_opt_mvco
                
            case 1
                savepath=fullfile(datapath,['model/onecomp_output_' folder_tag '/']); %Change path to sosiknas here!
                if ~exist(savepath,'dir'), mkdir(savepath), end % make directory if doesn't exist yet
                
                addpath(model_path)
                addpath(fullfile(model_path,'model_src/one_population'))
                
                call_to_opt_mvco_onecomp                
        end
        
    end

    %make model result figures:
    if do_modelfit_movie
        modelres_path=fullfile(datapath,'model',['output_' folder_tag],'/');  %path to model results
        setupdays_path=fullfile(datapath,'model',['input_beadmean_' folder_tag],'/'); %path to setup_days input
        addpath(fullfile('~/NES-LTER/phyto-division-rate-model/setup_and_postprocessing/','plotting_tools/subplot_tight/')); %path for making 'tighter' subplots :)

        modelfit_movies
    end

    %rerun specified additional days
    if redo_model
        disp(['Redoing days in ' num2str(year2do)])
        input_path=fullfile(datapath,'/model/input_beadmean_July2016/'); %path to setup_days input
        savepath=fullfile(datapath,'/model/output_July2016/');  %path to model results
        if ~exist(savepath,'dir'), mkdir(savepath), end %make directory if doesn't exist yet
        [days2redo]=exclude_modeldata(year2do); % get the list of days to redo

        if ~isempty(days2redo)
            days2redo=str2num(cell2mat(days2redo(:,1))); %turn cell array into day list
            call_to_opt_redo120
        end
    end

    %clearvars -except year2do
    % close all
end