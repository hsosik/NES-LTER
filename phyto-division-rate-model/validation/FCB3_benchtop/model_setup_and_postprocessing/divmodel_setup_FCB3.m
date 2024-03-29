
clear
close all

%FCB3 benchtop validation experiments:

valtype='benchtop'; %'dilseries

if ~isempty(strfind(computer,'WIN'))
    rootpath='\\sosiknas1\lab_data\MVCO\FCB\';
else
%     rootpath='/Volumes/Lab_data/MVCO/FCB/';
    datapath='/Users/kristenhunter-cevera/Documents/MATLAB/FCB3_benchtop/SynLab2006/model/';
    codepath='/Users/kristenhunter-cevera/NES-LTER/phyto-division-rate-model/';
end

addpath(fullfile(codepath, '../fcb_processing/miscellaneous/')) %contains helpful scripts :)

do_Solar = 0;
solarplotflag=0;

do_setupdays = 0;
folder_tag='dmn_14par';
do_setupdays_movie = 0;

do_model = 1;
do_modelfit_movie = 0;
redo_model=0;

% 
% if do_Solar
%     solarsavepath=fullfile(datapath, '/model/');
%     if ~exist(solarsavepath, 'dir'), mkdir(solarsavepath), end
%     getSolar
% end
% 
% if do_setupdays
%     mergedpath0 = fullfile([datapath,'data/processed/grouped/merged/']); %the '0' designation is important - mergedpath is a save variable name in the files that will be downloaded....
%     groupedpath = fullfile(datapath,'data/processed/grouped/');
%     beadpath=fullfile(datapath, 'data/processed/beads/');
%     modelinputpath = fullfile(datapath,[ 'model/input_beadmean_' folder_tag '/']);
%     if ~exist(modelinputpath, 'dir'), mkdir(modelinputpath), end %where daily input will go....
%     plotflag=1;
%     setup_days_all
% end
% 
% % make movies of input days:
% if do_setupdays_movie
%     modelinputpath = fullfile(datapath,[ 'model/input_beadmean_' folder_tag '/']);
%     MVCO_setup_days_QC
% end

%run the model:
if do_model
    addpath(fullfile(codepath,'model_src/')) %assumes you are in the setup_and_postprocessing folder - should fix this!
    addpath(fullfile(codepath,'validation/FCB3_benchtop/model_setup_and_postprocessing/'))
    pathname=fullfile(datapath, 'input/');
    savepath=fullfile(codepath,['validation/FCB3_benchtop/output_' folder_tag '/']); %Change path to sosiknas here!
    %savepath=fullfile(datapath,['output_' folder_tag '\']); %Change path to sosiknas here!
    if ~exist(savepath,'dir'), mkdir(savepath), end % make directory if doesn't exist yet
    call_to_opt_mvco
end

%make model result figures:
% if do_modelfit_movie
%     modelres_path=fullfile(datapath,['\model\output_' folder_tag '\']);  %path to model results
%     setupdays_path=fullfile(datapath,['model/input_beadmean_' folder_tag '/']); %path to setup_days input
%     addpath(fullfile(rootpath,'Syn_divrate_model/model_code/')); %path to access code to project model day forward
%     addpath(fullfile(rootpath,'Syn_divrate_model/subplot_tight/')); %path for making 'tighter' subplots :)
%     
%     modelfit_movies
% end
% 
% %rerun specified additional days
% if redo_model
%     disp(['Redoing days in ' num2str(year2do)])
%     input_path=fullfile(datapath,'/model/input_beadmean_July2016/'); %path to setup_days input
%     savepath=fullfile(datapath,'/model/output_July2016/');  %path to model results
%     if ~exist(savepath,'dir'), mkdir(savepath), end %make directory if doesn't exist yet
%     [days2redo]=exclude_modeldata(year2do); % get the list of days to redo
%     
%     if ~isempty(days2redo)
%         days2redo=str2num(cell2mat(days2redo(:,1))); %turn cell array into day list
%         call_to_opt_redo120
%     end
% end

