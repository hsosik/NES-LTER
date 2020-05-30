%Taylor Crockfokrd July 2018
%compile long single list of all files in cruise dir to paste into excel
%doc for manual tagging
clear all

%% VARIABLES TO EDIT
ifcb = 'IFCB101';
% start = '23 Mar 2018';
start = '11 Feb 2017';
 stop = '23 Feb 2017';
%stop = now + + (5/24); %account for UTC time else won't get most recent 4 or 5 hours of data

%% Choose dashboard to use
dashboards2choose = {'NESLTER_transect'; 'NESLTER_broadscale'; 'SPIROPA'; 'OTZ'; 'WHOI_Dock'; 'EXPORTS'; 'other'}; 
for n=1:length(dashboards2choose),fprintf('%2d    %s\n',n,char(dashboards2choose(n))),end
fprintf('\n')
pick = input('Pick the dashboard to tag by entering the number of the count listed above:');
dashboard2tag = char(dashboards2choose(pick));

%% If not a dashboard pre-listed, carefully manually enter dashboard
if regexp(dashboard2tag, 'other')
fprintf('\n')
fprintf('You chose other.')
fprintf('\n')
dashboard2tag = input('Enter name of dashboard to tag exactly as dashboard name appears in sosiknas IFCB_data dir: ','s'); 
end
fprintf('\n')
fprintf(['Dashboard to use is: ' dashboard2tag])
fprintf('\n')

%%
%%If a NESLTER_broadscale dashboard, give opportunity to make additional
%%tag.
if regexp(dashboard2tag, 'NESLTER_broadscale')
fprintf('\n')
fprintf('You chose NESLTER_broadscale. Often times, a cruise type such as ECOMON is associated with each cruise.')
fprintf('\n')
answer = input('Would you like to add a tag for the cruise type? (Answer y or n using.)','s');
if regexp(lower(answer),'y')
broadscale_type = {'ecomon'; 'cyst'; 'amaps'; 'other'}; 
for n=1:length(broadscale_type),fprintf('%2d    %s\n',n,char(broadscale_type(n))),end
fprintf('\n')
pick = input('Pick the dashboard to tag by entering the number of the count listed above: ');
tag1 = char(broadscale_type(pick));
if regexp(tag1, 'other')
fprintf('\n')
fprintf('You chose other.')
fprintf('\n')
tag1 = input('Enter name of the tag exactly as you would like it to appear: ','s'); 
end
fprintf('\n')
fprintf(['Tag to apply to this broadscale TO_TAG file is: ' tag1])
fprintf('\n')
end
end

%% Directory to get data from and where to output excel file created here
dirpath = ['\\sosiknas1\IFCB_data\' dashboard2tag '\data\' start(end-3:end) '\'];
dir2save = ['\\sosiknas1\IFCB_data\' dashboard2tag '\to_tag\'];
fprintf('\n')
fprintf(['For the ' dashboard2tag, ' dashboard with cruise dates:  ' start '  to  ' stop]);
fprintf('\n')

%% Enter cruise to use. 
cruise = input('What cruise label should be applied? (example: AR24A):   ','s');
    cruise = upper(cruise);
    fprintf('\n')
    fprintf('\n')
    disp(['The cruise tag will be: ' cruise])
    fprintf('\n')
    answer = input('Is this the correct tag? (y/n)','s');
    if regexp(lower(answer),'n')
        error('Run this file again. You messed up. You have not created a to tag excel doc.')
    elseif ~regexp(lower(answer),'y')
        error('YOU DID NOT ANSWER Y OR N DUMMY!!');
    end
%
%% Show what output file is named and where it will be saved
excelfile2save = [dir2save dashboard2tag '_' cruise '_to_tag_incomplete'];
fprintf('\n')
disp('Tag file being created is named: ');
disp(excelfile2save);

%% Get directories pertaining to given start and stop dates
d = dir([dirpath 'D*']);
isub = [d(:).isdir]; %# returns logical vector
d = d(isub); % only folders

startdate =  datenum(start);
stopdate =  datenum(stop);
dirname = struct2cell(d);
dirname = dirname(1,:); %names of dirs as cell array

tempdir = char(dirname');
tempdir = tempdir(:,2:end);
tempdir = datenum(tempdir,'yyyymmdd');
ind = find(tempdir >= startdate & tempdir <= stopdate);
dirlist = dirname(ind)';
clear ind temp* d isub

%% Compile all filenames. Only of particular IFCB entered at top of file
files = {};
file_size = [];
for j=1:length(dirlist)
    roidir      = [dirpath cell2mat(dirlist(j))];
    roifiles    = dir([roidir '\*' ifcb '.roi']);
    temp={roifiles.name}';
    if ~isempty(temp) 
        temp = char(temp);
        temp=temp(:,1:end-4);
        files = [files; cellstr(temp)];
        file_size = [file_size; [roifiles.bytes]'];
    end
    clear temp
end

%% Make cruise field length of file list for when outputting csv
cruise_field = repmat({cruise},length(files),1);

%% Find which files are empty and should be skipped
ind= find(file_size == 0);
toskip = repmat({''},length(files),1);
toskip(ind) = {'y'};

%% Label what sample_type is for entire list. Might need spot edits in outputted fiel later
sampletypes2choose = {'underway';'cast';'other'};
fprintf('\n'); fprintf('\n');
disp('Now we will apply a SAMPLE TYPE to the ENTIRE LIST OF FILES.');
fprintf('\n'); fprintf('\n');
for n=1:length(sampletypes2choose),fprintf('%2d    %s\n',n,char(sampletypes2choose(n))),end
fprintf('\n')
pick = input('Pick the sample type by entering the number of the count listed above:');
sampletype2use = char(sampletypes2choose(pick));
clear pick n sampletypes2choose

if regexp(sampletype2use,'other')
    sampletype2use = input('You chose other. What would you like the sample type to be? Type exactly as you would like it to appear. (example: underway)','s');
    fprintf('\n'); fprintf('\n');
    disp(['The sample type will be: ' sampletype2use]);
end
fprintf('\n');

%% Make final matrix of file list with cruise, sample_type, skips, and tag1 if exist. Will be final output to excel
if exist('tag1','var')
        A = {'Filename' 'cruise' 'skip' 'sample_type' 'tag1'};
        sampletype2use = repmat({sampletype2use},length(files),1);
        tag1 = repmat({tag1},length(files),1);
        A = [A; files cruise_field toskip sampletype2use tag1];
else
        A = {'Filename' 'cruise' 'skip' 'sample_type'};
        sampletype2use = repmat({sampletype2use},length(files),1);
        A = [A; files cruise_field toskip sampletype2use];
end
    fprintf('\n');
    disp(['NOTE: THE ONLY LABELS THAT HAVE BEEN APPLIED TO THIS FILE LIST ARE: cruise field = ' cruise ' and sample type = ' char(sampletype2use(1))]);
    fprintf('\n');
    disp('YOU STILL NEED TO REVIEW THIS FILE AND ADD ANY ADDITIONAL TAGS MANUALLY.');
    fprintf('\n');
    disp('THE EXCEL FILE TITLE HAS INCOMPLETE AT THE END FOR A REASON');

%% Make an excel file with the compiled filelist and 
xlswrite(excelfile2save,A)



%% Give list of sometimes incorrect empty files are made when instr screws up. should skip
disp('Empty files that need to be skipped: ')
disp(files(ind))

