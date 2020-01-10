%Taylor Crockfokrd July 2018
%compile long single list of all files in cruise dir to paste into excel
%doc for manual tagging
clear all

% start = '23 Mar 2018';
start = '1 Feb 2019';
 stop = '6 Feb 2019';
%stop = now + + (5/24); %account for UTC time else won't get most recent 4 or 5 hours of data

% dirpath = ['\\sosiknas1\Backup\SPIROPA\20180414_AR29\IFCB_data\data_underway\' start(end-3:end) '\'];
% dirpath = ['\\sosiknas1\Backup\LTER\20180404_AR28\IFCB_data\data_discrete\'];
% dirpath = '\\sosiknas1\IFCB_data\NESLTER_transect\data\2019\';
dirpath = '\\sosiknas1\IFCB_data\NESLTER_transect\data\2019\';


dashboards2choose = {'NESLTER_transect'; 'NESLTER_broadscale'; 'OTZ'; 'WHOI_Dock'; 'EXPORTS'; 'other'}; 
for n=1:length(dashboards2choose),fprintf('%2d    %s\n',n,char(dashboards2choose(n))),end
fprintf('\n')
pick = input('Pick the dashboard to tag by entering the number of the count listed above:');
dashboard2tag = char(dashboards2choose(pick));
dir2save = ['\\sosiknas1\IFCB_data\' dashboard2tag '\to_tag\'];
fprintf('\n')
cruise = input('What cruise tag should be applied? (example: cruise_AR24A):   ','s');
    fprintf('\n')
    fprintf('\n')
    disp(['The cruise tag will be: ' cruise])
    fprintf('\n')
    answer = input('Is this the correct tag? (y/n)','s');
    if strcmp(lower(answer),'n')
        error('Run this file again. You messed up. You have not created a to tag excel doc.')
    elseif ~strcmp(lower(answer),'y')
        error('YOU DID NOT ANSWER Y OR N DUMMY!!');
    end

excelfile2save = [dir2save dashboard2tag '_' cruise(8:end) '_to_tag_incomplete'];
fprintf('\n')
disp('Tag file being created is named: ');
disp(excelfile2save);


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

files = {};
file_size = [];
for j=1:length(dirlist)
    roidir      = [dirpath cell2mat(dirlist(j))];
    roifiles    = dir([roidir '\*IFCB127.roi']);
    temp={roifiles.name}';
    temp = char(temp);
    temp=temp(:,1:end-4);
    files = [files; cellstr(temp)];
    file_size = [file_size; [roifiles.bytes]'];
    clear temp
end
% files = char(files);
cruise_tag = repmat({cruise},length(files),1);
going2dashboard = repmat({dashboard2tag},length(files),1);

fprintf('\n')
fprintf('\n')
answer = input('Would you like to apply a 2nd tag (such as "underway") TO THE ENTIRE LIST OF FILES besides the cruise tag? (y/n)','s');
if strcmp(lower(answer),'y')
    tag2 = input('What would you like the 2nd tag to be? Type exactly as you would like it to appear. (example: underway)','s');
    fprintf('\n');
    fprintf('\n');
    disp(['The 2nd tag will be: ' tag2]);
    fprintf('\n');
    answer = input('Are you positive you want this tag applied to every single file on the list? (y/n)','s');
    if strcmp(lower(answer),'y')
        A = {'going_to_dashboard' 'Filename' 'Tag1' 'Tag2'};
        tag2 = repmat({tag2},length(files),1);
        A = [A; going2dashboard files cruise_tag tag2];
    elseif strcmp(lower(answer),'n')
        error('Run this file again. You messed up. You have not created a to tag excel doc.');
    end
    fprintf('\n');
    disp(['NOTE: THE ONLY TAGS THAT HAVE BEEN APPLIED TO THIS FILE LIST ARE: ' cruise ' and ' char(tag2(1))]);
    fprintf('\n');
    disp('YOU STILL NEED TO REVIEW THIS FILE AND ADD ANY ADDITIONAL TAGS MANUALLY.');
    fprintf('\n');
    disp('THE EXCEL FILE TITLE HAS INCOMPLETE AT THE END FOR A REASON');
elseif strcmp(lower(answer),'n')
    A = {'going_to_dashboard' 'Filename' 'Tag1'};
    A = [A; going2dashboard files cruise_tag];
    fprintf('\n')
    disp(['NOTE: THE ONLY TAG THAT HAVE BEEN APPLIED TO THIS FILE LIST IS: ' cruise ])
    fprintf('\n')
    disp('YOU STILL NEED TO REVIEW THIS FILE AND ADD ANY ADDITIONAL TAGS MANUALLY.')
    fprintf('\n')
    disp('THE EXCEL FILE TITLE HAS INCOMPLETE AT THE END FOR A REASON')
else error('YOU DID NOT ANSWER Y OR N DUMMY!!');
end

%make an excel file with the compiled filelist and 
xlswrite(excelfile2save,A)



%sometimes incorrect empty files are made when instr screws up. should skip
%them
ind= find(file_size == 0);
display('Empty files that need to be skipped: ')
display(files(ind))

