%Taylor Crockfokrd July 2018
%compile long single list of all files in cruise dir to paste into excel
%doc for manual tagging
clear all

ifcb = 'IFCB109';
% start = '23 Mar 2018';
start = '16 Apr 2018';
 stop = '30 Apr 2018';
%stop = now + + (5/24); %account for UTC time else won't get most recent 4 or 5 hours of data

% dirpath = '\\sosiknas1\IFCB_data\NESLTER_transect\data\2019\';
dirpath = '\\sosiknas1\IFCB_data\NESLTER_transect\data\2019\';
dirpath = '\\sosiknas1\IFCB_data\SPIROPA\data\2018\';

dashboards2choose = {'NESLTER_transect'; 'NESLTER_broadscale'; 'SPIROPA'; 'OTZ'; 'WHOI_Dock'; 'EXPORTS'; 'other'}; 
for n=1:length(dashboards2choose),fprintf('%2d    %s\n',n,char(dashboards2choose(n))),end
fprintf('\n')
pick = input('Pick the dashboard to tag by entering the number of the count listed above:');
dashboard2tag = char(dashboards2choose(pick));
dir2save = ['\\sosiknas1\IFCB_data\' dashboard2tag '\to_tag\'];
fprintf('\n')
fprintf(['For the ' dashboard2tag, ' dashboard with cruise dates:  ' start '  to  ' stop]);
fprintf('\n')
cruise = input('What cruise label should be applied? (example: AR24A):   ','s');
    cruise = upper(cruise);
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

excelfile2save = [dir2save dashboard2tag '_' cruise '_to_tag_incomplete'];
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
% files = char(files);
cruise_field = repmat({cruise},length(files),1);
%find which files are empty and should be skipped
ind= find(file_size == 0);
toskip = repmat({''},length(files),1);
toskip(ind) = {'y'};

fprintf('\n')
fprintf('\n')
answer = input('Would you like to apply a tag (such as "underway") TO THE ENTIRE LIST OF FILES besides the cruise tag? (y/n)','s');
if strcmp(lower(answer),'y')
    tag1 = input('What would you like the tag to be? Type exactly as you would like it to appear. (example: underway)','s');
    fprintf('\n');
    fprintf('\n');
    disp(['The 2nd tag will be: ' tag1]);
    fprintf('\n');
    answer = input('Are you positive you want this tag applied to every single file on the list? (y/n)','s');
    if strcmp(lower(answer),'y')
        A = {'Filename' 'cruise' 'skip' 'Tag1'};
        tag1 = repmat({tag1},length(files),1);
        A = [A; files cruise_field toskip tag1];
    elseif strcmp(lower(answer),'n')
        error('Run this file again. You messed up. You have not created a to tag excel doc.');
    end
    fprintf('\n');
    disp(['NOTE: THE ONLY LABELS THAT HAVE BEEN APPLIED TO THIS FILE LIST ARE: cruise field = ' cruise ' and tag = ' char(tag1(1))]);
    fprintf('\n');
    disp('YOU STILL NEED TO REVIEW THIS FILE AND ADD ANY ADDITIONAL TAGS MANUALLY.');
    fprintf('\n');
    disp('THE EXCEL FILE TITLE HAS INCOMPLETE AT THE END FOR A REASON');
elseif strcmp(lower(answer),'n')
    A = {'Filename' 'cruise' 'skip'};
    A = [A; files cruise_field toskip];
    fprintf('\n')
    disp(['NOTE: THE ONLY TAG THAT HAVE BEEN APPLIED TO THIS FILE LIST IS: ' cruise ]);
    fprintf('\n');
    disp('YOU STILL NEED TO REVIEW THIS FILE AND ADD ANY ADDITIONAL TAGS MANUALLY.');
    fprintf('\n');
    disp('THE EXCEL FILE TITLE HAS INCOMPLETE AT THE END FOR A REASON');
else error('YOU DID NOT ANSWER Y OR N DUMMY!!');
end

%make an excel file with the compiled filelist and 
xlswrite(excelfile2save,A)



%sometimes incorrect empty files are made when instr screws up. should skip
%them
display('Empty files that need to be skipped: ')
display(files(ind))

