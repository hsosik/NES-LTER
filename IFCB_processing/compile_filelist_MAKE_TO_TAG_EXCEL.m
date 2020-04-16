%Taylor Crockfokrd July 2018
%compile long single list of all files in cruise dir to paste into excel
%doc for manual tagging
clear all

ifcb = 'IFCB011';
% start = '23 Mar 2018';
start = '14 Aug 2018';
 stop = '29 Aug 2018';
%stop = now + + (5/24); %account for UTC time else won't get most recent 4 or 5 hours of data

% dirpath = '\\sosiknas1\IFCB_data\NESLTER_transect\data\2019\';
dirpath = '\\sosiknas1\IFCB_data\NESLTER_transect\data\2019\';
% dirpath = '\\sosiknas1\IFCB_data\EXPORTS\data\2028\';

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

sampletypes2choose = {'underway';'cast';'other'};
fprintf('\n')
fprintf('\n')
disp('Now we will apply a SAMPLE TYPE to the ENTIRE LIST OF FILES.');
fprintf('\n')
fprintf('\n')
for n=1:length(sampletypes2choose),fprintf('%2d    %s\n',n,char(sampletypes2choose(n))),end
fprintf('\n')
pick = input('Pick the sample type by entering the number of the count listed above:');

if pick == length(sampletypes2choose)
    sampletype2use = input('What would you like the sample type to be? Type exactly as you would like it to appear. (example: underway)','s');
    fprintf('\n');
    fprintf('\n');
    disp(['The sample type will be: ' sampletype2use]);
    fprintf('\n');
    answer = input('Are you positive you want this tag applied to every single file on the list? (y/n)','s');
    if strcmp(lower(answer),'y')
        A = {'Filename' 'cruise' 'skip' 'sample_type'};
        sampletype2use = repmat({sampletype2use},length(files),1);
        A = [A; files cruise_field toskip sampletype2use];
    elseif strcmp(lower(answer),'n')
        error('Run this file again. You messed up. You have NOT created a to tag excel doc.');
    end
    fprintf('\n');
    disp(['NOTE: THE ONLY LABELS THAT HAVE BEEN APPLIED TO THIS FILE LIST ARE: cruise field = ' cruise ' and sample type = ' char(sampletype2use(1))]);
    fprintf('\n');
    disp('YOU STILL NEED TO REVIEW THIS FILE AND ADD ANY ADDITIONAL TAGS MANUALLY.');
    fprintf('\n');
    disp('THE EXCEL FILE TITLE HAS INCOMPLETE AT THE END FOR A REASON');
elseif pick > length(sampletypes2choose)
     error('THAT NUMBER IS NOT AN OPTION DUMMY!!');
else
    sampletype2use = char(sampletypes2choose(pick));
    A = {'Filename' 'cruise' 'skip' 'sample_type'};
    sampletype2use = repmat({sampletype2use},length(files),1);
    A = [A; files cruise_field toskip sampletype2use];
    fprintf('\n')
    disp(['NOTE: THE ONLY TAG THAT HAVE BEEN APPLIED TO THIS FILE LIST IS: ' cruise ]);
    fprintf('\n');
    disp('YOU STILL NEED TO REVIEW THIS FILE AND ADD ANY ADDITIONAL TAGS MANUALLY.');
    fprintf('\n');
    disp('THE EXCEL FILE TITLE HAS INCOMPLETE AT THE END FOR A REASON');
end

%make an excel file with the compiled filelist and 
xlswrite(excelfile2save,A)



%sometimes incorrect empty files are made when instr screws up. should skip
%them
display('Empty files that need to be skipped: ')
display(files(ind))

