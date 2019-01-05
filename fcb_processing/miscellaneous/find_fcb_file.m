function [time_window,timepath,grouped_window,groupedpath,merged_window,mergedpath] = find_fcb_file(time1,time2)
%given anytime window, find the corresponding files in either time, merged
%or grouped steps....
%should be able to handle either matlab time entry or string entry

%extract year:
[yy1,~,~,~,~,~]=datevec(time1);
[yy2,~,~,~,~,~]=datevec(time2);

if yy1==yy2
    %find the correct year label:
    switch yy1
        case 2003
            yearlabel='May';
        case 2004
            yearlabel='Apr';
        case 2005
            yearlabel='Apr';
        case 2006
            yearlabel='May';
        case 2007
            yearlabel='Mar';
        otherwise
            yearlabel='Jan';
    end
else
    disp('hmmm - have not coded the piece to deal with timespan more than 1 year...')
    keyboard
end

%and the correct partial path:
if ~isempty(strfind(computer,'WIN'))
    rootpath='\\sosiknas1\lab_data\MVCO\FCB\';
else
    rootpath='/Volumes/Lab_data/MVCO/FCB/';
end

%generate time list:
timepath=fullfile(rootpath,['MVCO_' yearlabel num2str(yy1)],'data/processed/time/');
timefiles=dir([timepath '*.mat']);
time_file_record=cell(length(timefiles),3);

for j=1:length(timefiles)
    temp=load([timepath timefiles(j).name],'-regexp','\d'); %only load in the variable that contains the time, loaded as a structure
    tempname=fieldnames(temp); tempname=char(tempname); %extract structure field name
    eval(['timemat=temp.' tempname ';'])
    time_file_record{j,1}=timefiles(j).name;
    time_file_record{j,2}=min(timemat(:,2)); %second col is start time
    time_file_record{j,3}=max(timemat(:,3)); %third col is start time
    clearvars timemat temp tempname
end

%make sure in correct temporal order:
[~,is]=sort(cell2mat(time_file_record(:,2)));
time_file_record=time_file_record(is,:);

%generate corresponding lists for merged and grouped files:
temp=regexp(time_file_record(:,1),'(?<part1>[A-Z0-9\_]*)time(?<part2>[A-Z0-9\_]*.mat)','names');

mergedpath=fullfile(rootpath,['MVCO_' yearlabel num2str(yy1)],'data/processed/grouped/merged/');
mergedfiles=cell(length(timefiles),1);
for j=1:length(timefiles)
    mergedfiles{j}=[temp{j}.part1 'merged' temp{j}.part2];
end

groupedpath=fullfile(rootpath,['MVCO_' yearlabel num2str(yy1)],'data/processed/grouped/');
groupedfiles=cell(length(timefiles),1);
for j=1:length(timefiles)
    groupedfiles{j}=[temp{j}.part1 temp{j}.part2];
end




%okay, so now have the list of time covered in each file - note
%extract desired file names:
jj1=find(cell2mat(time_file_record(:,2)) <= datenum(time1));
jj2=find(cell2mat(time_file_record(:,3)) >= datenum(time2));

if isempty(jj1)
    jj1=1;
elseif isempty(jj2) %rare cases were at end of the record
    jj2=size(time_file_record,1); %chose ending index
end

si=jj1(end); %starting index
ei=jj2(1); %ending index

time_window=time_file_record(si:ei,1);
grouped_window=groupedfiles(si:ei,1);
merged_window=mergedfiles(si:ei,1);






end