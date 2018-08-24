%script to import all the Tioga processed CTD data (.cnv files) from the
%processed_CTD_data folder:

%if processing on Kristen's machine:
addpath /Users/kristenhunter-cevera/MVCO_light_at_depth/seawater_ver3_2/
pathname='\\sosiknas1\Lab_data\MVCO\processed_CTD_casts_20171206\';
% pathname='/Volumes/Lab_data/MVCO/processed_CTD_casts/';
%pathname = '/Volumes/TaylorF/from_Samwise/data/MVCO/SurveyCruises/'; %once maddie is mounted
%pathname = '/Volumes/J_data/MVCO/SurveyCruises/';

filelist=dir(pathname);
filenames=extractfield(filelist,'name'); %cell array of folder names
temp=regexp(filenames,'\.cnv'); %find only the .cnv files
ii=find(cellfun('isempty',temp)==0);
filenames=filenames(ii)';
temp=regexp(filenames,'\w*docktest\w*'); %find "docktest" files that are junk
ii=find(cellfun('isempty',maybe)==1); %remove docktest files from filenames list
filenames=filenames(ii);

%find the Tioga folders:
%ii=find(cellfun('isempty',strfind(folder_names','Tioga'))==0);
load(fullfile(pathname,'list_and_location_of_raw_ctd_files.mat'))

CTD=struct('cast_name',{},'file_location',{},'lat',{},'lon',{},'UTC',{},'upload_time',{},'data_hdr',{},'data',{},'notes',{});

for q=1:length(filenames)

    disp(filenames{q})

    CTD(q).cast_name=filenames{q};

    jj=find(cellfun('isempty',regexp(datafiles,filenames{q}(1:end-4)))==0);
    CTD(q).file_location=datafiles{jj};

    if regexp(filenames{q},'\w*handcast\w*')
        lat = 41.345; lon =

    [lat,lon, UTC_time,upload_time,header,data]=import_cnvfile([pathname filenames{q}]);

    if ~isempty(lat) || ~isempty(lon)
        CTD(q).lat=[str2num(lat.deg)+str2num(lat.min)/60];
        CTD(q).lon=-[str2num(lon.deg)+str2num(lon.min)/60];

    elseif ~isnan(upload_time) %if is empty, we can do a double check based on time stamp

        [mm, im]=min(abs(MVCO_event_time-upload_time));
        disp(abs(MVCO_event_time(im)-upload_time))

        if abs(MVCO_event_time(im)-upload_time) > 0.5 %greater than half a day...
            disp('Is this the correct time match for the file to the MVCO log?')
            keyboard
        end

        %Maybe just print out whether or not a multi-day event?
        tt=find(floor(upload_time)==floor(MVCO_event_time));
        if length(tt)==1
            disp('Part of multi-day event at MVCO')
        end

        %if found a match, use that lat and lon:
        CTD(q).lat=MVCO_bottle_events_unq{im,5};
        CTD(q).lon=MVCO_bottle_events_unq{im,6};
        CTD(q).notes=['lat/lon matched from log; ' MVCO_bottle_events_unq{im,1}];

    end

    CTD(q).UTC=UTC_time;
    CTD(q).upload_time=upload_time;

    %calculate and potential density:
    if any(~cellfun('isempty',data)) %if non-empty cell exists:

        %make sure you've got the right headers:
        s=find(cellfun('isempty',regexp(header,'Salinity'))==0);
        t=find(cellfun('isempty',regexp(header,'Temperature'))==0);
        p=find(cellfun('isempty',regexp(header,'Pressure'))==0);

        pdens=sw_pden(data{s},data{t},data{p},0); %0 refers to reference pressure
        data{end+1}=pdens;
        header{end+1}='Potential Density';
    end

    CTD(q).data_hdr=header;
    CTD(q).data=data;

    if length(header) ~= length(data)
        keyboard
    end
end

%% and to save...

save(fullfile(pathname,'CTD_30Mar2018'),'CTD')

<<<<<<< HEAD
%% check to see what is missing:
=======
ii=find(cellfun('isempty',{CTD(:).lat}')==1);
{CTD(ii).cast_name}';
%Okay, only about 8 or so actual cruises where we are missing this data...
>>>>>>> 093bc179247bcef8b709d5fa0291403ba4fbf747

ii=find(cellfun(@(x) isnan(x),{CTD(:).lat}')==1);
{CTD(ii).cast_name}'
%should be empty!

%% find only trips to tower or node:

templat=cell2mat({CTD(:).lat}');
templon=cell2mat({CTD(:).lon}');

box_tower=[-70.58 -70.53 41.315 41.33];
box_node=[-70.58 -70.53 41.33 41.345];

mvco_ind=find(templat > 41.3 & templat < 41.35 & templon < -70.53 & templon > -70.60 & cellfun('isempty',regexp({CTD(:).cast_name}','(deck)|(test)'))==1);
%this should be about ~99 casts...

%% and if curious, see where all the casts come from...

plot(templon,templat,'.')
patch([-70.53 -70.53 -70.60 -70.60],[41.3 41.35 41.35 41.3],'k','facecolor','none')
