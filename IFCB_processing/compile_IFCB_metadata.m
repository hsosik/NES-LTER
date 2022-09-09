function [] = compile_IFCB_metadata(ToTag_xlsFile)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[~,f] = fileparts(ToTag_xlsFile);
cruise = strsplit(f, '_');
cruise = cruise{3};
apibase = 'https://nes-lter-data.whoi.edu/api/';
myreadtable = @(filename)readtable(filename,'Delimiter','comma');
options = weboptions('ContentReader',myreadtable);
%options = weboptions('ContentType', 'table');

totag = readtable(ToTag_xlsFile);
%avoid case mis-matches
totag.Properties.VariableNames = lower(totag.Properties.VariableNames);
%initialize the new numbers
t = NaN(size(totag(:,1)));
totag = addvars(totag,t,t,t,t, 'NewVariableNames', {'lat' 'lon' 'depth' 'niskin'});
%t = repmat(' ', size(totag(:,1),1),1);
%totag = addvars(totag,t, 'NewVariableNames', {'datetime'})

%load the ship's underway data
uw = webread([apibase 'underway/' cruise '.csv'], options);

%load the event log
%event_log = webread([apibase 'events/' cruise '.csv'], options);

%% find the underway matchups
%find the underway rows
%tagstr = 'tag1';
%uwind = strmatch('underway', totag.(tagstr), 'exact');
%if isempty(uwind)
%    tagstr = 'tag2';
%    uwind = strmatch('underway', totag.(tagstr));
%end
%if isempty(uwind)
tagstr = 'sample_type';
uwind = strmatch('underway', totag.(tagstr));
%end

%get the underway match up
IFCB_mdate = IFCB_file2date(cellstr(totag.filename(uwind)));
IFCB_match_uw_results = IFCB_match_uw(totag.filename(uwind), IFCB_mdate, uw);
IFCB_match_uw_results.cruise = repmat(cellstr(cruise),size(IFCB_match_uw_results,1),1);

totag.lat(uwind) = IFCB_match_uw_results.lat;
totag.lon(uwind) = IFCB_match_uw_results.lon;
totag.depth(uwind) = NaN;
totag.datetime(uwind) = {''};
totag.cast(uwind) = {''};
totag.niskin(uwind) = NaN;

%% Find and assign the Cast and Niskin numbers from the log file
%find the cast rows in totag
castind = strmatch('cast', totag.(tagstr));
if ~isempty(castind)
    %load the ship's CTD btl data
    bottle_data = webread([apibase 'ctd/' cruise '/bottles.csv'], options);
    
    %find and read the cruise-specific IFCB logfile
    logfilelist = readtable('\\sosiknas1\IFCB_data\NESLTER_transect\to_tag\NESLTER_transect_IFCB_log_tag_filelist.xlsx');
    cruise_ind = strmatch(cruise, logfilelist.cruise);
    if ~isempty(logfilelist.log_file{cruise_ind})
        % set opts so make sure 'Cast' columne is read a char to handle option
        % for underway_discretes as 'UW1', 'UW2', etc.
        if ~isequal('NA', logfilelist.log_file{cruise_ind})
            opts = detectImportOptions(logfilelist.log_file{cruise_ind});
            opts = setvartype(opts, 'Cast', 'char');
            opts.DataRange = 'A2';
            IFCBlog = readtable(logfilelist.log_file{cruise_ind}, opts);
            %avoid case mis-matches
            IFCBlog.Properties.VariableNames = lower(IFCBlog.Properties.VariableNames);
            %check if any files are missing from tag file
            [~,ia] = setdiff(totag.filename(castind), IFCBlog.filename);
            if ~isempty(ia)
                disp('Missing from log file: ')
                disp(totag.filename(castind(ia)))
                keyboard
            end
        end
        %check if any files missing from log file
        %[~,ia] = setdiff( IFCBlog.filename, totag.filename(castind));
        %if ~isempty(ia)
        %    disp('Missing from log file: ')
        %    disp(IFCBlog.filename(ia))
        %    keyboard
        %end
        
        %now do the match up and assign cast and niskin in totag
        if ~isempty(castind)
            [~,ia,ib] = intersect(totag.filename(castind), IFCBlog.filename);
            if ~iscell(IFCBlog.cast(ib(1)))
                totag.cast(castind(ia)) = cellstr(num2str(IFCBlog.cast(ib)));
            else
                totag.cast(castind(ia)) = IFCBlog.cast(ib);
            end
            if iscell(IFCBlog.niskin(ib(1)))
                totag.niskin(castind(ia)) = str2double(IFCBlog.niskin(ib));
            else
                totag.niskin(castind(ia)) = IFCBlog.niskin(ib);
            end
            %% Find the cast matchup data
            IFCB_match_btl_results = IFCB_match_btl(totag.filename(castind),str2num(char(totag.cast(castind))), totag.niskin(castind), bottle_data);
            IFCB_match_btl_results.cruise = repmat(cellstr(cruise),size(IFCB_match_btl_results,1),1);
            totag.lat(castind) = IFCB_match_btl_results.lat;
            totag.lon(castind) = IFCB_match_btl_results.lon;
            totag.datetime(castind) = IFCB_match_btl_results.datetime;
            totag.depth(castind) = IFCB_match_btl_results.depth;
        end
        %check for casts with no info in bottle file
        ind = find(isnan(totag.lat(castind)));
        if ~isempty(ind)
            ctd_meta = webread([apibase 'ctd/' cruise '/metadata.csv'], options);
            unqcast = unique(totag.cast(castind(ind)));
            for count = 1:length(unqcast)
                %iii = find(totag.cast(castind(ind)) == unqcast(count));
                %cind = find(ctd_meta.cast == unqcast(count)));
                iii = strmatch(unqcast(count), totag.cast(castind(ind)));
                cind = find(ctd_meta.cast == str2num(unqcast{count}));
                totag.lat(castind(ind(iii))) = ctd_meta.latitude(cind);
                totag.lon(castind(ind(iii))) = ctd_meta.longitude(cind);
                totag.datetime(castind(ind(iii))) = ctd_meta.date(cind);
            end
            if ismember('depth_override', IFCBlog.Properties.VariableNames)
                totag.depth(castind(ind)) = IFCBlog.depth_override(ib(ind));
            else
                disp(totag(castind(ind),:))
                disp('Seems like some bottle files are missing. Need to add depth_override column to IFCBlog file for above cases.')
            end
            disp('No match up with bottle file info for samples listed below. Lat/Lon from CTD metadata, depth from IFCB_log')
            disp(totag(castind(ind),:))
            disp('No match up with bottle file info for samples listed above. Lat/Lon from CTD metadata, depth from IFCB_log')
            disp('Hit any key to continue')
            pause
        end
        %now check if any IFCB_match_btl results need to be filled with basic meta data
        if ~isempty(castind)
            ind = find(isnan(IFCB_match_btl_results.lat));
            [~,a,b] = intersect(IFCB_match_btl_results.pid(ind), totag.filename);
            IFCB_match_btl_results.datetime(ind(a)) = totag.datetime(b);
            IFCB_match_btl_results.depth(ind(a)) = totag.depth(b);
            IFCB_match_btl_results.lat(ind(a)) = totag.lat(b);
            IFCB_match_btl_results.lon(ind(a)) = totag.lon(b);
            IFCB_match_btl_restuls.mdate = datenum(IFCB_match_btl_results.datetime, 'yyyy-mm-dd hh:MM:SS+00:00');
        end
    end %end if log file exists
end %end if ~isempty(castind)

%% include any tags and comments from totag file and from IFCB log file
if exist(IFCBlog, 'var')
    [~,ia,ib] = intersect(totag.filename, IFCBlog.filename);
    tagindlog = find(strncmp('tag', IFCBlog.Properties.VariableNames,3));
    tagindtag = find(strncmp('tag', totag.Properties.VariableNames,3));
    for c = 1:length(ia)
        %[IFCBlog{ib,tagindlog} totag{ia,tagindtag}]
        tags_temp = setdiff(unique([IFCBlog{ib(c),tagindlog} totag{ia(c),tagindtag}]),char([]));
        for tn = 1:length(tags_temp)
            totag.(['tag' num2str(tn)])(ia(c)) = tags_temp(tn);
        end
        comments_temp = setdiff(unique([IFCBlog.comments(ib(c)) totag.comments(ia(c))]),char([]));
        if ~isempty(comments_temp)
            totag.comments(ia(c)) = join(comments_temp, '; ');
        end
    end   
end

%% find the underway discrete matchups
%uwdind = strmatch('underway_discrete', totag.tag2);
uwdind = find(ismember(totag.(tagstr), {'underway_discrete' 'bucket'}));

if ~isempty(uwdind)
    [~,ia,ib] = intersect(totag.filename(uwdind), IFCBlog.filename);
    if numel(ia) ~= numel(uwdind)
        disp('Missing underway_discrete match up')
        keyboard
    end
    temp_mdate = IFCB_file2date(cellstr(totag.filename(uwdind)));
    IFCB_mdate = NaN(size(temp_mdate));
    %%case for using the event log -- ABANDONED in favor of date/time of first IFCB file
    %%ic = strmatch('IFCB discrete', event_log.Action);
    %%[~,id] = ismember(IFCBlog.cast(ib), event_log.Cast(ic));
    %%ie = find(id);
    %%IFCB_mdate = totag.datetime(uwdind);
    %%IFCB_mdate(ia(ie)) = cellstr(event_log.dateTime8601(ic(id(ie))));
    
    %assume time of first file in set is approximately sample collection time
        unquw = unique(IFCBlog.cast(ib));
    for ii = 1:length(unquw)
        ind = strmatch(unquw(ii), IFCBlog.cast(ib), 'exact');
        IFCB_mdate(ia(ind)) = min(temp_mdate(ia(ind)));
    end
    %these are the non-blank datetime entries from the original log file
    if ismember('datetime', IFCBlog.Properties.VariableNames)
        ind = setdiff(1:length(uwdind), strmatch(' ', IFCBlog.datetime(ib)));
        IFCB_mdate(ia(ind)) = datenum(IFCBlog.datetime(ib(ind)), 'yyyy-mm-dd hh:MM:ss+00:00');
    end
    IFCB_match_uwdiscrete_results = IFCB_match_uw(totag.filename(uwdind), IFCB_mdate, uw);
    IFCB_match_uwdiscrete_results.cruise = repmat(cellstr(cruise),size(IFCB_match_uwdiscrete_results,1),1);
    totag.lat(uwdind) = IFCB_match_uwdiscrete_results.lat;
    totag.lon(uwdind) = IFCB_match_uwdiscrete_results.lon;
    totag.depth(uwdind) = NaN;
    totag.datetime(uwdind) = cellstr(datestr(IFCB_mdate,'yyyy-mm-dd hh:MM:ss+00:00'));
    %totag.cast(uwdind) = NaN;
    totag.cast(uwdind(ia)) = IFCBlog.cast(ib);
    totag.niskin(uwdind) = NaN;
end
totag.depth(find(ismember(totag.(tagstr), {'bucket'}))) = 0;

%% save results
%totag.Properties.VariableNames(strmatch('Tag1', totag.Properties.VariableNames)) = {'Cruise'}
if strmatch(tagstr, 'tag2') %old case
    totag.cruise = repmat(cellstr(cruise),size(totag,1),1);
end
f = strsplit(ToTag_xlsFile, '.');
writetable(totag, [f{1} '_meta.csv']);
disp(['CSV file for dashboard upload: ' f{1} '_meta.csv'])
[p f] = fileparts(f{1});
p = regexprep(p, 'to_tag', 'match_up\');
f = regexprep(f, 'to_tag', '');
if ~exist(p, 'dir'), mkdir(p), end
disp('Match-up ancillary data files: ')
save([p f 'uw_match'], 'IFCB_match_uw_results')
disp([p f 'uw_match.mat'])
if ~isempty(castind)
    save([p f 'cast_match'], 'IFCB_match_btl_results')
    disp([p f 'cast_match.mat'])
end
if ~isempty(uwdind)
    save([p f 'uwdiscrete_match'], 'IFCB_match_uwdiscrete_results')
    disp([p f 'uwdiscrete_match.mat'])
end
end

