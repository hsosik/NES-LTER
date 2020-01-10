function [] = compile_IFCB_metadata(ToTag_xlsFile)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[~,f] = fileparts(ToTag_xlsFile);
cruise = strsplit(f, '_');
cruise = cruise{3};
apibase = 'https://nes-lter-data.whoi.edu/api/';
options = weboptions('ContentType', 'table');

totag = readtable(ToTag_xlsFile);
%avoid case mis-matches
totag.Properties.VariableNames = lower(totag.Properties.VariableNames);
%load the ship's underway data
uw = webread([apibase 'underway/' cruise '.csv'], options);

%load the ship's underway data
bottle_data = webread([apibase 'ctd/' cruise '/bottles.csv'], options);

%load the event log
event_log = webread([apibase 'events/' cruise '.csv'], options);

%% find the underway matchups
%find the underway rows
uwind = strmatch('underway', totag.tag2); 
%get the underway match up
IFCB_mdate = IFCB_file2date(cellstr(totag.filename(uwind)));
IFCB_match_uw_results = IFCB_match_uw(totag.filename(uwind), IFCB_mdate, uw);
totag.lat(uwind) = IFCB_match_uw_results.lat;
totag.lon(uwind) = IFCB_match_uw_results.lon;
totag.depth(uwind) = NaN;
%totag.datetime(uwind) = {''};
totag.cast(uwind) = NaN;
totag.niskin(uwind) = NaN;

%% Find and assign the Cast and Niskin numbers from the log file
%find the cast rows in totag
castind = strmatch('cast', totag.tag2);
%find and read the cruise-specific IFCB logfile
logfilelist = readtable('\\sosiknas1\IFCB_data\NESLTER_transect\to_tag\NESLTER_transect_IFCB_log_tag_filelist.xlsx');
cruise_ind = strmatch(cruise, logfilelist.cruise);
IFCBlog = readtable(logfilelist.log_file{cruise_ind});
%avoid case mis-matches
IFCBlog.Properties.VariableNames = lower(IFCBlog.Properties.VariableNames); 
%check if any files are missing from tag file
[~,ia] = setdiff(totag.filename(castind), IFCBlog.filename);
if ~isempty(ia)
    disp('Missing from tag file: ')
    disp(totag.filename(castind(ia)))
    keyboard
end
%check if any files missing from log file
%[~,ia] = setdiff( IFCBlog.filename, totag.filename(castind));
%if ~isempty(ia)
%    disp('Missing from log file: ')
%    disp(IFCBlog.filename(ia))
%    keyboard
%end

%now do the match up and assign cast and niskin in totag
[~,ia,ib] = intersect(totag.filename(castind), IFCBlog.filename);
if iscell(IFCBlog.cast(ib(1)))
    totag.cast(castind(ia)) = str2double(IFCBlog.cast(ib))
else
    totag.cast(castind(ia)) = IFCBlog.cast(ib);
end    
if iscell(IFCBlog.niskin(ib(1)))
    totag.niskin(castind(ia)) = str2double(IFCBlog.niskin(ib));
else
    totag.niskin(castind(ia)) = IFCBlog.niskin(ib);
end
%% Find the cast matchup data
IFCB_match_btl_results = IFCB_match_btl(totag.filename(castind),totag.cast(castind), totag.niskin(castind), bottle_data);
totag.lat(castind) = IFCB_match_btl_results.lat;
totag.lon(castind) = IFCB_match_btl_results.lon;
totag.datetime(castind) = IFCB_match_btl_results.datetime;
totag.depth(castind) = IFCB_match_btl_results.depth;

%% find the underway discrete matchups
uwdind = strmatch('underway_discrete', totag.tag2);
[~,ia,ib] = intersect(totag.filename(uwdind), IFCBlog.filename);
%[~,ic,id] = intersect(event.Action, IFCBlog.cast(ib));
ic = strmatch('IFCB discrete', event_log.Action);
[~,id] = ismember(IFCBlog.cast(ib), event_log.Cast(ic));
ie = find(id);
%id(~id) = [];
IFCB_mdate = totag.datetime(uwdind);
IFCB_mdate(ia(ie)) = cellstr(event_log.dateTime8601(ic(id(ie))));

IFCB_match_uwdiscrete_results = IFCB_match_uw(totag.filename(uwdind), IFCB_mdate, uw);
totag.lat(uwdind) = IFCB_match_uwdiscrete_results.lat;
totag.lon(uwdind) = IFCB_match_uwdiscrete_results.lon;
totag.depth(uwdind) = NaN;
totag.datetime(uwdind) = {''};
totag.cast(uwdind) = NaN;
totag.niskin(uwdind) = NaN;

%% save results
%totag.Properties.VariableNames(strmatch('Tag1', totag.Properties.VariableNames)) = {'Cruise'}
totag.cruise = repmat(cellstr(cruise),size(totag,1),1);
f = strsplit(ToTag_xlsFile, '.');
writetable(totag, [f{1} '_meta.csv']);
disp(['CSV file for dashboard upload: ' f{1} '_meta.csv']) 
[p f] = fileparts(f{1});
p = regexprep(p, 'to_tag', 'match_up\');
f = regexprep(f, 'to_tag', '');
if ~exist(p, 'dir'), mkdir(p), end
save([p f 'cast_match'], 'IFCB_match_btl_results')
save([p f 'uw_match'], 'IFCB_match_uw_results')
disp('Match-up ancillary data files: ')
disp([p f 'cast_match.mat'])
disp([p f 'uw_match.mat'])

end

