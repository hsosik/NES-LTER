function [] = compile_IFCB_metadata_broadscale(ToTag_xlsFile, uw_compiled_file)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[~,f] = fileparts(ToTag_xlsFile);
cruise = strsplit(f, '_');
cruise = cruise{3};

totag = readtable(ToTag_xlsFile);
%avoid case mis-matches
totag.Properties.VariableNames = lower(totag.Properties.VariableNames);

%load the ship's underway data
load(uw_compiled_file)

%% find the underway matchups
%find the underway rows
uwind = strmatch('underway', totag.sample_type); 
%get the underway match up
IFCB_mdate = IFCB_file2date(cellstr(totag.filename(uwind)));
IFCB_match_uw_results = IFCB_match_uw(totag.filename(uwind), IFCB_mdate, uw);
totag.lat(uwind) = IFCB_match_uw_results.lat;
totag.lon(uwind) = IFCB_match_uw_results.lon;

%% save results
totag.cruise = repmat(cellstr(cruise),size(totag,1),1);
f = strsplit(ToTag_xlsFile, '.');
%writetable(totag, [f{1} '_meta.csv']);
disp(['CSV file for dashboard upload: ' f{1} '_meta.csv']) 
[p f] = fileparts(f{1});
p = regexprep(p, 'to_tag', 'match_up\');
f = regexprep(f, 'to_tag', '');
if ~exist(p, 'dir'), mkdir(p), end
save([p f 'uw_match'], 'IFCB_match_uw_results')
disp('Match-up ancillary data file: ')
disp([p f 'uw_match.mat'])

end

