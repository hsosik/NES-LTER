function Attune_uw_match = match_Attune_underway_LTER(AttuneTable_fullname,uw_fullname)
% Load an AttuneTable summary mat file (output from process_Attune...); load
% ship's underway data from CSV, add closest match in time underway data
% as new columns to AttuneTable; save and output the resulting table
% Handles underway data input either as ship-specific CSV file read from
% disk or webread of result from NES-LTER API call 
%
% For API example:
%   uw_fullname = 'https://nes-lter-data.whoi.edu/api/underway/en608.csv';
%   AttuneTable_fullname = '\\sosiknas1\Lab_data\Attune\cruise_data\20180131_EN608\Summary\AttuneTable';
%   Attune_uw_match = match_Attune_underway_LTER(AttuneTable_fullname,uw_fullname);
%
% For ship's CSV from disk example:
%   uw_fullname = '\\sosiknas1\Lab_data\LTER\20200201_EN649\scs\proc\cruise\Data60Sec_Cruise_20200201-004500.csv';
%   AttuneTable_fullname = '\\sosiknas1\Lab_data\Attune\cruise_data\20200201_EN649\Summary\AttuneTable';
%   Attune_uw_match = match_Attune_underway_LTER(AttuneTable_fullname,uw_fullname);

load( AttuneTable_fullname )

if strncmp ('http', uw_fullname, 4) %case for API 
    uw = webread(uw_fullname);
    dt = datetime(uw.date, 'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');
    uw_mdate = datenum(dt);
else % works for R/V Endeavor cruise files
    fid = fopen(uw_fullname);
    t = fgetl(fid); fclose(fid);
    numHdrLines = str2num(regexprep(t, '#DataStartLine:',''));
    opts = delimitedTextImportOptions('VariableNamesLine', numHdrLines, 'DataLines', numHdrLines+1);
    uw = readtable(uw_fullname, opts);    
    dt = datetime(uw.DateTime_ISO8601, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS''Z');
    uw_mdate = datenum(dt);
% new case needed for Armstrong--easiest to use API if possible!
end

tdiff = NaN(size(AttuneTable,1),1);
match_ind = tdiff;
for ii = 1:length(tdiff)
    [tdiff(ii), match_ind(ii)] = min(abs(datenum(AttuneTable.StartDate(ii)-uw_mdate)));
end

uw_match = uw(match_ind,:);
Attune_uw_match = [AttuneTable uw_match];
%good = find(AttuneTable.QC_flowrate_std<2 & AttuneTable.QC_flowrate_median<1.5);

[AttuneTable_path, AttuneTable_file] = fileparts(AttuneTable_fullname); 

outFullFileName = fullfile(AttuneTable_path,[AttuneTable_file '_uw_match']);
save(outFullFileName, 'Attune_uw_match')%, 'good')
disp('results saved: ') 
disp(outFullFileName)
disp(['Max time difference in minutes: ' num2str(max(tdiff)*24*60)])

end

