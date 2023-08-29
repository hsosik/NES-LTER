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

% uw_fullname = '\\sosiknas1\Lab_data\Attune\cruise_data\20220216_AT46\at46_underway.csv'; 
load( AttuneTable_fullname )

if contains (uw_fullname, 'underway/ar') %Armstrong cruises don't import correctly using normal strat
    websave([AttuneTable_fullname(1:end-15) '\underway.csv'], uw_fullname)
    opts = delimitedTextImportOptions("NumVariables", 38);
    opts.DataLines = [1, Inf];
    opts.Delimiter = [","];
    opts.VariableNames = ["date", "dec_lat", "dec_lon", "spd", "hdt", "cog", "sog", "wxtp_ta", "wxts_ta", "wxtp_pa", "wxts_pa", "wxtp_ri", "wxts_ri", "wxtp_rc", "wxts_rc", "wxtp_dm", "wxts_dm", "wxtp_sm", "wxts_sm", "wxtp_ua", "wxts_ua", "wxtp_ts", "wxts_ts", "wxtp_td", "wxts_td", "barom_p", "barom_s", "rad_sw", "rad_lw", "par", "sbe45s", "sbe48t", "flr", "flow", "ssvdslog", "depth12", "depth35", "em122"];
    opts.VariableTypes = ["string", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "string", "string", "string"];
    uw = readtable([AttuneTable_fullname(1:end-15) '\underway.csv'], opts);
    dt = uw.date;
    for i = 1:length(dt)
        stupiddate = dt{i}; 
        uw_mdate(i) = datetime(stupiddate(1:19), 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        %uw_mdate(i) = datetime(stupiddate, 'InputFormat', 'MM/dd/yyyy HH:mm');
    end
     uw_mdate = datenum(uw_mdate); 
%elseif contains (uw_fullname, 'underway/hrs') %neither does HRS cruise
elseif strncmp ('http', uw_fullname, 4) %case for API this line is for standard case where all is working, also works for hrs cruise
    websave([AttuneTable_fullname(1:end-15) '\underway.csv'], uw_fullname)
    uw = readtable([AttuneTable_fullname(1:end-15) '\underway.csv'], 'Delimiter',',');
    dt = uw.date;
    for i = 1:length(dt)
        stupiddate = dt{i}; 
        uw_mdate(i) = datetime(stupiddate(1:19), 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        %uw_mdate(i) = datetime(stupiddate, 'InputFormat', 'MM/dd/yyyy HH:mm');
    end
     uw_mdate = datenum(uw_mdate); 
elseif contains (uw_fullname, 'at46')
    %atlantis cruise
    uw = readtable(uw_fullname,'Delimiter',',');
    dt = datetime(uw.date, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSS+00:00');
    uw_mdate = datenum(dt);

elseif contains(uw_fullname, '.mat') 
    uw = load(uw_fullname);
    uw = uw.uw; 
    uw_mdate = uw.matdate;
    %hrs = string(uw.TIME_GMT); 
    %for i = 1:length(hrs)
    %    hr = hrs{i}; 
    %    stupiddate = [uw.DATE_GMT{i} ' ' hr(1:8)]; 
    %    uw_mdate(i) = datetime(stupiddate, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
    %end
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


[AttuneTable_path, AttuneTable_file] = fileparts(AttuneTable_fullname); 

outFullFileName = fullfile(AttuneTable_path,[AttuneTable_file '_uw_match']);
save(outFullFileName, 'Attune_uw_match')
disp('results saved: ') 
disp(outFullFileName)
disp(['Max time difference in minutes: ' num2str(max(tdiff)*24*60)])

end
