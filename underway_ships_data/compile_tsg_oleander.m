function [tsg] = compile_tsg_oleander %(year2do)
% function [tsg] = compile_tsg_oleander(year2do)
% function to compile underway TSG info for Oleander from CSV files for each transit
%
% Heidi M. Sosik, Woods Hole Oceanographic Institution, April 2025
% used to compile metadata for underway IFCB sampling

basepath = '\\sosiknas1\Lab_data\Oleander\tsg_csv\';
flist = dir([basepath '*.csv']);

outfile = [basepath '\compiled_underway\'];
if ~exist(outfile, 'dir')
    mkdir(outfile)
end
outfile = [outfile 'uw_compiled'];

uw = readtable([basepath flist(1).name]);
uw.Properties.VariableNames = lower(uw.Properties.VariableNames);
for ii = 2:length(flist)
    disp(ii)
    %uw = [uw; readtable([basepath flist(ii).name])];
    temp = readtable([basepath flist(ii).name], 'MissingRule','fill', 'TreatAsMissing',"NODATA");
    temp.Properties.VariableNames = lower(temp.Properties.VariableNames);
    uw = outerjoin(uw, temp, 'MergeKeys', true);
end
uw.Properties.VariableNames = lower(uw.Properties.VariableNames);
uw.datetime = datetime(uw.date_gmt, 'InputFormat', 'uuuu/MM/dd')+uw.time_gmt;
uw.matdate = datenum(uw.datetime); %for input to IFCB_match_uw

notes = {'Heidi Sosik, WHOI, produced with compile_tsg_oleander.m from '};
save(outfile, 'uw', 'notes')
disp('results saved: ') 
disp(outfile)

%%
%now save the IFCB_match table
basepath = '\\sosiknas1\IFCB_data\Oleander\match_up\';
metaT =  webread('https://ifcb-data.whoi.edu/api/export_metadata/Oleander', weboptions('timeout', 30));
iso8601format = 'yyyy-MM-dd HH:mm:ss+00:00';

metaT = metaT(~metaT.skip,:);
metaT = metaT(~strcmp(metaT.sample_type, 'beads'),:);

IFCB_mdate = datenum(datetime(metaT.sample_time, 'InputFormat', iso8601format));

IFCB_match_uw_results = IFCB_match_uw(metaT.pid, IFCB_mdate, uw);

save([basepath 'uw_match'], 'IFCB_match_uw_results')
disp('Match-up ancillary data file: ')
disp([basepath 'uw_match.mat'])

%%
match_uw = table;
slist = {'pid' 'lat' 'lon' 'temperature' 'salinity' 'mdate'};
orig_var = {'sbe45_temp' 'sbe45_sal'};
IFCB_match_uw_results.Properties.VariableNames(orig_var) = {'temperature' 'salinity'};
[~,ia] = ismember(slist, IFCB_match_uw_results.Properties.VariableNames);
match_uw = IFCB_match_uw_results(:,ia);

save([basepath 'compiledTS_tables'], 'match*')
disp('results saved:')
disp([basepath 'compiledTS_tables'])

end