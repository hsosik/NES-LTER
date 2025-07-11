% HB2503 SAMOS data is missing all of 07 June 2025, fill the long gap with
% files from the pCO2 system
% Note, these don't match precisely wrt time/lat/lon and temperature and
% salinty where both systems have data -- the differences are minor in the
% big scheme...but not as close as I expected 
%
% read this without the header line because it seems messed up...
T = readtable("\\sosiknas1\Lab_data\LTER\NESLTER_broadscale\20250528_HB2503\PCO2_underway\pCO2-15sec.tmp.ELG", 'FileType', 'text');

T = renamevars(T, {'Var1' 'Var2' 'Var5' 'Var6' 'Var9' 'Var10'}, ["date" "time" "lattemp" "lontemp" "SSPS" "TS"]);

t = char(strrep(T.lattemp,'N',''));
T.latitude_fullres = double(string(t(:,1:2))) + double(string(t(:,3:end-1)))/60;

t = char(strrep(T.lontemp,'W',''));
T.longitude_fullres = -1*(double(string(t(:,1:2))) + double(string(t(:,3:end-1)))/60);

T.datetime = datetime(T.date)+T.time; T.datetime2 = datetime(T.date)+T.Var4;
T.matdate = datenum(T.datetime); 
T.mdate_fullres = T.matdate;

Tind = T.datetime >= datetime(2025,6,7,0,0,0) & T.datetime <= datetime(2025,6,7,24,0,0);
S = T(Tind,{'matdate', 'mdate_fullres', 'latitude_fullres', 'longitude_fullres', 'TS', 'SSPS'});

load('\\sosiknas1\Lab_data\LTER\NESLTER_broadscale\20250528_HB2503\compiled_underway\HB2503uw_compiled.mat')
uw(end+1:end+height(S),S.Properties.VariableNames)=S;
uw{end-height(S)+1:end,~ismember(uw.Properties.VariableNames,S.Properties.VariableNames)} = NaN;

save('\\sosiknas1\Lab_data\LTER\NESLTER_broadscale\20250528_HB2503\compiled_underway\HB2503uw_compiled.mat', 'uw', '-append')


