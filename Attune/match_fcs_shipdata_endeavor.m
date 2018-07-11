basepath = '\\sosiknas1\Lab_data\Attune\EN608\'
f = [basepath 'Summary\FCSfileinfo.mat']
if exist(f,'file')
    load(f)
else
   [ FCSfileinfo ] = FCS_DateTimeList( [basepath 'ExportedFCS\']);
   save(f, 'FCSfileinfo')
end

load([basepath 'Summary\compiled_stats'])
[~,a,b] = intersect(FCSfileinfo.filelist, fcsfile_syn);
fcsmatch.mdate_start(b) = FCSfileinfo.matdate_start(a);

t = importdata('\\sosiknas1\Backup\LTER\20180131_EN608\underway\scs_GPS_met\proc\cruise\Data60Sec_Cruise_20180117-160100_time change.csv', ',', 134);
yd = t.textdata(135:end, 3);
mdate = datenum('1-1-2018') + str2num(char(yd));
lat = str2num(char(t.textdata(135:end,5)));
lon = str2num(char(t.textdata(135:end,6)));
temperature = t.textdata(135:end,72);%Met-RMY-Trans-SST5
s = char(temperature);
s = strmatch(' ', s(:,1));
temperature(s) = {'NaN'};
temperature = str2num(char(temperature));

for ii = 1:length(fcsmatch.mdate_start)
    [~,a] = min(abs(mdate - fcsmatch.mdate_start(ii)));
    fcsmatch.lat(ii) = lat(a);
    fcsmatch.lon(ii) = lon(a);
    fcsmatch.temperature(ii) = temperature(a);
end

clear ii s lat lon temperature mdate yd t a b f
