load C:\work\LTER\Attune\code\FCSfileinfo  %from FCS_DateTimeList
load C:\work\LTER\Attune\code\compiled_stats
[~,a,b] = intersect(FCSfileinfo.filelist, fcsfile_syn);
fcsmatch.mdate_start(b) = FCSfileinfo.matdate_start(a);

t = importdata('C:\work\LTER\Attune\EN608\Data60Sec_Cruise_20180117-160100.csv', ',', 134);
yd = t.textdata(135:end,3);
mdate = datenum('1-1-2018') + str2num(char(yd));
lat = str2num(char(t.textdata(135:end,5)));
lon = str2num(char(t.textdata(135:end,6)));
temperature = t.textdata(135:end,72);%Met-RMY-Trans-SST5
s = char(temperature);
s = strmatch(' ', s(:,1));
temperature(s) = {'NaN'};
temperature = str2num(char(temperature));

for ii = 1:length(fcsmatch.mdate_start)
    [~,a] = min(abs(mdate-fcsmatch.mdate_start(ii)));
    fcsmatch.lat(ii) = lat(a);
    fcsmatch.lon(ii) = lon(a);
    fcsmatch.temperature(ii) = temperature(a);
end


