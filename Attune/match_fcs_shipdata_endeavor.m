%function [fcsmatch]=match_fcs_shipdata_endeavor(basepath)

% basepath = '\\sosiknas1\Lab_data\Attune\EN608\'
basepath = '\\sosiknas1\Backup\LTER\20180720_EN617\Attune';
%basepath = 'E:\Attune_Data\EN608\';
%basepath = '\\sosiknas1\Backup\LTER\20180404_AR28\'

fpath = [basepath '\ExportedFCS\'];
f = [basepath '\Summary\Attune'];
outpath = [basepath '\Summary'];

load(f)

% Extracting files out of the directory sorts NES out from SFD
%first it will populate with NES titled files but if empty will go for SFD
%file string
filelist = dir([fpath 'NES*']);
filelist = {filelist.name}';
% flistchar = char(filelist);
% dstr = flistchar(:,15:end-27);
% mdate = datenum(dstr);
% [~,s] = sort(mdate);
% sortedlist = filelist(s);

[~,a,b] = intersect(Attune.FCSfileinfo.filelist, filelist);
Attune.fcsmatch.mdate_start(b) = Attune.FCSfileinfo.matdate_start(a);

% filename = '\\sosiknas1\Backup\LTER\20180131_EN608\underway\scs_GPS_met\proc\cruise\Data60Sec_Cruise_20180117-160100_time change.csv';
filename =  '\\sosiknas1\Backup\LTER\20180720_EN617\scs\proc\cruise\Data60Sec_Cruise_20180716-191500.csv';
% t = importdata(filename,',', 134);
t = importdata(filename,',', 136);

yd = t.textdata(137:end,3);
mdate = datenum('1-1-2018') + str2num(char(yd));
lat = str2num(char(t.textdata(137:end,5)));
lon = str2num(char(t.textdata(137:end,6)));
temperature = t.textdata(137:end,72);%Met-RMY-Trans-SST5
s = char(temperature);
s = strmatch(' ', s(:,1));
temperature(s) = {'NaN'};
temperature = str2num(char(temperature));
Attune.fcsmatch.lat = [];
Attune.fcsmatch.lon = [];
Attune.fcsmatch.temperature = [];
Attune.fcsmatch.salinity = [];
for ii = 1:length(Attune.fcsmatch.mdate_start)
    [~,a] = min(abs(mdate - Attune.fcsmatch.mdate_start(ii)));
    Attune.fcsmatch.lat(ii) = lat(a);
    Attune.fcsmatch.lon(ii) = lon(a);
    Attune.fcsmatch.temperature(ii) = temperature(a);
end

clear ii s lat lon temperature mdate yd t a b f
save([basepath '\Summary\Attune'],'Attune')
% save([outpath '\compiled_stats'],'fcsmatch')
