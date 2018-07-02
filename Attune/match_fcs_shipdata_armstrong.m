basepath = '\\sosiknas1\Backup\SPIROPA\20180414_AR29\';
fpath = [basepath '\underway\proc\'];
%fpath = '\\sosiknas1\Backup\LTER\20180404_AR28\underway\proc\';

f = [basepath 'Attune\Summary\FCSfileinfo.mat'];  %AR29 Cruise
if exist(f,'file')
    load(f)
else
   [FCSfileinfo] = FCS_DateTimeList([basepath '\Attune\FCSexport\']); %AR
    save(f, 'FCSfileinfo')
end

load([basepath '\Attune\Summary\compiled_stats.mat']) %AR29
[~,a,b] = intersect(FCSfileinfo.filelist, fcsfile_syn);
fcsmatch.mdate_start(b) = FCSfileinfo.matdate_start(a);


flist = dir([fpath 'AR*.csv']);
temp = flist(1).name;
if temp(1:4) == 'AR29' %change to is equal or strmatch if it doesn't work
    flist = flist(3:end); %remove two days before cruise
end;
mdate = [];
lat = [];
lon = [];
flr = [];
sbe45S = [];
sbe45T = [];
for ii = 1:length(flist)
    t = importdata([fpath char(flist(ii).name)]);
    s = char(t.textdata(3:end,1));
    s2 = char(t.textdata(3:end,2));
    s3 = repmat(' ', size(s,1), 1);
    mdate = [mdate; datenum([s s2 s3])];
    lat = [lat; t.data(:,1)];
    lon = [lon; t.data(:,2)];
    sbe45S = [sbe45S; t.data(:,30)];
    sbe45T = [sbe45T; t.data(:,31)];
    flr = [flr; t.data(:,32)];
end
clear s s2 s3 t ii fpath flist

for ii = 1:length(fcsmatch.mdate_start)
    [~,a] = min(abs(mdate-fcsmatch.mdate_start(ii)));
    fcsmatch.lat(ii) = lat(a);
    fcsmatch.lon(ii) = lon(a);
    fcsmatch.sbe45T(ii) = sbe45T(a);
    fcsmatch.sbe45S(ii) = sbe45S(a);
    fcsmatch.flr(ii) = flr(a);
end


%save AR29_underway