fpath = '\\sosiknas1\Backup\SPIROPA\20180414_AR29\underway\proc\';
%fpath = '\\sosiknas1\Backup\LTER\20180404_AR28\underway\proc\';

%f = '\\sosiknas1\Lab_data\Attune\EN608\Summary\FCSfileinfo.mat';  %FIX
f = '\\sosiknas1\Backup\SPIROPA\20180414_AR29\Attune\FCSfileinfo.mat';  %AR29 Cruise
if exist(f,'file')
    load(f)
else
   [ FCSfileinfo ] = FCS_DateTimeList( '\\sosiknas1\Backup\SPIROPA\20180414_AR29\Attune\FCSexport' ); %FIX
   %[ FCSfileinfo ] = FCS_DateTimeList( '\\sosiknas1\Lab_data\Attune\EN608\ExportedFCS' ); %FIX
    save(f, 'FCSfileinfo')
end

load \\sosiknas1\Backup\SPIROPA\20180414_AR29\Attune\Summary\compiled_stats.mat %FIX
[~,a,b] = intersect(FCSfileinfo.filelist, fcsfile_syn);
fcsmatch.mdate_start(b) = FCSfileinfo.matdate_start(a);


flist = dir([fpath 'AR18*.csv']);
flist = flist(3:end); %remove two days before cruise
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