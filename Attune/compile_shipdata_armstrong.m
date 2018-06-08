fpath = '\\multiproc\data_on_memory\underway\proc\';

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
    mdate = [mdate; datenum([s s3 s2])];
    lat = [lat; t.data(:,1)];
    lon = [lon; t.data(:,2)];
    sbe45S = [sbe45S; t.data(:,30)];
    sbe45T = [sbe45T; t.data(:,31)];
    flr = [flr; t.data(:,32)];
end
clear s s2 s3 t ii fpath flist

save AR29_underway