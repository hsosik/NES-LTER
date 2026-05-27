function [header, dat, headertitles, dattitles] = fcbreadraw2(datapath, filename);
%version to reflect extra header byte added March, 13 2006; modified from
%fcbreadraw.m
%special case for set of files in 2024 with only 3 records per set instead
%of 100 (data transmission troubles)

fid = fopen([datapath filename],'r','ieee-be');            
datraw=fread(fid,'uint16');
fclose(fid);

 year = datenum(str2num(filename(6:9)),0,0,0,0,0);

%setbytes = 1812; %9*200 (2-byte) + 12 (10x2-byte and 4x1-byte)
%setbytes = 1212; %12*100 (2-byte) + 12 (10x2-byte and 4x1-byte)
dtfile = datetime(str2num(filename(6:9)),1,str2num(filename(11:13)), str2num(filename(15:16)), str2num(filename(17:18)), str2num(filename(19:20)));
if dtfile > datetime(2024,9,12,0,0,0) & dtfile < datetime(2024,10,9,12,0,0)
    setbytes = 49; % 12*3 + 13
else
    setbytes = 1213; %12*100 (2-byte) + 13 (11x2-byte and 4x1-byte) % after adding 2 bytes to header for comparator value for solenoids (4 V = operating, i.e., open to ocean)
end
headerlength = 15;  %number of data values in header, 1 and 2 byte counted separately

%if length(datraw) < 1000*setbytes; 
    setnum = floor(length(datraw)/setbytes);
    dat = reshape(datraw(1:setnum*setbytes),setbytes,setnum);
%else
    %dat = reshape(datraw,length(datraw)/500,500);  %split into 500 records (500*200 = 100000 events per file)
%    dat = reshape(datraw,setbytes,1000);  %split into 500 records (500*200 = 100000 events per file)
%end;

%this section parses out the 1-byte values into separate columns, assumes
%their location in header
% temp = [floor(dat(9,:)/256); rem(dat(9,:),256); floor(dat(10,:)/256); rem(dat(10,:),256)];
% dat = [dat(1:8,:); temp; dat(11:end,:)];
temp = [floor(dat(9,:)/256); rem(dat(9,:),256); floor(dat(10,:)/256); rem(dat(10,:),256)];
dat = [dat(1:8,:); temp; dat(11:end,:)];

%read elements 1:? according to header length, 15 data values as of Mar 2006
header = dat(1:headerlength,:)';
% dat = dat(15:end,:);
% dat = reshape(dat(:),9,length(dat(:))/9)';         
dat = dat(headerlength+1:end,:);
dat = reshape(dat(:),12,length(dat(:))/12)';         

headertitles = {'start hr' 'start sec' 'start msec' 'end hr' 'end sec' 'end msec' 'temperature' 'humidity' 'start port' 'end port' 'start syr#' 'end syr#' 'start syr pos' 'end syr pos.'};
% dattitles = {'PE,lo gn' 'PE, hi gn' 'FLS,lo gn' 'FLS, hi gn' 'CHL,lo gn' 'CHL, hi gn' 'SSC,lo gn' 'SSC, hi gn' 'CHL peak' };
dattitles = {'PE,lo gn' 'PE, hi gn' 'FLS,lo gn' 'FLS, hi gn' 'CHL,lo gn' 'CHL, hi gn' 'SSC,lo gn' 'SSC, hi gn' 'CHL peak' };