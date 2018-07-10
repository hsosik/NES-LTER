function [header, dat, headertitles, dattitles] = fcbreadraw(datapath, filename);

fid = fopen([datapath filename],'r','ieee-be');            
datraw=fread(fid,'uint16');
fclose(fid);

setbytes = 1812; %9*200 (2-byte) + 12 (10x2-byte and 4x1-byte)
if length(datraw) < 500*setbytes; 
    setnum = floor(length(datraw)/setbytes);
    dat = reshape(datraw(1:setnum*setbytes),setbytes,setnum);
else
    %dat = reshape(datraw,length(datraw)/500,500);  %split into 500 records (500*200 = 100000 events per file)
    dat = reshape(datraw,setbytes,500);  %split into 500 records (500*200 = 100000 events per file)
end;
temp = [floor(dat(9,:)/256); rem(dat(9,:),256); floor(dat(10,:)/256); rem(dat(10,:),256)];
dat = [dat(1:8,:); temp; dat(11:end,:)];

header = dat(1:14,:)';
dat = dat(15:end,:);
dat = reshape(dat(:),9,length(dat(:))/9)';         

headertitles = {'start hr' 'start sec' 'start msec' 'end hr' 'end sec' 'end msec' 'temperature' 'humidity' 'start port' 'end port' 'start syr#' 'end syr#' 'start syr pos' 'end syr pos.'};
dattitles = {'PE,lo gn' 'PE, hi gn' 'FLS,lo gn' 'FLS, hi gn' 'CHL,lo gn' 'CHL, hi gn' 'SSC,lo gn' 'SSC, hi gn' 'CHL peak' };