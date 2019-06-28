function [ FCSfileinfo ] = FCS_DateTimeList(fcs_path, FCSfileinfo_name)
%function [ FCSfileinfo ] = FCS_DateTimeList(fcs_path, FCSfileinfo_name)
%input: 
%   fcspath - path to a directory of fcs files
%   FCSfileinfo - optional input of existing compiled info (file with full path) if already available for append of new files
%output: a structure called fcsfileinfo with matlab date start and stop
%times and other fcs hdr information
%
%Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2019

if ~exist('fcs_path', 'var')
    fcs_path = uigetdir(pwd, 'Pick a Directory of FCS files')
else
    if ~exist(fcs_path, 'dir')
        disp('WARNING: Directory not found. Check input path.')
        FCSfileinfo = [];
        return
    end
end

fcslist = dir(fullfile(fcs_path, '*.fcs'));
fcslist = {fcslist.name}';

if exist('FCSfileinfo_name', 'var')
    load(FCSfileinfo_name)
    fcslist = setdiff(fcslist, FCSfileinfo.filelist);
    a = length(FCSfileinfo.filelist);
    b = length(fcslist);
    FCSfileinfo.filelist(a+1:a+b) = fcslist;
    FCSfileinfo.matdate_start(a+1:a+b) = NaN;
    FCSfileinfo.matdate_stop(a+1:a+b) = NaN;
    FCSfileinfo.vol_analyzed(a+1:a+b) = NaN;
else  %initialize
    FCSfileinfo.filelist = fcslist;
    FCSfileinfo.matdate_start = NaN(size(fcslist));
    FCSfileinfo.matdate_stop = FCSfileinfo.matdate_start;
    FCSfileinfo.vol_analyzed = FCSfileinfo.matdate_start;
    a = 0;
    b = length(fcslist);
end

if b > 0
    for ii = a+1:a+b
        if ~rem(ii,10)
            disp([num2str(ii) ' of ' num2str(a+b)])
        end
        [~,fcshdr] = fca_readfcs(fullfile(fcs_path, fcslist{ii-a}));
        FCSfileinfo.matdate_start(ii) = datenum([fcshdr.date ', ' fcshdr.starttime]);
        FCSfileinfo.matdate_stop(ii) = datenum([fcshdr.date ', ' fcshdr.stoptime]);
        FCSfileinfo.vol_analyzed(ii) = fcshdr.VOL;
    end
end

[~,ind] = sort(FCSfileinfo.matdate_start);
f = fields(FCSfileinfo)
for ii = 1:length(f)
    FCSfileinfo.(f{ii}) = FCSfileinfo.(f{ii})(ind,:);
end
