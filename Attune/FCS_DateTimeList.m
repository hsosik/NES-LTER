function [ FCSfileinfo ] = FCS_DateTimeList( fcs_path )
%UNTITLED2
%input: path to a directory of fcs files
%output: a structure called fcsfileinfo with matlab date starts and stops
%   Detailed explanation goes here

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
FCSfileinfo.matdate_start = NaN(size(fcslist));
FCSfileinfo.matdate_stop = FCSfileinfo.matdate_start;
for ii = 1:length(fcslist)
    if ~rem(ii,10)
        disp([num2str(ii) ' of ' num2str(length(fcslist))])
    end
    [~,fcshdr] = fca_readfcs(fullfile(fcs_path, fcslist{ii}));
    FCSfileinfo.matdate_start(ii) = datenum([fcshdr.date ', ' fcshdr.starttime]);
    FCSfileinfo.matdate_stop(ii) = datenum([fcshdr.date ', ' fcshdr.stoptime]);
end

FCSfileinfo.filelist = fcslist;
