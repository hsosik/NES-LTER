function [ FCSfileinfo ] = FCS_DateTimeList(fcs_path, FCSfileinfo_name)
% modified 4/26/2021 to make output into a table, and to include quality
% control data

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

%first check if we are starting from scratch 
if exist('FCSfileinfo_name', 'var') 
    load(FCSfileinfo_name)
        if istable('FCSfileinfo_name') %if the file is a table 
            fcslist = setdiff(fcslist, FCSfileinfo.filelist);
            a = length(FCSfileinfo.filelist);
            b = length(fcslist);
            FCSfileinfo.filelist(a+1:a+b) = fcslist;
            FCSfileinfo.matdate_start(a+1:a+b) = NaN;
            FCSfileinfo.matdate_stop(a+1:a+b) = NaN;
            FCSfileinfo.vol_analyzed(a+1:a+b) = NaN;
        else % the variable is old and still a structure. Needs to be converted. 
            FCSfileinfo = struct2table(FCSfileinfo); 
            fcslist = setdiff(fcslist, FCSfileinfo.filelist);
            a = length(FCSfileinfo.filelist);
            b = length(fcslist);
            FCSfileinfo.filelist(a+1:a+b) = fcslist;
            FCSfileinfo.matdate_start(a+1:a+b) = NaN;
            FCSfileinfo.matdate_stop(a+1:a+b) = NaN;
            FCSfileinfo.vol_analyzed(a+1:a+b) = NaN;
        end
else  %initialize
            FCSfileinfo = table(fcslist);
            FCSfileinfo.matdate_start = NaN(size(fcslist));
            FCSfileinfo.matdate_stop = FCSfileinfo.matdate_start;
            FCSfileinfo.vol_analyzed = FCSfileinfo.matdate_start;
            FCSfileinfo.QC_flag = FCSfileinfo.matdate_start;
            a = 0;
             b = length(fcslist);
end

if b > 0 %now add to existing table
    for ii = a+1:a+b
        if ii == 20; %~rem(ii,10)
            disp([num2str(ii) ' of ' num2str(a+b)])
        end
        [fcsdat,fcshdr] = fca_readfcs(fullfile(fcs_path, fcslist{ii-a}));
        if ~(fcshdr.TotalEvents==0)
            FCSfileinfo.matdate_start(ii) = datenum([fcshdr.date ', ' fcshdr.starttime]);
            FCSfileinfo.matdate_stop(ii) = datenum([fcshdr.date ', ' fcshdr.stoptime]);
            FCSfileinfo.vol_analyzed(ii) = fcshdr.VOL;
            
            %adding quality control flags 
            t = find(fcsdat(:,12)>200 & fcsdat(:,3)>200);
            QC_flowrate(1) = (median(fcsdat(t,3)./fcsdat(t,12)));
            QC_flowrate(2) = (std(fcsdat(t,3)./fcsdat(t,12)));
            QC_flag = 0; %default bad
            if (QC_flowrate(2)<2 && QC_flowrate(1)<1.5)
                QC_flag = 1; %set to good
            end
            FCSfileinfo.QC_flag(ii) = QC_flag; 
            clear QC_flag QC_flowrate t
                        
        end
    end
end

%put table in chronological order
[~,ind] = sort(FCSfileinfo.matdate_start);
FCSfileinfo = FCSfileinfo(ind, :); 


end