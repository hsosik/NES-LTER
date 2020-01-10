
% dirpath = '\\sosiknas1\Backup\LTER\20171022_AR24\AR24A\IFCB_data\data9\'; %AR24A - Leg1 Pioneer 9
% dirpath = '\\sosiknas1\Backup\LTER\20171022_AR24\AR24B\IFCB_data\data9\'; %AR24B - Leg2 Pioneer 9
dirpath = '\\sosiknas1\Backup\LTER\20171022_AR24\AR24C\IFCB_data\data9\'; %AR24C - Leg3 Pioneer 9
dirpath = 'D:\LTER\20190809_EN644\IFCB_data\discrete\';
d = dir([dirpath 'D*']);
isub = [d(:).isdir]; %# returns logical vector
d = d(isub); % only folders
dirname = struct2cell(d);
dirlist = dirname(1,:)'; %names of dirs as cell array
filename={};
for count=1:length(dirlist)
    filedir      = [dirpath cell2mat(dirlist(count))];
    tempfiles    = dir([filedir '\*.hdr']);
    filename=[filename; {tempfiles.name}'];
end
temp=char(filename);
temp=temp(:,1:end-4);
open temp
