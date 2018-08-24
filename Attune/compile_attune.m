function [compiled_stats] =compile_attune(basepath) % Input: path to directory of exported stats files
%output 
fpath = [basepath '\ExportedStats\'];
outpath = [basepath '\Summary\'];

% Extracting files out of the directory sorts NES out from SFD
%first it will populate with NES titled files but if empty will go for SFD
%PROBLEM: ONLY 2565 files are analyzed when there are more
filelist = dir([fpath 'NES*']); 
if isempty(filelist) == 1
     filelist = dir([fpath 'SFD*']);
end

filelist = {filelist.name}';
flistchar = char(filelist);
% 
% dstr = flistchar(:,15:end-5);
% dstr = flistchar(:,10:end-5)
% mdate = datenum(dstr);
% [~,s] = sort(mdate);

% filelist = filelist(s);

SynConc = [];
SynCount = [];
SynYcv = [];
EukConc = SynConc;
EukCount = [];
EukYcv = [];
fcsfile_syn = SynConc;
fcsfile_euk = SynConc;

for count = 1:length(filelist)
    disp(filelist(count))
    disp(count)
    itable = importfile([fpath filelist{count}]);
    ii = strmatch( 'Syn', itable.Gate);
    sample = itable.Sample(ii); 
    a = ~isnan(sample);clc
    ii = ii(a);
    SynConc = [SynConc; itable.Concentration(ii)];
    SynCount = [SynCount; itable.Count(ii)];
    SynYcv = [SynYcv; itable.YCV(ii)];
    exp = itable.Experiment(ii);
    sample = itable.Sample(ii); 
    temp = cellstr([char(exp) repmat('_Group_day0_Sample(', length(exp),1) num2str(sample) repmat(').fcs', length(exp),1)]);
    temp = regexprep(temp, '( ', '(');
    temp = regexprep(temp, '( ', '(');
    fcsfile_syn = [fcsfile_syn; temp]; clear temp 
    ii = strmatch( 'Euk', itable.Gate);
    sample = itable.Sample(ii); 
    a = ~isnan(sample);
    ii = ii(a);
    EukConc = [EukConc; itable.Concentration(ii)]; 
    EukCount = [EukCount; itable.Count(ii)];
    EukYcv = [EukYcv; itable.YCV(ii)];
    exp = itable.Experiment(ii);
    sample = itable.Sample(ii); 
    temp = cellstr([char(exp) repmat('_Group_day0_Sample(', length(exp),1) num2str(sample) repmat(').fcs', length(exp),1)]);
    temp = regexprep(temp, '( ', '(');
    temp = regexprep(temp, '( ', '(');
    temp = regexprep(temp, '(NaN)', '');
    fcsfile_euk = [fcsfile_euk; temp]; clear temp 
    FileSampleCount(count) = length(ii);
end

compiled_stats.SynConc = SynConc;
compiled_stats.EukConc = EukConc;
compiled_stats.SynCount = SynCount;
compiled_stats.SynYcv = SynYcv;
compiled_stats.EukCount = EukCount;
compiled_stats.EukYcv = EukYcv;
save([outpath 'compiled_stats'], 'fcsfile*', 'SynConc', 'EukConc','SynCount','SynYcv', 'EukCount','EukYcv');
end