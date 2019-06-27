function [] = process_attune_v2(basepath, assign_class_function, plot_flag, filetype2exclude)
%function [] = process_attune_v2(basepath, assign_class_function, plot_flag)
%
%input: basepath for cruise or project with Attune NxT data
%       assign_class_function: name of function specifying gating details for class boundaries
%       plot_flag: true to plot cytograms with classes by color
%output: results stored to data file in summary folder under basepath
%
%Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2019

fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'Summary' filesep];
classpath = [outpath filesep 'class' filesep];
if ~exist(outpath, 'dir')
    mkdir(outpath)
end
if ~exist(classpath, 'dir')
    mkdir(classpath)
end
if plot_flag
    warning off
end
if exist([outpath 'FCSfileinfo.mat'], 'file')
    Attune = load([outpath 'FCSfileinfo']);
    [ FCSfileinfo ] = FCS_DateTimeList(fpath, [outpath 'FCSfileinfo']); %check if any new files to append  
else
    [ FCSfileinfo ] = FCS_DateTimeList(fpath);
    Attune.FCSfileinfo = FCSfileinfo; clear FCSfileinfo
end
save([outpath 'FCSfileinfo'], 'FCSfileinfo')

for iii = 1:length(filetype2exclude)    
    t = strmatch(filetype2exclude{iii}, Attune.FCSfileinfo.filelist);
    if ~isempty(t)
        f = fieldnames(Attune.FCSfileinfo);
        for ii = 1:length(f)
            Attune.FCSfileinfo.(f{ii})(t) = [];
        end
    end
end
    
filelist = Attune.FCSfileinfo.filelist;
% Creating the variables

temp = NaN(length(filelist),3);
SynCount = temp;
SynBiovol = temp;
SynCarbon = temp;
temp = NaN(length(filelist),6);
EukCount = temp;
EukBiovol = temp;
EukCarbon = temp;
QC_flowrate = NaN(length(filelist),2);

for count = 10:length(filelist)
    if ~rem(count,10)
        disp([num2str(count) ' of ' num2str(length(filelist))])
    end
    filename = [fpath filelist{count}];
    disp(filename)
    %reading in each FCS file with fca_readfcs
    [fcsdat,fcshdr] =fca_readfcs(filename);
    %t = find(fcsdat(:,12)>500 & fcsdat(:,3)>100);  %April 2018
    t = find(fcsdat(:,12)>200 & fcsdat(:,3)>200);  %RBH??
    QC_flowrate(count,1) = (median(fcsdat(t,3)./fcsdat(t,12)));
    QC_flowrate(count,2) = (std(fcsdat(t,3)./fcsdat(t,12)));

    QC_flag = 0; %default bad
    if (QC_flowrate(count,2)<2 & QC_flowrate(count,1)<1.5)
        QC_flag = 1; %set to good 
    end
    
    [~,fname] = fileparts(filename);
    [ class ] = eval([assign_class_function '( fcsdat, fcshdr, plot_flag, fname, QC_flag );']); clear fname
    
    %compute the cell volume from side scattering (Area, integrated signal)
    temp = fcsdat(:,3); %SSC-A
    temp(temp<0) = NaN; %FIGURE OUT LATER -- why are there SSC-A < 0 !!!???
    volume = 10.^(1.3.*log10(temp) - 2.9); %cubic microns, APPROX calibration from Summer 2018
    carbon = biovol2carbon( volume, 0 ); % carbon, picograms per cell
    notes = ['Class 1= Euk, Class 2 = Syn, Class 0 = junk; Cell volume in cubic microns; assign_class_function = ' assign_class_function];
    
    save([classpath regexprep(filelist{count}, '.fcs', '')], 'class', 'volume', 'notes')
    
    SynCount(count,1) = sum(class==2);
    EukCount(count,1) = sum(class==1);
    SynBiovol(count,1) = sum(volume(class==2));
    EukBiovol(count,1) = nansum(volume(class==1));  %%CHECK PROBLEM WHERE SOME SSC-A values are negative
    SynCarbon(count,1) = sum(carbon(class==2));
    EukCarbon(count,1) = nansum(carbon(class==1));
    
    diam = (volume*3/4/pi).^(1/3)*2; %equivalent spherical diam, micrometers
    ind = find(diam<=2 & class==2);
    SynCount(count,2) = length(ind);
    SynBiovol(count,2) = sum(volume(ind));
    SynCarbon(count,2) = sum(carbon(ind));
    ind = find(diam>2 & diam<=5 & class==2);
    SynCount(count,3) = length(ind);
    SynBiovol(count,3) = sum(volume(ind));
    SynCarbon(count,3) = sum(carbon(ind));
    
    ind = find(diam<=2 & class==1);
    EukCount(count,2) = length(ind);
    EukBiovol(count,2) = sum(volume(ind));
    EukCarbon(count,2) = sum(carbon(ind));
    ind = find(diam>2 & diam<=5 & class==1);
    EukCount(count,3) = length(ind);
    EukBiovol(count,3) = sum(volume(ind));
    EukCarbon(count,3) = sum(carbon(ind));
    ind = find(diam>5 & diam<=10 & class==1);
    EukCount(count,4) = length(ind);
    EukBiovol(count,4) = sum(volume(ind));
    EukCarbon(count,4) = sum(carbon(ind));
    ind = find(diam>10 & diam<=20 & class==1);
    EukCount(count,5) = length(ind);
    EukBiovol(count,5) = sum(volume(ind));
    EukCarbon(count,5) = sum(carbon(ind));
    ind = find(diam>20 & diam<=50 & class==1);
    EukCount(count,6) = length(ind);
    EukBiovol(count,6) = sum(volume(ind));
    EukCarbon(count,6) = sum(carbon(ind));
    ind = find(diam>50 & class==1);
    EukCount(count,7) = length(ind);
    EukBiovol(count,7) = sum(volume(ind));
    EukCarbon(count,7) = sum(carbon(ind));
    
    
    AttuneTable = table(Attune.FCSfileinfo.filelist, datetime(Attune.FCSfileinfo.matdate_start, 'ConvertFrom', 'datenum'), datetime(Attune.FCSfileinfo.matdate_stop, 'ConvertFrom', 'datenum'), Attune.FCSfileinfo.vol_analyzed/1e6, 'VariableNames', {'Filename' 'StartDate' 'StopDate' 'VolAnalyzed_ml'});
    AttuneTable = [AttuneTable array2table(SynCount, 'VariableNames', {'SynCountTotal' 'SynCountlt2' 'SynCount2_5'})];
    AttuneTable = [AttuneTable array2table(EukCount, 'VariableNames', {'EukCountTotal' 'EukCountlt2' 'EukCount2_5' 'EukCount5_10' 'EukCount10_20' 'EukCount20_50' 'EukCountgt50'})];
    AttuneTable = [AttuneTable array2table(SynBiovol, 'VariableNames', {'SynBiovolTotal' 'SynBiovollt2' 'SynBiovol2_5'})];
    AttuneTable = [AttuneTable array2table(EukBiovol, 'VariableNames', {'EukBiovolTotal' 'EukBiovollt2' 'EukBiovol2_5' 'EukBiovol5_10' 'EukBiovol10_20' 'EukBiovol20_50' 'EukBiovolgt50'})];
    AttuneTable = [AttuneTable array2table(SynCarbon, 'VariableNames', {'SynCarbonTotal' 'SynCarbonlt2' 'SynCarbon2_5'})];
    AttuneTable = [AttuneTable array2table(EukCarbon, 'VariableNames', {'EukCarbonTotal' 'EukCarbonlt2' 'EukCarbon2_5' 'EukCarbon5_10' 'EukCarbon10_20' 'EukCarbon20_50' 'EukCarbongt50'})];
    AttuneTable.QC_flowrate_median = QC_flowrate(:,1);
    AttuneTable.QC_flowrate_std = QC_flowrate(:,2);
    
end
AttuneTable = sortrows(AttuneTable, 'StartDate');

save([outpath '\AttuneTable'],'AttuneTable')
disp(['Result file saved:'])
disp([outpath '\AttuneTable'])

end