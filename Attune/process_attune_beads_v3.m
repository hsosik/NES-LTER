function [] = process_attune_v3(basepath, assign_class_function, plot_flag, filetype2include)
%function [] = process_attune_v3(basepath, assign_class_function, plot_flag)
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

for iii = 1:length(filetype2include)    
    t = strmatch(filetype2include{iii}, Attune.FCSfileinfo.filelist);
    %if ~isempty(t)
        f = fieldnames(Attune.FCSfileinfo);
        for ii = 1:length(f)
            Attune.FCSfileinfo.(f{ii}) = Attune.FCSfileinfo.(f{ii})(t);
        end
    %end
end
    
filelist = Attune.FCSfileinfo.filelist;
AttuneTable = table(Attune.FCSfileinfo.filelist, datetime(Attune.FCSfileinfo.matdate_start, 'ConvertFrom', 'datenum'), datetime(Attune.FCSfileinfo.matdate_stop, 'ConvertFrom', 'datenum'), Attune.FCSfileinfo.vol_analyzed/1e6, 'VariableNames', {'Filename' 'StartDate' 'StopDate' 'VolAnalyzed_ml'});

% Creating the variables
numClass = 3; %3 beads
numpar = 27;
Count = NaN(length(filelist),numClass);
Mean = NaN(length(filelist),numpar,numClass);
Mode = Mean;
Std = Mean;
bins = logspace(1,log10(1048576),256);
binc = [0 bins(1:end-1)+diff(bins)]; binc(end) = 1048575; 

QC_flowrate = NaN(length(filelist),2);

for filecount = 1:length(filelist)
    if ~rem(filecount,10)
        disp([num2str(filecount) ' of ' num2str(length(filelist))])
    end
    filename = [fpath filelist{filecount}];
    disp(filename)
    %reading in each FCS file with fca_readfcs
    [fcsdat,fcshdr] =fca_readfcs(filename);
    %t = find(fcsdat(:,12)>500 & fcsdat(:,3)>100);  %April 2018
    t = find(fcsdat(:,12)>200 & fcsdat(:,3)>200);  %RBH??
    QC_flowrate(filecount,1) = (median(fcsdat(t,3)./fcsdat(t,12)));
    QC_flowrate(filecount,2) = (std(fcsdat(t,3)./fcsdat(t,12)));

    QC_flag = 0; %default bad
    if (QC_flowrate(filecount,2)<2 & QC_flowrate(filecount,1)<1.5)
        QC_flag = 1; %set to good 
    end
    
    [~,fname] = fileparts(filename);
    [ class ] = eval([assign_class_function '( fcsdat, fcshdr, plot_flag, fname, QC_flag );']); clear fname
    
    %notes = ['Class 1= Euk, Class 2 = Syn, Class 0 = junk; Cell volume in cubic microns; assign_class_function = ' assign_class_function];
    notes = ['Class 1= Euk, Class 2 = Syn, Class 3 = lowPEeuks, Class 4 = hiPEeuks, Class 5 = oddPE_SSC_Chl, Class 0 = junk; Cell volume in cubic microns; assign_class_function = ' assign_class_function];
    class_labels = {'0.5micron_' '1micron_' 'McLane_'};
    
    for ii = 1:numClass
        Count(filecount,ii) = sum(class==ii);
        Mean(filecount,:,ii) = mean(fcsdat(class==ii,2:end));
        temp = mode(discretize(fcsdat(class==ii,2:end),[0 bins]));
        ind = (~isnan(temp));
        Mode(filecount,ind,ii) = binc(temp(ind));
        Std(filecount,:,ii) = std(fcsdat(class==ii,2:end));
    end
    
end
   
ParNames = matlab.lang.makeValidName({fcshdr.par(2:end).name}); 
AttuneTable.QC_flowrate_median = QC_flowrate(:,1);
AttuneTable.QC_flowrate_std = QC_flowrate(:,2);

AttuneTable_bead1 = [AttuneTable array2table(Count(:,1), 'VariableNames', {'Count'})];
AttuneTable_bead1 = [AttuneTable_bead1 array2table(squeeze(Mean(:,:,1)), 'VariableNames', append('Mean_', ParNames))];
AttuneTable_bead1 = [AttuneTable_bead1 array2table(squeeze(Mode(:,:,1)), 'VariableNames', append('Mode_', ParNames))];
AttuneTable_bead1 = [AttuneTable_bead1 array2table(squeeze(Std(:,:,1)), 'VariableNames', append('StdDev_', ParNames))];
AttuneTable_bead2 = [AttuneTable array2table(Count(:,2), 'VariableNames', {'Count'})];
AttuneTable_bead2 = [AttuneTable_bead2 array2table(squeeze(Mean(:,:,2)), 'VariableNames', append('Mean_', ParNames))];
AttuneTable_bead2 = [AttuneTable_bead2 array2table(squeeze(Mode(:,:,2)), 'VariableNames', append('Mode_', ParNames))];
AttuneTable_bead2 = [AttuneTable_bead2 array2table(squeeze(Std(:,:,2)), 'VariableNames', append('StdDev_', ParNames))];
AttuneTable_bead3 = [AttuneTable array2table(Count(:,3), 'VariableNames', {'Count'})];
AttuneTable_bead3 = [AttuneTable_bead3 array2table(squeeze(Mean(:,:,3)), 'VariableNames', append('Mean_', ParNames))];
AttuneTable_bead3 = [AttuneTable_bead3 array2table(squeeze(Mode(:,:,3)), 'VariableNames', append('Mode_', ParNames))];
AttuneTable_bead3 = [AttuneTable_bead3 array2table(squeeze(Std(:,:,3)), 'VariableNames', append('StdDev_', ParNames))];
    
AttuneTable_bead1 = sortrows(AttuneTable_bead1, 'StartDate');
AttuneTable_bead2 = sortrows(AttuneTable_bead2, 'StartDate');
AttuneTable_bead3 = sortrows(AttuneTable_bead3, 'StartDate');

f = fullfile(outpath, 'AttuneTable_bead')
save([outpath '\AttuneTable_bead'],'AttuneTable_bead*')
disp(['Result file saved:'])
disp(f)

end