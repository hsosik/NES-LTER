function [] = process_attune_lter(basepath, assign_class_function, startdate, plot_flag, filetype2exclude)

% make output directories and create metadata structure (FCSfileinfo)

fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'bead_calibrated' filesep];
classpath = [outpath 'class' filesep];

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
    FCSfileinfo = FCS_DateTimeList(fpath, [outpath 'FCSfileinfo']); %check if any new files to append  
else
    FCSfileinfo = FCS_DateTimeList(fpath);
end

Attune.FCSfileinfo = FCSfileinfo;
save([outpath 'FCSfileinfo'], 'FCSfileinfo')
clear FCSfileinfo

for iii = 1:length(filetype2exclude)    
    t = strmatch(filetype2exclude{iii}, Attune.FCSfileinfo.filelist);
    if ~isempty(t)
        f = fieldnames(Attune.FCSfileinfo);
        for ii = 1:length(f)
            Attune.FCSfileinfo.(f{ii})(t) = [];
        end
    end
end

% identify bead run files and determine mean bead SSC-H for cruise
% may need to adjust epsilon and/or minpts depending on cruise

bead_files = dir([fpath '\FCB_bead_check*']);
bead_ssch = NaN(length(bead_files), 1);
bead_qc = bead_ssch;
hv = bead_ssch;
bead_time = NaT(length(bead_files), 1);
for i = 1:length(bead_files)
    [fcsdat, fcshdr] = fca_readfcs([fpath bead_files(i).name]);
    [~, ~, m1, ~, QC_flag] = assign_class_beads_algorithm_v3(fcsdat, 0.2, 100, startdate, 0);
    bead_ssch(i) = m1;
    bead_qc(i) = QC_flag;
    bead_time(i) = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    if datetime(startdate) < datetime('1-Aug-2019')
        hv(i) = fcshdr.par(12).hv;
    else
        hv(i) = fcshdr.par(19).hv;
    end
end
bead_qc = logical(bead_qc);
bead_mean = NaN(length(unique(hv)), 2);
bead_mean(:,1) = unique(hv);
for i = 1:length(bead_mean)
    bead_mean(i,2) = mean(bead_ssch(bead_qc & hv==bead_mean(i,1)));
end

% plot bead data for user verification

figure; hold on;
hv_list = unique(hv);
for i=1:length(hv_list)
    plot(bead_time(hv==hv_list(i) & bead_qc), bead_ssch(hv==hv_list(i) & bead_qc), '.')
end
xlabel('Time')
ylabel('Bead SSC-H')
hold off

% read in and process Attune data files

filelist = Attune.FCSfileinfo.filelist;
AttuneTable = table(Attune.FCSfileinfo.filelist, datetime(Attune.FCSfileinfo.matdate_start, 'ConvertFrom', 'datenum'), datetime(Attune.FCSfileinfo.matdate_stop, 'ConvertFrom', 'datenum'), Attune.FCSfileinfo.vol_analyzed/1e6, 'VariableNames', {'Filename' 'StartDate' 'StopDate' 'VolAnalyzed_ml'});

% Creating the variables
numClass = 6;
diamEdges = [0 2 5 10 20 50 inf];
numBins = length(diamEdges)-1;
Count = NaN(length(filelist),numClass);
Biovol = Count;
Carbon = Count;
CountBin = NaN(length(filelist),numBins);
BiovolBin = CountBin;
CarbonBin = CountBin;
QC_flowrate = NaN(length(filelist),2);

for count = 1:length(filelist)
    if ~rem(count,10)
        disp([num2str(count) ' of ' num2str(length(filelist))])
    end
    filename = [fpath filelist{count}];
    disp(filename)
    [fcsdat,fcshdr] = fca_readfcs(filename);
    if datetime(startdate) < datetime('1-Aug-2019')
        file_hv(i) = fcshdr.par(12).hv;
    else
        file_hv(i) = fcshdr.par(19).hv;
    end
    t = find(fcsdat(:,12)>200 & fcsdat(:,3)>200);
    QC_flowrate(count,1) = (median(fcsdat(t,3)./fcsdat(t,12)));
    QC_flowrate(count,2) = (std(fcsdat(t,3)./fcsdat(t,12)));

    QC_flag = 0; %default bad
    if (QC_flowrate(count,2)<2 & QC_flowrate(count,1)<1.5)
        QC_flag = 1; %set to good 
    end
    
    [~,fname] = fileparts(filename);
    class = eval([assign_class_function '( fcsdat, fcshdr, plot_flag, fname, QC_flag, Attune.FCSfileinfo.matdate_start(count) );']); clear fname
    volume = 10.^(1.4225*log10(fcsdat(:,12)./bead_mean(bead_mean(:,1)==file_hv,2)) + 1.1432);
    carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
    notes = ['Class 1= Euk, Class 2 = Syn, Class 5 = lowPEeuks, Class 4 = hiPEeuks, Class 5 = Syn_euk_coincident1, Class 6 = Syn_euk_coincident2, Class 7 = noise; Class 0 = junk; Cell volume in cubic microns; assign_class_function = ' assign_class_function];
    class_labels = {'Euk_' 'Syn_' 'lowPEeuk_' 'hiPEeuk_' 'SynEuk1_' 'SynEuk2_'};
    save([classpath regexprep(filelist{count}, '.fcs', '')], 'class', 'volume', 'notes')
    
    for ii = 1:numClass
        Count(count,ii) = sum(class==ii);
        Biovol(count,ii) = nansum(volume(class==ii));
        Carbon(count,ii) = nansum(carbon(class==ii));
    end
    
    diam = (volume*3/4/pi).^(1/3)*2; %equivalent spherical diam, micrometers
    for ii = 1:length(diamEdges)-1
        ind = find(diam>=diamEdges(ii) & diam<diamEdges(ii+1) & class~=0);
        CountBin(count,ii) = size(ind,1);
        BiovolBin(count,ii) = sum(volume(ind));
        CarbonBin(count,ii) = nansum(carbon(ind));
    end
end

for ii = 1:numBins
    binlabel{ii} = ['X' num2str(diamEdges(ii)) 'to' num2str(diamEdges(ii+1))];
end
AttuneTable = [AttuneTable array2table(Count, 'VariableNames', regexprep(class_labels, '_', '_count'))];
AttuneTable = [AttuneTable array2table(Biovol, 'VariableNames', regexprep(class_labels, '_', '_biovolume'))];
AttuneTable = [AttuneTable array2table(Carbon, 'VariableNames', regexprep(class_labels, '_', '_carbon'))];
AttuneTable = [AttuneTable array2table(CountBin, 'VariableNames', regexprep(binlabel, 'X', 'count_'))];
AttuneTable = [AttuneTable array2table(BiovolBin, 'VariableNames', regexprep(binlabel, 'X', 'biovolume'))];
AttuneTable = [AttuneTable array2table(CarbonBin, 'VariableNames', regexprep(binlabel, 'X', 'carbon'))];

AttuneTable.QC_flowrate_median = QC_flowrate(:,1);
AttuneTable.QC_flowrate_std = QC_flowrate(:,2);
    
AttuneTable = sortrows(AttuneTable, 'StartDate');

save([outpath '\AttuneTable'],'AttuneTable', 'assign_class_function')
disp(['Result file saved:'])
disp([outpath '\AttuneTable'])

end
