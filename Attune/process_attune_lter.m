function [] = process_attune_lter(basepath, assign_class_function, plot_flag, filetype2exclude)

% make output directories and create metadata structure (FCSfileinfo)

fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'bead_calibrated_test' filesep];
beadfigpath = [outpath filesep 'bead_plots'];
classpath = [outpath 'class' filesep];

if ~exist(outpath, 'dir')
    mkdir(outpath)
end
if ~exist(beadfigpath, 'dir')
    mkdir(beadfigpath)
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

startdate = min(FCSfileinfo.matdate_start);
Attune.FCSfileinfo = FCSfileinfo;
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

% identify bead run files and determine mean bead SSC-H for cruise
% may need to adjust epsilon and/or minpts depending on cruise

bead_files = dir([fpath '\FCB_bead_check*']);
bead_files = {bead_files.name};
[~,ia,ib] = intersect(FCSfileinfo.filelist, bead_files);
[~,is] = sort(FCSfileinfo.matdate_start(ia)); 
bead_files = bead_files(ib(is));
f = fieldnames(FCSfileinfo);
for ii = 1:length(f)
      bead_FCSfileinfo.(f{ii}) = FCSfileinfo.(f{ii})(ia(is));
end
clear FCSfileinfo

if isempty(bead_files)
    bead_files = dir([fpath '\Daily bead check*']);
    disp('STOP--need new code for cases with PT beads')
    keyboard
end
bead_ssch = NaN(length(bead_files), 1);
bead_qc = bead_ssch;
hv = bead_ssch;
bead_time = NaT(length(bead_files), 1);
ssc_ch = 19; %GL1-H, low sensitivity SSC, with OD2
if startdate < datenum('1-Aug-2019')
    ssc_ch = 12;
end
   
beadstat = table;
beadstat_temp = table;
for ii = 1:length(bead_files)
%    disp(ii)
    [fcsdat, fcshdr] = fca_readfcs([fpath bead_files{ii}]);
    [~, ~, m1, ~, temp_table, QC_flag] = assign_class_beads_algorithm_v4(fcsdat, fcshdr, 1);
    beadstat(ii,:) = temp_table;
    beadstat_temp.hv(ii,:) = {fcshdr.par.hv};
    bead_ssch(ii) = m1;
    bead_qc(ii) = QC_flag;
    bead_time(ii) = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    hv(ii) = fcshdr.par(ssc_ch).hv;
    figure(99)
    subplot(2,2,1), title(bead_files(ii), 'interpreter', 'none')
    subplot(2,2,2), title([datestr(bead_time(ii)) '; SSC hv = ' num2str(hv(ii)) ' (' fcshdr.par(ssc_ch).name ')'])
   % print(figure(99), fullfile(beadfigpath, regexprep(bead_files{ii}, '.fcs', '.png')), '-dpng')
 %   pause
end
for ii = 1:length(beadstat_temp.hv(:)), if length(beadstat_temp.hv{ii})==0, beadstat_temp.hv{ii} = [NaN]; end; end;
beadstat.hv = cell2mat(beadstat_temp.hv); clear beadstat_temp
bead_qc = logical(bead_qc);
bead_mean = NaN(length(unique(hv)), 2);
bead_mean(:,1) = unique(hv);
%for i = 1:length(bead_mean)
for ii = 1:size(bead_mean,1) %heidi
    bead_mean(ii,2) = mean(bead_ssch(bead_qc & hv==bead_mean(ii,1)));
end

% plot bead data for user verification

figure; hold on;
hv_list = unique(hv);
for ii=1:length(hv_list)
    plot(bead_time(hv==hv_list(ii) & bead_qc), bead_ssch(hv==hv_list(ii) & bead_qc), '.')
end
legend(num2str(hv_list(:)))
xlabel('Time')
ylabel('Bead SSC-H')
hold off
parname = {fcshdr.par.name};
save([outpath 'beadstat'],'bead*', 'parname', 'ssc_ch') 
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
    %file_hv(i) = fcshdr.par(19).hv; 
    file_hv = fcshdr.par(ssc_ch).hv; %heidi
    t = find(fcsdat(:,12)>200 & fcsdat(:,3)>200);
    QC_flowrate(count,1) = (median(fcsdat(t,3)./fcsdat(t,12)));  %Heidi: DOUBLE CHECK IF THIS SHOULD BE ssc_ch instead of 12
    QC_flowrate(count,2) = (std(fcsdat(t,3)./fcsdat(t,12)));

    QC_flag = 0; %default bad
    if (QC_flowrate(count,2)<2 & QC_flowrate(count,1)<1.5)
        QC_flag = 1; %set to good 
    end
    
    [~,fname] = fileparts(filename);
    class = eval([assign_class_function '( fcsdat, fcshdr, plot_flag, fname, QC_flag, Attune.FCSfileinfo.matdate_start(count) );']); clear fname
%Heidi: WELL, THIS SEEMS WRONG...if beads are ch 19, then should get cells ch 19??
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
