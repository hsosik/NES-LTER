% process attune files for AR29

% make output directories and create metadata structure (FCSfileinfo)

filetype2exclude = {'SFD_AR29_Dilution', 'SFD_AR29_Grazer'};
plot_flag = 0;

fpath = '\\sosiknas1\Lab_data\Attune\cruise_data\20180414_AR29\FCS\';
outpath = '\\sosiknas1\Lab_data\Attune\cruise_data\20180414_AR29\bead_calibrated\';
classpath = '\\sosiknas1\Lab_data\Attune\cruise_data\20180414_AR29\bead_calibrated\class\';

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

% assign scattering channels
ssc_ch = 3;
ssch = 12;


% calculate PT bead mean and convert to FCB equivalent
pt_mean = 1.7486e4; % this is the "bead settings" file in AR29 FCS folder
fcb_mean = 0.0835*pt_mean;

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
    file_hv = fcshdr.par(ssc_ch).hv;
    t = find(fcsdat(:,12)>200 & fcsdat(:,3)>200);
    QC_flowrate(count,1) = (median(fcsdat(t,3)./fcsdat(t,12))); 
    QC_flowrate(count,2) = (std(fcsdat(t,3)./fcsdat(t,12)));

    QC_flag = 0; %default bad
    if (QC_flowrate(count,2)<2 & QC_flowrate(count,1)<1.5)
        QC_flag = 1; %set to good 
    end
    
    % correct for negative SSC-A values
    cf = fitlm(fcsdat(fcsdat(:,ssch)<1000,ssch), fcsdat(fcsdat(:,ssch)<1000, ssc_ch), 'Intercept', false);
    cf = cf.Coefficients.Estimate;
    fcsdat(fcsdat(:,ssc_ch)<0, ssc_ch) = cf*fcsdat(fcsdat(:,ssc_ch)<0, ssch);
   
    [~,fname] = fileparts(filename);
    class = assign_class_spiropa(fcsdat, fcshdr, plot_flag, fname, QC_flag, Attune.FCSfileinfo.matdate_start(count)); clear fname;
    ssca = fcsdat(:,ssc_ch);
    volume = 10.^(1.2232*log10(ssca./fcb_mean) + 1.0868);
    carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
    notes = 'Class 1= Euk, Class 2 = Syn, Class 5 = lowPEeuks, Class 4 = hiPEeuks, Class 5 = Syn_euk_coincident1, Class 6 = Syn_euk_coincident2, Class 7 = noise; Class 0 = junk; Cell volume in cubic microns';
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

save([outpath '\AttuneTable'],'AttuneTable')
disp(['Result file saved:'])
disp([outpath '\AttuneTable'])
