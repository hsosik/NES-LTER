function [] = process_attune_beads(basepath, assign_class_function, plot_flag, filetype2include)

% make output directories and create metadata structure (FCSfileinfo)

outpath = [basepath filesep 'bead_results' filesep];
beadfigpath = [outpath filesep 'bead_plots'];
%classpath = [outpath 'class' filesep];

if ~exist(outpath, 'dir')
    mkdir(outpath)
end
if ~exist(beadfigpath, 'dir')
    mkdir(beadfigpath)
end

% if ~exist(classpath, 'dir')
%     mkdir(classpath)
% end

if plot_flag
    warning off
end

if exist([outpath 'FCSfileinfo.mat'], 'file')
    Attune = load([outpath 'FCSfileinfo']);
    FCSfileinfo = FCS_DateTimeList(basepath, [outpath 'FCSfileinfo']); %check if any new files to append  
else
    FCSfileinfo = FCS_DateTimeList(basepath);
end

startdate = min(FCSfileinfo.matdate_start);
Attune.FCSfileinfo = FCSfileinfo;
save([outpath 'FCSfileinfo'], 'FCSfileinfo')

% for iii = 1:length(filetype2exclude)    
%     t = strmatch(filetype2exclude{iii}, Attune.FCSfileinfo.filelist);
%     if ~isempty(t)
%         f = fieldnames(Attune.FCSfileinfo);
%         for ii = 1:length(f)
%             Attune.FCSfileinfo.(f{ii})(t) = [];
%         end
%     end
% end

% identify bead run files and determine mean bead SSC-H for cruise
% may need to adjust epsilon and/or minpts depending on cruise

bead_files = dir([basepath '\FCB_bead_check*']);
bead_files = {bead_files.name};
[~,ia,ib] = intersect(FCSfileinfo.filelist, bead_files);
[~,is] = sort(FCSfileinfo.matdate_start(ia)); 
bead_files = bead_files(ib(is));
f = fieldnames(FCSfileinfo);
for ii = 1:length(f)
      bead_FCSfileinfo.(f{ii}) = FCSfileinfo.(f{ii})(ia(is));
end

if isempty(bead_files)
    bead_files = dir([basepath '\Daily bead check*']);
    disp('STOP--need new code for cases with PT beads')
    keyboard
end
bead_ssch = NaN(length(bead_files), 1);
bead_qc = bead_ssch;
hv = bead_ssch;
bead_time = NaT(length(bead_files), 1);
   
beadstat = table;
beadstat_temp = table;
for ii = 1:length(bead_files)
    disp(ii)
    [fcsdat, fcshdr] = fca_readfcs([basepath bead_files{ii}]);
    [~, ~, m1, ~, temp_table, QC_flag] = assign_class_beads_algorithm_v4(fcsdat, fcshdr, 1);
    beadstat(ii,:) = temp_table;
    beadstat_temp.hv(ii,:) = {fcshdr.par.hv};
    bead_ssch(ii) = m1;
    bead_qc(ii) = QC_flag;
    bead_time(ii) = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    ssc_ch(ii) = 19; %GL1-H, low sensitivity SSC, with OD2
    if startdate < datenum('1-Aug-2019')
        ssc_ch(ii) = 12;
    end
    hv(ii) = fcshdr.par(ssc_ch).hv;
    figure(99)
    subplot(2,2,1), title(bead_files(ii), 'interpreter', 'none')
    subplot(2,2,2), title([datestr(bead_time(ii)) '; SSC hv = ' num2str(hv(ii)) ' (' fcshdr.par(ssc_ch).name ')'])
    print(figure(99), fullfile(beadfigpath, regexprep(bead_files{ii}, '.fcs', '.png')), '-dpng')
    %pause
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
legend(num2str(hv_list'))
xlabel('Time')
ylabel('Bead SSC-H')
hold off
parname = {fcshdr.par.name};
save([outpath 'beadstat'],'bead*', 'parname', 'ssc_ch', 'hv') 

disp(['Result file saved:'])
disp([outpath 'beadstat'])

end
