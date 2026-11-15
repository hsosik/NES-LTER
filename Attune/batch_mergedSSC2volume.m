function batch_mergedSSC2volume(p)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%get the SSC merge stats
%load([p.outpath 'mergeTable'])
load([p.outpath 'FCSfileinfo'])

fcslist = dir([p.classpath '*.mat']); %just get the ones in the class folder (not the beads, etc.)
fcslist = regexprep({fcslist.name}', '.mat', '.fcs');
mergeT = table;
mergeT.filename = fcslist;
[~,ia] = ismember(fcslist, FCSfileinfo.fcslist);
mergeT.filetime = datetime(FCSfileinfo.matdate_start(ia), 'ConvertFrom', 'datenum');

addpath(fullfile(fileparts(mfilename('fullpath')), 'classifier_dev'))
beadstat = load_beadstat(p.outpath); % from process_beads_PT_adjust_3

% cruises before Aug 2019: one scatter channel (SSC 488/10; GL1-H is PE
% fluorescence in that layout, not side scatter), so no merge
if isfield(p, 'OD2setting') && ~strcmp(p.OD2setting, 'GL1')
    volume_488_2019(p, mergeT, beadstat)
    return
end

%get the 1 micron bead peaks (KDE mode; was centend2) for SSC and GL1, QC'd runs
%starting on cruise day and later
ind = beadstat.QC1 & beadstat.time>dateshift(min(mergeT.filetime),'start', 'day');
B1 = beadstat.peak1(ind, :);
beadT = table(B1.SSC_A, B1.SSC_H, beadstat.hv.SSC(ind), 'VariableNames', {'SSCA_1micron' 'SSCH_1micron' 'SSC_hv'});
beadSSCmean = groupsummary(beadT,'SSC_hv','mean'); %compute mean grouped by HV
beadT = table(B1.GL1_A, B1.GL1_H, beadstat.hv.GL1(ind), 'VariableNames', {'GL1A_1micron' 'GL1H_1micron' 'GL1_hv'});
beadGL1mean = groupsummary(beadT,'GL1_hv','mean'); %compute mean grouped by HV

sscstr = strcat("SSC-", p.SSCDIM);
gl1str = strcat("GL1-", p.SSCDIM);

SSCAmin = 500; %use SSC-H below this value
% SSC/GL1 merge window in 1 um bead units of GL1 (gl1_bdnorm), shared with the
% classifier (classifier_config cfg.gl1_merge_bead, default [0.5 2]): SSC below,
% GL1 above, mean of the two between. Was raw GL1 counts 500-1000, which moves in
% bead units when the GL1 voltage changes (0.45-0.9x bead at GL1 220 V, 2.6-5.2x
% at 180 V on AR104, where SSC is already compressed). Nov 2026,
% classifier_dev/pt_transfer/README.md
cfg = classifier_config();
GL1_merge_bead = cfg.gl1_merge_bead;

for filecount = 1:height(mergeT)
   if ~rem(filecount,50)
        disp(['volume calculating ' num2str(filecount) ' of ' num2str(height(mergeT))])
   end
 %   if ~isnan(mergeT.slope_median(filecount))
        filename = [p.fpath, mergeT.filename{filecount}];
        [fcsdat,fcshdr] = fca_readfcs(filename);
        fcsdat = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});

        file_hv = array2table([NaN [fcshdr.par.hv]], 'VariableNames', {fcshdr.par.name});
        mergeT.SSC_hv(filecount) = file_hv.(sscstr);
        mergeT.GL1_hv(filecount) = file_hv.(gl1str);
        %mergeT.filetime(filecount) = datetime([fcshdr.date, ' ', fcshdr.starttime]);
        
        bdrow = find(mergeT.SSC_hv(filecount)==beadSSCmean.SSC_hv);  %FIX THIS LATER IF KEEPING
        if p.SSCDIM=='A'
            ssc_bdnorm = fcsdat.(sscstr)./beadSSCmean.mean_SSCA_1micron(bdrow); % *.8; %
            gl1_bdnorm = fcsdat.(gl1str)./beadGL1mean.mean_GL1A_1micron(bdrow); % *.4;  %.42;  %HACK FIX!!!
        else
            fcsdat.(sscstr)(fcsdat.(sscstr)==0) = 1;
            ssc_bdnorm = fcsdat.(sscstr)./beadSSCmean.mean_SSCH_1micron(bdrow); % *.8; %
            gl1_bdnorm = fcsdat.(gl1str)./beadGL1mean.mean_GL1H_1micron(bdrow); % *.4;  %.42;  %HACK FIX!!!
        end
        sscmerge_bdnorm = ssc_bdnorm;
        itemp = gl1_bdnorm > GL1_merge_bead(2);
%     itemp = ssc_bdnorm>2 & gl1_bdnorm>0;
     sscmerge_bdnorm(itemp) = gl1_bdnorm(itemp);
     %%%%%itemp = ssc_bdnorm>0.5 & ssc_bdnorm<2 & gl1_bdnorm>0;
  %   itemp = ssc_bdnorm>1 & ssc_bdnorm<2 & gl1_bdnorm>0;
     itemp = gl1_bdnorm > GL1_merge_bead(1) & gl1_bdnorm <= GL1_merge_bead(2);
     sscmerge_bdnorm(itemp) = mean([ssc_bdnorm(itemp) gl1_bdnorm(itemp)],2);
     %c.merge_info = array2table([GL1_overlap_min GL1_overlap_max mean(fcsdat.(sscstr)(itemp)./fcsdat.(gl1str)(itemp)) median(fcsdat.(sscstr)(itemp)./fcsdat.(gl1str)(itemp))], 'variablenames', {'GL1min' 'GL1max' 'mean_SSC-to-GL1_overlap' 'median_SSC-to-GL1_overlap'});
     itemp = itemp & fcsdat.("SSC-H")./fcsdat.("BL3-H")<100; %get stats without non-chl triggers with odd SSC-GL1 relationship
    % mergeT.mean_SSC2GL1_overlap(filecount) = mean(fcsdat.(sscstr)(itemp)./fcsdat.(gl1str)(itemp));
    % mergeT.median_SSC2GL1_overlap(filecount) = median(fcsdat.(sscstr)(itemp)./fcsdat.(gl1str)(itemp));
     mergeT.mean_SSCA2GL1A_overlap(filecount) = mean(fcsdat.("SSC-A")(itemp)./fcsdat.("GL1-A")(itemp));
     mergeT.median_SSCA2GL1A_overlap(filecount) = median(fcsdat.("SSC-A")(itemp)./fcsdat.("GL1-A")(itemp));
     mergeT.mean_SSCH2GL1H_overlap(filecount) = mean(fcsdat.("SSC-H")(itemp)./fcsdat.("GL1-H")(itemp));
     mergeT.median_SSCH2GL1H_overlap(filecount) = median(fcsdat.("SSC-H")(itemp)./fcsdat.("GL1-H")(itemp));
     mergeT.mean_SSCA2SSCH_overlap(filecount) = mean(fcsdat.("SSC-A")(itemp)./fcsdat.("SSC-H")(itemp));
     mergeT.median_SSCA2SSCH_overlap(filecount) = median(fcsdat.("SSC-A")(itemp)./fcsdat.("SSC-H")(itemp));
     mergeT.mean_GL1A2GL1H_overlap(filecount) = mean(fcsdat.("GL1-A")(itemp)./fcsdat.("GL1-H")(itemp));
     mergeT.median_GL1A2GL1H_overlap(filecount) = median(fcsdat.("GL1-A")(itemp)./fcsdat.("GL1-H")(itemp));
     mergeT.n_overlap(filecount) = numel(itemp);
     volume_cubic_micron = NaN(size(sscmerge_bdnorm)); 
        volume_cubic_micronH = volume_cubic_micron;
        ii = sscmerge_bdnorm>0;
        prelim_flag = 0;
        [volume_cubic_micron(ii) vol_func_string] = Attune_SC2vol(sscmerge_bdnorm(ii),strcat('SSC-',p.SSCDIM), prelim_flag);
        if p.SSCDIM=='A'
            ii = fcsdat.(sscstr)<SSCAmin; %for small signals use SSC-H instead (A goes negative)
            volume_cubic_micron(ii) = Attune_SC2vol(fcsdat.("SSC-H")(ii)./beadSSCmean.mean_SSCH_1micron(bdrow),strcat('SSC-H'),0);
        end
        % sscmerge = fcsdat.(sscstr);
        % itemp = fcsdat.(gl1str)>gl1max | fcsdat.(sscstr)>sscmax; %use merged GL1 for large signals
        % sscmerge(itemp) = fcsdat.(gl1str)(itemp)*mergeT.slope_median(filecount); 
        % %now use merged GL1 for every other point the overlap
        % itemp = find(fcsdat.(gl1str)>gl1min & fcsdat.(gl1str)<gl1max & fcsdat.(sscstr)<sscmax);
        % sscmerge(itemp(1:2:end)) = fcsdat.(gl1str)(itemp(1:2:end)) *mergeT.slope_median(filecount);
        % bdrow = find(mergeT.SSC_hv(filecount)==beadSSCmean.SSC_hv);
        % sscmerge_bdnorm = sscmerge/beadSSCmean.mean_SSCA_1micron(bdrow);
        % volume_cubic_micron = NaN(size(sscmerge_bdnorm)); 
        % volume_cubic_micronH = volume_cubic_micron;
        % ii = sscmerge>0;
        % [volume_cubic_micron(ii) vol_func_string] = Attune_SC2vol(sscmerge_bdnorm(ii),strcat('SSC-',p.SSCDIM));
        % ii = sscmerge<SSCAmin; %for small signals use SSC-H instead (A goes negative)
        % volume_cubic_micron(ii) = Attune_SC2vol(fcsdat.("SSC-H")(ii)./beadSSCmean.mean_SSCH_1micron(bdrow),strcat('SSC-H'));
        %%ctemp = load([p.classpath regexprep(mergeT.filename{filecount},'.fcs', '.mat')]);
        %%c.class = ctemp.class;  clear ctemp %TEMP get rid of old stuff
        c = load([p.classpath regexprep(mergeT.filename{filecount},'.fcs', '.mat')]);
        c.vol_notes = {strcat('calibrated: ', string(datetime())); strcat(' using SSC-', p.SSCDIM); vol_func_string; ...
            sprintf('SSC/GL1 merge window %g-%g x 1 um bead (GL1)', GL1_merge_bead)};
        c.ssc_merge_bdnorm = sscmerge_bdnorm; c.beadSSCmean = beadSSCmean; c.beadGL1mean = beadGL1mean; c.file_hv.SSC = mergeT.SSC_hv(filecount); c.file_hv.GL1 = mergeT.GL1_hv(filecount);
        c.volume_cubic_microns = volume_cubic_micron;
%        c.merge_info = array2table([mergeT.slope_median(filecount) gl1min gl1max sscmax SSCAmin], 'variablenames', {'slope' 'GL1min' 'GL1max' 'SSCmax' 'SSCAmin'}); 
        save([p.classpath regexprep(mergeT.filename{filecount},'.fcs', '.mat')],'-struct', "c")
%    end
end
save([p.outpath 'mergeTable'], "mergeT", "GL1_merge_bead")
end

function volume_488_2019(p, mergeT, beadstat)
% Volumes for cruises before the Aug 2019 layout change (Nov 2026,
% classifier_dev/pt_transfer/README.md):
%   - SSC-H relative to the 1 um bead at the file's time and voltage
%     (bead_reference, as the classifier; beadstat may be FCB runs or the PT
%     performance-test transfer, whose voltages differ from the samples')
%   - p.OD2setting 'None' (AR34B, RB1904, EN627): decompress large-particle
%     SSC first (ssc_decompress, at the file's SSC voltage). Beyond the
%     invertible range the volume is a lower bound and vol_flag is true.
%   - p.OD2setting 'SSC': SSC with OD2 is linear for these signals
%   - volume: Attune_SC2vol(..., 'SSC-H_488_2019') (March/April 2019 cultures)
if p.SSCDIM ~= 'H'
    warning('batch_mergedSSC2volume:SSCDIM', 'pre-Aug 2019 calibration is for SSC-H; using SSC-H')
end
cfg = classifier_config();
warning('off', 'bead_reference:noQC') % no GL1 scatter channel before Aug 2019
decomp = strcmp(p.OD2setting, 'None');
cruise_note = "";
if contains(p.basepath, 'EN627')
    cruise_note = "EN627: scatter optics degraded (Feb 2019); SSC-based volumes uncertain, flag in analyses";
end
for filecount = 1:height(mergeT)
    if ~rem(filecount, 50)
        disp(['volume calculating ' num2str(filecount) ' of ' num2str(height(mergeT))])
    end
    [fcsdat, fcshdr] = fca_readfcs([p.fpath, mergeT.filename{filecount}]);
    j = strcmp({fcshdr.par.name}, 'SSC-H');
    hv = fcshdr.par(j).hv;
    ref = bead_reference(beadstat, datenum(mergeT.filetime(filecount)), cfg);
    beadpos = ref.bead_log.SSC + ref.k.SSC * log10(hv / ref.bead_hv.SSC);
    raw = fcsdat(:, j);
    s = log10(max(raw, 1)) - beadpos; % log10 SSC-H / 1 um bead
    flag = false(size(s));
    t = s;
    if decomp
        [t, flag, ~, tmax] = ssc_decompress(s, hv);
        t(flag) = tmax; % lower bound
    end
    volume_cubic_micron = NaN(size(s));
    ii = raw > 0;
    [volume_cubic_micron(ii), vol_func_string] = Attune_SC2vol(10.^t(ii), 'SSC-H_488_2019', 0);
    mergeT.SSC_hv(filecount) = hv;
    mergeT.n_vol_flag(filecount) = sum(flag);
    c = load([p.classpath regexprep(mergeT.filename{filecount}, '.fcs', '.mat')]);
    c.vol_notes = {strcat('calibrated: ', string(datetime())); ' using SSC-H, pre-Aug 2019 488 nm calibration'; vol_func_string; ...
        sprintf('OD2 %s; decompression %s', p.OD2setting, string(decomp)); ...
        sprintf('1 um bead (bead_reference): log10 SSC-H %.3f at %g V', ref.bead_log.SSC, ref.bead_hv.SSC); cruise_note};
    c.ssc_bdnorm = 10.^s;
    c.ssc_linear_bdnorm = 10.^t;
    c.vol_flag = flag; % true: SSC beyond the decompression range, volume is a lower bound
    c.file_hv = struct('SSC', hv);
    c.volume_cubic_microns = volume_cubic_micron;
    save([p.classpath regexprep(mergeT.filename{filecount}, '.fcs', '.mat')], '-struct', "c")
end
OD2setting = p.OD2setting;
save([p.outpath 'mergeTable'], "mergeT", "OD2setting")
end
