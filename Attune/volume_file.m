function [c, row] = volume_file(c, fcsdat, fcshdr, V, filetime)
% function [c, row] = volume_file(c, fcsdat, fcshdr, V, filetime)
% Cell volumes for one FCS file, added to its class-file struct c
% (volume_cubic_microns, vol_notes and the scatter values used). row: the file's
% fields for mergeTable (preserved samples: the scatter voltages only).
%
% fcsdat, fcshdr  from fca_readfcs
% V               volume_setup
% filetime        datetime of the file start (underway, before Aug 2019: bead_reference)
%
% Errors when no bead position can be found (e.g. no bead run at the file's SSC
% voltage in 'merge' mode); the caller decides whether to skip the file.
%
% Nov 2026, split from batch_mergedSSC2volume and preserved_volume

row = struct();
switch V.mode
    case 'merge'
        [c, row] = merge_volume(c, fcsdat, fcshdr, V);
    case '488_2019'
        [c, row] = volume_488_2019(c, fcsdat, fcshdr, V, filetime);
    case 'preserved'
        c = preserved_file(c, fcsdat, fcshdr, V);
        row.SSC_hv = c.file_hv.SSC; row.GL1_hv = c.file_hv.GL1;
    otherwise
        error('volume_file:mode', 'unknown volume mode "%s"', V.mode)
end
end


function [c, row] = merge_volume(c, fcsdat, fcshdr, V)
% underway, Aug 2019 on: SSC merged with GL1 in 1 um bead units, Attune_SC2vol
SSCAmin = 500; %use SSC-H below this value
GL1_merge_bead = V.GL1_merge_bead; % SSC below, GL1 above, mean of the two between (classifier_config cfg.gl1_merge_bead)
sscstr = strcat("SSC-", V.SSCDIM);
gl1str = strcat("GL1-", V.SSCDIM);
beadSSCmean = V.beadSSCmean; beadGL1mean = V.beadGL1mean;
fcsdat = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
file_hv = array2table([NaN [fcshdr.par.hv]], 'VariableNames', {fcshdr.par.name});
row.SSC_hv = file_hv.(sscstr);
row.GL1_hv = file_hv.(gl1str);
bdrow = find(row.SSC_hv == beadSSCmean.SSC_hv);
if isempty(bdrow)
    error('volume_file:nobeads', 'no QC''d 1 um bead run at SSC %g V', row.SSC_hv)
end
if V.SSCDIM == 'A'
    ssc_bdnorm = fcsdat.(sscstr) ./ beadSSCmean.mean_SSCA_1micron(bdrow);
    gl1_bdnorm = fcsdat.(gl1str) ./ beadGL1mean.mean_GL1A_1micron(bdrow);
else
    fcsdat.(sscstr)(fcsdat.(sscstr) == 0) = 1;
    ssc_bdnorm = fcsdat.(sscstr) ./ beadSSCmean.mean_SSCH_1micron(bdrow);
    gl1_bdnorm = fcsdat.(gl1str) ./ beadGL1mean.mean_GL1H_1micron(bdrow);
end
sscmerge_bdnorm = ssc_bdnorm;
itemp = gl1_bdnorm > GL1_merge_bead(2);
sscmerge_bdnorm(itemp) = gl1_bdnorm(itemp);
itemp = gl1_bdnorm > GL1_merge_bead(1) & gl1_bdnorm <= GL1_merge_bead(2);
sscmerge_bdnorm(itemp) = mean([ssc_bdnorm(itemp) gl1_bdnorm(itemp)], 2);
itemp = itemp & fcsdat.("SSC-H") ./ fcsdat.("BL3-H") < 100; %stats without non-chl triggers with odd SSC-GL1 relationship
row.mean_SSCA2GL1A_overlap = mean(fcsdat.("SSC-A")(itemp) ./ fcsdat.("GL1-A")(itemp));
row.median_SSCA2GL1A_overlap = median(fcsdat.("SSC-A")(itemp) ./ fcsdat.("GL1-A")(itemp));
row.mean_SSCH2GL1H_overlap = mean(fcsdat.("SSC-H")(itemp) ./ fcsdat.("GL1-H")(itemp));
row.median_SSCH2GL1H_overlap = median(fcsdat.("SSC-H")(itemp) ./ fcsdat.("GL1-H")(itemp));
row.mean_SSCA2SSCH_overlap = mean(fcsdat.("SSC-A")(itemp) ./ fcsdat.("SSC-H")(itemp));
row.median_SSCA2SSCH_overlap = median(fcsdat.("SSC-A")(itemp) ./ fcsdat.("SSC-H")(itemp));
row.mean_GL1A2GL1H_overlap = mean(fcsdat.("GL1-A")(itemp) ./ fcsdat.("GL1-H")(itemp));
row.median_GL1A2GL1H_overlap = median(fcsdat.("GL1-A")(itemp) ./ fcsdat.("GL1-H")(itemp));
row.n_overlap = numel(itemp);
volume_cubic_micron = NaN(size(sscmerge_bdnorm));
ii = sscmerge_bdnorm > 0;
[volume_cubic_micron(ii), vol_func_string] = Attune_SC2vol(sscmerge_bdnorm(ii), strcat('SSC-', V.SSCDIM), 0);
if V.SSCDIM == 'A'
    ii = fcsdat.(sscstr) < SSCAmin; %for small signals use SSC-H instead (A goes negative)
    volume_cubic_micron(ii) = Attune_SC2vol(fcsdat.("SSC-H")(ii) ./ beadSSCmean.mean_SSCH_1micron(bdrow), strcat('SSC-H'), 0);
end
c.vol_notes = {strcat('calibrated: ', string(datetime())); strcat(' using SSC-', V.SSCDIM); vol_func_string; ...
    sprintf('SSC/GL1 merge window %g-%g x 1 um bead (GL1)', GL1_merge_bead)};
c.ssc_merge_bdnorm = sscmerge_bdnorm; c.beadSSCmean = beadSSCmean; c.beadGL1mean = beadGL1mean;
c.file_hv.SSC = row.SSC_hv; c.file_hv.GL1 = row.GL1_hv;
c.volume_cubic_microns = volume_cubic_micron;
end


function [c, row] = volume_488_2019(c, fcsdat, fcshdr, V, filetime)
% Volumes for cruises before the Aug 2019 layout change (Nov 2026,
% classifier_dev/pt_transfer/README.md):
%   - SSC-H relative to the 1 um bead at the file's time and voltage
%     (bead_reference, as the classifier; beadstat may be FCB runs or the PT
%     performance-test transfer, whose voltages differ from the samples')
%   - OD2setting 'None' (AR34B, RB1904, EN627): decompress large-particle
%     SSC first (ssc_decompress, at the file's SSC voltage). Beyond the
%     invertible range the volume is a lower bound and vol_flag is true.
%   - OD2setting 'SSC': SSC with OD2 is linear for these signals
%   - volume: Attune_SC2vol(..., 'SSC-H_488_2019') (March/April 2019 cultures)
if V.SSCDIM ~= 'H'
    warning('batch_mergedSSC2volume:SSCDIM', 'pre-Aug 2019 calibration is for SSC-H; using SSC-H')
end
ws = warning('off', 'bead_reference:noQC'); % no GL1 scatter channel before Aug 2019
j = strcmp({fcshdr.par.name}, 'SSC-H');
hv = fcshdr.par(j).hv;
ref = bead_reference(V.beadstat, datenum(filetime), V.cfg);
warning(ws);
beadpos = ref.bead_log.SSC + ref.k.SSC * log10(hv / ref.bead_hv.SSC);
raw = fcsdat(:, j);
s = log10(max(raw, 1)) - beadpos; % log10 SSC-H / 1 um bead
flag = false(size(s));
t = s;
if V.decomp
    [t, flag, ~, tmax] = ssc_decompress(s, hv);
    t(flag) = tmax; % lower bound
end
volume_cubic_micron = NaN(size(s));
ii = raw > 0;
[volume_cubic_micron(ii), vol_func_string] = Attune_SC2vol(10.^t(ii), 'SSC-H_488_2019', 0);
row.SSC_hv = hv;
row.n_vol_flag = sum(flag);
c.vol_notes = {strcat('calibrated: ', string(datetime())); ' using SSC-H, pre-Aug 2019 488 nm calibration'; vol_func_string; ...
    sprintf('OD2 %s; decompression %s', V.OD2setting, string(V.decomp)); ...
    sprintf('1 um bead (bead_reference): log10 SSC-H %.3f at %g V', ref.bead_log.SSC, ref.bead_hv.SSC); V.cruise_note};
c.ssc_bdnorm = 10.^s;
c.ssc_linear_bdnorm = 10.^t;
c.vol_flag = flag; % true: SSC beyond the decompression range, volume is a lower bound
c.file_hv = struct('SSC', hv);
c.volume_cubic_microns = volume_cubic_micron;
end


function c = preserved_file(c, fcsdat, fcshdr, V)
% preserved samples: same calibration as underway (bead-relative SSC-H merged
% with GL1-H, Attune_SC2vol 'SSC-H'), 1 um bead from the class file's bead
% reference (class_meta.ref, same-day lab beads) moved to the file's voltages
if ~isfield(c, 'class_meta') || ~isfield(c.class_meta, 'ref') || isempty(c.class_meta.ref)
    error('volume_file:noref', 'no bead reference')
end
ref = c.class_meta.ref;
W = V.GL1_merge_bead;
pn = {fcshdr.par.name};
is = strcmp(pn, 'SSC-H'); ig = strcmp(pn, 'GL1-H');
hvS = fcshdr.par(is).hv; hvG = fcshdr.par(ig).hv;
bS = 10.^(ref.bead_log.SSC + ref.k.SSC * log10(hvS / ref.bead_hv.SSC)); % 1 um bead SSC-H at hvS
bG = 10.^(ref.bead_log.GL1 + ref.k.GL1 * log10(hvG / ref.bead_hv.GL1));
ssc = fcsdat(:, is); ssc(ssc == 0) = 1;
s = ssc / bS; g = fcsdat(:, ig) / bG;
m = s;
hi = g > W(2); mid = g > W(1) & ~hi;
m(hi) = g(hi);
m(mid) = (s(mid) + g(mid)) / 2;
v = NaN(size(m)); ii = m > 0;
[v(ii), vf] = Attune_SC2vol(m(ii), 'SSC-H', 0);
c.volume_cubic_microns = v;
c.ssc_merge_bdnorm = m;
c.file_hv = struct('SSC', hvS, 'GL1', hvG);
c.vol_notes = {"calibrated: " + string(datetime('now')); "using SSC-H merged with GL1-H"; string(vf); ...
    sprintf('SSC/GL1 merge window %g-%g x 1 um bead (GL1)', W); ...
    sprintf('1 um bead from class_meta.ref (lab beads): log10 SSC-H %.3f at %g V, GL1-H %.3f at %g V', ...
    ref.bead_log.SSC, ref.bead_hv.SSC, ref.bead_log.GL1, ref.bead_hv.GL1); ...
    "preserved sample: live size calibration applied as is (no preservation correction)"};
end
