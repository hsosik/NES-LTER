function V = volume_setup(p, first_time)
% function V = volume_setup(p, first_time)
% Settings and bead positions for cell volumes from SSC, shared by the
% process_cruise volume step (batch_mergedSSC2volume, preserved_volume) and the
% classify step, which computes volumes as it classifies when both steps are run
% (each FCS file read once). Per-file volumes: volume_file(c, fcsdat, fcshdr, V, t).
%
% p           process_cruise settings (outpath, SSCDIM, OD2setting, sample)
% first_time  datetime of the cruise's first sample file (underway, Aug 2019 on:
%             bead runs from that day on are used)
%
% V.mode  'merge'     underway, Aug 2019 on: SSC merged with GL1, 1 um bead means
%                     per SSC voltage over the cruise (beadSSCmean, beadGL1mean)
%         '488_2019'  underway before Aug 2019: SSC-H only, bead_reference
%         'preserved' preserved samples: bead reference saved with the classes
%                     (class_meta.ref, same-day lab beads)
%
% Nov 2026, split from batch_mergedSSC2volume

addpath(fullfile(fileparts(mfilename('fullpath')), 'classifier_dev'))
cfg = classifier_config();
V = struct('mode', '', 'cfg', cfg, 'GL1_merge_bead', cfg.gl1_merge_bead, 'SSCDIM', 'H', 'OD2setting', '');
if isfield(p, 'SSCDIM'), V.SSCDIM = p.SSCDIM; end
if isfield(p, 'OD2setting'), V.OD2setting = p.OD2setting; end
if isfield(p, 'sample') && p.sample == "preserved"
    V.mode = 'preserved';
    return
end
beadstat = load_beadstat(p.outpath); % from process_beads_PT_adjust_3
if ~strcmp(V.OD2setting, 'GL1')
    % cruises before Aug 2019: one scatter channel (SSC 488/10; GL1-H is PE
    % fluorescence in that layout, not side scatter), so no merge
    V.mode = '488_2019';
    V.beadstat = beadstat;
    V.decomp = strcmp(V.OD2setting, 'None');
    V.cruise_note = "";
    if contains(p.basepath, 'EN627')
        V.cruise_note = "EN627: scatter optics degraded (Feb 2019); SSC-based volumes uncertain, flag in analyses";
    end
    return
end
V.mode = 'merge';
%1 micron bead peaks (KDE mode) for SSC and GL1, QC'd runs starting on cruise day and later
ind = beadstat.QC1 & beadstat.time > dateshift(first_time, 'start', 'day');
B1 = beadstat.peak1(ind, :);
beadT = table(B1.SSC_A, B1.SSC_H, beadstat.hv.SSC(ind), 'VariableNames', {'SSCA_1micron' 'SSCH_1micron' 'SSC_hv'});
V.beadSSCmean = groupsummary(beadT, 'SSC_hv', 'mean'); %mean grouped by HV
beadT = table(B1.GL1_A, B1.GL1_H, beadstat.hv.GL1(ind), 'VariableNames', {'GL1A_1micron' 'GL1H_1micron' 'GL1_hv'});
V.beadGL1mean = groupsummary(beadT, 'GL1_hv', 'mean');
end
