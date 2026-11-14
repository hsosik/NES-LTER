function [beadstat, meta] = load_beadstat(bp)
% function [beadstat, meta] = load_beadstat(bp)
% Bead statistics written by process_beads_PT_adjust_3 (Nov 2026 redesign).
% bp - the cruise's bead_calibrated folder (or the cruise folder)
%
% beadstat - one row per bead run: filename, time, hv, and per population
%            P = 05, 1, 6: nP, QCP (true = good), whyP, peakP, medP, rsdP
%            (see measure_beads)
% meta     - cruise configuration, bead signature, processing notes
%
% Nov 2026, bead stats redesign

fn = fullfile(bp, 'beadstat.mat');
if ~exist(fn, 'file'), fn = fullfile(bp, 'bead_calibrated', 'beadstat.mat'); end
if ~exist(fn, 'file')
    error('load_beadstat:missing', 'no beadstat.mat in %s: run process_beads_PT_adjust_3', bp)
end
s = load(fn, 'beadstat', 'meta');
beadstat = s.beadstat;
meta = s.meta;
end
