function sig = bead_signature()
% function sig = bead_signature()
% Expected log10 peak heights of the three FCB bead populations at the
% standard bead-check voltages (sig.hvref), and the channels used to find each
% population. Positions are the medians over 2020-2026 cruises of the bead
% medians in beadstat_2026 (output/bead_runs.mat), already corrected to sig.hvref.
% (All those runs are at sig.hvref, so the positions do not depend on k.)
% At other voltages: expected = pos + k * log10(hv / hvref), k = cfg.k from
% classifier_config.
%
% Spectral signatures that make the populations easy to separate:
%   0.5 um - bright on GL2 (4.3), dim on BL1 (1.9)
%   1 um   - bright on BL1 (4.3), dim on GL3 (1.9)
%   6 um   - saturates BL3, GL1 (OD2 scatter) ~4.3, GL3 ~3.4
%
% Nov 2026, bead stats redesign prototype

sig.chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
sig.hvref = struct('FSC', 280, 'SSC', 220, 'BL1', 280, 'BL2', 360, 'BL3', 340, ...
    'GL1', 220, 'GL2', 360, 'GL3', 320, 'GL4', 260);
% voltage exponents: one source, shared with the classifier
cfg = classifier_config();
sig.k = cfg.k;
sig.saturation = log10(1048575) - 0.02; % peaks above this are clipped
sig.floor = 2.8; % a bead peak expected below this is too close to noise to gate on

% gate: gating channels; alt: used instead when a gating channel would be
% saturated, below sig.floor (e.g. 0.5 um on SSC behind the OD2 filter), or not
% installed at this file's settings
p = struct('name', {}, 'size_um', {}, 'gate', {}, 'alt', {}, 'pos', {}, 'max_rsd', {}, 'shift', {}, 'limit', {}, 'sigch', {});
%                         FSC   SSC   BL1   BL2   BL3   GL1   GL2   GL3   GL4
p(1).name = 'b05'; p(1).size_um = 0.5; p(1).gate = {'SSC' 'GL2'}; p(1).alt = {'BL2' 'GL2'};
p(1).pos = mk(sig.chans, [3.26  3.86  1.90  3.81  2.99  2.20  4.30  3.24  2.59]);
p(1).max_rsd = [0.09 0.05]; p(1).limit = {'BL1' 0.6}; p(1).sigch = {'SSC' 'BL2' 'GL2' 'GL3'};
p(2).name = 'b1';  p(2).size_um = 1;   p(2).gate = {'SSC' 'BL1'};
p(2).pos = mk(sig.chans, [4.21  4.70  4.27  4.59  3.13  3.00  3.20  1.90  1.72]);
p(2).max_rsd = [0.06 0.06]; p(2).limit = {'GL3' 0.6}; p(2).sigch = {'FSC' 'SSC' 'BL1' 'BL2' 'GL1'};
p(3).name = 'b6';  p(3).size_um = 6;   p(3).gate = {'GL1' 'GL3'}; p(3).alt = {'FSC' 'GL3'};
p(3).pos = mk(sig.chans, [5.45  5.63  1.90  2.70  6.02  4.32  2.55  3.39  5.91]);
p(3).max_rsd = [0.10 0.06]; p(3).sigch = {'FSC' 'SSC' 'GL1' 'GL3'}; p(3).limit = {'BL1' 1.0}; % 6 um BL1 varies 1.8-2.5 between bead lots
% max_rsd: limit on the robust SD (log10) on each gating channel, from the
% 95-99th percentile over all runs. shift: struct of per-channel offsets of
% this cruise's beads from pos (empty = none), learned in a first pass.
% limit: {channel, d} drops events more than d above the expected position on
% that channel: aggregates of a 1 um bead with a 0.5 um bead match single
% 1 um beads on SSC and BL1 but carry the 0.5 um bead's GL3 (~1.3 log above).
% Likewise 6 um beads stuck to a 1 um bead carry its BL1 (up to 68% of the
% 6 um gate on EN649; GL3, FSC and SSC are unchanged, but the gate stays clean).
for i = 1:3, p(i).shift = struct(); end
sig.pop = p;

% search window half-width around the expected position (log10), in each
% gating channel; covers between-cruise drift (up to ~0.15) and hv model error
sig.window = 0.35;
% QC limits per population
sig.qc.min_n = 100;          % events in the final gate
sig.qc.max_offset = 0.25;    % |found - expected| on each gating channel
% whole-signature check: on each of the population's signature channels
% (pop.sigch: bright and stable ones; dim channels such as 1 um on GL2 or
% 6 um on BL1/BL2 drift by up to 0.4 within a cruise) that is on scale at
% this file's voltages (sig.sig_range), its peak must be within sig.sig_tol
% of the expected position. This is what rejects a wrong cluster.
sig.sig_range = [2.3 5.8];
sig.sig_tol = 0.4;
% pass 1 (learn the cruise's per-channel offsets): wide window, no position
% or signature checks. pass 2 (with pop.shift set): narrow window, tight checks
sig.pass1 = struct('window', 0.6, 'max_offset', Inf, 'sig_tol', Inf);
sig.pass2 = struct('window', 0.25, 'max_offset', 0.2, 'sig_tol', 0.2);
% (max_offset stays loose: a real instrument change within a cruise moves
% every population together, e.g. SSC +0.12 in the first EN644 runs)
end

function s = mk(chans, v)
s = struct();
for i = 1:numel(chans), s.(chans{i}) = v(i); end
end
