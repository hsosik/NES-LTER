function [beadstat, meta] = load_beadstat(bp, use_decisions)
% function [beadstat, meta] = load_beadstat(bp, use_decisions)
% Bead statistics written by process_beads_PT_adjust_3 (Nov 2026 redesign).
% bp - the cruise's bead_calibrated folder (or the cruise folder)
% use_decisions - true (default): apply the bead review (bead_review_decisions.csv
%            beside beadstat.mat, written by the bead_review page): runs marked
%            exclude get QC05/QC1/QC6 false and "excluded in bead review" in whyP,
%            so bead_reference (classification) and the volume step skip them.
%            false: the bead stats as measured (bead_review itself); the marks
%            are still listed in meta.excluded
%
% beadstat - one row per bead run: filename, time, hv, and per population
%            P = 05, 1, 6: nP, QCP (true = good), whyP, peakP, medP, rsdP
%            (see measure_beads); excluded (true = excluded by the bead review;
%            all false when use_decisions is false)
% meta     - cruise configuration, bead signature, processing notes;
%            excluded (file names marked exclude), decisions (the file read)
%
% Nov 2026, bead stats redesign; bead review exclusions

if nargin < 2 || isempty(use_decisions), use_decisions = true; end
fn = fullfile(bp, 'beadstat.mat');
if ~exist(fn, 'file'), fn = fullfile(bp, 'bead_calibrated', 'beadstat.mat'); end
if ~exist(fn, 'file')
    error('load_beadstat:missing', 'no beadstat.mat in %s: run process_beads_PT_adjust_3', bp)
end
s = load(fn, 'beadstat', 'meta');
beadstat = s.beadstat;
meta = s.meta;

% bead review decisions
dec = fullfile(fileparts(fn), 'bead_review_decisions.csv');
meta.decisions = string(dec);
meta.excluded = strings(0, 1);
beadstat.excluded = false(height(beadstat), 1);
if exist(dec, 'file')
    if ~exist('review_decisions', 'file'), addpath(fullfile(fileparts(mfilename('fullpath')), 'review')), end
    T = review_decisions(dec);
    [~, f1, e1] = fileparts(cellstr(beadstat.filename)); % compare names without folders
    [~, f2, e2] = fileparts(T.file(strcmp(T.status, 'exclude')));
    x = ismember(strcat(f1, e1), strcat(f2, e2));
    meta.excluded = string(beadstat.filename(x));
    if use_decisions, beadstat.excluded = x(:); end
end
if any(beadstat.excluded)
    x = beadstat.excluded;
    for P = {'05' '1' '6'}
        beadstat.(['QC' P{1}])(x) = false;
        w = beadstat.(['why' P{1}]);
        w(x) = strtrim(w(x) + " excluded in bead review");
        beadstat.(['why' P{1}]) = w;
    end
    fprintf('load_beadstat: %d bead run(s) excluded in the bead review: %s\n', sum(x), strjoin(meta.excluded, ', '))
end
end
