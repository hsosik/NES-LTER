function S = generate_beadstats(cruises)
% function S = generate_beadstats(cruises)
% Batch bead processing: run process_beads_PT_adjust_3 for every cruise with
% FCB bead-check files (or the ones listed) and write each cruise's
% bead_calibrated\beadstat.mat, overwriting any earlier one (older
% beadstat_2026/2021 files are left in place). Use it after a change to the
% bead method; for a single new cruise, wrapper step 4 does the same.
% Returns a per-cruise summary.
%
% cruises - optional cellstr of cruise folder names (default: all)
%
% Nov 2026, bead stats redesign

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
warning('off', 'stats:robustcov:SingularCovariance')
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
if nargin < 1 || isempty(cruises)
    D = dir([root '20*_*']); D = D([D.isdir]); cruises = {D.name};
end
beadfiles2include = {'FCB_bead' 'fcb_bead' 'FCB_beadcheck'};
S = table;
for i = 1:numel(cruises)
    bp = fullfile(root, cruises{i});
    F = dir(fullfile(bp, 'FCS', '*.fcs'));
    if ~any(contains({F.name}, beadfiles2include)), continue, end
    conf = attune_config(extractAfter(cruises{i}, '_'));
    if ~conf.fcb_beads, continue, end
    FCSfileinfo = [];
    fi = fullfile(bp, 'bead_calibrated', 'FCSfileinfo.mat');
    if exist(fi, 'file'), s = load(fi, 'FCSfileinfo'); FCSfileinfo = s.FCSfileinfo; end
    if ~exist(fullfile(bp, 'bead_calibrated'), 'dir'), mkdir(fullfile(bp, 'bead_calibrated')), end
    try
        B = process_beads_PT_adjust_3(bp, FCSfileinfo, beadfiles2include, 'FCB', conf.OD2);
        S = [S; table(string(cruises{i}), height(B), sum(B.QC05), sum(B.QC1), sum(B.QC6), ...
            'VariableNames', {'cruise' 'runs' 'QC05' 'QC1' 'QC6'})]; %#ok<AGROW>
    catch err
        fprintf('%s FAILED: %s\n', cruises{i}, err.message)
    end
end
disp(S)
end
