function T = review_prototype(cruise, nfiles, first, outdir, use_old)
% function T = review_prototype(cruise, nfiles, first, outdir, use_old)
% Prototype of the cruise review: classify nfiles consecutive samples (in time order,
% starting at the first-th) with assign_class_v2, write one review frame per sample
% and the HTML viewer into outdir. Class files are not written.
% use_old (default true): add old-class and disagreement rows from bead_calibrated\class.
% Classifications are cached in outdir\cache so frames can be redrawn quickly.
% Returns a table of per-file timings (read, classify, frame) in seconds.
% e.g. review_prototype('EN727', 20, 1000, 'C:\...\output\review_EN727')
%
% Nov 2026, classification review prototype

if nargin < 5 || isempty(use_old), use_old = true; end
here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here), fileparts(fileparts(here)))
if review_is_cruise_dir(outdir) % a cruise's own review must match its class files
    error('review_prototype:cruiseReview', 'Not into a cruise''s own review folder; build that with process_cruise(cruise, Steps="classify").');
end
s = attune_cruise_settings(cruise);
c = s.config;
bp = fullfile(s.basepath, 'bead_calibrated');
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
Fi = FCSfileinfo(~contains(FCSfileinfo.fcslist, s.exclude), :);
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
Fi = sortrows(Fi, 'matdate_start');
Fi = Fi(first:min(height(Fi), first + nfiles - 1), :);
n = height(Fi);

bs = load_beadstat(bp);
cfg = classifier_config();
cfg.gl1_is_od2_scatter = strcmp(c.OD2, 'GL1');
cfg.cruise = s.cruise;
R = review_config();
for sd = {'frames' 'cache'}, if ~exist(fullfile(outdir, sd{1}), 'dir'), mkdir(fullfile(outdir, sd{1})); end, end

Q = cell(n, 1); tm = nan(n, 3);
for k = 1:n
    fn = Fi.fcslist{k};
    t0 = tic;
    [d, h] = fca_readfcs(fullfile(s.basepath, 'FCS', fn));
    tm(k, 1) = toc(t0); t0 = tic;
    cache = fullfile(outdir, 'cache', regexprep(fn, '\.fcs$', '.mat'));
    if exist(cache, 'file') %classification kept from an earlier run (delete cache\ to redo)
        load(cache, 'cls', 'info')
    else
        ref = bead_reference(bs, Fi.matdate_start(k), cfg);
        [cls, info] = assign_class_v2(d, h, ref, cfg, Fi.pro_measured(k));
        save(cache, 'cls', 'info')
    end
    tm(k, 2) = toc(t0); t0 = tic;
    old = [];
    cf = fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat'));
    if use_old && exist(cf, 'file')
        o = load(cf, 'class'); old = o.class(:, 1);
        if numel(old) ~= numel(cls), old = []; end
    end
    Q{k} = review_sample(d, h, cls, info, old, outdir, fn, Fi.matdate_start(k), Fi.pro_measured(k), R);
    tm(k, 3) = toc(t0);
    fprintf('%2d/%d %s  read %.2f  classify %.2f  frame %.2f s\n', k, n, fn, tm(k, :))
end
review_viewer(outdir, Q, R, s.cruise);
T = array2table(tm, 'VariableNames', {'read' 'classify' 'frame'});
T.file = Fi.fcslist;
fprintf('median s/file: read %.2f, classify %.2f, frame+QC %.2f\n', median(tm))
end
