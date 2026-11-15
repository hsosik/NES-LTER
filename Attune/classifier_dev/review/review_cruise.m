function Q = review_cruise(cruise, opts)
% function Q = review_cruise(cruise, opts)
% Classification review for a whole cruise, without writing class files: classify
% every sample with assign_class_v2 (in parallel), write one review frame per sample
% and the HTML viewer (index.html) into the review folder. Open index.html in a browser.
% process_cruise can produce the same output during classification by calling
% review_sample in its loop and review_viewer after it.
%
% opts (all optional):
%   .outdir    review folder (default <cruise>\bead_calibrated\review)
%   .root      cruise_data root passed to attune_cruise_settings
%   .files     [] all; a number n: n files spread over the cruise; or a list of names
%   .old       true (default): compare with existing class files in bead_calibrated\class
%              (agreement metric, old counts) where they match the event count
%   .parallel  true (default when a pool can start): parfor over files
%   .decisions review decisions file the viewer reads and saves to (default
%              <cruise>ead_calibratedeview_decisions.csv; see review_decisions)
% Returns the QC structs (also saved as review_qc.mat in outdir); files that fail are
% listed in review_errors.txt.
% e.g. review_cruise('EN727')
%
% Nov 2026, classification review

if nargin < 2, opts = struct(); end
def = struct('outdir', '', 'root', '', 'files', [], 'old', true, 'parallel', [], 'decisions', '');
f = fieldnames(def);
for i = 1:numel(f), if ~isfield(opts, f{i}), opts.(f{i}) = def.(f{i}); end, end
here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here), fileparts(fileparts(here)))

s = attune_cruise_settings(cruise, opts.root);
c = s.config;
bp = fullfile(s.basepath, 'bead_calibrated');
outdir = opts.outdir; if isempty(outdir), outdir = fullfile(bp, 'review'); end
decfile = opts.decisions; if isempty(decfile), decfile = fullfile(bp, 'review_decisions.csv'); end
build = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')); % this classification run
if ~exist(fullfile(outdir, 'frames'), 'dir'), mkdir(fullfile(outdir, 'frames')); end

load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
Fi = FCSfileinfo(~contains(FCSfileinfo.fcslist, s.exclude), :);
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
Fi = sortrows(Fi, 'matdate_start');
if iscell(opts.files) || isstring(opts.files)
    Fi = Fi(ismember(Fi.fcslist, cellstr(opts.files)), :);
elseif ~isempty(opts.files) && opts.files < height(Fi)
    Fi = Fi(unique(round(linspace(1, height(Fi), opts.files))), :);
end
n = height(Fi);
fprintf('review_cruise %s: %d files -> %s\n', s.cruise, n, outdir)

bs = load_beadstat(bp);
cfg = classifier_config();
cfg.gl1_is_od2_scatter = strcmp(c.OD2, 'GL1');
cfg.cruise = s.cruise;
R = review_config();
fcsfile = fullfile(s.basepath, 'FCS', Fi.fcslist);
classfile = fullfile(bp, 'class', regexprep(Fi.fcslist, '\.fcs$', '.mat'));
names = Fi.fcslist; t = Fi.matdate_start; pro = Fi.pro_measured; use_old = opts.old;

use_par = opts.parallel;
if isempty(use_par), use_par = license('test', 'Distrib_Computing_Toolbox') && n > 4; end
if use_par && isempty(gcp('nocreate')), parpool; end
M = 0; if use_par, M = Inf; end
Q = cell(n, 1); err = strings(n, 1);
t0 = tic;
parfor (k = 1:n, M)
    try
        warning('off', 'stats:statrobustfit:IterationLimit')
        [d, h] = fca_readfcs(fcsfile{k});
        ref = bead_reference(bs, t(k), cfg);
        [cls, info] = assign_class_v2(d, h, ref, cfg, pro(k));
        old = [];
        if use_old && exist(classfile{k}, 'file')
            o = load(classfile{k}, 'class'); old = o.class(:, 1);
            if numel(old) ~= numel(cls), old = []; end
        end
        Q{k} = review_sample(d, h, cls, info, old, outdir, names{k}, t(k), pro(k), R);
    catch ME
        err(k) = string(names{k}) + ": " + string(ME.message);
    end
end
fprintf('review_cruise: %d files in %.1f min, %d failed\n', n, toc(t0) / 60, sum(err ~= ""))
ok = ~cellfun(@isempty, Q);
Q = Q(ok);
save(fullfile(outdir, 'review_qc.mat'), 'Q', 'err', 'build', 'decfile')
if any(err ~= "")
    writelines(err(err ~= ""), fullfile(outdir, 'review_errors.txt'));
end
review_viewer(outdir, Q, R, s.cruise, decfile, build);
fprintf('viewer: %s\n', fullfile(outdir, 'index.html'))
end
