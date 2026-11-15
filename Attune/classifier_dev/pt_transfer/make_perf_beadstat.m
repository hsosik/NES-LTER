function info = make_perf_beadstat(cruisedir, opts)
% function info = make_perf_beadstat(cruisedir, opts)
% Build and install the bead table (beadstat.mat) for a cruise WITHOUT FCB
% bead runs, from the Attune performance tests (PT transfer). One call per
% cruise; this is the entry point for the processing wrapper's bead step when
% attune_config(cruise).fcb_beads is false.
%
% Steps:
%   1. performance-test table (pt_perf_table): read from
%      classifier_dev/output/pt_perf.mat, or rebuilt from
%      \\sosiknas1\Lab_data\Attune\performance_test_fcs when missing or
%      opts.rebuild
%   2. FCB-vs-PT offsets (validate_perf_transfer, from FCB bead runs paired
%      with performance tests on AR34B, RB1904, TN368, EN644): read from
%      classifier_dev/output/perf_transfer_validation.mat, or rebuilt likewise
%   3. perf_beadstats for this cruise -> output/pt_transfer/<cruise>_beadstat.mat
%   4. install as <cruisedir>/bead_calibrated/beadstat.mat (opts.install):
%      - no table there: write it
%      - an earlier PT-transfer table: back it up, then replace
%      - any other table (FCB pipeline, old PT conversion, ...): refuse
%        unless opts.force; with force, back it up, then replace
%      Backups: beadstat_backup_<yyyymmddTHHMMSS>.mat in the same folder.
%
% cruisedir - cruise folder, full path or name under cruise_data
%             (e.g. '20180131_EN608')
% opts      - struct, all optional:
%   .install  (true)  copy into bead_calibrated
%   .force    (false) replace a table that is not a PT transfer; also allow a
%                     cruise that attune_config says has FCB beads
%   .rebuild  (false) rebuild the performance-test table and offsets
%   .bead_dir (cruisedir/bead_calibrated) install folder (for testing)
%
% info - struct: cruise, ntests, built (file in output/pt_transfer),
%        installed (path or ""), backup (path or ""), previous_source
%
% Cruises this applies to (attune_config fcb_beads false): EN608, AR28B,
% AR29, EN617, AR31A (installed Nov 2026) and EN627 (its existing PT
% conversion is kept unless forced; its SSC is unreliable anyway).
%
% Nov 2026, PT transfer for early cruises

if nargin < 2, opts = struct(); end
def = struct('install', true, 'force', false, 'rebuild', false, 'bead_dir', '');
f = fieldnames(def);
for i = 1:numel(f), if ~isfield(opts, f{i}), opts.(f{i}) = def.(f{i}); end, end

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(here, dev, fileparts(dev))
root = ['\' '\sosiknas1\lab_data\Attune\cruise_data\'];
if ~contains(cruisedir, filesep), cruisedir = fullfile(root, cruisedir); end
cruisedir = regexprep(cruisedir, '[\\/]+$', '');
[~, cdir] = fileparts(cruisedir);
cruise = extractAfter(cdir, '_');
if isempty(cruise), error('make_perf_beadstat:name', 'cruise folder should be <yyyymmdd>_<CRUISE>: %s', cdir); end

conf = attune_config(cruise);
if conf.fcb_beads && ~opts.force
    error('make_perf_beadstat:fcb', ['%s has FCB bead runs (attune_config); use the FCB bead ' ...
        'pipeline, or opts.force = true to build a PT-transfer table anyway'], cruise)
end

% ---- 1-2: performance tests and offsets
fR = fullfile(dev, 'output', 'pt_perf.mat');
fV = fullfile(dev, 'output', 'perf_transfer_validation.mat');
if opts.rebuild || ~isfile(fR)
    fprintf('building performance-test table ...\n')
    R = pt_perf_table([], fR); %#ok<NASGU>
end
s = load(fR, 'R'); R = s.R;
if opts.rebuild || ~isfile(fV)
    fprintf('building FCB-vs-PT offsets ...\n')
    V = validate_perf_transfer([], [], R);
    save(fV, 'V')
end
s = load(fV, 'V'); V = s.V;

% ---- 3: table for this cruise
outdir = fullfile(dev, 'output', 'pt_transfer');
out = perf_beadstats({cdir}, V, R, outdir);
built = fullfile(outdir, [cruise '_beadstat.mat']);
k = matlab.lang.makeValidName(cruise);
info = struct('cruise', string(cruise), 'ntests', 0, 'built', "", 'installed', "", ...
    'backup', "", 'previous_source', "");
if ~isfield(out, k) || height(out.(k)) == 0
    error('make_perf_beadstat:none', '%s: no good performance tests in the cruise window and optical state (%s)', cruise, conf.OD2)
end
info.ntests = height(out.(k));
info.built = string(built);
if ~opts.install, return, end

% ---- 4: install
bdir = opts.bead_dir;
if isempty(bdir), bdir = fullfile(cruisedir, 'bead_calibrated'); end
if ~exist(bdir, 'dir'), mkdir(bdir); end
dest = fullfile(bdir, 'beadstat.mat');
if isfile(dest)
    prev = whos('-file', dest);
    src = "unknown (no meta)";
    if ismember('meta', {prev.name})
        m = load(dest, 'meta');
        if isfield(m.meta, 'source'), src = string(m.meta.source); end
    end
    info.previous_source = src;
    if src ~= "PT performance test transfer" && ~opts.force
        error('make_perf_beadstat:exists', ['%s already has a bead table from "%s"; not replaced ' ...
            '(opts.force = true replaces it after a backup). New table: %s'], dest, src, built)
    end
    bk = fullfile(bdir, ['beadstat_backup_' char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmss')) '.mat']);
    copyfile(dest, bk);
    info.backup = string(bk);
end
copyfile(built, dest);
info.installed = string(dest);
fprintf('%s: %d performance tests -> %s', cruise, info.ntests, dest)
if info.backup ~= "", fprintf(' (previous "%s" table backed up as %s)', info.previous_source, info.backup), end
fprintf('\n')
end
