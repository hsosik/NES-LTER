function p = process_cruise(cruise, opt)
% function p = process_cruise(cruise, Name=Value, ...)
% Process one cruise of Attune data end to end, or any subset of the steps.
% Replaces the copied process_wrapper_* scripts: every per-cruise setting
% comes from attune_cruise_settings (attune_cruise_settings.csv plus
% classifier_dev/attune_config), and every cruise uses the same classifier,
% assign_class_v2.
%
% e.g. process_cruise("EN727")                                   % core steps
%      process_cruise("EN727", Steps=["classify" "volume" "table"])
%      process_cruise("EN727", Steps="classify", Files=8, OutPath="C:\scratch\EN727")
%      process_cruise("EN727", Steps="classify", NewOnly=true)   % at sea: new files only
%
% Steps (run in this order, however they are listed)
%   files     FCS_DateTimeList -> bead_calibrated\FCSfileinfo.mat
%             An existing FCSfileinfo.mat is extended with new files (Force: rebuild);
%             pro_measured comes from attune_pro_measured (trigger settings and voltages).
%   beads     bead table bead_calibrated\beadstat.mat
%             - FCB bead cruises (attune_config fcb_beads): process_beads_PT_adjust_3.
%               It will not overwrite a table made another way (meta.source set)
%               unless Force.
%             - no FCB beads (EN608 to EN627): pt_transfer\make_perf_beadstat.
%               EN627 keeps its PT conversion unless Force.
%             A newly written table (FCB or PT transfer) also gets the bead review
%             (classifier_dev\review\bead_review, page in bead_calibrated\bead_review\;
%             for PT transfer it shows the performance tests) unless Review is false.
%             Runs or tests marked exclude there (bead_review_decisions.csv) are
%             skipped by classify and volume (load_beadstat); after changing
%             exclusions rerun from classify (without NewOnly).
%   classify  assign_class_v2 on every sample file -> bead_calibrated\class\<file>.mat
%             holding class, syn_type, flag, notes, class_meta (see below), and the
%             classification review in bead_calibrated\review\ (unless Review is
%             false): one cytogram frame per file, per-file QC and index.html, a viewer
%             that steps through the cruise, flags unusual files and saves review
%             decisions to bead_calibrated\review_decisions.csv (classifier_dev\review)
%   volume    batch_mergedSSC2volume: adds volume_cubic_microns etc. to the class files
%   table     generate_attune_table -> bead_calibrated\AttuneTable.mat
%   metadata  match_Attune_underway_LTER -> bead_calibrated\AttuneTable_uw_match.mat
%   movie     attune_lter_moviemaker (optional)
%   voldists  get_cruise_voldists_fromEDItable2, Plot_Voldists_function (optional)
% "core" = files beads classify volume table metadata; "all" adds movie and voldists.
%
% Options (Name=Value, all optional)
%   Steps     step labels, "core" (default) or "all"
%   Sample    "underway" (default) or "preserved": discrete samples from <cruise>\preserved\FCS,
%             output in <cruise>\preserved\bead_calibrated\. Built so far: files (adds
%             cast, niskin, runtype: preserved_fileinfo), beads (lab bead runs:
%             preserved_bead_step, and the bead review page in bead_calibrated\bead_review\
%             unless Review is false; exclusions apply to this cruise), metadata (CTD bottles / event log:
%             preserved_metadata, writes BottleMetadata.mat), classify (by run type:
%             assign_class_discrete for CHL/PE/pro runs, assign_class_hbac for SYBR runs;
%             same-day lab beads, discrete_bead_reference; per-event time_gate mask),
%             volume (preserved_volume: same calibration, beads from class_meta.ref),
%             table (preserved_table: per-Niskin DiscreteTable.mat, each group from one
%             run type; joined with BottleMetadata by the metadata step).
%   Files     classify only these: a number n (n files spread through the cruise) or a
%             list of FCS names (default: all sample files)
%   NewOnly   true: classify only files without a class file (default false)
%   Force     true: rebuild what steps would otherwise extend or protect: a fresh
%             FCSfileinfo, a bead table replacing one not made from FCB beads
%             (default false)
%   OutPath   output folder instead of <cruise>\bead_calibrated\ (testing). Steps after
%             files read FCSfileinfo.mat and beadstat.mat from here, so copy them in first.
%   Parallel  true/false: parfor over files in classify (default: true when a pool can start)
%   Review    true (default): build the classification and bead reviews. A run on some
%             files (Files, NewOnly) updates their frames and keeps the rest.
%   Root      folder of cruise folders (default \\sosiknas1\Lab_data\Attune\cruise_data\)
% Options kept in a struct s can be passed as c = namedargs2cell(s); process_cruise(cruise, c{:})
%
% Class file variables written by classify (volume loads and re-saves the file, keeping them)
%   class      per event, double: 0 unassigned, 1 Euk, 2 Syn, 3 dim crypto, 4 bright crypto,
%              5 Euk with Syn, 7 Pro, 8 Syn on particle
%   syn_type   per event, uint8: 1/2/3 Syn type A/B/C (PUB:PEB low/mid/high),
%              9 Syn-bearing but spectrum not usable, 0 not Syn-bearing (classes 2, 5, 8 only)
%   bact_type  SYBR (hbac) preserved runs only: per event, uint8, 1 LNA, 2 HNA bacteria,
%              0 not bacteria (class 9 only; method and cut in class_meta.hbac.hna)
%   flag       per event, uint8, only on class 0: 1 Syn-pigment particle with low chl
%              (degraded Syn/Euk, grazer, detritus), 2 below the Euk chl gap, 3 trigger-
%              noise pile at the chl edge (pre-2019), 4 Pro-like event in a file without
%              Pro settings (Pro lower bound)
%   notes      string: (1) class codes, (2) classifier name, (3) date, (4) syn_type codes,
%              (5) flag codes
%   class_meta struct: classifier, date, code version, cruise configuration, bead
%              reference used, classifier settings, mixture components,
%              pro_near_trigger_frac (Pro-setting files: share of Pro within 0.1 log of
%              the SSC or chl trigger edge)
%
% Each run appends its settings to bead_calibrated\process_cruise_log.mat.
%
% Nov 2026, processing wrapper redesign (replaces process_wrapper_2024b);
% 19 Nov 2026: Name=Value options and step labels replace step numbers and an opts struct

arguments
    cruise {mustBeTextScalar}
    opt.Steps string = "core"
    opt.Sample (1,1) string {mustBeMember(opt.Sample, ["underway" "preserved"])} = "underway"
    opt.Files = []
    opt.NewOnly (1,1) logical = false
    opt.Force (1,1) logical = false
    opt.OutPath (1,1) string = ""
    opt.Parallel {mustBeScalarOrEmpty} = []
    opt.Review (1,1) logical = true
    opt.Root (1,1) string = ""
end
steps = step_list(opt.Steps);
pres = opt.Sample == "preserved";
if pres && any(ismember(steps, ["movie" "voldists"]))
    error('process_cruise:preserved', 'movie and voldists are not used for preserved samples')
end
if ~(isempty(opt.Files) || (isnumeric(opt.Files) && isscalar(opt.Files)) || iscellstr(opt.Files) || isstring(opt.Files))
    error('process_cruise:files', 'Files must be a number of files or a list of FCS names')
end
opts = struct('root', char(opt.Root), 'outpath', char(opt.OutPath), 'files', {opt.Files}, ...
    'newonly', opt.NewOnly, 'parallel', opt.Parallel, 'force', opt.Force, 'review', opt.Review);
do = @(name) any(steps == name);

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'preserved'), fullfile(here, 'classifier_dev'), fullfile(here, 'classifier_dev', 'pt_transfer'), ...
    fullfile(here, 'classifier_dev', 'review'), fullfile(here, 'classifier_dev', 'discrete'))

%% settings
s = attune_cruise_settings(cruise, opts.root);
c = s.config;
if pres
    % preserved samples are run in the lab after the cruise: current layout, OD2 on GL1
    c = attune_config('LAB');
end
p.cruise = s.cruise;
p.basepath = [s.basepath filesep];
p.fpath = fullfile(s.basepath, 'FCS', filesep);
if pres, p.fpath = fullfile(s.preserved_path, 'FCS', filesep); end
if isempty(opts.outpath) && pres
    p.outpath = fullfile(s.preserved_path, 'bead_calibrated', filesep);
elseif isempty(opts.outpath)
    p.outpath = fullfile(s.basepath, 'bead_calibrated', filesep);
else
    p.outpath = [regexprep(char(opts.outpath), '[\\/]+$', '') filesep];
end
p.classpath = fullfile(p.outpath, 'class', filesep);
p.OD2setting = c.OD2;      % read by batch_mergedSSC2volume and the movie maker
p.SSCDIM = 'H';            % batch_mergedSSC2volume: height channels
p.beadfiles2include = {'FCB_bead' 'fcb_bead'};
p.filetype2exclude = s.exclude;
if pres, p.filetype2exclude = s.exclude_preserved; end
p.uw_fullname = s.underway;
p.assign_class_function = 'assign_class_v2';
p.config = c;
p.steps = steps;          % step labels
p.sample = opt.Sample;
p.opts = opts;
if strcmp(c.layout, 'v2019'), p.moviechannels = 'early'; else, p.moviechannels = 'late'; end
p.framemaker = 'make_movieframe_density';
p.stepsize = 1;

if pres
    fprintf('%s preserved samples: %s\n', p.cruise, p.fpath)
else
    fprintf('%s: %s\n  OD2 on %s, layout %s, FCB beads %d\n', p.cruise, p.basepath, c.OD2, c.layout, c.fcb_beads)
end
if ~isempty(s.notes), fprintf('  note: %s\n', s.notes); end
if ~exist(p.classpath, 'dir'), mkdir(p.classpath); end
write_log(p)

%% files
if do("files")
    fi = [p.outpath 'FCSfileinfo.mat'];
    if exist(fi, 'file') && ~opts.force
        FCSfileinfo = FCS_DateTimeList(p.fpath, fi); % reads only new files; fills voltages for old rows
    else
        FCSfileinfo = FCS_DateTimeList(p.fpath);
    end
    if pres, FCSfileinfo = preserved_fileinfo(FCSfileinfo); end % cast, niskin, runtype...
    save(fi, 'FCSfileinfo')
    fprintf('files: %d FCS files listed\n', height(FCSfileinfo))
end

%% beads
if do("beads") && pres
    made = preserved_bead_step(p, opts.force); % lab bead runs; exclusions from this cruise's bead review
    if made && opts.review
        try
            bead_review(p.cruise, struct('root', opts.root, 'beaddir', p.outpath));
        catch ME
            warning('process_cruise:bead_review', 'bead review not built: %s', ME.message)
        end
    end
elseif do("beads")
    made = bead_step(p, c, opts.force);
    if made && opts.review
        try
            bead_review(p.cruise, struct('root', opts.root, 'beaddir', p.outpath));
        catch ME
            warning('process_cruise:bead_review', 'bead review not built: %s', ME.message)
        end
    end
end

%% classify
if do("classify")
    classify_step(p, c, opts)
end

%% volume
if do("volume") && pres
    use_par = opts.parallel; if isempty(use_par), use_par = license('test', 'Distrib_Computing_Toolbox'); end
    preserved_volume(p, use_par); % bead reference from class_meta (same-day lab beads)
elseif do("volume")
    batch_mergedSSC2volume(p)
end

%% table
if do("table") && pres
    preserved_table(p); % per-Niskin DiscreteTable.mat
elseif do("table")
    generate_attune_table(p.classpath, [p.outpath 'FCSfileinfo.mat'])
end

%% metadata (underway match; bottle match for preserved samples)
if do("metadata") && pres
    preserved_metadata(p, s);
elseif do("metadata")
    match_Attune_underway_LTER([p.outpath 'AttuneTable.mat'], p.uw_fullname);
end

%% movie
if do("movie")
    close all
    figure('WindowState', 'maximized')
    attune_lter_moviemaker(p.fpath, p.classpath, p.OD2setting, p.framemaker, p.moviechannels, p.stepsize)
end

%% volume distributions
if do("voldists")
    get_cruise_voldists_fromEDItable2(p.basepath)
    Plot_Voldists_function(p.basepath, p.outpath)
end
end


function made = bead_step(p, c, force)
% made: true when a bead table was written or installed (the bead review is built from it)
made = false;
fn = [p.outpath 'beadstat.mat'];
if c.fcb_beads
    if exist(fn, 'file')
        m = load(fn, 'meta');
        if isfield(m, 'meta') && isfield(m.meta, 'source') && ~force
            fprintf('beads: keeping %s (meta.source "%s"); Force=true to replace it with an FCB table\n', fn, m.meta.source)
            return
        end
    end
    load([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo')
    B = process_beads_PT_adjust_3(p.basepath, FCSfileinfo, p.beadfiles2include, 'FCB', c.OD2);
    made = height(B) > 0;
else
    try
        info = make_perf_beadstat(regexprep(p.basepath, '[\\/]+$', ''), struct('force', force, 'bead_dir', p.outpath));
        fprintf('beads: PT-transfer bead table from %d performance tests: %s\n', info.ntests, info.installed)
        made = strlength(string(info.installed)) > 0;
    catch ME
        if strcmp(ME.identifier, 'make_perf_beadstat:exists')
            m = load(fn, 'meta');
            fprintf('beads: keeping %s (meta.source "%s")\n', fn, m.meta.source)
        else
            rethrow(ME)
        end
    end
end
end


function classify_step(p, c, opts)
load([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo')
Fi = FCSfileinfo;
Fi = Fi(~contains(Fi.fcslist, p.filetype2exclude), :);
pres = isfield(p, 'sample') && p.sample == "preserved";
if pres, Fi = Fi(~strcmp(Fi.runtype, 'other'), :); end
Fi = sortrows(Fi, 'matdate_start');
if iscell(opts.files) || isstring(opts.files)
    Fi = Fi(ismember(Fi.fcslist, cellstr(opts.files)), :);
elseif ~isempty(opts.files) && opts.files < height(Fi)
    Fi = Fi(unique(round(linspace(1, height(Fi), opts.files))), :);
end
classfile = strcat(p.classpath, regexprep(Fi.fcslist, '\.fcs$', '.mat'));
if opts.newonly
    have = cellfun(@(x) exist(x, 'file') == 2, classfile);
    Fi = Fi(~have, :); classfile = classfile(~have);
end
if ~pres && ~ismember('trig_BL3_eq340', Fi.Properties.VariableNames) && height(Fi) > 0
    % FCSfileinfo from before Nov 2026: old pro_measured rule, no voltages
    fprintf('classify: FCSfileinfo predates the voltage-based pro_measured rule; recomputing it for these files (rerun the files step to save it)\n')
    Fi = attune_pro_measured(Fi, p.fpath);
end
n = height(Fi);
fprintf('classify: %d files\n', n)
if n == 0, return, end

[bs, bmeta] = load_beadstat(p.outpath); % skips runs marked exclude in the bead review
bexcl = bmeta.excluded;
if pres
    % preserved samples: lab runs (current layout, OD2 on GL1), preserved-sample settings
    cfg = discrete_config();
    cfg.gl1_is_od2_scatter = true;
    cfg.cruise = 'LAB'; %bead_reference: attune_config('LAB') = OD2 on GL1, current layout
else
    cfg = classifier_config();
    cfg.gl1_is_od2_scatter = strcmp(c.OD2, 'GL1');
    cfg.cruise = p.cruise; %bead_reference: OD2 placement (ref.OD2) and channel map from attune_config
end
[~, h] = system(['git -C "' fileparts(mfilename('fullpath')) '" rev-parse --short HEAD']);
[~, dirty] = system(['git -C "' fileparts(mfilename('fullpath')) '" status --porcelain -- classifier_dev']);
version = strtrim(h);
if ~isempty(strtrim(dirty)), version = [version ' + uncommitted classifier_dev changes']; end
src = "FCB beads";
if isfield(bmeta, 'source'), src = string(bmeta.source); end

codes = class_codes();
fcsfile = strcat(p.fpath, Fi.fcslist);
names = Fi.fcslist; t = Fi.matdate_start; pro = Fi.pro_measured;
nev = nan(n, 1); ncomp = nan(n, 1); counts = nan(n, 10); err = strings(n, 1);
tgf = 0.8; % preserved: count only the last 80% of acquisition time (as the hand gates)
if pres
    rt = Fi.runtype; pmeta = table2struct(Fi(:, {'cast' 'niskin' 'sample_label' 'repeat' 'sizefrac' 'vol_analyzed'}));
    codes.class = codes.class + "; 9 = heterotrophic bacteria (SYBR-stained runs; Pro there are 7)";
else
    rt = repmat({''}, n, 1); pmeta = repmat(struct(), n, 1);
end
% classification review: frame and QC per file in the loop, viewer after it
do_review = opts.review; revdir = fullfile(p.outpath, 'review'); R = []; Q = cell(n, 1);
if do_review
    if ~exist(fullfile(revdir, 'frames'), 'dir'), mkdir(fullfile(revdir, 'frames')); end
    R = review_config();
end
use_par = opts.parallel;
if isempty(use_par), use_par = license('test', 'Distrib_Computing_Toolbox') && n > 4; end
if use_par && isempty(gcp('nocreate')), parpool; end
M = 0; if use_par, M = Inf; end
parfor (k = 1:n, M)
    try
        [fcsdat, fcshdr] = fca_readfcs(fcsfile{k});
        if pres
            ref = discrete_bead_reference(bs, t(k), cfg);
            if isempty(ref), error('process_cruise:nobeads', 'no lab bead run on %s', datestr(t(k), 'dd-mmm-yyyy')), end
            if strcmp(rt{k}, 'hbac')
                classifier = "assign_class_hbac"; [cls, info] = assign_class_hbac(fcsdat, fcshdr, ref, cfg);
            else
                classifier = "assign_class_discrete"; [cls, info] = assign_class_discrete(fcsdat, fcshdr, ref, cfg, rt{k});
            end
        else
            ref = bead_reference(bs, t(k), cfg);
            classifier = "assign_class_v2"; [cls, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, pro(k));
        end
        o = struct();
        o.class = double(cls(:));
        o.syn_type = field_or(info, 'syn_type', zeros(numel(cls), 1, 'uint8'));
        o.flag = field_or(info, 'flag', zeros(numel(cls), 1, 'uint8'));
        o.notes = [codes.class; classifier; string(datetime('now')); codes.syn_type; codes.flag];
        o.class_meta = struct('classifier', classifier, 'date', datetime('now'), 'version', version, ...
            'cruise', p.cruise, 'config', c, 'bead_source', src, 'bead_excluded', bexcl, 'ref', ref, 'cfg', cfg, ...
            'pro_measured', pro(k), 'pro_near_trigger_frac', pro_near_trigger(cls, info, pro(k)), ...
            'ncomp', field_or(info, 'ncomp', NaN), 'comp', field_or(info, 'comp', table), ...
            'trigger_edges', field_or(info, 'trigger_edges', field_or(info, 'edges', struct()))); % hbac: info.edges
        if pres
            o.time_gate = time_gate(fcsdat, fcshdr, tgf);
            o.notes(6) = "time_gate: events in the last " + 100*tgf + "% of acquisition time; " + ...
                "counts per ml use time_gate events and class_meta.vol_analyzed_ml";
            m = pmeta(k);
            o.class_meta.sample = "preserved";
            o.class_meta.runtype = string(rt{k});
            o.class_meta.cast = m.cast; o.class_meta.niskin = m.niskin;
            o.class_meta.sample_label = string(m.sample_label); o.class_meta.repeat = m.repeat;
            o.class_meta.sizefrac = string(m.sizefrac);
            o.class_meta.time_gate_fraction = tgf;
            o.class_meta.vol_analyzed_ml = tgf * m.vol_analyzed / 1e6;
            o.class_meta.bead_fallback_days = field_or(ref, 'fallback', struct());
            o.class_meta.pro = field_or(info, 'pro', struct()); % lower_bound used from CHL runs only (n: Pro count)
            o.class_meta.large_pro = field_or(info, 'large_pro', struct());
            o.class_meta.euk_autofluor = field_or(info, 'euk_autofluor', struct());
            o.class_meta.discrete_cfg = field_or(info, 'discrete_cfg', struct()); % settings used (version tag)
            if strcmp(rt{k}, 'hbac')
                o.class_meta.hbac = struct('n_core', field_or(info, 'n_core', NaN), ...
                    'chlx', field_or(info, 'chlx', []), 'pex', field_or(info, 'pex', []), 'hna', field_or(info, 'hna', struct()));
                o.bact_type = field_or(info, 'bact_type', zeros(numel(cls), 1, 'uint8'));
                o.notes(7) = codes.bact_type;
            end
        end
        save_classfile(classfile{k}, o)
        if do_review && ~strcmp(rt{k}, 'hbac') % SYBR runs: no review frame yet
            try
                warning('off', 'stats:statrobustfit:IterationLimit')
                Q{k} = review_sample(fcsdat, fcshdr, cls, info, [], revdir, names{k}, t(k), pro(k), R);
            catch ME
                warning('process_cruise:review', '%s: review frame failed: %s', names{k}, ME.message)
            end
        end
        nev(k) = numel(cls); ncomp(k) = field_or(info, 'ncomp', NaN);
        counts(k, :) = [sum(cls == 0) sum(cls == 1) sum(cls == 2) sum(cls == 3) sum(cls == 4) ...
            sum(cls == 5) sum(cls == 7) sum(cls == 8) sum(cls == 9) sum(o.flag > 0)];
    catch ME
        err(k) = string(ME.identifier) + ": " + string(ME.message);
        warning('process_cruise:classify', '%s: %s', Fi.fcslist{k}, ME.message)
    end
end
L = [table(string(Fi.fcslist), t, pro, nev, ncomp, 'VariableNames', {'file' 'time' 'pro_measured' 'n_events' 'ncomp'}) ...
    array2table(counts, 'VariableNames', {'n0' 'n1_euk' 'n2_syn' 'n3_crypto_dim' 'n4_crypto_bright' 'n5_euk_syn' 'n7_pro' 'n8_syn_particle' 'n9_bact' 'n_flagged'}) ...
    table(err, 'VariableNames', {'error'})];
if pres, L.runtype = string(rt); end
save([p.outpath 'class_log.mat'], 'L')
fprintf('classify: %d files classified, %d failed; log in %sclass_log.mat\n', sum(err == ""), sum(err ~= ""), p.outpath)
if do_review
    review_step(revdir, Q, R, p, version)
end
end


function review_step(revdir, Q, R, p, version)
% merge this run's review QC into the cruise's (files classified earlier keep theirs if
% their frame is still there) and rewrite the viewer
Q = Q(~cellfun(@isempty, Q));
qf = fullfile(revdir, 'review_qc.mat');
if exist(qf, 'file')
    o = load(qf, 'Q');
    if isfield(o, 'Q') && ~isempty(o.Q)
        oldQ = o.Q(:);
        newf = cellfun(@(q) q.file, Q, 'UniformOutput', false);
        keep = ~ismember(cellfun(@(q) q.file, oldQ, 'UniformOutput', false), newf) ...
            & cellfun(@(q) exist(fullfile(revdir, q.frame), 'file') == 2, oldQ);
        Q = [oldQ(keep); Q(:)];
    end
end
if isempty(Q), return, end
build = [char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')) ' (' version ')'];
decfile = fullfile(p.outpath, 'review_decisions.csv');
err = strings(0, 1); %#ok<NASGU> same variables as review_cruise writes
save(qf, 'Q', 'err', 'build', 'decfile')
review_viewer(revdir, Q, R, p.cruise, decfile, build);
fprintf('classify: review of %d files: %s\n', numel(Q), fullfile(revdir, 'index.html'))
end


function codes = class_codes()
codes.class = "class: 0 = unassigned, 1 = Euk (no PE), 2 = Syn, 3 = dim crypto (PE euk, red-shifted PE), " + ...
    "4 = bright crypto (PE euk), 5 = Euk with Syn (Syn-like PE, possible grazer), 7 = Pro, " + ...
    "8 = Syn on particle (count with Syn abundance, not Syn size or carbon); cell volume in cubic microns";
codes.syn_type = "syn_type (classes 2, 5, 8): 1 = type A (low PUB:PEB), 2 = type B (mid), 3 = type C (high PUB:PEB), " + ...
    "9 = spectrum not usable, 0 = not Syn-bearing";
codes.bact_type = "bact_type (SYBR runs, class 9 only): 1 = LNA, 2 = HNA (low / high nucleic acid; cut on BL1 - slope*(SSC+2), " + ...
    "class_meta.hbac.hna: method 'fixed' = no SYBR valley, split nominal), 0 = not bacteria";
codes.flag = "flag (class 0 only): 1 = Syn-pigment particle with low chl (degraded Syn/Euk, grazer, detritus), " + ...
    "2 = below the Euk chl gap, 3 = trigger-noise pile at the chl edge (pre-2019 files), " + ...
    "4 = Pro-like (no PE, no larger than Syn) in a file run without Pro settings: Pro lower bound, 0 = none";
end


function g = time_gate(fcsdat, fcshdr, frac)
% preserved samples: events in the last frac of the acquisition time (the hand-gated
% analysis left out the first 20% of each run; ApplyAWSgates_hierarchical)
it = find(strcmpi({fcshdr.par.name}, 'Time'), 1);
if isempty(it), it = 1; end
full = 1000 * etime(datevec([fcshdr.date ' ' fcshdr.stoptime]), datevec([fcshdr.date ' ' fcshdr.starttime]));
g = fcsdat(:, it) <= full & fcsdat(:, it) > full * (1 - frac);
end


function f = pro_near_trigger(cls, info, pro)
% Pro-setting files: fraction of Pro (class 7) within 0.1 log of the SSC or chl (BL3)
% trigger edge, i.e. where the trigger may be cutting off small, dim Pro (a quality
% flag, not a correction). NaN for other files or without Pro
f = NaN;
p = cls == 7;
if ~pro || ~any(p) || ~isfield(info, 'features') || ~isfield(info, 'trigger_edges'), return, end
near = false(sum(p), 1);
for ch = {'SSC' 'BL3'}
    if isfield(info.trigger_edges, ch{1})
        near = near | info.features.(ch{1})(p) - info.trigger_edges.(ch{1}) < 0.1;
    end
end
f = mean(near);
end


function v = field_or(s, f, default)
if isfield(s, f), v = s.(f); else, v = default; end
end


function save_classfile(fn, o)
save(fn, '-struct', 'o')
end


function write_log(p)
fn = [p.outpath 'process_cruise_log.mat'];
run = struct('date', datetime('now'), 'p', p);
runs = run;
if exist(fn, 'file')
    r = load(fn, 'runs');
    if isfield(r, 'runs'), runs = [r.runs(:); run]; end
end
save(fn, 'runs')
end


function steps = step_list(req)
% step labels in run order; "core" and "all" expand
order = ["files" "beads" "classify" "volume" "table" "metadata" "movie" "voldists"];
req = lower(strtrim(req(:)'));
if any(req == "all"), req = [req order]; end
if any(req == "core"), req = [req order(1:6)]; end
req = req(~ismember(req, ["core" "all"]));
bad = setdiff(req, order);
if ~isempty(bad)
    error('process_cruise:steps', 'Unknown step "%s". Steps are: %s, or "core" / "all"', bad(1), strjoin(order, ', '))
end
steps = order(ismember(order, req));
end
