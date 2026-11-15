function p = process_cruise(cruise, steps, opts)
% function p = process_cruise(cruise, steps, opts)
% Process one cruise of Attune data end to end, or any subset of the steps.
% Replaces the copied process_wrapper_* scripts: every per-cruise setting
% comes from attune_cruise_settings (attune_cruise_settings.csv plus
% classifier_dev/attune_config), and every cruise uses the same classifier,
% assign_class_v2.
%
% e.g. process_cruise('EN727')            % all main steps (1:6)
%      process_cruise('EN727', [3 4 5])   % reclassify, volumes, table
%      process_cruise('EN727', 3, struct('files', 5, 'outpath', 'C:\scratch\EN727\'))
%
% Steps
%   1 files     FCS_DateTimeList -> bead_calibrated\FCSfileinfo.mat
%   2 beads     bead table bead_calibrated\beadstat.mat
%               - FCB bead cruises (attune_config fcb_beads): process_beads_PT_adjust_3.
%                 It will not overwrite a table made another way (meta.source set)
%                 unless opts.force.
%               - no FCB beads (EN608 to EN627): pt_transfer\make_perf_beadstat.
%                 EN627 keeps its PT conversion unless opts.force.
%   3 classify  assign_class_v2 on every sample file -> bead_calibrated\class\<file>.mat
%               holding class, syn_type, flag, notes, class_meta (see below)
%   4 volume    batch_mergedSSC2volume: adds volume_cubic_microns etc. to the class files
%   5 table     generate_attune_table -> bead_calibrated\AttuneTable.mat
%   6 underway  match_Attune_underway_LTER -> bead_calibrated\AttuneTable_uw_match.mat
%   7 movie     attune_lter_moviemaker (optional)
%   8 voldists  get_cruise_voldists_fromEDItable2, Plot_Voldists_function (optional)
%
% opts (all optional)
%   .root      folder of cruise folders (default \\sosiknas1\Lab_data\Attune\cruise_data\)
%   .outpath   output folder instead of <cruise>\bead_calibrated\ (testing). Steps 3 to 6
%              read FCSfileinfo.mat and beadstat.mat from here, so copy them in first.
%   .files     classify only these: a number n (n files spread through the cruise) or a
%              cellstr of FCS names
%   .append    true: skip files that already have a class file (default false)
%   .parallel  true: parfor over files in step 3 (default true when a pool can start)
%   .force     passed to the bead step (default false)
%
% Class file variables written by step 3 (step 4 loads and re-saves the file, keeping them)
%   class      per event, double: 0 unassigned, 1 Euk, 2 Syn, 3 dim crypto, 4 bright crypto,
%              5 Euk with Syn, 7 Pro, 8 Syn on particle
%   syn_type   per event, uint8: 1/2/3 Syn type A/B/C (PUB:PEB low/mid/high),
%              9 Syn-bearing but spectrum not usable, 0 not Syn-bearing (classes 2, 5, 8 only)
%   flag       per event, uint8, only on class 0: 1 Syn-pigment particle with low chl
%              (degraded Syn/Euk, grazer, detritus), 2 below the Euk chl gap, 3 Pro tail
%              at the chl trigger
%   notes      string: (1) class codes, (2) classifier name, (3) date, (4) syn_type codes,
%              (5) flag codes
%   class_meta struct: classifier, date, code version, cruise configuration, bead
%              reference used, classifier settings, mixture components
%
% Each run appends its settings to bead_calibrated\process_cruise_log.mat.
%
% Nov 2026, processing wrapper redesign (replaces process_wrapper_2024b)

if nargin < 2 || isempty(steps), steps = 1:6; end
if nargin < 3, opts = struct(); end
def = struct('root', '', 'outpath', '', 'files', [], 'append', false, 'parallel', [], 'force', false);
f = fieldnames(def);
for i = 1:numel(f), if ~isfield(opts, f{i}), opts.(f{i}) = def.(f{i}); end, end

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'classifier_dev'), fullfile(here, 'classifier_dev', 'pt_transfer'))

%% settings
s = attune_cruise_settings(cruise, opts.root);
c = s.config;
p.cruise = s.cruise;
p.basepath = [s.basepath filesep];
p.fpath = fullfile(s.basepath, 'FCS', filesep);
if isempty(opts.outpath)
    p.outpath = fullfile(s.basepath, 'bead_calibrated', filesep);
else
    p.outpath = [regexprep(char(opts.outpath), '[\\/]+$', '') filesep];
end
p.classpath = fullfile(p.outpath, 'class', filesep);
p.OD2setting = c.OD2;      % read by batch_mergedSSC2volume and the movie maker
p.SSCDIM = 'H';            % batch_mergedSSC2volume: height channels
p.beadfiles2include = {'FCB_bead' 'fcb_bead'};
p.filetype2exclude = s.exclude;
p.uw_fullname = s.underway;
p.assign_class_function = 'assign_class_v2';
p.config = c;
p.steps = steps;
p.opts = opts;
if strcmp(c.layout, 'v2019'), p.moviechannels = 'early'; else, p.moviechannels = 'late'; end
p.framemaker = 'make_movieframe_density';
p.stepsize = 1;

fprintf('%s: %s\n  OD2 on %s, layout %s, FCB beads %d\n', p.cruise, p.basepath, c.OD2, c.layout, c.fcb_beads)
if ~isempty(s.notes), fprintf('  note: %s\n', s.notes); end
if ~exist(p.classpath, 'dir'), mkdir(p.classpath); end
write_log(p)

%% 1 files
if ismember(1, steps)
    FCSfileinfo = FCS_DateTimeList(p.fpath);
    save([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo')
    fprintf('step 1: %d FCS files listed\n', height(FCSfileinfo))
end

%% 2 beads
if ismember(2, steps)
    bead_step(p, c, opts.force)
end

%% 3 classify
if ismember(3, steps)
    classify_step(p, c, opts)
end

%% 4 volume
if ismember(4, steps)
    batch_mergedSSC2volume(p)
end

%% 5 table
if ismember(5, steps)
    generate_attune_table(p.classpath, [p.outpath 'FCSfileinfo.mat'])
end

%% 6 underway
if ismember(6, steps)
    match_Attune_underway_LTER([p.outpath 'AttuneTable.mat'], p.uw_fullname);
end

%% 7 movie
if ismember(7, steps)
    close all
    figure('WindowState', 'maximized')
    attune_lter_moviemaker(p.fpath, p.classpath, p.OD2setting, p.framemaker, p.moviechannels, p.stepsize)
end

%% 8 volume distributions
if ismember(8, steps)
    get_cruise_voldists_fromEDItable2(p.basepath)
    Plot_Voldists_function(p.basepath, p.outpath)
end
end


function bead_step(p, c, force)
fn = [p.outpath 'beadstat.mat'];
if c.fcb_beads
    if exist(fn, 'file')
        m = load(fn, 'meta');
        if isfield(m, 'meta') && isfield(m.meta, 'source') && ~force
            fprintf('step 2: keeping %s (meta.source "%s"); opts.force = true to replace it with an FCB table\n', fn, m.meta.source)
            return
        end
    end
    load([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo')
    process_beads_PT_adjust_3(p.basepath, FCSfileinfo, p.beadfiles2include, 'FCB', c.OD2);
else
    try
        info = make_perf_beadstat(regexprep(p.basepath, '[\\/]+$', ''), struct('force', force, 'bead_dir', p.outpath));
        fprintf('step 2: PT-transfer bead table from %d performance tests: %s\n', info.ntests, info.installed)
    catch ME
        if strcmp(ME.identifier, 'make_perf_beadstat:exists')
            m = load(fn, 'meta');
            fprintf('step 2: keeping %s (meta.source "%s")\n', fn, m.meta.source)
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
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
Fi = sortrows(Fi, 'matdate_start');
if iscell(opts.files) || isstring(opts.files)
    Fi = Fi(ismember(Fi.fcslist, cellstr(opts.files)), :);
elseif ~isempty(opts.files) && opts.files < height(Fi)
    Fi = Fi(unique(round(linspace(1, height(Fi), opts.files))), :);
end
classfile = strcat(p.classpath, regexprep(Fi.fcslist, '\.fcs$', '.mat'));
if opts.append
    have = cellfun(@(x) exist(x, 'file') == 2, classfile);
    Fi = Fi(~have, :); classfile = classfile(~have);
end
n = height(Fi);
fprintf('step 3: classifying %d files\n', n)
if n == 0, return, end

[bs, bmeta] = load_beadstat(p.outpath);
cfg = classifier_config();
cfg.gl1_is_od2_scatter = strcmp(c.OD2, 'GL1');
cfg.cruise = p.cruise; %bead_reference: OD2 placement (ref.OD2) and channel map from attune_config
[~, h] = system(['git -C "' fileparts(mfilename('fullpath')) '" rev-parse --short HEAD']);
[~, dirty] = system(['git -C "' fileparts(mfilename('fullpath')) '" status --porcelain -- classifier_dev']);
version = strtrim(h);
if ~isempty(strtrim(dirty)), version = [version ' + uncommitted classifier_dev changes']; end
src = "FCB beads";
if isfield(bmeta, 'source'), src = string(bmeta.source); end

codes = class_codes();
fcsfile = strcat(p.fpath, Fi.fcslist);
t = Fi.matdate_start; pro = Fi.pro_measured;
nev = nan(n, 1); ncomp = nan(n, 1); counts = nan(n, 9); err = strings(n, 1);
use_par = opts.parallel;
if isempty(use_par), use_par = license('test', 'Distrib_Computing_Toolbox') && n > 4; end
if use_par && isempty(gcp('nocreate')), parpool; end
M = 0; if use_par, M = Inf; end
parfor (k = 1:n, M)
    try
        [fcsdat, fcshdr] = fca_readfcs(fcsfile{k});
        ref = bead_reference(bs, t(k), cfg);
        [cls, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, pro(k));
        o = struct();
        o.class = double(cls);
        o.syn_type = field_or(info, 'syn_type', zeros(size(cls), 'uint8'));
        o.flag = field_or(info, 'flag', zeros(size(cls), 'uint8'));
        o.notes = [codes.class; "assign_class_v2"; string(datetime('now')); codes.syn_type; codes.flag];
        o.class_meta = struct('classifier', "assign_class_v2", 'date', datetime('now'), 'version', version, ...
            'cruise', p.cruise, 'config', c, 'bead_source', src, 'ref', ref, 'cfg', cfg, ...
            'pro_measured', pro(k), 'ncomp', info.ncomp, 'comp', info.comp, ...
            'trigger_edges', field_or(info, 'trigger_edges', struct()));
        save_classfile(classfile{k}, o)
        nev(k) = numel(cls); ncomp(k) = info.ncomp;
        counts(k, :) = [sum(cls == 0) sum(cls == 1) sum(cls == 2) sum(cls == 3) sum(cls == 4) ...
            sum(cls == 5) sum(cls == 7) sum(cls == 8) sum(o.flag > 0)];
    catch ME
        err(k) = string(ME.identifier) + ": " + string(ME.message);
        warning('process_cruise:classify', '%s: %s', Fi.fcslist{k}, ME.message)
    end
end
L = [table(string(Fi.fcslist), t, pro, nev, ncomp, 'VariableNames', {'file' 'time' 'pro_measured' 'n_events' 'ncomp'}) ...
    array2table(counts, 'VariableNames', {'n0' 'n1_euk' 'n2_syn' 'n3_crypto_dim' 'n4_crypto_bright' 'n5_euk_syn' 'n7_pro' 'n8_syn_particle' 'n_flagged'}) ...
    table(err, 'VariableNames', {'error'})];
save([p.outpath 'class_log.mat'], 'L')
fprintf('step 3: %d files classified, %d failed; log in %sclass_log.mat\n', sum(err == ""), sum(err ~= ""), p.outpath)
end


function codes = class_codes()
codes.class = "class: 0 = unassigned, 1 = Euk (no PE), 2 = Syn, 3 = dim crypto (PE euk, red-shifted PE), " + ...
    "4 = bright crypto (PE euk), 5 = Euk with Syn (Syn-like PE, possible grazer), 7 = Pro, " + ...
    "8 = Syn on particle (count with Syn abundance, not Syn size or carbon); cell volume in cubic microns";
codes.syn_type = "syn_type (classes 2, 5, 8): 1 = type A (low PUB:PEB), 2 = type B (mid), 3 = type C (high PUB:PEB), " + ...
    "9 = spectrum not usable, 0 = not Syn-bearing";
codes.flag = "flag (class 0 only): 1 = Syn-pigment particle with low chl (degraded Syn/Euk, grazer, detritus), " + ...
    "2 = below the Euk chl gap, 3 = Pro tail at the chl trigger, 0 = none";
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
