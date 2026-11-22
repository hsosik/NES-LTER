function made = preserved_bead_step(p, force)
% function made = preserved_bead_step(p, force)
% Bead table for a cruise's preserved (discrete) samples, which are run in the
% lab: process_cruise beads step with Sample="preserved".
%
% 1. The shared lab bead table (lab FCB bead checks between cruises,
%    classifier_dev\discrete\lab_beadstat, same method as cruise bead tables)
%    is kept on the NAS next to its runs (labfile below) and rebuilt when the
%    folder has bead runs it has not seen (or force).
% 2. The cruise gets <outpath>\beadstat.mat: the lab runs from the days its
%    samples were run, plus for each such day and channel without a good
%    same-day bead population the nearest good runs before and after that day
%    (discrete_bead_reference's per-channel fallback), so a cruise's bead
%    reference does not change when the lab table grows. meta.days lists each
%    run day: samples, same-day bead runs, channels needing a fallback; days
%    with no lab bead run are reported, and their samples are skipped by classify.
% 3. Bead review exclusions: runs marked exclude in this cruise's
%    <outpath>\bead_review_decisions.csv (bead review page) do not count as good
%    when fallback runs are chosen, and load_beadstat drops them for classify
%    (discrete_bead_reference). They apply to this cruise only, not to the shared
%    lab table. After changing them, rerun process_cruise from the beads step.
%
% made - true when the cruise table was written
%
% Nov 2026, process_cruise preserved mode

if nargin < 2, force = false; end
labdir = [char([92 92]) 'sosiknas1\Lab_data\Attune\cruise_data\beads\FCB_bead_mix_experiment_settings\between_cruises'];
labfile = fullfile(labdir, 'outputs', 'lab_beadstat.mat');
[lab, labmeta] = update_lab_table(labdir, labfile, force);
% this cruise's bead review exclusions (not the shared lab table's)
excl = false(height(lab), 1);
dec = [p.outpath 'bead_review_decisions.csv'];
if exist(dec, 'file')
    T = review_decisions(dec);
    [~, f1, e1] = fileparts(cellstr(lab.filename)); [~, f2, e2] = fileparts(T.file(strcmp(T.status, 'exclude')));
    excl = ismember(strcat(f1, e1), strcat(f2, e2)); excl = excl(:);
    if any(excl), fprintf('beads: %d lab bead run(s) excluded in this cruise''s bead review\n', sum(excl)); end
end

load([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo')
F = FCSfileinfo(~contains(FCSfileinfo.fcslist, p.filetype2exclude) & ~strcmp(FCSfileinfo.runtype, 'other') ...
    & isfinite(FCSfileinfo.matdate_start), :);
days = unique(floor(F.matdate_start));
lday = floor(datenum(lab.time));
cfg = classifier_config();
P = {'05' '1' '6'};
keep = false(height(lab), 1);
D = table;
for d = days'
    same = lday == d;
    keep = keep | same;
    fb = strings(0, 1);
    if any(same)
        for c = cfg.chans
            ch = c{1}; pp = P{cfg.bead_pop.(ch)};
            good = lab.(['QC' pp]) & lab.(['peak' pp]).([ch '_H']) > 0 & lab.hv.(ch) > 0;
            good = good & ~excl;
            if any(good & same), continue, end
            fb(end+1) = ch; %#ok<AGROW>
            b = find(good & lday < d, 1, 'last'); a = find(good & lday > d, 1, 'first');
            keep([b; a]) = true;
        end
    end
    D = [D; table(datetime(d, 'ConvertFrom', 'datenum'), sum(floor(F.matdate_start) == d), sum(same), ...
        strjoin(fb, ' '), 'VariableNames', {'day' 'n_samples' 'n_bead_runs' 'fallback_channels'})]; %#ok<AGROW>
end
beadstat = lab(keep, :);
meta.source = "lab FCB bead checks";
meta.lab_table = labfile;
meta.lab_dir = labdir; % bead run FCS files (bead review)
meta.config = labmeta.config; meta.signature = labmeta.signature; meta.nshift = labmeta.nshift;
meta.processed = labmeta.processed;
meta.made = datetime('now');
meta.cruise = p.cruise;
meta.days = D;
meta.notes = ["Lab FCB bead checks (between cruises) from the run days of this cruise's preserved samples, " + ...
    "plus the nearest good runs on other days where a channel had no good same-day run."; ...
    "References: same-day runs only (discrete_bead_reference); a channel without a good same-day " + ...
    "population uses the nearest good run (orange ticks)."; ...
    "Exclusions here apply to this cruise only; rerun process_cruise from the beads step after changing them."];
meta.code = "preserved\preserved_bead_step (lab table: measure_beads on the lab bead checks)";
fn = [p.outpath 'beadstat.mat'];
if exist(fn, 'file') && ~force
    old = load(fn, 'meta');
    if ~isfield(old, 'meta') || ~isfield(old.meta, 'source') || old.meta.source ~= meta.source
        fprintf('beads: keeping %s (not a lab bead table); Force=true to replace it\n', fn)
        made = false; return
    end
end
save(fn, 'beadstat', 'meta')
made = true;
fprintf('beads: %d lab bead runs for %d run days -> %s\n', height(beadstat), numel(days), fn)
nob = D.n_bead_runs == 0;
if any(nob)
    fprintf('beads: no lab bead run on %s (%d samples will be skipped)\n', ...
        strjoin(string(D.day(nob), 'dd-MMM-yyyy'), ', '), sum(D.n_samples(nob)))
end
fbd = strlength(D.fallback_channels) > 0 & ~nob;
if any(fbd)
    fprintf('beads: fallback to nearest good run on %d days (%s)\n', sum(fbd), strjoin(unique(split(strjoin(D.fallback_channels(fbd), ' '))), ' '))
end
end


function [lab, meta] = update_lab_table(labdir, labfile, force)
% the shared lab bead table: measure_beads on every lab bead check since 16 Aug 2019
% (current filter layout, OD2 on GL1), as classifier_dev\discrete\lab_beadstat but
% keeping measure_beads' meta (bead signature with the learned offsets: the bead
% review re-gates runs with it). Rebuilt when the folder has runs it has not
% scanned, when the table predates that meta, or with force.
f = dir(fullfile(labdir, '*.fcs')); names = {f.name}';
if exist(labfile, 'file') && ~force
    s = load(labfile, 'beadstat', 'meta');
    if isfield(s.meta, 'files_scanned') && isempty(setdiff(names, s.meta.files_scanned)) && isfield(s.meta, 'signature')
        lab = s.beadstat; meta = s.meta; return
    end
end
fprintf('beads: building the lab bead table from %d runs in %s\n', numel(names), labdir)
if ~exist(fileparts(labfile), 'dir'), mkdir(fileparts(labfile)); end
dat = cell(numel(names), 1); hdr = dat; keep = false(numel(names), 1);
for k = 1:numel(names)
    try
        [d, h] = fca_readfcs(fullfile(labdir, names{k}));
        if datenum([h.date ' ' h.starttime]) >= datenum(2019, 8, 16)
            dat{k} = single(d); hdr{k} = h; keep(k) = true;
        end
    catch err
        warning('preserved_bead_step:read', 'could not read %s: %s', names{k}, err.message)
    end
end
conf = attune_config('LAB'); % current layout, OD2 on GL1
[beadstat, meta] = measure_beads(dat(keep), hdr(keep), names(keep), conf);
beadstat = sortrows(beadstat, 'time');
meta.source = labdir;
meta.processed = datetime('now');
meta.code = "preserved\preserved_bead_step (measure_beads)";
meta.files_scanned = names;
save(labfile, 'beadstat', 'meta')
fprintf('beads: %d lab bead runs; QC good 0.5 um %d, 1 um %d, 6 um %d -> %s\n', height(beadstat), ...
    sum(beadstat.QC05), sum(beadstat.QC1), sum(beadstat.QC6), labfile)
lab = beadstat;
end
