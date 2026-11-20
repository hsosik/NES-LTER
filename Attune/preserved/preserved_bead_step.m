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
%
% made - true when the cruise table was written
%
% Nov 2026, process_cruise preserved mode

if nargin < 2, force = false; end
labdir = [char([92 92]) 'sosiknas1\Lab_data\Attune\cruise_data\beads\FCB_bead_mix_experiment_settings\between_cruises'];
labfile = fullfile(labdir, 'outputs', 'lab_beadstat.mat');
lab = update_lab_table(labdir, labfile, force);

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
            if ismember('excluded', lab.Properties.VariableNames), good = good & ~lab.excluded; end
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
meta.made = datetime('now');
meta.cruise = p.cruise;
meta.days = D;
meta.code = "preserved\preserved_bead_step (lab table: classifier_dev\discrete\lab_beadstat)";
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


function lab = update_lab_table(labdir, labfile, force)
% rebuild the lab table when the folder has bead runs it has not scanned
f = dir(fullfile(labdir, '*.fcs')); names = {f.name}';
if exist(labfile, 'file') && ~force
    s = load(labfile, 'beadstat', 'meta');
    if isfield(s.meta, 'files_scanned') && isempty(setdiff(names, s.meta.files_scanned))
        lab = s.beadstat; return
    end
end
fprintf('beads: building the lab bead table from %d runs in %s\n', numel(names), labdir)
if ~exist(fileparts(labfile), 'dir'), mkdir(fileparts(labfile)); end
lab_beadstat(labfile);
s = load(labfile, 'beadstat', 'meta');
meta = s.meta; meta.files_scanned = names; beadstat = s.beadstat; %#ok<NASGU>
save(labfile, 'beadstat', 'meta')
lab = s.beadstat;
end
