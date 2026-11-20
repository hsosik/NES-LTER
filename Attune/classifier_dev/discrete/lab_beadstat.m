function beadstat = lab_beadstat(outfile)
% function beadstat = lab_beadstat(outfile)
% Bead statistics for the lab FCB bead checks run between cruises, which is
% when preserved (discrete) samples are run. Same method as the cruise bead
% tables (measure_beads, process_beads_PT_adjust_3), so bead_reference can
% normalize preserved runs the same way as underway runs.
% Only runs from 16 Aug 2019 on (current filter layout, OD2 on GL1).
%
% outfile - default classifier_dev/output/discrete/lab_beadstat.mat
% Read it with s = load(outfile); bead_reference(s.beadstat, t, cfg)
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
if nargin < 1 || isempty(outfile)
    outfile = fullfile(fileparts(here), 'output', 'discrete', 'lab_beadstat.mat');
end
if ~exist(fileparts(outfile), 'dir'), mkdir(fileparts(outfile)), end
fpath = '\\sosiknas1\lab_data\Attune\cruise_data\beads\FCB_bead_mix_experiment_settings\between_cruises';
F = dir(fullfile(fpath, '*.fcs'));
files = {F.name}';
dat = cell(numel(files), 1); hdr = cell(numel(files), 1); keep = false(numel(files), 1);
for f = 1:numel(files)
    try
        [d, h] = fca_readfcs(fullfile(fpath, files{f}));
        t = datenum([h.date ' ' h.starttime]);
        if t >= datenum(2019, 8, 16)
            dat{f} = single(d); hdr{f} = h; keep(f) = true;
        end
    catch err
        warning('could not read %s: %s', files{f}, err.message)
    end
end
files = files(keep); dat = dat(keep); hdr = hdr(keep);
conf = attune_config('EN720'); %any post-Aug 2019 cruise: current layout, OD2 on GL1
beadstat = measure_beads(dat, hdr, files, conf);
beadstat = sortrows(beadstat, 'time');
meta.source = fpath;
meta.processed = datetime('now');
meta.code = "classifier_dev/discrete/lab_beadstat (measure_beads)";
save(outfile, 'beadstat', 'meta')
fprintf('%d lab bead runs; QC good 0.5 um %d, 1 um %d, 6 um %d\nsaved %s\n', height(beadstat), ...
    sum(beadstat.QC05), sum(beadstat.QC1), sum(beadstat.QC6), outfile)
end
