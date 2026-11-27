function preserved_volume(p, parallel)
% function preserved_volume(p, parallel)
% Cell volumes for preserved (discrete) samples: process_cruise volume step
% with Sample="preserved". Same calibration as underway samples (bead-relative
% SSC-H merged with GL1-H, Attune_SC2vol 'SSC-H'), but the 1 um bead positions
% come from the bead reference saved with the classes (class_meta.ref: same-day
% lab bead runs, discrete_bead_reference), moved to the file's SSC and GL1
% voltages with the bead voltage exponent (ref.k). Volume and classes therefore
% use the same beads, and no bead run at the file's exact voltage is needed.
% Per-file calculation: volume_file (mode 'preserved'), also used by the classify
% step when the classify and volume steps are run together.
%
% The live size calibration is applied as is: preserved cells are smaller than
% live ones from the same water (classifier_dev\output\discrete\
% preservation_correction_report.md); no correction yet, noted in vol_notes.
%
% Adds to each class file: volume_cubic_microns, ssc_merge_bdnorm, file_hv,
% vol_notes. Files without a bead reference (no lab beads that day) are skipped.
%
% Nov 2026, process_cruise preserved mode

if nargin < 2 || isempty(parallel), parallel = license('test', 'Distrib_Computing_Toolbox'); end
d = dir(fullfile(p.classpath, '*.mat'));
V = volume_setup(p);
n = numel(d); skipped = strings(n, 1); rows = cell(n, 1);
names = regexprep({d.name}', '\.mat$', '.fcs');
load(fullfile(p.outpath, 'FCSfileinfo.mat'), 'FCSfileinfo')
[in, ia] = ismember(names, FCSfileinfo.fcslist);
filetime = NaT(n, 1); filetime(in) = datetime(FCSfileinfo.matdate_start(ia(in)), 'ConvertFrom', 'datenum');
if parallel && isempty(gcp('nocreate')), parpool; end
M = 0; if parallel, M = Inf; end
fpath = p.fpath;
parfor (k = 1:n, M)
    fn = fullfile(d(k).folder, d(k).name);
    c = load(fn);
    if ~isfield(c, 'class_meta') || ~isfield(c.class_meta, 'ref') || isempty(c.class_meta.ref)
        skipped(k) = "no bead reference"; continue
    end
    try
        [fcsdat, fcshdr] = fca_readfcs(fullfile(fpath, regexprep(d(k).name, '\.mat$', '.fcs')));
        [c, rows{k}] = volume_file(c, fcsdat, fcshdr, V, []);
    catch ME
        skipped(k) = string(ME.message); continue
    end
    save_struct(fn, c)
end
volume_merge_table(p, V, names, filetime, rows, false); % which files have volumes (cruise_status)
fprintf('volume: %d preserved class files calibrated, %d skipped\n', sum(skipped == ""), sum(skipped ~= ""))
if any(skipped ~= "")
    q = find(skipped ~= "");
    for k = q(1:min(5, end))', fprintf('  %s: %s\n', d(k).name, skipped(k)); end
end
end


function save_struct(fn, c)
save(fn, '-struct', 'c')
end
