function preserved_volume(p, parallel)
% function preserved_volume(p, parallel)
% Cell volumes for preserved (discrete) samples: process_cruise volume step
% with Sample="preserved". Same calibration as underway samples (bead-relative
% SSC-H merged with GL1-H, Attune_SC2vol 'SSC-H'), but the 1 um bead positions
% come from the bead reference saved with the classes (class_meta.ref: same-day
% lab bead runs, discrete_bead_reference), moved to the file's SSC and GL1
% voltages with the bead voltage exponent (ref.k). Volume and classes therefore
% use the same beads, and no bead run at the file's exact voltage is needed.
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
cfg = classifier_config();
W = cfg.gl1_merge_bead; % SSC below, GL1 above, mean of the two between (1 um bead units of GL1)
n = numel(d); skipped = strings(n, 1);
if parallel && isempty(gcp('nocreate')), parpool; end
M = 0; if parallel, M = Inf; end
fpath = p.fpath;
parfor (k = 1:n, M)
    fn = fullfile(d(k).folder, d(k).name);
    c = load(fn);
    if ~isfield(c, 'class_meta') || ~isfield(c.class_meta, 'ref') || isempty(c.class_meta.ref)
        skipped(k) = "no bead reference"; continue
    end
    ref = c.class_meta.ref;
    try
        [fcsdat, fcshdr] = fca_readfcs(fullfile(fpath, regexprep(d(k).name, '\.mat$', '.fcs')));
    catch ME
        skipped(k) = string(ME.message); continue
    end
    pn = {fcshdr.par.name};
    is = strcmp(pn, 'SSC-H'); ig = strcmp(pn, 'GL1-H');
    hvS = fcshdr.par(is).hv; hvG = fcshdr.par(ig).hv;
    bS = 10.^(ref.bead_log.SSC + ref.k.SSC * log10(hvS / ref.bead_hv.SSC)); % 1 um bead SSC-H at hvS
    bG = 10.^(ref.bead_log.GL1 + ref.k.GL1 * log10(hvG / ref.bead_hv.GL1));
    ssc = fcsdat(:, is); ssc(ssc == 0) = 1;
    s = ssc / bS; g = fcsdat(:, ig) / bG;
    m = s;
    hi = g > W(2); mid = g > W(1) & ~hi;
    m(hi) = g(hi);
    m(mid) = (s(mid) + g(mid)) / 2;
    v = NaN(size(m)); ii = m > 0;
    [v(ii), vf] = Attune_SC2vol(m(ii), 'SSC-H', 0);
    c.volume_cubic_microns = v;
    c.ssc_merge_bdnorm = m;
    c.file_hv = struct('SSC', hvS, 'GL1', hvG);
    c.vol_notes = {"calibrated: " + string(datetime('now')); "using SSC-H merged with GL1-H"; string(vf); ...
        sprintf('SSC/GL1 merge window %g-%g x 1 um bead (GL1)', W); ...
        sprintf('1 um bead from class_meta.ref (lab beads): log10 SSC-H %.3f at %g V, GL1-H %.3f at %g V', ...
        ref.bead_log.SSC, ref.bead_hv.SSC, ref.bead_log.GL1, ref.bead_hv.GL1); ...
        "preserved sample: live size calibration applied as is (no preservation correction)"};
    save_struct(fn, c)
end
fprintf('volume: %d preserved class files calibrated, %d skipped\n', sum(skipped == ""), sum(skipped ~= ""))
if any(skipped ~= "")
    q = find(skipped ~= "");
    for k = q(1:min(5, end))', fprintf('  %s: %s\n', d(k).name, skipped(k)); end
end
end


function save_struct(fn, c)
save(fn, '-struct', 'c')
end
