% debug_syn_checks: which new Syn check removes old-Syn events in a summer file?
addpath('..')
root = '\\sosiknas1\lab_data\Attune\cruise_data\20250814_AR92';
bp = fullfile(root, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
cfg = classifier_config();
Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & ~contains(FCSfileinfo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
Fi = sortrows(Fi, 'matdate_start');
for k = round(linspace(1, height(Fi), 4))
    fn = Fi.fcslist{k};
    if exist(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'file') ~= 2, continue, end
    [d, h] = fca_readfcs(fullfile(root, 'FCS', fn));
    ref = bead_reference(bs, Fi.matdate_start(k), cfg);
    [new, info] = assign_class_v2(d, h, ref, cfg, Fi.pro_measured(k));
    c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class'); old = c.class(:, 1);
    F = info.features; r = info.syn_ref; lost = old == 2 & new == 0; kept = old == 2 & new == 2;
    fprintf('%s pro=%d: old syn %d, lost %d\n', fn, Fi.pro_measured(k), sum(old == 2), sum(lost))
    [~, fi] = attune_features(d, h, ref);
    fprintf('  GL3 floor %.2f, lost events with GL3 > floor+1: %d\n', fi.floor.GL3, sum(lost & F.GL3 > fi.floor.GL3 + 1))
    fprintf('  syn ref size %.2f GL2GL3 %.2f PEchl %.2f\n', r.size, r.GL2GL3, r.PEchl)
    fprintf('  lost: too small %d, GL2GL3 off %d, PEchl off %d, too big %d, chlSSC<det %d\n', sum(lost & F.size < r.size - 0.5), ...
        sum(lost & abs(F.GL2GL3 - r.GL2GL3) > 0.15), sum(lost & abs(F.PEchl - r.PEchl) > 0.7), sum(lost & F.size > r.size + 1), sum(lost & F.chlSSC < -0.8))
    fprintf('  kept GL2GL3 median %.2f IQR %.2f; lost GL2GL3 median %.2f, size median %.2f, GL3 median %.2f, GL2 %.2f, BL3 %.2f\n', ...
        median(F.GL2GL3(kept)), iqr(F.GL2GL3(kept)), median(F.GL2GL3(lost)), median(F.size(lost)), median(F.GL3(lost)), median(F.GL2(lost)), median(F.BL3(lost)))
    fprintf('  kept GL3 median %.2f, GL2 %.2f, size %.2f, BL3 %.2f\n', median(F.GL3(kept)), median(F.GL2(kept)), median(F.size(kept)), median(F.BL3(kept)))
end
