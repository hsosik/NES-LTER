% debug_example1: which rule drops the bright Syn-like events in example 1 (AR92 15Aug(3) Sample(18)),
% and do they look like elongated / dividing Syn (longer pulses, same pigment ratios)?
addpath('..')
root = '\\sosiknas1\lab_data\Attune\cruise_data\20250814_AR92';
bp = fullfile(root, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
cfg = classifier_config(); sig = cfg.sig;
fn = 'NESLTER_AR92_15Aug2025(23)_Group_day0_Pro(5).fcs';
k = find(strcmp(FCSfileinfo.fcslist, fn));
[d, h] = fca_readfcs(fullfile(root, 'FCS', fn));
ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
[new, info] = assign_class_v2(d, h, ref, cfg, FCSfileinfo.pro_measured(k));
c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class'); old = c.class(:, 1);
F = info.features; r = info.syn_ref; t = info.trace; r0 = t.syn_ref0;
disp(info.comp)
fprintf('syn ref before check: size %.2f BL3 %.2f PEchl %.2f GL2GL3 %.2f\n', r0.size, r0.BL3, r0.PEchl, r0.GL2GL3)
fprintf('syn ref after check:  size %.2f BL3 %.2f PEchl %.2f GL2GL3 %.2f\n', r.size, r.BL3, r.PEchl, r.GL2GL3)
lost = old == 2 & new == 0;
hc = t.has_chl;
gl = zeros(size(old)); gl(hc) = t.gmm_label; go = false(size(old)); go(hc) = t.gmm_outlier;
for grp = {'bright (BL3 > ref+0.5)', 'other'}
    if grp{1}(1) == 'b', g = lost & F.BL3 > r.BL3 + 0.5; else, g = lost & ~(F.BL3 > r.BL3 + 0.5); end
    fprintf('\n%s lost Syn: %d\n', grp{1}, sum(g))
    fprintf('  median BL3 %.2f GL2 %.2f size %.2f PEchl %.2f GL2GL3 %.2f chlSSC %.2f\n', median(F.BL3(g)), median(F.GL2(g)), median(F.size(g)), median(F.PEchl(g)), median(F.GL2GL3(g)), median(F.chlSSC(g)))
    fprintf('  not has_chl %d\n', sum(g & ~hc))
    fprintf('  GMM: outlier %d; component label 0/1/2/3/7: %s\n', sum(g & go), mat2str(arrayfun(@(x) sum(g & hc & ~go & gl == x), [0 1 2 3 7])))
    fprintf('  after GMM class 2: %d; detritus rule removes %d\n', sum(g & t.after_gmm == 2), sum(g & t.after_gmm ~= 0 & t.after_detritus == 0))
    sy = g & t.after_detritus == 2;
    fprintf('  Syn check removes %d: PEchl %d, too big %d, too small %d, GL2GL3 %d\n', sum(sy & t.after_syncheck == 0), ...
        sum(sy & abs(F.PEchl - r0.PEchl) > sig.syn_PEchl_tol), sum(sy & F.size > r0.size + sig.syn_size_tol), sum(sy & F.size < r0.size - sig.syn_size_tol_below), ...
        sum(sy & F.GL3 > info.floor.GL3 + sig.ratio_min_above_floor & abs(F.GL2GL3 - r0.GL2GL3) > sig.syn_GL2GL3_tol))
    fprintf('  PE-euk step removes %d\n', sum(g & t.after_syncheck ~= 0 & t.after_peeuk == 0))
end

% Pro dropped by v2: where do they fall out?
e = info.trigger_edges; disp(e)
chl_min = max(sig.chl_min, e.BL3 + sig.trigger_margin);
fprintf('chl_min used %.2f (sig.chl_min %.2f, BL3 edge %.2f + margin %.2f)\n', chl_min, sig.chl_min, e.BL3, sig.trigger_margin)
lp = old == 7 & new == 0; kp = old == 7 & new == 7;
fprintf('old Pro %d: kept %d, dropped %d\n', sum(old == 7), sum(kp), sum(lp))
fprintf('dropped: BL3 <= chl_min %d; not has_chl %d; GMM outlier %d; GMM label 0/1/2/3/7 %s; detritus %d; euk-size %d\n', ...
    sum(lp & F.BL3 <= chl_min), sum(lp & ~hc), sum(lp & go), mat2str(arrayfun(@(x) sum(lp & hc & ~go & gl == x), [0 1 2 3 7])), ...
    sum(lp & t.after_gmm ~= 0 & t.after_detritus == 0), sum(lp & t.after_syncheck == 1 & t.after_peeuk == 0))
fprintf('BL3 of dropped: p5 %.2f median %.2f p95 %.2f; kept Pro BL3 p5 %.2f median %.2f\n', prctile(F.BL3(lp), [5 50 95]), prctile(F.BL3(kp), [5 50]))
fprintf('SSC of dropped median %.2f, kept %.2f; pro_max_SSC %.2f pro_max_BL3 %.2f\n', median(F.size(lp)), median(F.size(kp)), sig.pro_max_SSC, sig.pro_max_BL3)
fprintf('old Pro BL3 min %.2f\n', min(F.BL3(old == 7)))
ap = old == 0 & new == 7; fprintf('added Pro %d: BL3 median %.2f size median %.2f\n', sum(ap), median(F.BL3(ap)), median(F.size(ap)))
