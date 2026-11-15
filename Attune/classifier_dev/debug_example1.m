% debug_example1: which rule drops the bright Syn-like events in example 1 (AR92 15Aug(3) Sample(18)),
% and do they look like elongated / dividing Syn (longer pulses, same pigment ratios)?
addpath('..')
root = '\\sosiknas1\lab_data\Attune\cruise_data\20250814_AR92';
bp = fullfile(root, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
cfg = classifier_config(); sig = cfg.sig;
fn = 'NESLTER_AR92_15Aug2025(3)_Group_day0_Sample(18).fcs';
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

% elongated or dividing Syn? pulse width and area/height vs core Syn
pn = {h.par.name};
col = @(n) d(:, strcmp(pn, n));
FSCW = log10(max(col('FSC-W'), 1));
AH = @(ch) log10(max(col([ch '-A']), 1) ./ max(col([ch '-H']), 1));
BL3AH = AH('BL3'); SSCAH = AH('SSC'); GL2AH = AH('GL2');
core = new == 2 & abs(F.BL3 - r.BL3) < 0.2;
brightkept = new == 2 & F.BL3 > r.BL3 + 0.5;
bl = lost & F.BL3 > r.BL3 + 0.5;
euk = new == 1 & abs(F.BL3 - median(F.BL3(bl))) < 0.3;
fprintf('\n%-26s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n', 'group', 'n', 'BL3', 'GL2', 'size', 'FSC', 'FSC-W', 'BL3AH', 'SSCAH', 'GL2AH')
grps = {core, brightkept, bl, euk}; gn = {'core Syn (v2)', 'bright Syn kept (v2)', 'bright Syn lost', 'Euks of same chl'};
for i = 1:numel(grps)
    g = grps{i};
    fprintf('%-26s %6d %6.2f %6.2f %6.2f %6.2f %6.2f %6.3f %6.3f %6.3f\n', gn{i}, sum(g), median(F.BL3(g)), median(F.GL2(g)), median(F.size(g)), ...
        median(F.FSC(g)), median(FSCW(g)), median(BL3AH(g)), median(SSCAH(g)), median(GL2AH(g)));
end
% does the lost group continue the Syn size-chl line, or sit off it?
p = polyfit(F.BL3(new == 2), F.size(new == 2), 1);
fprintf('Syn size-on-chl slope %.2f; size residual: lost bright %.2f, kept bright %.2f\n', p(1), ...
    median(F.size(bl) - polyval(p, F.BL3(bl))), median(F.size(brightkept) - polyval(p, F.BL3(brightkept))))
% does the PE spectrum (GL2:GL3) drift with brightness, or with pulse length?
s2 = old == 2;
edges = -1:0.25:1.25;
fprintf('\nold Syn by chl bin: BL3 | n | GL2GL3 median | BL3AH median | frac GL2GL3 > ref+0.15\n')
for b = 1:numel(edges) - 1
    g = s2 & F.BL3 >= edges(b) & F.BL3 < edges(b + 1);
    fprintf('%5.2f-%5.2f %6d %6.2f %6.3f %6.3f\n', edges(b), edges(b + 1), sum(g), median(F.GL2GL3(g)), median(BL3AH(g)), mean(F.GL2GL3(g) > r0.GL2GL3 + sig.syn_GL2GL3_tol))
end
g = s2 & F.BL3 > 0.25;
fprintf('bright old Syn: corr(GL2GL3, BL3AH) %.2f, corr(GL2GL3, GL2AH) %.2f, corr(GL2GL3, BL3) %.2f\n', corr(F.GL2GL3(g), BL3AH(g)), corr(F.GL2GL3(g), GL2AH(g)), corr(F.GL2GL3(g), F.BL3(g)))
GL3AH = AH('GL3');
fprintf('GL3AH: core %.3f, bright kept %.3f, bright lost %.3f\n', median(GL3AH(core)), median(GL3AH(brightkept)), median(GL3AH(bl)))
