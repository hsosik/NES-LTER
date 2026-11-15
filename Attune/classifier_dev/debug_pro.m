function debug_pro(cruisedir, fn)
% function debug_pro(cruisedir, fn)
% Trace Pro calls in one file: which step drops old Pro, and the PE (GL2) of v2 Pro
% compared with Euks and with old Pro.
addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = fullfile('\\sosiknas1\lab_data\Attune\cruise_data', cruisedir);
bp = fullfile(root, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
cfg = classifier_config(); sig = cfg.sig;
k = find(strcmp(FCSfileinfo.fcslist, fn));
[d, h] = fca_readfcs(fullfile(root, 'FCS', fn));
ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
[new, info] = assign_class_v2(d, h, ref, cfg, FCSfileinfo.pro_measured(k));
c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class'); old = c.class(:, 1);
F = info.features; t = info.trace; r = info.syn_ref;
fprintf('%s  pro_measured %d  GL2 floor %.2f  GL2 edge %s\n', fn, FCSfileinfo.pro_measured(k), info.floor.GL2, mat2str(struct2array(info.trigger_edges), 3))
m = info.comp; disp(m(:, {'SSC' 'BL3' 'GL2' 'PEchl' 'GL2GL3' 'chlSSC' 'label' 'weight'}))
fprintf('syn ref size %.2f BL3 %.2f\n', r.size, r.BL3)
hc = t.has_chl; gl = nan(size(old)); gl(hc) = t.gmm_label; go = false(size(old)); go(hc) = t.gmm_outlier;
lp = old == 7 & new == 0;
fprintf('\nold Pro %d, v2 Pro %d, old Pro dropped %d, old Pro to Syn %d\n', sum(old == 7), sum(new == 7), sum(lp), sum(old == 7 & new == 2))
fprintf(' dropped: not has_chl %d, GMM outlier %d, GMM label 0/1/2/3/7 %s\n', sum(lp & ~hc), sum(lp & go), mat2str(arrayfun(@(x) sum(lp & hc & ~go & gl == x), [0 1 2 3 7])))
fprintf(' dropped: detritus %d, Syn check %d, Euk min size %d, Euk chl:size %d, PE-euk step %d\n', ...
    sum(lp & t.after_gmm ~= 0 & t.after_detritus == 0), sum(lp & t.after_detritus == 2 & t.after_syncheck == 0), ...
    sum(lp & t.after_syncheck == 1 & F.size < r.size - sig.euk_min_size_below_syn), ...
    sum(lp & t.after_syncheck == 1 & F.size >= r.size - sig.euk_min_size_below_syn & t.after_peeuk == 0), sum(lp & t.after_syncheck == 3 & t.after_peeuk == 0))
fprintf(' dropped size p5/50/95 %s, kept Pro %s; euk_min size cut %.2f; pro_max_SSC %.2f\n', mat2str(prctile(F.size(lp), [5 50 95]), 3), ...
    mat2str(prctile(F.size(old == 7 & new == 7), [5 50 95]), 3), r.size - sig.euk_min_size_below_syn, sig.pro_max_SSC)
smalleuk = new == 1 & F.size < prctile(F.size(new == 1), 25);
fprintf('\nGL2 p5/50/95: v2 Pro %s, old Pro %s, small Euks %s, all Euks %s\n', mat2str(prctile(F.GL2(new == 7), [5 50 95]), 3), ...
    mat2str(prctile(F.GL2(old == 7), [5 50 95]), 3), mat2str(prctile(F.GL2(smalleuk), [5 50 95]), 3), mat2str(prctile(F.GL2(new == 1), [5 50 95]), 3))
fprintf('PEchl p5/50/95: v2 Pro %s, old Pro %s, Euks %s\n', mat2str(prctile(F.PEchl(new == 7), [5 50 95]), 3), ...
    mat2str(prctile(F.PEchl(old == 7), [5 50 95]), 3), mat2str(prctile(F.PEchl(new == 1), [5 50 95]), 3))
end
