function debug_euk(cruisedir, fn, min_size)
% function debug_euk(cruisedir, fn, min_size)
% Trace old-class Euks that v2 leaves unassigned (optionally only those with
% size > min_size): which step removes them.
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
F = info.features; t = info.trace;
if nargin < 3, min_size = -Inf; end
lost = old == 1 & new == 0 & F.size > min_size;
kept = old == 1 & new == 1 & F.size > min_size;
m = info.comp; disp(m(:, {'SSC' 'BL3' 'GL2' 'PEchl' 'chlSSC' 'label' 'weight'}))
hc = t.has_chl; gl = nan(size(old)); gl(hc) = t.gmm_label; go = false(size(old)); go(hc) = t.gmm_outlier;
fprintf('old Euk lost %d (kept %d) with size > %.1f\n', sum(lost), sum(kept), min_size)
fprintf(' not has_chl %d; GMM outlier %d; GMM label 0/1/2/3/7 %s\n', sum(lost & ~hc), sum(lost & go), mat2str(arrayfun(@(x) sum(lost & hc & ~go & gl == x), [0 1 2 3 7])))
fprintf(' detritus chl:size rule %d; Syn check %d; Euk chl:size band %d; PE-euk step %d\n', ...
    sum(lost & t.after_gmm ~= 0 & t.after_detritus == 0), sum(lost & t.after_detritus == 2 & t.after_syncheck == 0), ...
    sum(lost & t.after_syncheck == 1 & t.after_peeuk == 0), sum(lost & t.after_syncheck == 3 & t.after_peeuk == 0))
fprintf(' lost: size %s BL3 %s chlSSC %s\n kept: size %s BL3 %s chlSSC %s\n', mat2str(prctile(F.size(lost), [10 50 90]), 3), ...
    mat2str(prctile(F.BL3(lost), [10 50 90]), 3), mat2str(prctile(F.chlSSC(lost), [10 50 90]), 3), mat2str(prctile(F.size(kept), [10 50 90]), 3), ...
    mat2str(prctile(F.BL3(kept), [10 50 90]), 3), mat2str(prctile(F.chlSSC(kept), [10 50 90]), 3))
e1 = t.after_syncheck == 1;
fprintf(' Euk chlSSC median (cut = median - %.2f): %.2f -> cut %.2f\n', sig.euk_chlSSC_tol, median(F.chlSSC(e1)), median(F.chlSSC(e1)) - sig.euk_chlSSC_tol)
% Mahalanobis outliers: which dimension
if any(lost & go)
    X = [F.size F.BL3 F.PEchl F.GL2GL3];
    fprintf(' outliers: PEchl %s GL2GL3 %s\n', mat2str(prctile(F.PEchl(lost & go), [10 50 90]), 3), mat2str(prctile(F.GL2GL3(lost & go), [10 50 90]), 3))
    fprintf(' kept Euks: PEchl %s GL2GL3 %s\n', mat2str(prctile(F.PEchl(kept), [10 50 90]), 3), mat2str(prctile(F.GL2GL3(kept), [10 50 90]), 3))
end
end
