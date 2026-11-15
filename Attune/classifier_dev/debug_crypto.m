function debug_crypto(cruisedir, fn)
% function debug_crypto(cruisedir, fn)
% Trace old-class PE euks (3, 4) that v2 does not call crypto (3/4): which step
% and which variables (PE, spectrum, saturation) put them elsewhere.
addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = fullfile('\\sosiknas1\lab_data\Attune\cruise_data', cruisedir);
bp = fullfile(root, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
cfg = classifier_config(); sig = cfg.sig;
k = find(strcmp(FCSfileinfo.fcslist, fn));
[d, h] = fca_readfcs(fullfile(root, 'FCS', fn));
ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
pm = 0; if ismember('pro_measured', FCSfileinfo.Properties.VariableNames), pm = FCSfileinfo.pro_measured(k); end
[new, info] = assign_class_v2(d, h, ref, cfg, pm);
c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class'); old = c.class(:, 1);
F = info.features; t = info.trace; pn = {h.par.name};
sat = @(ch) d(:, strcmp(pn, [ch '-H'])) >= sig.sat_frac * max(d(:, strcmp(pn, [ch '-H'])));
satG2 = sat('GL2'); satB2 = sat('BL2'); satG3 = sat('GL3');
disp(info.crypto_split)
for o = [3 4]
    g = old == o;
    fprintf('\nold %d (%d): v2 class 0/1/2/3/4/5/8: %s\n', o, sum(g), mat2str(arrayfun(@(x) sum(g & new == x), [0 1 2 3 4 5 8])))
    lost = g & ~ismember(new, [3 4]);
    fprintf('  not crypto in v2: %d; after GMM class 0/1/2/3: %s; after PE-euk step 0/1/2/3: %s\n', sum(lost), ...
        mat2str(arrayfun(@(x) sum(lost & t.after_gmm == x), [0 1 2 3])), mat2str(arrayfun(@(x) sum(lost & t.after_peeuk == x), [0 1 2 3])))
    fprintf('  lost: GL2sat %d, BL2sat %d, GL3sat %d; size %s; GL2 %s; BL2 %s; GL3 %s; chlSSC %s; PEchl %s\n', sum(lost & satG2), sum(lost & satB2), sum(lost & satG3), ...
        mat2str(prctile(F.size(lost), [10 50 90]), 2), mat2str(prctile(F.GL2(lost), [10 50 90]), 2), mat2str(prctile(F.BL2(lost), [10 50 90]), 2), ...
        mat2str(prctile(F.GL3(lost), [10 50 90]), 2), mat2str(prctile(F.chlSSC(lost), [10 50 90]), 2), mat2str(prctile(F.PEchl(lost), [10 50 90]), 2))
end
k4 = new == 4;
fprintf('\nv2 bright crypto: GL2sat %d BL2sat %d GL3sat %d; GL2 %s BL2 %s GL3 %s\n', sum(k4 & satG2), sum(k4 & satB2), sum(k4 & satG3), ...
    mat2str(prctile(F.GL2(k4), [10 50 90]), 2), mat2str(prctile(F.BL2(k4), [10 50 90]), 2), mat2str(prctile(F.GL3(k4), [10 50 90]), 2))
fprintf('pe_min_GL2 %.2f, peeuk_min_GL2 %.2f, pe_min_PEchl %.2f, syn size %.2f\n', sig.pe_min_GL2, sig.peeuk_min_GL2, sig.pe_min_PEchl, info.syn_ref.size)
end
