function q = review_qc_hbac(cls, info, cfg)
% function q = review_qc_hbac(cls, info, cfg)
% Per-file QC for one SYBR-stained bacteria run (assign_class_hbac), for the
% review viewer in mode 'hbac': counts per class with bacteria split into LNA (n91)
% and HNA (n92), the HNA fraction and how the split was made, the bacteria core
% and spillover lines, and q.lines, the class boundaries the viewer draws on each
% panel of review_config_hbac (struct array: v panel, a intercept, b slope, x for
% a vertical line, col, dash, label).
%
% cfg - discrete_config (the cfg.hbac thresholds used)
%
% Nov 2026, preserved-sample review (SYBR runs)

h = cfg.hbac;
cls = double(cls(:));
bt = zeros(size(cls)); if isfield(info, 'bact_type'), bt = double(info.bact_type(:)); end
q = struct();
q.n_events = numel(cls);
for k = [0 1 2 7 9], q.(sprintf('n%d', k)) = sum(cls == k); end
q.n91 = sum(cls == 9 & bt ~= 2); q.n92 = sum(cls == 9 & bt == 2);
hn = struct('cut', NaN, 'slope', NaN, 'method', '', 'dip', NaN, 'valley', NaN);
if isfield(info, 'hna'), hn = info.hna; end
q.hna_frac = q.n92 / max(1, q.n9);
q.hna_method = char(fo(hn, 'method', ''));
q.hna_valley = double(strcmp(q.hna_method, 'valley'));
q.hna_cut = fo(hn, 'cut', NaN); q.hna_slope = fo(hn, 'slope', NaN); q.hna_dip = fo(hn, 'dip', NaN);
q.n_core = fo(info, 'n_core', NaN);
cx = [NaN NaN NaN]; px = cx;
if isfield(info, 'chlx'), cx = info.chlx; end
if isfield(info, 'pex'), px = info.pex; end
q.chlx_off = cx(1); q.chlx_sd = cx(2); q.chlx_cut = cx(3);
q.pex_off = px(1); q.pex_sd = px(2); q.pex_cut = px(3);
e = fo(info, 'edges', struct('BL1', -Inf, 'SSC', -Inf));
q.trig = struct();
for f = {'BL1' 'SSC'}
    if isfield(e, f{1}) && isfinite(e.(f{1})), q.trig.(f{1}) = e.(f{1}); end
end

% class boundaries per panel (review_config_hbac order)
L = struct('v', {}, 'a', {}, 'b', {}, 'x', {}, 'col', {}, 'dash', {}, 'label', {});
add = @(L, v, a, b, x, col, dash, lab) [L struct('v', v, 'a', a, 'b', b, 'x', x, 'col', col, 'dash', dash, 'label', lab)];
bc = [.45 .2 .6]; hc = [.75 0 .2]; ec = [0 .45 0]; sc = [.85 .45 0]; pc = [0 .25 .85]; k0 = [0 0 0];
if isfield(e, 'BL1') && isfinite(e.BL1)
    L = add(L, 1, e.BL1 + h.bl1_edge_margin, 0, NaN, k0, '1 2', 'SYBR edge + margin: bacteria above');
end
L = add(L, 1, h.stain_min, 1, NaN, bc, '5 3', 'stain line BL1 = SSC + stain_min: bacteria above');
L = add(L, 1, NaN, NaN, h.bact_max_SSC, bc, '2 2', 'bacteria size limit');
if isfinite(q.hna_cut)
    dash = '6 2'; if ~q.hna_valley, dash = '1 3'; end
    L = add(L, 1, q.hna_cut + q.hna_slope * 2, q.hna_slope, NaN, hc, dash, ['HNA / LNA cut (' q.hna_method ')']);
    L = add(L, 2, NaN, NaN, q.hna_cut, hc, dash, ['HNA / LNA cut (' q.hna_method ')']);
end
if isfinite(q.chlx_cut)
    L = add(L, 3, q.chlx_off, 1, NaN, bc, '1 2', 'bacteria chl spillover line');
    L = add(L, 3, q.chlx_cut, 1, NaN, ec, '5 3', 'chl excess cut: phytoplankton above');
end
if isfinite(q.pex_cut)
    L = add(L, 4, q.pex_off, 1, NaN, bc, '1 2', 'bacteria PE spillover line');
    L = add(L, 4, q.pex_cut, 1, NaN, sc, '5 3', 'PE excess cut');
end
L = add(L, 5, NaN, NaN, h.pro_max_SSC, pc, '5 3', 'Pro (left) / Euk (right) size split');
L = add(L, 5, h.min_chlSSC, 1, NaN, ec, '2 2', 'chl:size minimum (Euk, Syn)');
L = add(L, 5, NaN, NaN, h.bact_max_SSC, bc, '2 2', 'bacteria size limit');
L = add(L, 6, h.syn_min_PEchl, 1, NaN, sc, '5 3', 'Syn PE:chl minimum');
% panel 7: PE:chl relative to the bacteria's (0 = the SYBR spectrum, where the unstained-particle band sits)
sp = q.pex_off - q.chlx_off;
if isfinite(sp)
    L = add(L, 7, 0, 0, NaN, bc, '1 2', 'bacteria PE:chl (SYBR spectrum)');
    L = add(L, 7, max(h.syn_min_PEchl - sp, fo(h, 'syn_above_sybr', -Inf)), 0, NaN, sc, '5 3', 'Syn PE:chl minimum (absolute and above the bacteria''s)');
    if h.euk_below_sybr == h.euk_stained_below_sybr
        L = add(L, 7, -h.euk_below_sybr, 0, NaN, ec, '5 3', 'Euk, stained Euk and large Pro: below this');
    else
        L = add(L, 7, -h.euk_below_sybr, 0, NaN, ec, '5 3', 'Euk: below this');
        L = add(L, 7, -h.euk_stained_below_sybr, 0, NaN, ec, '1 3', 'stained Euk and large Pro: below this');
    end
end
L = add(L, 7, NaN, NaN, h.pro_max_SSC, pc, '5 3', 'Pro (left) / Euk (right) size split');
L = add(L, 7, NaN, NaN, h.bact_max_SSC, bc, '2 2', 'bacteria size limit');
% panel 8: large deep Pro limits (Pro-like "Euks" become Pro where they dominate)
if isfield(h, 'lp')
    L = add(L, 8, NaN, NaN, h.lp.max_chlSSC, pc, '5 3', 'large Pro: chl:size below this');
    L = add(L, 8, h.lp.max_BL3, 0, NaN, pc, '5 3', 'large Pro: chl below this');
end
L = add(L, 8, NaN, NaN, h.min_chlSSC, ec, '2 2', 'chl:size minimum (Euk, Syn)');
q.lines = L;
end

function v = fo(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
    if isnumeric(v) || islogical(v), v = double(v(1)); end
else
    v = d;
end
end
