function [class, info] = assign_class_hbac(fcsdat, fcshdr, ref, cfg)
% function [class, info] = assign_class_hbac(fcsdat, fcshdr, ref, cfg)
% Classify one SYBR-stained bacteria run (hbac_SYBR_SSC, preserved samples).
% Shared class codes: 0 unassigned, 1 Euk, 2 Syn, 7 Pro, 9 heterotrophic bacteria.
% (Pro is 7 here too, so bacteria are heterotrophs only; "bacteria including
% Pro" is 7 + 9 in the table step.)
%
% Bead-normalized features (attune_features, same-day lab beads). SYBR (BL1)
% spills into the chl (BL3) and PE (GL2) channels along fixed lines, so per
% file:
%   chlx = BL3 - BL1, pex = GL2 - BL1, centred on the stained bacteria
%   (median and robust SD of the core: stained, small, above the BL1 edge)
%   chl excess = chlx above the bacteria by more than cfg.hbac.phyto_sd SDs
%   (and at least cfg.hbac.phyto_min), and chl:size (BL3 - SSC) above
%   cfg.hbac.min_chlSSC (Pro-sized: cfg.hbac.min_chlSSC_pro); PE excess likewise for pex
% Then:
%   Syn (2)  PE excess with Syn-like PE:chl (GL2 - BL3 > cfg.hbac.syn_min_PEchl)
%            and chl:size above cfg.hbac.min_chlSSC, or chl excess with that
%            PE:chl; SSC below cfg.hbac.syn_max_SSC (large PE-signal particles
%            with little chl, EN695 C03N06, are not Syn)
%   Euk (1)  chl excess, not Syn, SSC above cfg.hbac.pro_max_SSC
%   Pro (7)  chl excess, not Syn, SSC at or below cfg.hbac.pro_max_SSC
%   bacteria (9) the rest that is stained (BL1 - SSC > cfg.hbac.stain_min),
%            above the BL1 trigger edge + cfg.hbac.bl1_edge_margin and no larger
%            than cfg.hbac.bact_max_SSC
% Hand-gated SYBR runs, 24 runs on 6 cruises (examine_hbac, Nov 2026):
%   bacteria chlx 0.35 +- 0.14, pex -0.52 +- 0.18 in every cruise;
%   BL1 - SSC 1st percentile ~0; SSC 5-95% -2.6 to -1.5
%
% HNA / LNA (hbac_hna_split, cfg.hbac.hna): per-event info.bact_type, uint8,
%   1 LNA, 2 HNA, 0 not bacteria. Cut on y = BL1 - slope * (SSC + 2) at the
%   run's SYBR density valley when it is clear (info.hna.method 'valley'),
%   else a fixed cut ('fixed': HNA and LNA merged, mostly summer surface
%   samples; the split there is nominal).
%
% info - core offsets and SDs, thresholds used, trigger edges, features,
%        bact_type, bact_type_names, hna (cut, slope, method, dip, counts)
%
% Nov 2026, discrete-sample classifier development

h = cfg.hbac;
[F, finfo] = attune_features(fcsdat, fcshdr, ref);
n = height(F);
class = zeros(n, 1);
e = struct('BL1', -Inf, 'SSC', -Inf);
trigs = [finfo.trigger finfo.trigger2]; levels = [finfo.trigger_level finfo.trigger2_level];
for t = 1:2
    ch = char(regexp(trigs(t), '(?<=_)\w+$', 'match', 'once'));
    if isfield(e, ch), e.(ch) = finfo.floor.(ch) + log10(levels(t)) - log10(cfg.floor_raw); end
end
chlx = F.BL3 - F.BL1; pex = F.GL2 - F.BL1; stain = F.BL1 - F.SSC; PEchl = F.GL2 - F.BL3;
core = stain > h.stain_min & F.SSC < h.core_max_SSC & F.BL1 > e.BL1 + h.bl1_edge_margin & isfinite(chlx) & isfinite(pex);
info = struct('n_core', sum(core), 'edges', e, 'features', F);
if sum(core) < h.min_core
    oc = h.default_chlx; sc = h.default_sd; op = h.default_pex; sp = h.default_sd;
else
    oc = median(chlx(core)); sc = max(1.4826 * mad(chlx(core), 1), h.min_sd);
    op = median(pex(core)); sp = max(1.4826 * mad(pex(core), 1), h.min_sd);
end
cut_c = oc + max(h.phyto_min, h.phyto_sd * sc);
cut_p = op + max(h.phyto_min, h.phyto_sd * sp);
chl_ex = chlx > cut_c; pe_ex = pex > cut_p;
% chl excess must also be real chl for the particle's size: large unstained
% particles have low BL1, so chlx is high without chl (v1: Syn 2.2x, Euk 3.3x hand)
chlSSC = F.BL3 - F.SSC;
chl_ex = chl_ex & (chlSSC > h.min_chlSSC | (F.SSC <= h.pro_max_SSC & chlSSC > h.min_chlSSC_pro));
synlike = PEchl > h.syn_min_PEchl;
syn = ((pe_ex & chlSSC > h.min_chlSSC) | chl_ex) & synlike & F.SSC < h.syn_max_SSC & (pe_ex | F.SSC > h.pro_max_SSC);
phyto = chl_ex | pe_ex;
class(syn) = 2;
class(chl_ex & ~syn & F.SSC > h.pro_max_SSC) = 1;
class(chl_ex & ~syn & F.SSC <= h.pro_max_SSC) = 7;
bact = ~phyto & stain > h.stain_min & F.BL1 > e.BL1 + h.bl1_edge_margin & F.SSC <= h.bact_max_SSC;
class(bact) = 9;
info.chlx = [oc sc cut_c]; info.pex = [op sp cut_p];
% HNA / LNA subtype of the bacteria (like info.syn_type; class 9 unchanged)
[info.bact_type, info.hna] = hbac_hna_split(F.BL1, F.SSC, class == 9, h.hna);
info.bact_type_names = info.hna.names;
info.counts = struct('syn', sum(class == 2), 'euk', sum(class == 1), 'pro', sum(class == 7), 'bact', sum(class == 9), ...
    'lna', info.hna.n_lna, 'hna', info.hna.n_hna);
end
