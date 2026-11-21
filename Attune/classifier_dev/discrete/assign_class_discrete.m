function [class, info] = assign_class_discrete(fcsdat, fcshdr, ref, cfg, runtype)
% function [class, info] = assign_class_discrete(fcsdat, fcshdr, ref, cfg, runtype)
% Classify one preserved (discrete) phytoplankton run: assign_class_v2 with the
% preserved-sample settings (discrete_config) and run-type adjustments.
% Class codes are the shared scheme (0 unassigned, 1 Euk, 2 Syn, 3 dim crypto,
% 4 bright crypto, 5 Euk with Syn, 7 Pro, 8 Syn on particle).
%
% ref     - discrete_bead_reference (same-day lab beads)
% cfg     - discrete_config (cfg.cruise, cfg.gl1_is_od2_scatter set by the caller)
% runtype - 'pro' (phyto_CHL_SSC_pro), 'CHL' (phyto_CHL_SSC) or 'PE' (phyto_PE_SSC);
%           logical true/false is taken as 'pro'/'CHL'
%
% Run-type adjustments:
%  - Pro runs keep the live detritus chl:size cut (cfg.sig.detritus_max_chlSSC_pro):
%    Pro chl:size is 0.27-1.32 (5-95%) in preserved samples, so the preserved cut
%    (0.6) named whole Pro components detritus (EN706 Pro runs: 26,000 hand Pro -> 0)
%  - CHL runs also assign Pro (cfg.sig.pro_in_chl_runs). Their chl trigger usually
%    cuts off most of the Pro cluster (pro_in_chl_runs: median 5% of a Pro run's
%    Pro recordable), but bright (e.g. deep) Pro can be complete. info.pro
%    reports how close the Pro cluster sits to the chl trigger edge, and
%    info.pro.lower_bound is true where the count is only a lower bound.
%
%  - Large, bright Pro (deep samples; larger than Syn) are named Euk by v2's
%    Pro size rule; a per-file 2-group fit on the small no-PE Euks in
%    [chl:size, GL2:GL3] moves a clearly separate lower group to Pro (large_pro)
%
% info.large_pro - accepted, n moved to Pro, group means, separation
% info.pro - n, peak_above_edge (KDE mode of Pro BL3 minus the chl trigger
%            edge), frac_near_edge (Pro within cfg.sig.pro_edge_band of the edge),
%            lower_bound (false only when the cluster is shown complete: at least
%            cfg.sig.pro_edge_min_n Pro clear of the edge; meaningful for CHL runs only,
%            since dim Pro sit at the edge in Pro runs too)
%
% Nov 2026, discrete-sample classifier development

if islogical(runtype) || isnumeric(runtype)
    if runtype, runtype = 'pro'; else, runtype = 'CHL'; end
end
runtype = char(runtype);
pro_run = strcmp(runtype, 'pro');
assign_pro = pro_run || (strcmp(runtype, 'CHL') && cfg.sig.pro_in_chl_runs);
if pro_run && isfield(cfg.sig, 'detritus_max_chlSSC_pro')
    cfg.sig.detritus_max_chlSSC = cfg.sig.detritus_max_chlSSC_pro;
end
[cfg.sig.pe_min_PEchl, info_pe] = euk_autofluor_cut(fcsdat, fcshdr, ref, cfg);
if ~pro_run, cfg.sig.pro_pick_min_n = Inf; end %v2's event-level Pro pick: Pro runs only
[class, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, assign_pro);
info.euk_autofluor = info_pe;
info.pick_trim = struct('applied', false, 'n_before', 0, 'n_after', 0);
if pro_run && isfield(info, 'pro_pick') && info.pro_pick.used
    [class, info.pick_trim] = pro_pick_trim(class, info, cfg);
end
info.runtype = runtype;
info.large_pro = struct('accepted', false, 'n', 0);
info.pro_swap = struct('applied', false, 'n_to_pro', 0, 'n_to_0', 0);
if pro_run
    [class, info.pro_swap] = pro_component_swap(class, info, cfg);
end
if assign_pro && ~info.pro_swap.applied %(the swap has already found the file's Pro cluster)
    [class, info.large_pro] = large_pro(class, info, cfg, pro_run);
end
% v2's wide-pulse rule (Pro-setting files, Nov 2026): events whose BL3 pulse is
% wider than the file's Syn-referenced A/H (info.trace.wide_pulse) are not Pro; the
% steps above can make Pro of events v2 left out, so the same exclusion applies here
info.n_pro_wide_wrapper = 0;
if isfield(info, 'trace') && isfield(info.trace, 'pro_wide_pulse')
    w = class == 7 & info.trace.pro_wide_pulse; %v2's adaptive Pro cut (info.pro_AH_max)
    class(w) = 0; info.n_pro_wide_wrapper = sum(w);
elseif isfield(info, 'trace') && isfield(info.trace, 'wide_pulse')
    w = class == 7 & info.trace.wide_pulse;
    class(w) = 0; info.n_pro_wide_wrapper = sum(w);
end
info.pro = pro_completeness(info, class, cfg);
info.discrete_cfg = cfg.sig;
end

function [class, s] = large_pro(class, info, cfg, pro_run)
% Large, bright Pro (deep samples, e.g. EN688 C1N18: size ~ -1.0, larger than Syn)
% fall outside v2's Pro size rule (pro_max_SSC) and are named Euk. Among the
% no-PE "Euks" (class 1) fit 1-3 groups in [chl:size, GL2:GL3] (BIC); the
% lowest chl:size group becomes Pro only if it is lower than the other Euk
% groups by at least lp_min_d_chlSSC (chl:size) and lp_min_d_GL2GL3, lp_min_sep
% of its SDs away, small (median size < Syn + lp_size_above_syn) and dim
% (median BL3 < lp_max_BL3); its events under those size and BL3 limits move.
% Hand-gated preserved Pro vs Euks at size -1.3 to -0.6: chl:size 1.09 vs 1.70,
% GL2:GL3 -0.10 vs 0.34 (single events overlap; GL3 near its floor).
g = cfg.sig;
s = struct('accepted', false, 'n', 0, 'method', '', 'mu', [], 'sep', NaN, 'ncand', 0, 'lda_n', 0, 'lda_median', NaN, 'lda_frac', NaN);
F = info.features;
if sum(class == 2) >= 20, syn_size = median(F.size(class == 2)); else, syn_size = g.syn_size_default; end
% the fixed-boundary fallback (large_pro_lda; its candidates need chl:size >= lp_min_chlSSC,
% since the boundary favours low chl:size and took detritus in EN655 C15N09): Pro runs, and CHL runs where v2 has
% already named a Pro component (v2 judges Pro against Syn since Nov 2026 and can name
% a slice of a large Pro cluster Pro, the rest Euk; the pool is then almost all Pro
% and the group fit finds no split: AR38 C11N13, EN688 C8N19 CHL runs)
use_lda = pro_run || (isfield(info, 'comp') && istable(info.comp) && any(strcmp(info.comp.Properties.VariableNames, 'label')) && any(info.comp.label == 7));
% v2 now unassigns small, dim no-PE "Euks" as Pro-like (euk_prolike, Nov 2026);
% those are exactly the large-Pro candidates, so they rejoin the Euk pool here
euk = class == 1;
if isfield(info, 'trace') && isfield(info.trace, 'after_syncheck')
    euk = euk | (class == 0 & info.trace.after_syncheck == 1); %(any later v2 step that unassigned a no-PE Euk: AR92 C20N07 CHL, 399 of 995 hand Pro)
    % ... and unassigned no-PE chlorophyll events that v2 never made Euk: since Nov 2026
    % v2 leaves components larger than pro_max_size_abs (-1.15) unnamed rather than Euk
    % (HB1907 C15N03 Pro run: 2415 of 4544 hand Pro, size ~ -1.07)
    euk = euk | (class == 0 & info.trace.has_chl & F.PEchl < g.pe_min_PEchl & F.chlSSC >= g.lp_min_chlSSC);
end
if sum(euk) < g.lp_min_n && use_lda, [class, s] = large_pro_lda(class, F, syn_size, g, s, euk); return, end
% Fit on all the file's Euks (the larger Euks are the reference group: in deep
% samples the small dim "Euks" can be almost all Pro); K = 1-3 by BIC
cand = find(euk & isfinite(F.GL2GL3) & isfinite(F.chlSSC));
s.ncand = numel(cand);
if numel(cand) < g.lp_min_n, return, end
X = [F.chlSSC(cand) F.GL2GL3(cand)];
if isfield(cfg, 'seed'), rng(cfg.seed, 'twister'), else, rng(1, 'twister'), end
opts = statset('MaxIter', 300);
best = []; bic = Inf;
for K = 1:3
    try
        gk = fitgmdist(X, K, 'RegularizationValue', 1e-3, 'Replicates', 3, 'Options', opts);
    catch
        continue
    end
    if gk.BIC < bic - 10 * (K > 1), best = gk; bic = gk.BIC; end
end
if isempty(best) || best.NumComponents < 2
    if use_lda, [class, s] = large_pro_lda(class, F, syn_size, g, s, euk); end
    return
end
p = posterior(best, X);
[~, lo] = min(best.mu(:, 1));
others = setdiff(1:best.NumComponents, lo);
w = best.ComponentProportion(others) / sum(best.ComponentProportion(others));
mo = w * best.mu(others, :); %weighted mean of the other Euk groups
d = mo - best.mu(lo, :);
S = best.Sigma(:, :, lo);
s.sep = sqrt(d / S * d');
inlo = p(:, lo) > 0.5;
s.mu = [best.mu(lo, :) median(F.size(cand(inlo))) median(F.BL3(cand(inlo))); mo NaN NaN];
ok = d(1) >= g.lp_min_d_chlSSC && d(2) >= g.lp_min_d_GL2GL3 && s.sep >= g.lp_min_sep && best.mu(lo, 1) >= g.lp_min_chlSSC ...
    && s.mu(1, 3) < syn_size + g.lp_size_above_syn && s.mu(1, 4) < g.lp_max_BL3;
pro = cand(inlo & F.size(cand) < syn_size + g.lp_size_above_syn & F.BL3(cand) < g.lp_max_BL3);
if ~ok || numel(pro) < g.lp_min_group_n
    if use_lda, [class, s] = large_pro_lda(class, F, syn_size, g, s, euk); end
    return
end
class(pro) = 7;
s.accepted = true; s.n = numel(pro); s.method = 'gmm';
end

function [class, s] = large_pro_lda(class, F, syn_size, g, s, euk)
% Pro runs: Pro can outnumber the Euks so much that the group fit splits the Pro
% instead of separating them from the Euks (AR38 C11N13, EN688 C11N18: 0 Pro vs
% hand 4232 and 1942). Fixed discriminant from hand-gated preserved Pro vs Euks
% (size -1.6 to -0.6, BL3 < 1; 90% of Pro, 6% of Euks on the Pro side). The
% small, dim no-PE "Euks" on the Pro side move to Pro if there are at least
% lp_min_group_n and their median score is at least lp_lda_min_median: files
% with hand Pro 1.3-9.5 (24 of 27 above 2), Euk-only files up to 1.8; or, where
% at least lp_lda_dominant(1) of the candidates are on the Pro side (Euk-only
% files up to 0.69; AR38 C11N13, EN688 C11N18 pro runs 0.84, 0.83), median at
% least lp_lda_dominant(2)
cand = find(euk & isfinite(F.GL2GL3) & isfinite(F.chlSSC) ...
    & F.size < syn_size + g.lp_size_above_syn & F.BL3 < g.lp_max_BL3 & F.chlSSC >= g.lp_min_chlSSC);
sc = g.lp_lda(1) * F.chlSSC(cand) + g.lp_lda(2) * F.GL2GL3(cand) + g.lp_lda(3);
side = sc > 0;
s.lda_n = sum(side); s.lda_median = median(sc(side)); s.lda_frac = mean(side);
dominant = s.lda_frac >= g.lp_lda_dominant(1) && s.lda_median >= g.lp_lda_dominant(2);
if s.lda_n < g.lp_min_group_n || ~(s.lda_median >= g.lp_lda_min_median || dominant), return, end
class(cand(side)) = 7;
s.accepted = true; s.n = s.lda_n; s.method = 'lda';
end

function [class, s] = pro_component_swap(class, info, cfg)
% Pro runs: v2 names a Pro-like component Euk when it is as large as Syn (size >=
% Syn - 0.1). Preserved Pro can be (EN655 C15N09: Pro component size -1.25, Syn
% -1.24), so the real, compact Pro cluster is named Euk (then unassigned as
% Pro-like), while dimmer groups below it are named Pro: the pile of small faint
% particles on the chl trigger floor (BL3 -1.03, edge -1.15) and a broad faint
% group between floor and Pro (BL3 -0.55). Where a compact no-PE component, at
% least pc_min_above_edge above the edge, no larger than Syn + pc_max_size_above_syn,
% with chl:size below Syn's by pc_chlSSC_below_syn(1)-(2) (Pro 0.28-0.86 below)
% and on the Pro side of the Pro-vs-Euk spectral boundary, is named Euk and every
% component named Pro is at least pc_min_below_cand dimmer, that component is the
% Pro and the dimmer Pro components are unassigned (events within its spread stay).
g = cfg.sig;
s = struct('applied', false, 'n_to_pro', 0, 'n_to_0', 0, 'pro_comp', [], 'floor_comp', []);
if ~isfield(info, 'gmm') || isempty(info.gmm) || ~isfield(info.trigger_edges, 'BL3'), return, end
mu = info.comp; gm = info.gmm; edge = info.trigger_edges.BL3;
if ~any(mu.label == 2), return, end
w = mu.weight(mu.label == 2) / sum(mu.weight(mu.label == 2));
syn_size = w' * mu.SSC(mu.label == 2); syn_chl = w' * mu.chlSSC(mu.label == 2);
sdBL3 = sqrt(squeeze(gm.Sigma(2, 2, :)));
if any(strcmp(mu.Properties.VariableNames, 'sd_BL3')), sdBL3 = mu.sd_BL3; end %v2's robust (narrow-pulse) SD where given
cand = mu.label ~= 2 & mu.PEchl < g.pe_min_PEchl & mu.BL3 > edge + g.pc_min_above_edge & mu.BL3 < g.pro_max_BL3 ...
    & mu.SSC < syn_size + g.pc_max_size_above_syn & sdBL3(:) <= g.pro_max_sd_BL3 ...
    & mu.chlSSC < syn_chl - g.pc_chlSSC_below_syn(1) & mu.chlSSC > syn_chl - g.pc_chlSSC_below_syn(2) ...
    & g.lp_lda(1) * mu.chlSSC + g.lp_lda(2) * mu.GL2GL3 + g.lp_lda(3) > 0; %Pro side of the Pro-vs-Euk boundary (large_pro_lda)
if ~any(cand), return, end
% every component named Pro must sit clearly below the candidate (otherwise v2's naming stands)
% (the candidate may itself be named Pro: since Nov 2026 v2 can name it Pro and still
% name the floor pile Pro too, EN655 C15N09). Every other component named Pro must sit
% clearly below the candidate, and at least one must exist (otherwise v2's naming stands)
on_floor = mu.label == 7 & ~cand; %(the components replaced; name kept for info)
if ~any(on_floor) || any(on_floor & mu.BL3 > max(mu.BL3(cand)) - g.pc_min_below_cand), return, end
F = info.features; hc = info.trace.has_chl;
comp = zeros(size(class));
comp(hc) = cluster(gm, [F.size(hc) F.BL3(hc) F.PEchl(hc) F.GL2GL3(hc)]);
to_pro = ismember(comp, find(cand)) & ismember(class, [0 1 7]) & F.PEchl < g.pe_min_PEchl;
% events v2 made Pro outside the new Pro component (the floor pile, and catch-all
% events given Pro by their closeness to it) stay Pro only within the new Pro
% component's own spread (cfg.outlier_p, as v2's outliers)
X = [F.size F.BL3 F.PEchl F.GL2GL3];
ic = find(cand); d2 = inf(size(class));
for k = ic(:)'
    dx = X - gm.mu(k, :);
    d2 = min(d2, sum((dx / gm.Sigma(:, :, k)) .* dx, 2));
end
to_0 = class == 7 & ~to_pro & ~(d2 < chi2inv(cfg.outlier_p, size(X, 2)));
% ... and with chl:size near the new Pro's (the spread test has no chl:size: EN655
% C15N09 kept 337 low-chl events, chl:size median 0.32 vs Pro 1.10)
cs = F.chlSSC(to_pro);
if numel(cs) >= 20
    to_0 = to_0 | (class == 7 & ~to_pro & F.chlSSC < median(cs) - g.pro_chlSSC_sd * 1.4826 * mad(cs, 1));
end
class(to_0) = 0; class(to_pro) = 7;
s = struct('applied', true, 'n_to_pro', sum(to_pro), 'n_to_0', sum(to_0), 'pro_comp', find(cand)', 'floor_comp', find(on_floor)');
end

function [class, s] = pro_pick_trim(class, info, cfg)
% v2's event-level Pro pick (Pro runs where the mixture names no Pro, e.g. winter
% Pro on the chl trigger floor, EN727 C19N13) is followed by v2's catch-all, which
% gives Pro to leftover cell-like events near the picked ones: 523 Pro in all vs
% 301 picked and a hand-gated compact cluster of 191 (late window). The catch-all
% additions were brighter (BL3 median -0.63 vs -0.99) and spread in chl:size and
% PE:chl, so with the pick only the picked events stay Pro
s = struct('applied', false, 'n_before', sum(class == 7), 'n_after', sum(class == 7));
if ~isfield(info.trace, 'catchall'), return, end
class(class == 7 & info.trace.catchall) = 0;
s = struct('applied', true, 'n_before', s.n_before, 'n_after', sum(class == 7));
end

function p = pro_completeness(info, class, cfg)
% where the Pro cluster sits relative to the chl (BL3) trigger edge
% lower_bound is false only where the cluster is shown complete; fewer than
% pro_edge_min_n Pro cannot show it (calibrate_pro_chl: CHL runs with 52-249
% Pro well above the edge had 1-2% of the Pro run's Pro)
p = struct('n', sum(class == 7), 'peak_above_edge', NaN, 'frac_near_edge', NaN, 'lower_bound', true);
if ~isfield(info.trigger_edges, 'BL3') || p.n < cfg.sig.pro_edge_min_n, return, end
x = info.features.BL3(class == 7) - info.trigger_edges.BL3;
g = linspace(0, 3, 301);
d = ksdensity(x, g, 'bandwidth', 0.05);
[~, k] = max(d);
p.peak_above_edge = g(k);
p.frac_near_edge = mean(x < cfg.sig.pro_edge_band);
p.lower_bound = p.peak_above_edge < cfg.sig.pro_lb_min_peak_above_edge | p.frac_near_edge > cfg.sig.pro_lb_max_frac_near_edge;
end

function [cut, s] = euk_autofluor_cut(fcsdat, fcshdr, ref, cfg)
% PE detected = PE:chl above the file's Euk autofluorescence: the Euks form the
% lowest major mode of PE:chl among cell-like events (chl:size above the
% detritus cut); the cut is that mode + cfg.sig.euk_af_above_mode. CHL runs: mode
% ~ -2.05, cut ~ -1.3; PE-triggered runs record only the brightest Euks (EN668
% C002N17: mode ~ -1.2, which the fixed cut called Syn). Without a clear Euk
% mode, cfg.sig.pe_min_PEchl.
s = struct('mode', NaN, 'n', 0);
cut = cfg.sig.pe_min_PEchl;
F = attune_features(fcsdat, fcshdr, ref);
pechl = F.GL2 - F.BL3;
cell_like = F.BL3 - F.SSC > cfg.sig.detritus_max_chlSSC & isfinite(pechl) & F.BL3 > prctile(F.BL3, 10);
x = pechl(cell_like);
s.n = numel(x);
if numel(x) < cfg.sig.euk_af_min_n, return, end
g = linspace(-3.5, 1, 226);
d = ksdensity(x, g, 'bandwidth', 0.08);
pk = find(d(2:end-1) > d(1:end-2) & d(2:end-1) >= d(3:end)) + 1;
pk = pk(d(pk) >= cfg.sig.euk_af_min_peak * max(d) & g(pk) < cfg.sig.euk_af_max_mode);
if isempty(pk), return, end
s.mode = g(pk(1));
cut = s.mode + cfg.sig.euk_af_above_mode;
end
