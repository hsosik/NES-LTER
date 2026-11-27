function [class, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, pro_measured)
% function [class, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, pro_measured)
% Prototype cruise-independent classifier.
%  1. put all channels on the bead-relative scale (attune_features)
%  2. drop events with no chlorophyll signal above the trigger edge
%  3. fit a Gaussian mixture to chlorophyll-containing events in
%     [SSC, BL3, PE:chl, GL2:GL3] and pick the number of components by BIC
%  4. name each component from its mean using pigment signatures that are
%     the same across cruises and seasons (cfg.sig)
%  5. events far from their component (Mahalanobis) are left as 0
% Class codes match the existing class files:
%  0 = not phytoplankton / unassigned, 1 = Euk without PE, 2 = Syn,
%  3 = dim cryptophyte (less PE, more red-shifted PE spectrum), 4 = bright cryptophyte,
%  5 = Euk with Syn (Syn-like PE, ~one Syn's worth; possible grazers), 7 = Pro
%  8 = particle-associated Syn (new code: Syn on or inside a larger particle,
%      ~1 Syn per event; count with Syn abundance, not Syn size or carbon)
%
% ref - from bead_reference; cfg - from classifier_config
%
% Nov 2026, classifier redesign prototype

% laser dropout, split Syn (Nov 2026): in EN661 5-6 Feb 2021 (while #LASER2DELAY was 305, with
% slow flow) and single files of EN649 and EN668 the instrument often recorded one Syn as two
% events with the same Time stamp: the blue-laser part (FSC, SSC, BL channels; no green-laser PE)
% and right after it the green-laser part (GL channels; no SSC/chl pulse). Unmerged, the green
% parts were a tiny, chl-poor cluster below the chl trigger (named Syn, or setting the Syn
% reference) and the blue parts were Euks (+34% old classes, +51% v2 over neighbouring files).
% A pair = a blue part (bright PE under blue with Syn-like PE:chl under blue, GL2 far below BL2,
% narrow GL2 pulse) followed within dropout_pair_dt Time units by an event with PE under green
% well above PE under blue (65% of blue parts in the EN661 episode). The pair is merged into the
% blue part's row before anything else: blue-laser channels from the blue part, green-laser
% channels from the green part; the green part's row is then left out (flag 6). Merged blue parts
% have SSC within ~7% of whole Syn in neighbouring files, chl ~23% lower (off-core passage).
% Only in files with at least dropout_min_split pairs: elsewhere the tests flagged real Syn under
% other PE gains (AR92, HRS2609) and preserved runs had blue-only events with no green parts
sig = cfg.sig; n = size(fcsdat, 1);
mergedg = false(n, 1); noblue = false(n, 1); nogreen = false(n, 1); info.n_split_pairs = 0;
pn = {fcshdr.par.name};
if sig.noblue_on && all(ismember({'GL2-A' 'BL2-A' 'BL3-A' 'SSC-W' 'BL3-W' 'GL2-W' 'Time'}, pn))
    lc = @(d, s) log10(max(d(:, find(strcmp(pn, s), 1)), 1)); wc = @(d, s) d(:, find(strcmp(pn, s), 1));
    blue_part = @(d) lc(d, 'BL2-A') > sig.nogreen_min_BL2 & lc(d, 'BL2-A') - lc(d, 'BL3-A') > sig.nogreen_min_BL2BL3 & ...
        lc(d, 'GL2-A') - lc(d, 'BL2-A') < -sig.noblue_min_GL2BL2 & wc(d, 'GL2-W') < sig.nogreen_max_GL2W;
    green_part = @(d) lc(d, 'GL2-A') - lc(d, 'BL2-A') > sig.noblue_min_GL2BL2 & wc(d, 'SSC-W') < sig.noblue_max_SSCW & ...
        wc(d, 'BL3-W') < sig.noblue_max_BL3W & lc(d, 'GL2-A') > sig.split_min_GL2;
    bp = blue_part(fcsdat); tm = wc(fcsdat, 'Time');
    rg = lc(fcsdat, 'GL2-A'); rb = lc(fcsdat, 'BL2-A');
    ib = find(bp(1:end - 1));
    pr = tm(ib + 1) - tm(ib) <= sig.dropout_pair_dt & rg(ib + 1) - rb(ib + 1) > sig.merge_min_GL2BL2 & rg(ib + 1) > sig.split_min_GL2 & ~bp(ib + 1);
    ib = ib(pr);
    info.n_split_pairs = numel(ib);
    if numel(ib) >= sig.dropout_min_split
        gcol = startsWith(pn, 'GL');
        fcsdat(ib, gcol) = fcsdat(ib + 1, gcol);
        mergedg(ib + 1) = true;
        noblue = green_part(fcsdat) & ~mergedg;  % green parts left unpaired
        nogreen = blue_part(fcsdat);             % blue parts left unpaired (merged rows now carry green PE)
    end
end
info.n_merged = sum(mergedg); info.n_unpaired_blue = sum(nogreen);
[F, finfo] = attune_features(fcsdat, fcshdr, ref);
n = height(F);
class = zeros(n, 1);
% size: SSC, switching to scatter through the OD2 filter (GL1) for large particles.
% SSC-H does not clip; it rolls off gradually from ~10x the 1 um bead while GL1
% stays linear. As in batch_mergedSSC2volume, SSC and GL1 each relative to their
% own 1 um bead peak are taken as equivalent (no fitted offset): GL1 above
% cfg.gl1_merge_bead(2) x the 1 um bead, the mean of the two between
% cfg.gl1_merge_bead(1) and (2) (bead units, so the window does not move with
% the GL1 voltage; was raw GL1-H 500-1000 counts)
% Only where GL1 is OD2 green scatter: cfg.gl1_is_od2_scatter (set per cruise from
% attune_config(cruise).OD2 == "GL1"), ref.OD2 when bead_reference knew the cruise,
% and GL1 actually measured (no green scatter channel before Aug 2019, when OD2 sat
% on SSC or nowhere: size is SSC alone)
F.size = F.SSC;
use_gl1 = cfg.gl1_is_od2_scatter && any(isfinite(F.GL1)) && (~isfield(ref, 'OD2') || strcmpi(ref.OD2, 'GL1'));
info.size_from_gl1 = use_gl1;
if use_gl1
    hi = F.GL1 > log10(cfg.gl1_merge_bead(2));
    mid = F.GL1 > log10(cfg.gl1_merge_bead(1)) & ~hi;
    F.size(hi) = F.GL1(hi);
    F.size(mid) = log10((10.^F.SSC(mid) + 10.^F.GL1(mid)) / 2);
end
F.PEchl = F.GL2 - F.BL3;
F.GL2GL3 = F.GL2 - F.GL3;
F.chlSSC = F.BL3 - F.size;
sig = cfg.sig;

% trigger edges in bead-relative units; events pile up just above them
edge = struct();
trigs = [finfo.trigger finfo.trigger2]; levels = [finfo.trigger_level finfo.trigger2_level];
for t = 1:2
    ch = char(regexp(trigs(t), '(?<=_)\w+$', 'match', 'once'));
    if ~isempty(ch) && isfield(finfo.floor, ch)
        edge.(ch) = finfo.floor.(ch) + log10(levels(t)) - log10(cfg.floor_raw);
    end
end
info.trigger_edges = edge;
info.floor = finfo.floor;

% pigment present: chlorophyll clear of its trigger edge, or PE clearly
% detected (OR_GL2 triggering picks up dim Syn by PE alone)
chl_min = sig.chl_min;
if isfield(edge, 'BL3'), chl_min = max(chl_min, edge.BL3 + sig.trigger_margin); end
pe_min = sig.pe_min_GL2;
if isfield(edge, 'GL2'), pe_min = max(pe_min, edge.GL2 + sig.trigger_margin); end
has_pe = F.GL2 > pe_min & F.PEchl > sig.pe_min_PEchl & F.BL3 > sig.chl_min;
% split-Syn parts (merged green parts, unpaired parts) are left out of every rule
has_pe(noblue | mergedg) = false;
has_chl = (F.BL3 > chl_min | has_pe) & F.SSC > finfo.floor.SSC + 0.1 & ~noblue & ~nogreen & ~mergedg;
X = [F.size F.BL3 F.PEchl F.GL2GL3];
X = X(has_chl, :);
info.n_has_chl = sum(has_chl);
% chl pulse shape (Nov 2026): small single cells give short BL3 pulses, like Syn
% (log10 A/H about 0). In AR95 Pro runs a second group of faint, small, no-PE
% signals with wide pulses (log A/H about 0.3-0.4, 150-550 per file, also present
% under rich Pro clusters) sat among the Pro and was counted as Pro (24% of v2 Pro
% in 21Oct AltPro(24)). Pulse shape is judged against small PE-bearing events
% (mostly Syn; 99th percentile 0.06-0.15 above their median). Events more than
% pro_max_AH_above_ref above that are 'wide' (NaN reference or column: none)
info.pulse_ref_AH = NaN;
wide = false(n, 1);
if ismember('BL3_AH', F.Properties.VariableNames)
    rs = has_pe & F.BL3 > chl_min + 0.5 & F.size < sig.syn_size_default + sig.pulse_ref_size_above & isfinite(F.BL3_AH);
    if sum(rs) >= 50, info.pulse_ref_AH = median(F.BL3_AH(rs)); end
    wide = F.BL3_AH > info.pulse_ref_AH + sig.pro_max_AH_above_ref; %false where the reference is NaN
end
info.trace.wide_pulse = wide;
% in Pro-setting files, small dim no-PE events with wide pulses are left out of the
% mixture fit (they are still assigned with it), so they neither pull Pro components
% into the smear nor stop a small compact Pro cluster from getting its own
% component (AR95 19Oct AltPro(30) Pro(10))
fitok = true(size(X, 1), 1);
if pro_measured
    wx = wide & ~has_pe & F.BL3 < sig.pro_max_BL3;
    fitok = ~wx(has_chl);
end
info.n_fit_wide_excluded = sum(~fitok);
if size(X, 1) < 50
    info.ncomp = 0; info.comp = table;
    info.syn_type = zeros(n, 1, 'uint8'); info.flag = zeros(n, 1, 'uint8'); info.features = F; %nothing assigned
    info.euk_chl_gap = -Inf; info.n_euk_below_gap = 0; info.syn_rescued = 0; info.n_euk_synpe = 0; info.n_noblue = sum(noblue); info.n_syn_noblue = 0; info.n_syn_nogreen = sum(nogreen); info.dropout_frac_syn = NaN;
    return
end

% fit on a subsample, choose K by BIC
if isfield(cfg, 'seed'), sd = cfg.seed; else, sd = 1; end
rng(sd, 'twister') %same fit in serial runs and on parallel workers (whose default generator differs)
nfit = min(size(X, 1), cfg.gmm_max_fit);
% trigger-noise pile (Nov 2026): in files run without Pro settings, no-PE events
% just above the chl trigger edge can be most of the file (AR34B, RB1904, EN627:
% 60-98% of chl events; post-2019 non-Pro files < 25%). Fitted at full weight they
% take nearly every component, real Syn are absorbed and a few bright cryptos get
% named Syn. Where the pile is over edge_pile_max_frac of the events, it enters the
% fit sample with at most edge_pile_fit_ratio x the other events (all events are
% still assigned with the fitted model; the Pro-tail rule below handles the pile)
pile = false(size(X, 1), 1);
if ~pro_measured && isfield(edge, 'BL3')
    pile = ~has_pe(has_chl) & X(:, 2) < edge.BL3 + sig.edge_pro_band;
end
info.edge_pile_frac = mean(pile);
info.edge_pile_capped = info.edge_pile_frac > sig.edge_pile_max_frac;
if info.edge_pile_capped
    ip = find(pile); in = find(~pile);
    % ... but enough events in all to fit the largest mixture (20 per component)
    npf = max(round(sig.edge_pile_fit_ratio * numel(in)), 20 * max(cfg.gmm_K) - numel(in));
    ip = ip(randperm(numel(ip), min(numel(ip), npf)));
    pool = [in; ip];
    nfit = min(numel(pool), cfg.gmm_max_fit);
    Xfit = X(pool(randperm(numel(pool), nfit)), :);
elseif all(fitok)
    Xfit = X(randperm(size(X, 1), nfit), :);
else
    idx = find(fitok);
    nfit = min(numel(idx), cfg.gmm_max_fit);
    Xfit = X(idx(randperm(numel(idx), nfit)), :);
end
best = []; bestbic = Inf;
opts = statset('MaxIter', 500);
warning('off', 'stats:gmdistribution:FailedToConvergeReps')
warning('off', 'stats:gmdistribution:FailedToConverge')
for K = cfg.gmm_K
    if nfit < 20 * K, break, end
    try
        g = fitgmdist(Xfit, K, 'RegularizationValue', 1e-3, 'Replicates', 2, 'Options', opts);
    catch
        continue
    end
    if g.BIC < bestbic, best = g; bestbic = g.BIC; end
end
if isempty(best) %no mixture could be fitted (too few events): nothing assigned
    info.ncomp = 0; info.comp = table;
    info.syn_type = zeros(n, 1, 'uint8'); info.flag = zeros(n, 1, 'uint8'); info.features = F; %nothing assigned
    info.euk_chl_gap = -Inf; info.n_euk_below_gap = 0; info.syn_rescued = 0; info.n_euk_synpe = 0; info.n_noblue = sum(noblue); info.n_syn_noblue = 0; info.n_syn_nogreen = sum(nogreen); info.dropout_frac_syn = NaN;
    return
end
info.ncomp = best.NumComponents;
info.gmm = best;

% assign events to components (used for naming below and for assignment after it)
[comp, ~, ~, ~, d2] = cluster(best, X);
% chl spread of each component for the Pro compactness test (Nov 2026): robust SD of
% the narrow-pulse events assigned to it. The fitted Gaussian SD is inflated by a
% component's tail toward the chl trigger (AR95 19Oct AltPro(29): 0.43 with a compact
% core); components with fewer than 20 such events keep the fitted SD
sdB = squeeze(sqrt(best.Sigma(2, 2, :)));
sdS = squeeze(sqrt(best.Sigma(1, 1, :))); %size spread, the same way (Pro size-spread cap below)
wx = wide(has_chl);
for cc = 1:best.NumComponents
    in_c = comp == cc & ~wx;
    if sum(in_c) >= 20, sdB(cc) = 1.4826 * mad(X(in_c, 2), 1); sdS(cc) = 1.4826 * mad(X(in_c, 1), 1); end
end
% name components from their means
mu = name_components(best, finfo, sig, pe_min, pro_measured, sdB);
mu.sd_BL3 = sdB(:);
mu.sd_size = sdS(:);
% Pro candidates (compact, dim, no PE) are judged against this file's Syn
% components (Nov 2026; survey of 108 Pro-setting files in 9 cruises,
% survey_pro / comp_detail in the README):
%  size: smaller than Syn by more than pro_min_size_below_syn. The fixed
%    pro_max_SSC (-1.44) named fall/spring Pro clusters only 0.15-0.5 smaller than
%    Syn (size about -1.40) as Euk (AR104 23Oct wPro(1), AR82A, AE2426, AR79),
%    leaving them to the catch-all; candidates that are not smaller go back to Euk.
%    And no more than pro_max_size_below_syn smaller: tiny high-chl artifacts
%    1.2-1.5 below Syn (AR104 23Oct wPro(9) Sample(1), AR79) are not Pro (real Pro
%    0.1-1.0 below)
%  chl per scatter: below Syn's (real Pro 0.1-1.0 below; the artifacts sit at or
%    above Syn), and within pro_max_chlSSC_below_syn of Syn; a compact group far
%    lower is faint detritus / small particles at the bottom of the background
%    smear (AR95 19Oct AltPro(30) Pro(10): 1.1-1.2 below), not Pro
% Files with no Syn component keep the fixed size window
info.pro_comp_rejected = zeros(0, 1);
if any(mu.label == 7)
    if any(mu.label == 2)
        ws = mu.weight(mu.label == 2);
        sync = sum(mu.chlSSC(mu.label == 2) .* ws) / sum(ws);
        syns = sum(mu.SSC(mu.label == 2) .* ws) / sum(ws);
        cand = mu.label == 7;
        % chl:scatter limit below Syn grows with the component's height above the chl
        % trigger edge (Nov 2026): near the edge (within pro_chlSSC_edge_start) it is
        % pro_max_chlSSC_below_syn, where faint detritus piles up (AR95 19Oct AltPro(30):
        % 1.13 below Syn, 0.20 above the edge; near-edge groups without old-class Pro in
        % AE2426, AR95, HRS2609, AR104: 1.03-1.16 below, 0.16-0.27 above). Compact clusters
        % well clear of the edge are real Pro at up to ~1.05 below (AR79 10Nov Pro(14):
        % 1.03 below, 0.55 above; AR79 09Nov Pro(2): ~1.0 below, ~0.35 above)
        chl_lim = repmat(sig.pro_max_chlSSC_below_syn, height(mu), 1);
        if isfield(edge, 'BL3')
            chl_lim = min(sig.pro_max_chlSSC_below_syn_far, sig.pro_max_chlSSC_below_syn + ...
                sig.pro_chlSSC_edge_slope * max(0, mu.BL3 - edge.BL3 - sig.pro_chlSSC_edge_start));
        end
        mu.pro_chl_lim = chl_lim;
        too_big = cand & mu.SSC >= syns - sig.pro_min_size_below_syn;
        % ... and no larger than any Pro cluster seen (bead-relative size, the same
        % scale in every cruise): where Syn are unusually large (AR104 24Oct shelf, Syn
        % at -0.41) a broad dim group at -0.88 was otherwise smaller enough than Syn
        bad = cand & ~too_big & (mu.SSC > sig.pro_max_size_abs | mu.SSC < syns - sig.pro_max_size_below_syn | ...
            mu.chlSSC > sync + sig.pro_max_chlSSC_above_syn | mu.chlSSC < sync - chl_lim);
        % ... and compact in size too (Nov 2026): real Pro components have a size spread
        % (robust SD) of 0.10-0.38; a broad column of faint small particles along the SSC
        % trigger was named Pro in preserved EN727 C19N13 (no real Pro; size SD ~0.5)
        broad = cand & ~too_big & mu.sd_size > sig.pro_max_sd_size;
        bad = bad | broad;
        % bottom slice of a Pro cluster (Nov 2026): where the chl trigger cuts a Pro
        % cluster, the mixture can give its lowest part a component of its own, at the
        % cluster's size but lower chl, so lower chl:scatter. A candidate that fails only
        % the chl:scatter limit (and is within pro_max_chlSSC_below_syn_far of Syn) is Pro
        % if an accepted Pro component lies directly above it: size within
        % pro_slice_size_tol, brighter in chl (HRS2609 18Aug(1) Pro(19): slice 1.04 below
        % Syn, 0.15 above the edge, 743 of its 842 events old-class Pro). The near-edge
        % false groups have no Pro cluster above them (AR95 19Oct AltPro(30), AE2426 09Nov(10))
        chl_only = bad & ~broad & mu.chlSSC < sync - chl_lim & mu.chlSSC >= sync - sig.pro_max_chlSSC_below_syn_far & ...
            ~(mu.SSC > sig.pro_max_size_abs | mu.SSC < syns - sig.pro_max_size_below_syn | mu.chlSSC > sync + sig.pro_max_chlSSC_above_syn);
        acc = find(cand & ~too_big & ~bad);
        info.pro_comp_slice = zeros(0, 1);
        for j = find(chl_only)'
            if any(abs(mu.SSC(acc) - mu.SSC(j)) <= sig.pro_slice_size_tol & mu.BL3(acc) > mu.BL3(j))
                bad(j) = false;
                info.pro_comp_slice(end + 1, 1) = j;
            end
        end
        % Syn-sized part of a Pro cluster (Nov 2026): the mixture can split a Pro cluster
        % so that one component sits within pro_min_size_below_syn of Syn's size; it was
        % named Euk, then unassigned as Pro-like (AR79 10Nov Pro(14): the cluster's centre,
        % 9,010 events, 0.09 smaller than Syn; the rest of the cluster accepted as Pro).
        % A too-large candidate is Pro if it is no larger than the Syn components, passes
        % the other Pro rules, and an accepted Pro component is within pro_slice_size_tol
        % in size. A lone Syn-sized group stays Euk
        adj = too_big & mu.SSC <= syns & mu.SSC <= sig.pro_max_size_abs & mu.sd_size <= sig.pro_max_sd_size & ...
            mu.chlSSC <= sync + sig.pro_max_chlSSC_above_syn & mu.chlSSC >= sync - chl_lim;
        info.pro_comp_adjoin = zeros(0, 1);
        for j = find(adj)'
            if any(abs(mu.SSC(acc) - mu.SSC(j)) <= sig.pro_slice_size_tol)
                too_big(j) = false;
                info.pro_comp_adjoin(end + 1, 1) = j;
            end
        end
        mu.label(too_big) = 1;
        mu.label(bad) = 0;
        info.pro_comp_rejected = find(bad);
    else
        too_big = mu.label == 7 & mu.SSC >= sig.pro_max_SSC;
        mu.label(too_big) = 1;
    end
end
info.comp = mu;

% assign events, leave outliers unassigned
d2own = d2(sub2ind(size(d2), (1:size(X, 1))', comp));
lab = mu.label(comp);
lab(d2own > chi2inv(cfg.outlier_p, size(X, 2))) = 0;
class(has_chl) = lab;
% with the pile under-weighted in the fit, its events can fall into Syn or PE-euk
% components; they have no PE, so they may only be Euks (cut later as the Pro tail
% where small) or unassigned, and stay out of the catch-all
pile_ev = false(n, 1);
if info.edge_pile_capped
    pile_ev(has_chl) = pile;
    class(pile_ev & class ~= 1) = 0;
end
% rare Syn (Nov 2026): where Syn are a tiny fraction of the file (~100-200 of 80-100k
% events in AR34B / RB1904 files) no component is named Syn and they end up as PE euks
% or unassigned. Without a Syn component, events with PE, a Syn-type PE spectrum, and
% PE:chl and size near the usual Syn (cfg defaults, Syn-check tolerances) start as Syn;
% the Syn checks below then work from them as usual
info.syn_rescued = 0;
if sig.syn_rescue && ~any(class == 2)
    gl3c = F.GL3 > finfo.floor.GL3 + sig.ratio_min_above_floor;
    rs = has_pe & ~pile_ev & (~gl3c | F.GL2GL3 > sig.syn_min_GL2GL3) & ...
        abs(F.PEchl - sig.syn_PEchl_default) < sig.syn_PEchl_tol & abs(F.size - sig.syn_size_default) < sig.syn_size_tol & ...
        F.chlSSC > sig.detritus_max_chlSSC;
    % ... and chl clearly above the chl trigger edge: in OR-trigger files faint events
    % that cleared only the PE trigger (chl below the edge) have Syn-like PE:chl and
    % spectrum, and seeding Syn from them (RB1904 12May Sample(1)) moved the file's Syn
    % reference onto that pile (Syn on particle, crypto PE scale and PE-euk test followed)
    if isfield(edge, 'BL3'), rs = rs & F.BL3 > edge.BL3 + sig.edge_pro_band; end
    if sum(rs) >= sig.syn_rescue_min_n
        class(rs) = 2;
        info.syn_rescued = sum(rs);
    end
end
info.trace.gmm_label = mu.label(comp); info.trace.gmm_outlier = d2own > chi2inv(cfg.outlier_p, size(X, 2)); info.trace.has_chl = has_chl;
info.trace.after_gmm = class;

% event-level checks: broad components pick up stray events
% too little chlorophyll for the particle size: detritus / non-phytoplankton
det = F.chlSSC < sig.detritus_max_chlSSC;
% Pro sit at the chl trigger edge, so the largest ones fall below the fixed cut
% without being detritus; drop Pro only if also far below the Pro cluster's own
% chl:size (in the sparse tail, not through the cluster's edge)
pro = class == 7;
if any(pro)
    det(pro) = det(pro) & F.chlSSC(pro) < median(F.chlSSC(pro)) - sig.pro_chlSSC_sd * 1.4826 * mad(F.chlSSC(pro), 1);
end
% same for Syn: dim-chl Syn near the chl trigger edge (or detected by PE alone)
synd = class == 2;
if any(synd)
    det(synd) = det(synd) & F.chlSSC(synd) < median(F.chlSSC(synd)) - sig.pro_chlSSC_sd * 1.4826 * mad(F.chlSSC(synd), 1);
end
class(det) = 0;
info.trace.after_detritus = class;
% Syn must sit near this file's Syn cluster in PE:chl, PE spectrum and size
syn = class == 2;
% the reference (centre and spread) is taken from Syn-named events no larger than
% syn_ref_max_size (Nov 2026). A broad PE-bearing component of large particles can be
% named Syn (preserved winter CHL runs EN627 C45N05, AR99 C02N16; underway AR88
% 26Apr(3) Sample(30), EN627 1Feb2019a Sample(1)); its events then outnumbered the real
% Syn, the reference moved to size +0.2 to +0.5, and "no Euk much smaller than Syn"
% removed the small Euks. Real Syn peak at -1.5 to -0.75 (winter EN627 the largest)
synr = syn;
if sum(syn & F.size <= sig.syn_ref_max_size) >= sig.syn_ref_min_n
    synr = syn & F.size <= sig.syn_ref_max_size;
end
% ... and, in files triggered on PE as well as chl, from Syn-named events clearly above
% the chl trigger edge when there are syn_ref_min_n of them (Nov 2026). In winter EN661
% 05Feb(20) Sample(27) a tight cluster of tiny PE-triggered events below the chl edge
% (size -2.5, PE:chl +1.4) outnumbered the real Syn (size -1.1, PE:chl 0); the reference
% went to it and the real Syn became Euk with Syn (the Syn rescue already has this rule).
% Only where those events have Syn-like PE:chl: where the Syn-named events above the edge
% are a broad particle component (EN661 06Feb(5) Sample(5), PE:chl -1.3) the rule is skipped
if sig.syn_ref_above_edge && isfield(edge, 'BL3')
    sa = synr & F.BL3 > edge.BL3 + sig.edge_pro_band;
    if sum(sa) >= sig.syn_ref_min_n && abs(median(F.PEchl(sa)) - sig.syn_PEchl_default) < sig.syn_PEchl_tol
        synr = sa;
    end
end
info.syn_ref_n_above = sum(syn & ~synr);
syn_PEchl = median(F.PEchl(synr)); syn_size = median(F.size(synr));
syn_GL2GL3 = median(F.GL2GL3(synr)); syn_BL3 = median(F.BL3(synr));
if ~any(syn), syn_PEchl = sig.syn_PEchl_default; syn_size = sig.syn_size_default; end
% GL2:GL3 is only meaningful where GL3 is clear of its floor (dim summer Syn are not)
gl3ok = F.GL3 > finfo.floor.GL3 + sig.ratio_min_above_floor;
% ... and the reference PE:chl from those of them with a measured Syn-type spectrum
% (GL3 clear of its floor, GL2:GL3 above syn_min_GL2GL3) when there are at least
% syn_ref_spec_min_n (Nov 2026). In MVCO MV337N04 (preserved) a cloud of dim PE
% particles (PE:chl ~ -0.7, GL3 at its floor so no spectrum test) was named Syn and
% outnumbered the 42 real Syn (PE:chl +0.35) ten to one, so the Syn check kept the
% cloud. Where dim Syn have GL3 at the floor (summer) the brighter Syn set the value
info.syn_ref_n_spec = sum(synr & gl3ok & F.GL2GL3 > sig.syn_min_GL2GL3);
if info.syn_ref_n_spec >= sig.syn_ref_spec_min_n
    syn_PEchl = median(F.PEchl(synr & gl3ok & F.GL2GL3 > sig.syn_min_GL2GL3));
end
% Syn pigment types differ in GL2:GL3 and a sample can hold a mixture (e.g. two
% types at ~0.26 and ~0.4 in AR92), so accept the range spanned by this file's Syn
% components (those with a Syn-type spectrum measured clear of the GL3 floor), not one centre
syn_range = syn_ratio_range(mu, syn_GL2GL3, sig);
% PE:chl below Syn (PE euks) uses the fixed tolerance; size uses it widened to the
% file's Syn spread; more PE per chl than Syn is allowed further (dim-chl Syn near the
% chl trigger edge or detected by PE alone), up to syn_PEchl_sd robust SDs
if any(syn)
    syn_PEchl_sdv = 1.4826 * mad(F.PEchl(synr), 1);
    % size centre and spread from the peak of the distribution, which neither an upward
    % smear of large Syn nor tiny signals below shift or stretch
    [syn_size, syn_size_sdv] = peak_centre(F.size(synr));
else
    syn_PEchl_sdv = 0; syn_size_sdv = 0;
end
% Syn PE before the check, for the PE gap below; events too large for Syn are left
% out (a bright-crypto component named Syn in AR34B files put the PE peak above the
% real Syn, and the gap below it removed every Syn). Only the upper size limit: the
% faint pile the gap is for is small
gs0 = F.GL2(syn & F.size <= syn_size + sig.syn_size_tol);
class(syn & (F.PEchl < syn_PEchl - sig.syn_PEchl_tol | F.PEchl > syn_PEchl + max(sig.syn_PEchl_tol, sig.syn_PEchl_sd * syn_PEchl_sdv) | ...
    F.size > syn_size + sig.syn_size_tol | F.size < syn_size - min(sig.syn_size_max_below, max(sig.syn_size_tol_below, sig.syn_size_sd * syn_size_sdv)) | ...
    (gl3ok & (F.GL2GL3 < syn_range(1) - sig.syn_GL2GL3_tol | F.GL2GL3 > syn_range(2) + sig.syn_GL2GL3_tol)))) = 0;
% boundaries of the Syn check, kept for plotting (plot_v2_rules)
info.rules.syn_check = struct('PEchl', syn_PEchl, 'PEchl_lo', syn_PEchl - sig.syn_PEchl_tol, ...
    'PEchl_hi', syn_PEchl + max(sig.syn_PEchl_tol, sig.syn_PEchl_sd * syn_PEchl_sdv), ...
    'size', syn_size, 'size_hi', syn_size + sig.syn_size_tol, ...
    'size_lo', syn_size - min(sig.syn_size_max_below, max(sig.syn_size_tol_below, sig.syn_size_sd * syn_size_sdv)), ...
    'GL2GL3_lo', syn_range(1) - sig.syn_GL2GL3_tol, 'GL2GL3_hi', syn_range(2) + sig.syn_GL2GL3_tol, ...
    'gl3ok_min', finfo.floor.GL3 + sig.ratio_min_above_floor);
% Syn must carry PE close to the Syn cluster's (dim-chl Syn keep normal PE; events
% low in both chl and PE are not Syn). Spreads are taken from the low side of the
% distribution, which an upward tail (winter smear) does not stretch
syn_PE_floor = -Inf;
if any(class == 2)
    gs = F.GL2(class == 2); gmed = median(gs);
    syn_PE_floor = gmed - max(sig.syn_PE_tol_below, sig.syn_PE_sd * (gmed - prctile(gs, 15.87)));
    % ... and above a PE gap. In OR-trigger files, faint events that cleared only the
    % GL2 trigger pile up at its edge, separated from the Syn cluster by a trough in
    % PE. The gap is found in the Syn PE distribution before the Syn check (which
    % removes most of the pile by size, leaving its Syn-sized part as a flat tail with
    % no visible trough): a local minimum below the Syn peak, under syn_PE_gap_frac of
    % the peak density, with density rising again below it by syn_PE_gap_rise into a
    % pile of at least syn_PE_gap_min_pile events. Syn whose PE runs smoothly down to
    % the trigger edge (summer) have no such gap and keep the statistical floor
    pe_gap = -Inf;
    if numel(gs0) >= sig.syn_PE_gap_min_n
        [fk, gk] = ksdensity(gs0, 'NumPoints', 400, 'Bandwidth', sig.syn_PE_gap_bw);
        [fm, ip] = max(fk); best = Inf;
        for j = ip - 1:-1:2
            if fk(j) <= fk(j - 1) && fk(j) <= fk(j + 1) && fk(j) < sig.syn_PE_gap_frac * fm && ...
                    max(fk(1:j - 1)) > sig.syn_PE_gap_rise * fk(j) && sum(gs0 < gk(j)) >= sig.syn_PE_gap_min_pile && fk(j) < best
                best = fk(j); pe_gap = gk(j);
            end
        end
    end
    syn_PE_floor = max(syn_PE_floor, pe_gap);
    info.syn_PE_gap = pe_gap;
    % no-PE noise floor (on when cfg sets syn_PE_noise_prc): upper percentile of GL2
    % in small Euks only (no larger than Syn + syn_PE_noise_size_above). Large Euks
    % show PE-channel signal from bright chl spillover, and the Euk set at this stage
    % still holds cryptos, so all Euks together overstate the noise
    info.syn_PE_noise = NaN;
    if isfield(sig, 'syn_PE_noise_prc')
        se = class == 1 & F.size < syn_size + sig.syn_PE_noise_size_above;
        if sum(se) >= sig.syn_PE_noise_min_n
            info.syn_PE_noise = prctile(F.GL2(se), sig.syn_PE_noise_prc);
            syn_PE_floor = max(syn_PE_floor, info.syn_PE_noise);
        end
        info.syn_PE_noise_n = sum(se);
    end
    class(class == 2 & F.GL2 < syn_PE_floor) = 0;
end
info.syn_PE_floor = syn_PE_floor;
syn_rejected = syn & class == 0; %re-checked as Euks below
info.trace.after_syncheck = class;
info.trace.syn_ref0 = struct('PEchl', syn_PEchl, 'GL2GL3', syn_GL2GL3, 'BL3', syn_BL3, 'size', syn_size, 'GL2GL3_range', syn_range);
% the tolerances above are broad, so take the file's Syn position from the checked events
syn = class == 2;
if any(syn)
    syn_PEchl = median(F.PEchl(syn)); syn_GL2GL3 = median(F.GL2GL3(syn)); syn_BL3 = median(F.BL3(syn));
    syn_size = peak_centre(F.size(syn));  %peak of the Syn size distribution (robust to a smear above or tiny signals below)
end
% no known Euks are much smaller than Syn (Ostreococcus is only slightly smaller);
% smaller chlorophyll-bearing particles without PE are Pro only in Pro-setting files
class(class == 1 & F.size < syn_size - sig.euk_min_size_below_syn) = 0;
% ... and Euks smaller than the Syn peak must be brighter in chl than Pro can be:
% small no-PE events dimmer than pro_max_BL3 are Pro-like (bright Pro in files run
% without Pro settings, or Pro the mixture did not separate), not Euks (Nov 2026)
euk_prolike = ~has_pe & F.size < syn_size & F.BL3 < sig.pro_max_BL3;
class(class == 1 & euk_prolike) = 0;
info.trace.euk_prolike = euk_prolike;
% Euks carry no PE: events named Euk that have PE with about the PE:chl of this file's
% Syn are not Euks (Nov 2026). In OR-trigger files faint Syn-like events detected by
% the PE trigger alone (chl below the chl edge) were put in Euk components (AR95,
% HRS2609, AR34B, RB1904: 40-90 per file); they are left for the Syn rules and the
% catch-all. (Replaces the Pro-tail cut, which removed no-PE events within 0.4 of
% the chl edge and no larger than Syn: in summer files that sliced the bottom off real
% small-Euk clusters, ~3,500 old-class Euks in 12 TN368 files; the Euk chl-gap rule
% handles the no-PE pile at the edge where there is one)
euk_synpe = class == 1 & has_pe & F.PEchl > syn_PEchl - sig.syn_PEchl_tol;
class(euk_synpe) = 0;
info.n_euk_synpe = sum(euk_synpe);
% where the trigger-noise pile dominates (capped fit, pre-2019), no-PE events in the
% edge band are left out at any size, as the old per-cruise classifiers did (raw BL3
% < 200 with a trigger of 100); flag 3
pile_cut = class == 1 & pile_ev;
class(pile_cut) = 0;
% events in a component named PE euk (broad large-particle components can hold
% bright cryptos and large Euks together) but without PE of their own are Euks
class(class == 3 & ~has_pe) = 1;
% Euks must have clearly more chl per scatter than non-phytoplankton particles.
% Euks are variable, so the cut is the looser of (Euk median - euk_chlSSC_tol) and
% the midpoint between the Euk median and the median of this file's large
% non-phytoplankton particles (scatter without matching chl)
euk = class == 1;
if any(euk)
    euk_med = median(F.chlSSC(euk));
    euk_cut = euk_med - sig.euk_chlSSC_tol;
    nonphyto = F.chlSSC < sig.detritus_max_chlSSC & F.size > syn_size + sig.euk_rescue_size_above_syn & F.SSC > finfo.floor.SSC + 0.1;
    if sum(nonphyto) >= sig.nonphyto_min_n
        euk_cut = min(euk_cut, (euk_med + median(F.chlSSC(nonphyto))) / 2);
    end
    euk_cut = max(euk_cut, sig.detritus_max_chlSSC);
    info.euk_chlSSC_cut = euk_cut;
    class(euk & F.chlSSC < euk_cut) = 0;
    % large Euks are sparse and can be absorbed by a broad detritus component (named
    % 0 from its mean); give back events in such components that individually sit
    % in this file's Euk chl:size band, have no PE, and are well larger than Syn
    in_det = false(n, 1); in_det(has_chl) = mu.label(comp) == 0 & d2own <= chi2inv(cfg.outlier_p, size(X, 2));
    euk_back = class == 0 & in_det & ~has_pe & F.chlSSC > euk_cut & F.size > syn_size + sig.euk_rescue_size_above_syn;
    class(euk_back) = 1;
    info.n_euk_rescued = sum(euk_back);
end
% events the mixture put with Syn but that failed the Syn check (too little PE
% for their chl, or too large) are re-checked as Euks: chl normal for a Euk of their
% size (within euksyn_max_chl_below_euk of the Euk chl-size line) and no smaller
% than Euks may be -> Euk with Syn (5) if they carry PE with a Syn-like spectrum,
% else Euk (1)
if any(syn_rejected) && sum(class == 1) >= 50
    be = robustfit(F.size(class == 1), F.BL3(class == 1));
    eukish = syn_rejected & F.BL3 - (be(1) + be(2) * F.size) > -sig.euksyn_max_chl_below_euk & ...
        F.size > syn_size - sig.euk_min_size_below_syn & F.chlSSC > sig.detritus_max_chlSSC;
    pe_syn = F.GL2 > pe_min & F.PEchl > median(F.PEchl(class == 1)) + sig.euksyn_min_PEchl_above_euk & ...
        (~gl3ok | (F.GL2GL3 > syn_range(1) - sig.syn_GL2GL3_tol & F.GL2GL3 < syn_range(2) + sig.syn_GL2GL3_tol));
    % events with Syn-level PE per chl (Euks have ~1000x less) that failed only the
    % spectrum test are not plain Euks: leave them for the catch-all
    synpe = has_pe & F.PEchl > syn_PEchl - sig.syn_PEchl_tol;
    % Euk with Syn here must meet the same rules as from the PE-euk step: clearly
    % larger than Syn and clearly less PE per chl than Syn (else they stay unassigned)
    euksyn_ok = F.size > syn_size + sig.peeuk_size_above_syn & F.PEchl < syn_PEchl - sig.peeuk_PEchl_below_syn;
    class(eukish & pe_syn & euksyn_ok) = 5;
    class(eukish & ~pe_syn & ~synpe) = 1;
    info.n_syn_rejected_rechecked = [sum(eukish & pe_syn & euksyn_ok) sum(eukish & ~pe_syn & ~synpe) sum(eukish & ~pe_syn & synpe)];
    info.trace.synpe_held = eukish & ~pe_syn & synpe;
end
% chl gap at the bottom of the Euks (Nov 2026): faint events that just clear the
% chl trigger pile up at its edge. Where the chl distribution of small no-PE
% particles (Euks and unassigned, smaller than Syn + euk_gap_size_above) has a deep
% minimum between that pile and the Euk peak, small Euks below the minimum are left
% unassigned (and the catch-all may not re-add them). Files whose smallest Euks run
% smoothly down to the trigger have no such minimum and keep the rules above.
% The population below the minimum must peak within euk_gap_pile_max_above_edge of the
% chl edge (Nov 2026): a cluster well above the edge is not trigger noise. In MVCO
% MV337N04 (preserved) the window took in the large Euks as the "Euk peak", and the
% whole small-Euk cluster (peak ~1.0 above the edge) fell below the minimum between
% the two; underway trigger-noise piles peak 0.37-0.60 above the edge
euk_below_gap = false(n, 1);
info.euk_chl_gap = -Inf; info.euk_gap_peak = NaN;
if isfield(edge, 'BL3')
    cand = ismember(class, [0 1]) & ~has_pe & F.size < syn_size + sig.euk_gap_size_above & F.BL3 > edge.BL3 - 0.05 & ...
        F.SSC > finfo.floor.SSC + 0.1;
    x = F.BL3(cand);
    if numel(x) >= sig.euk_gap_min_n
        [fk, gk] = ksdensity(x, 'NumPoints', 400, 'Bandwidth', sig.euk_gap_bw);
        % Euk peak: the highest density at least euk_gap_min_above_edge above the edge
        up = find(gk >= edge.BL3 + sig.euk_gap_min_above_edge);
        if ~isempty(up)
            [fm, ip] = max(fk(up)); ip = up(ip); best = Inf;
            % a large trigger pile can reach past euk_gap_min_above_edge: the window's
            % highest density is then the pile's falling slope at the window's lower
            % end, not a peak, and the gap above it is never searched (AR39B
            % 08Oct2019d Sample(113): pile 5x the Euk peak, gap at BL3 0-0.3 between
            % them, 670 pile events named Euk). Use the highest local maximum above it
            if sig.euk_gap_peak_local && ip == up(1) && ip > 1 && fk(ip - 1) > fk(ip)
                lm = up(2:end - 1); lm = lm(fk(lm) >= fk(lm - 1) & fk(lm) >= fk(lm + 1));
                if ~isempty(lm), [fm, k] = max(fk(lm)); ip = lm(k); end
            end
            info.euk_gap_peak = gk(ip);
            for j = ip - 1:-1:2
                if gk(j) <= edge.BL3, break, end
                if fk(j) <= fk(j - 1) && fk(j) <= fk(j + 1) && fk(j) < sig.euk_gap_frac * fm && ...
                        max(fk(1:j - 1)) > sig.euk_gap_rise * fk(j) && sum(x < gk(j)) >= sig.euk_gap_min_pile && fk(j) < best && ...
                        pile_near_edge(fk(1:j - 1), gk(1:j - 1), edge.BL3, sig.euk_gap_pile_max_above_edge)
                    best = fk(j); info.euk_chl_gap = gk(j);
                end
            end
        end
    end
    euk_below_gap = F.BL3 < info.euk_chl_gap & F.size < syn_size + sig.euk_gap_size_above & ~has_pe;
    class(class == 1 & euk_below_gap) = 0;
end
info.n_euk_below_gap = sum(euk_below_gap & info.trace.after_gmm == 1);
% Pro have no PE: GL2 must be at the noise level seen in this file's Euks
% (which also lack PE), else the file's GL2 floor plus a default
pro = class == 7;
if any(pro)
    if sum(class == 1) >= 20
        pro_GL2_max = median(F.GL2(class == 1)) + sig.pro_GL2_above_euk;
    else
        pro_GL2_max = finfo.floor.GL2 + sig.pro_GL2_default_above_floor;
    end
    % ... and well outside the Pro cluster's own GL2 spread, so the cut falls in the
    % sparse tail instead of slicing through the top of the cluster
    pro_GL2_max = max(pro_GL2_max, median(F.GL2(pro)) + sig.pro_GL2_sd * 1.4826 * mad(F.GL2(pro), 1));
    class(pro & F.GL2 > pro_GL2_max) = 0;
    info.pro_GL2_max = pro_GL2_max;
end
% Pro are small single cells: chl pulse no wider than pro_max_AH_above_ref above the
% pulse reference (see 'wide' above), and no wider than the Pro's own pulse-shape
% spread allows: median + pro_AH_sd x the spread on the narrow side (which a
% wide-pulse group above cannot stretch). Summer Pro near the chl trigger have a
% broad, single-peaked A/H distribution whose tail the fixed cut sliced (AR92,
% HRS2609: 280-430 old-class Pro per file); tight fall Pro (AR95, spread ~0.05) keep
% a cut near the fixed one, below their separate wide group at +0.3 to +0.45 (Nov 2026)
pro = class == 7;
pro_AH_max = info.pulse_ref_AH + sig.pro_max_AH_above_ref;
if any(pro) && isfinite(pro_AH_max)
    a = F.BL3_AH(pro); am = median(a);
    pro_AH_max = max(pro_AH_max, am + sig.pro_AH_sd * (am - prctile(a, 15.87)));
end
info.pro_AH_max = pro_AH_max;
pro_wide_pulse = F.BL3_AH > pro_AH_max; %false where pro_AH_max is NaN
info.trace.pro_wide_pulse = pro_wide_pulse;
pro_wide = pro & pro_wide_pulse;
class(pro_wide) = 0;
info.n_pro_wide = sum(pro_wide);
info.trace.after_pro_pe = class;

% PE-containing euks are often too rare to get their own component; pick them
% out event by event: PE present, larger than this file's Syn, and either less
% PE per chlorophyll than Syn or a cryptophyte-type PE spectrum
% PE spectrum relative to this file's Syn from BL2:GL3 (blue- vs green-excited PE;
% both stay on scale where GL2 saturates, and cryptophyte PE sits well below Syn).
% Bright cryptos near GL2 saturation can have Syn-like PE:chl and a GL2:GL3 that is
% compressed, so the BL2:GL3 test is needed to recognise them
% physical channel -> file column (attune_config map; the PE channel is GL1-H before Aug 2019)
satch = @(ch) sat_col(fcsdat, finfo.col.(ch), sig.sat_frac);
satG2 = satch('GL2');
okB = ~satch('BL2') & ~satch('GL3') & F.BL2 > finfo.floor.BL2 + 1 & F.GL3 > finfo.floor.GL3 + 1;
syB = class == 2 & okB & ~satG2;
% crypto-type spectrum: BL2:GL3 below this file's Syn by more than the larger of
% peeuk4_dBL2GL3 and crypto_dB_sd robust SDs of the Syn spectrum. Syn spectra are
% tight in winter (SD ~0.08) but wide in summer/fall when pigment types mix
% (SD up to ~0.3), where a fixed 0.3 offset is ordinary Syn spread
dB_cut = sig.peeuk4_dBL2GL3;
if sum(syB) >= 50
    dB = (F.BL2 - F.GL3) - median(F.BL2(syB) - F.GL3(syB));
    dB_cut = min(dB_cut, -sig.crypto_dB_sd * 1.4826 * mad(dB(syB), 1));
else
    dB = nan(n, 1);
end
% PE euks: larger than Syn with clearly less PE per chl than Syn, or with a
% crypto-type spectrum; the spectrum route also needs PE per chl below the Syn
% median (Syn-like PE:chl on a larger particle is Syn or Syn on a particle), except
% where GL2 saturates and PE:chl is compressed (bright cryptos)
spec = F.GL2GL3 < sig.peeuk_max_GL2GL3 | (okB & dB < dB_cut);
spec_ok = satG2 | F.PEchl < syn_PEchl - sig.peeuk_spec_PEchl_below_syn;
% PE present on a large particle: clearly above the PE-channel signal of this
% file's large Euks, which is chl spillover (no PE). Their PE:chl falls with chl,
% so the floor is the Euks' PE:chl-vs-chl trend plus peeuk_spill_sd robust SDs
% (fixed pe_min_PEchl / peeuk_min_GL2 cuts ran through a separate PE band in
% summer files). Fixed cuts are the fallback when there are too few large Euks
big = F.size > syn_size + sig.peeuk_size_above_syn;
eub = class == 1 & big & F.BL3 > chl_min;
pe_big = has_pe & F.GL2 > sig.peeuk_min_GL2;
info.peeuk_spill = []; info.trace.peeuk_fit = eub;
if sum(eub) >= sig.peeuk_spill_min_n
    bs = robustfit(F.BL3(eub), F.PEchl(eub));
    rs = F.PEchl(eub) - (bs(1) + bs(2) * F.BL3(eub));
    spill_sd = 1.4826 * mad(rs, 1);
    % the large-Euk set can hold PE euks (and particles) beside the Euks (AT46 C17N02:
    % 290 hand PE euks and 434 particles beside 736 Euks, SD 0.47 vs 0.20 for the
    % Euks; floor above every crypto). Where the residuals show a separate upper mode
    % above the dominant (spillover) one, the floor goes no higher than the valley
    % between them ('valley'). Two Euk modes without PE (EN655, EN668 summer) are
    % closer (0.4-0.6) with a shallow valley and leave the floor as it is
    spill_off = sig.peeuk_spill_sd * spill_sd; spill_cut = Inf;
    if strcmp(sig.peeuk_spill_method, 'valley')
        spill_cut = spill_valley(rs, sig);
        spill_off = min(spill_off, spill_cut);
    end
    pe_big = F.GL2 > pe_min & F.PEchl > bs(1) + bs(2) * F.BL3 + spill_off;
    info.peeuk_spill = [bs' spill_sd spill_cut];
end
pe_euk = pe_big & F.size > syn_size + sig.peeuk_size_above_syn & ...
    F.chlSSC > sig.detritus_max_chlSSC & ...
    (F.PEchl < syn_PEchl - sig.peeuk_PEchl_below_syn | (spec & spec_ok));
if sig.peeuk_sat_spectrum
    % chl saturated: chl is capped, so PE:chl is overstated, and the brightest Euks'
    % chl spillover alone saturates GL2 (EN627 C45N05: BL3 at its ceiling, BL2/GL3 on
    % scale with a spillover/Syn-like BL2:GL3). Spillover and Syn PE share that
    % spectrum, so only a crypto-type spectrum shows PE there
    pe_euk = pe_euk & ~(satch('BL3') & okB & ~(dB < dB_cut));
end
info.crypto_dB_cut = dB_cut;
class(pe_euk) = 3;
class(class == 3 & ~pe_euk) = 0; %component-level PE euk calls must pass the same test
info.trace.after_peeuk = class;

% split PE euks (Nov 2026) by PE spectrum and amount:
%  5 = Euk with Syn: Syn-like PE spectrum, Euk-level chl, PE ~ one Syn's worth
%      (grazers with ingested Syn or Syn attached; not chance coincidences);
%      those with much less chl than Euks of their size go to 8 (Syn on a particle)
%  3 = dim cryptophytes and 4 = bright cryptophytes: cryptophyte-type spectrum
%      (BL2:GL3 below Syn); in winter/fall these form two compact groups that differ
%      mainly in PE and spectrum (dim ones more red-shifted; also ~3x less scatter), split per file when the
%      two groups are clearly there; otherwise split by a PE floor (dimmer = 3)
% The spectrum is read from BL2:GL3 (blue- vs green-excited PE; both stay on scale
% where GL2 saturates) relative to this file's Syn; GL2:GL3 where BL2 or GL3 is
% unusable; events saturated in GL2 with no usable BL2/GL3 are crypto-type.
% PE amount: GL2, or BL2 scaled to GL2 (from this file's Syn) where GL2 saturates.
pe3 = class == 3;
info.crypto_split = struct('n', 0, 'accepted', false);
if any(pe3)
    satcr = satG2; %GL2-saturated events taken as crypto-type (see peeuk_sat_spectrum below)
    sy = class == 2 & okB & ~satG2;
    if sum(sy) >= 50
        dB = (F.BL2 - F.GL3) - median(F.BL2(sy) - F.GL3(sy));
        dG = F.GL2GL3 - median(F.GL2GL3(sy));
        % spectral cuts relative to this file's Syn spread (see dB_cut above)
        dB_cut = min(sig.peeuk4_dBL2GL3, -sig.crypto_dB_sd * 1.4826 * mad(dB(sy), 1));
        dG_cut = min(sig.peeuk4_dGL2GL3, -sig.crypto_dB_sd * 1.4826 * mad(dG(sy), 1));
        info.crypto_dB_cut = dB_cut;
        % Syn can be a mixture of pigment types (e.g. a second type ~0.7 below the
        % main one in BL2:GL3 in HRS2609 17Aug(11)), so a spectrum below the Syn
        % median is not enough: crypto-type means below every Syn type in the file
        % (the low end of the Syn spectrum distribution, minus a margin); otherwise
        % the PE is Syn-type (Euk with Syn). GL2-saturated PE euks are bright
        % cryptos regardless of spectrum (set after the split below)
        info.syn_ref_BL2GL3 = median(F.BL2(sy) - F.GL3(sy));
        PE = F.GL2; k = satG2 & okB;
        PE(k) = F.BL2(k) - median(F.BL2(sy) - F.GL2(sy));
        PEs = PE - median(F.GL2(sy));          %log10 Syn equivalents
        dB_low = min(dB_cut, prctile(dB(sy), sig.crypto_syn_low_prc) - sig.crypto_syn_low_margin);
        dG_low = min(dG_cut, prctile(dG(sy), sig.crypto_syn_low_prc) - sig.crypto_syn_low_margin);
        info.crypto_dB_low = dB_low;
        % events with far more PE than a Syn (> crypto_min_PE_above_syn, i.e. >10 Syn;
        % Euk with Syn holds ~one Syn's worth) keep the median-relative spectrum cut:
        % unsaturated bright cryptos in summer files sit within the Syn spectral spread
        loose = (okB & dB < dB_cut) | (~okB & ~satG2 & gl3ok & dG < dG_cut);
        % GL2-saturated events are crypto-type unless (peeuk_sat_spectrum) BL2:GL3 is
        % usable and not crypto-type (bright Euks whose chl spillover saturates GL2)
        satcr = satG2; if sig.peeuk_sat_spectrum, satcr = satG2 & ~(okB & ~(dB < dB_cut)); end
        crypto = (okB & dB < dB_low) | (~okB & ~satG2 & gl3ok & dG < dG_low) | satcr | (loose & PEs > sig.crypto_min_PE_above_syn);
    else
        crypto = F.GL2GL3 < sig.peeuk_max_GL2GL3 | satG2;
        dB = nan(n, 1); PEs = F.GL2 - syn_PEchl - syn_BL3;
    end
    cr = pe3 & crypto;
    class(pe3 & ~crypto) = 5;
    class(cr) = 4;
    % two crypto groups? fit two components on [PE, spectrum] and accept only if
    % both are populated, compact, well separated, and the dimmer one is the more
    % red-shifted (as in EN727 / AR95)
    % the dim/bright fit uses events with a crypto-type spectrum (saturated events with
    % compressed, Syn-like spectra would widen the bright group; they are set to 4 below)
    cspec = cr & okB & isfinite(dB) & dB < dB_cut;
    fitset = cspec;
    info.crypto_split.n = sum(cr);
    if sum(fitset) >= sig.crypto_split_min_n
        Xc = [PEs(fitset) dB(fitset)];
        try
            rng(sd + 7, 'twister')
            gc = fitgmdist(Xc, 2, 'Replicates', 3, 'RegularizationValue', 1e-3, 'Options', opts);
            [~, o] = sort(gc.mu(:, 1));            %o(1) dimmer, o(2) brighter
            S1 = gc.Sigma(:, :, o(1)); S2 = gc.Sigma(:, :, o(2));
            dm = gc.mu(o(2), :) - gc.mu(o(1), :);
            sep = sqrt(dm / ((S1 + S2) / 2) * dm');
            nmin = min(gc.ComponentProportion) * size(Xc, 1);
            sdc = sqrt([S1(1, 1) S2(1, 1); S1(2, 2) S2(2, 2)]); %rows PE, dB; columns dim, bright
            if sig.crypto_robust_compact
                % robust SD of each group's hard-assigned members: a few stragglers
                % inflate the Gaussian SD past the limit and veto a clean split
                % (EN727 27Jan2025(4) Sample(41): dim-group dB SD 0.217, robust ~0.1)
                pf = posterior(gc, Xc); [~, hard] = max(pf, [], 2);
                for j = 1:2
                    mj = hard == o(j);
                    if sum(mj) >= 5, sdc(:, j) = 1.4826 * mad(Xc(mj, :), 1)'; end
                end
            end
            compact = all(sdc(1, :) <= sig.crypto_max_sd_PE) && all(sdc(2, :) <= sig.crypto_max_sd_dB);
            redder = gc.mu(o(1), 2) < gc.mu(o(2), 2) - sig.crypto_min_dB_offset;
            ok = nmin >= sig.crypto_min_group_n && compact && sep >= sig.crypto_min_sep && redder;
            info.crypto_split = struct('n', sum(cr), 'accepted', ok, 'sep', sep, 'nmin', nmin, 'compact', compact, 'redder', redder, ...
                'mu_dim', gc.mu(o(1), :), 'mu_bright', gc.mu(o(2), :));
            info.crypto_gmm = gc; info.crypto_gmm_order = o;
            if ok
                % events with a usable spectrum by posterior; the rest by PE alone
                post = posterior(gc, [PEs dB]);
                small = cr & okB & isfinite(dB) & post(:, o(1)) > 0.5;
                small = small | (cr & ~(okB & isfinite(dB)) & abs(PEs - gc.mu(o(1), 1)) < abs(PEs - gc.mu(o(2), 1)));
                class(small) = 3;
            end
        catch
        end
    end
    if ~info.crypto_split.accepted
        % no two separate groups. If the crypto-type events form one compact group,
        % label the whole group by where its median PE sits (a PE floor would cut a
        % single group in two); if they are diffuse (summer), split by PE event by event
        one = cspec;
        single = sum(one) >= sig.crypto_min_group_n && 1.4826 * mad(PEs(one), 1) <= sig.crypto_max_sd_PE && ...
            1.4826 * mad(dB(one), 1) <= sig.crypto_max_sd_dB;
        info.crypto_split.single_group = single;
        if single
            % whole group by its median PE; events far outside the group's PE spread
            % (both sides) take the other label, so a few dim events are not called
            % bright with a bright group, nor bright ones dim with a dim group
            gmed = median(PEs(one)); gsd = 1.4826 * mad(PEs(one), 1);
            if gmed < sig.crypto_min_PE_above_syn, class(cr) = 3; end
            class(cr & PEs < gmed - sig.crypto_outlier_sd * gsd) = 3;
            class(cr & PEs > gmed + sig.crypto_outlier_sd * gsd) = 4;
        else
            class(cr & PEs < sig.crypto_min_PE_above_syn) = 3;
        end
    end
    % PE euks at or near GL2 saturation are bright cryptos: no other class carries
    % that much PE (Syn on a particle and Euk with Syn hold ~one Syn's worth), and
    % saturation compresses their PE:chl and GL2:GL3 toward Syn-like values
    % (with peeuk_sat_spectrum, only those without a usable non-crypto spectrum)
    class(pe3 & satcr) = 4;
    % Euk-with-Syn events far below the Euks' chl-size line are not Euks carrying
    % a Syn but aggregates / detritus with a Syn and a little chl (two groups in
    % EN727 and AR92, ~1.2 and ~0.2 below the line): count them as Syn on a particle
    if any(class == 5) && sum(class == 1) >= 50
        b = robustfit(F.size(class == 1), F.BL3(class == 1));
        euksyn_lowchl = class == 5 & F.BL3 - (b(1) + b(2) * F.size) < -sig.euksyn_max_chl_below_euk;
        class(euksyn_lowchl) = sig.euksyn_low_chl_class; %8 = Syn on particle, or 0 = leave unassigned
        info.trace.euksyn_lowchl = euksyn_lowchl;
        info.trace.from_euksyn = class == 8; %Euk with Syn moved to Syn on particle (trace only)
    end
    % class 6, PE euk cluster (Nov 2026): Euk-with-Syn events that are one compact
    % population with their own PE. AT46 C17N02: 270 class-5 events with a PE spectrum
    % outside nearly all the file's Syn (GL2:GL3 0.62 vs Syn 0.29, 99th pct 0.59), PE
    % rising with chl (slope 0.94, r 0.88) and PE:chl robust SD 0.13; ingested Syn
    % (Poisson numbers of the file's Syn) would give PE:chl SD 0.23-0.41 and a Syn
    % spectrum. The same spectrum (GL2:GL3 0.55-0.61, BL2-GL2 -0.34 to -0.42) recurs
    % in EN727 C20N02, AR43 C02N24, MVCO MV417N04 and MV432N04 (preserved); class 5
    % in live files is diffuse (PE:chl SD 0.4-1.6, PE unrelated to chl). Class-5
    % events with usable spectra: if their PE:chl and GL2:GL3 are tight (robust SD),
    % the events near the centre show PE proportional to chl, and their GL2:GL3 sits
    % clearly above the file's Syn, those events take cl.class
    % What they are (27 Nov 2026): Euks whose PE-channel fluorescence is raised by
    % preservation, not a PE population. The live file nearest in time (AT46 C17N02,
    % AR43 C02N24, EN727 C02N02 and C20N12, 1-75 min apart) has the same Euks (same
    % size and chl) with GL2 1.5-2.3 log10 lower and GL2:GL3 -0.3 to -0.5, and no events
    % where class 6 sits; all preserved Euks shift up (Syn ~0.3). Downstream (tables,
    % scoring) count class 6 as Euk
    info.peeuk_cluster = struct('n', 0, 'accepted', false);
    cl = sig.peeuk_cluster;
    c5 = class == 5 & okB & ~satG2 & isfinite(dB);
    if cl.on && sum(c5) >= cl.min_n && sum(sy) >= 50
        cP = median(F.PEchl(c5)); sP = 1.4826 * mad(F.PEchl(c5), 1);
        cG = median(F.GL2GL3(c5)); sG = 1.4826 * mad(F.GL2GL3(c5), 1);
        mem = c5 & abs(F.PEchl - cP) < cl.win_sd * sP & abs(F.GL2GL3 - cG) < cl.win_sd * sG;
        pc = struct('n', sum(mem), 'accepted', false, 'sd_PEchl', sP, 'sd_GL2GL3', sG, 'PEchl', cP, 'GL2GL3', cG, ...
            'dGL2GL3_syn', cG - median(F.GL2GL3(sy)), 'slope', NaN, 'r', NaN);
        if sum(mem) >= cl.min_n
            bp = robustfit(F.BL3(mem), F.GL2(mem)); pc.slope = bp(2); pc.r = corr(F.BL3(mem), F.GL2(mem));
            pc.accepted = sP <= cl.max_sd_PEchl && sG <= cl.max_sd_GL2GL3 && pc.slope >= cl.slope(1) && pc.slope <= cl.slope(2) && ...
                pc.r >= cl.min_r && pc.dGL2GL3_syn >= cl.min_dGL2GL3_syn;
            if pc.accepted, class(mem) = cl.class; end
        end
        info.peeuk_cluster = pc;
    end
end

% particle-associated Syn: Syn spectrum and about one Syn's pigment on a much
% larger particle (aggregate, detritus, or a grazer that has eaten Syn). These
% form a continuous tail in scatter rather than a cluster, so they are picked
% event by event; no upper size limit, and exempt from the chl:size detritus rule
if any(class == 2)
    % core Syn centre and spread from the peak of the size distribution
    [syn_core, syn_size_sd] = peak_centre(F.size(class == 2));
    synpart_cut = syn_core + max(sig.synpart_min_size_above, sig.synpart_size_sd * syn_size_sd);
    info.syn_core_size = [syn_core syn_size_sd];
    % Syn chl-size trend from the core cluster: large Syn follow it (chl rises with
    % scatter along the phyto diagonal); a Syn on a larger particle carries about one
    % Syn's chl on more scatter, so it falls clearly below the trend extended
    corez = class == 2 & abs(F.size - syn_core) < 2 * syn_size_sd;
    if sum(corez) >= 50
        bt = robustfit(F.size(corez), F.BL3(corez));
    else
        bt = [syn_BL3 - syn_core; 0];
    end
    below_trend = F.BL3 < bt(1) + bt(2) * F.size - sig.synpart_below_trend;
    info.syn_chl_trend = bt';
    synpart_rule = has_pe & ...
        abs(F.PEchl - syn_PEchl) < sig.synpart_PEchl_tol & (~gl3ok | (F.GL2GL3 > syn_range(1) - sig.synpart_GL2GL3_tol & F.GL2GL3 < syn_range(2) + sig.synpart_GL2GL3_tol)) & ...
        F.BL3 > syn_BL3 - sig.synpart_BL3_below & F.BL3 < syn_BL3 + sig.synpart_BL3_above & ...
        F.size > synpart_cut & below_trend;
    % ... with chl above the chl trigger edge: about one Syn's chl must be measured, not
    % cut by the trigger. Where the chl trigger cuts through the Syn cluster, PE-triggered
    % large particles with faint PE at the trigger limits formed a band of false Syn on
    % particle (EN617 20Jul2018A Sample(1): 2,201 of 2,269 below the chl edge) (Nov 2026)
    if isfield(edge, 'BL3'), synpart_rule = synpart_rule & F.BL3 > edge.BL3; end
    synpart = ismember(class, [0 2]) & synpart_rule;
    class(synpart) = 8;
    info.trace.synpart_main = synpart; %(trace only)
    % larger than the core Syn, below the Syn chl trend, but not a Syn on a particle
    % (e.g. too little chl even for one Syn): not Syn. Large Syn on the trend stay Syn
    class(class == 2 & F.size > synpart_cut & below_trend) = 0;
end
% catch-all (Nov 2026): no event that looks like a cell is left unassigned. An
% event is cell-like if it carries PE like a cell or chl like a cell for its size
% (above the detritus cut and, for particles well larger than Syn, above this
% file's Euk chl:size cut). Unassigned cell-like events take the majority class of
% their nearest assigned neighbours in [size, chl, PE, BL2, GL3], if they also meet
% that class's own rules (sizes relative to Syn; never Pro, see below); otherwise
% they stay unassigned. The Pro tail at the chl trigger in files run without Pro
% settings stays unassigned.
if exist('euk_cut', 'var'), ccut = euk_cut; else, ccut = sig.detritus_max_chlSSC; end
cell_chl = has_chl & F.chlSSC > sig.detritus_max_chlSSC & (F.size <= syn_size + sig.euk_rescue_size_above_syn | F.chlSSC > ccut);
% PE like a cell: at least ~1/3 of a Syn's PE (detected-but-faint PE near the GL2
% trigger is common on non-phyto particles)
if any(class == 2)
    cell_pe = has_pe & F.GL2 > median(F.GL2(class == 2)) - sig.cell_PE_below_syn;
else
    cell_pe = has_pe & F.GL2 > pe_min + sig.cell_PE_below_syn;
end
if ~isfield(info.trace, 'euksyn_lowchl'), info.trace.euksyn_lowchl = false(n, 1); end
left = class == 0 & (cell_pe | cell_chl) & ~euk_below_gap & ~pile_ev & ~(info.trace.euksyn_lowchl & sig.euksyn_low_chl_class == 0);
refset = ismember(class, [1 2 3 4 5 6 7 8]);
info.n_catchall = 0;
if any(left) && sum(refset) >= 50
    Z = [F.size F.BL3 F.GL2 F.BL2 F.GL3];
    zs = std(Z(refset, :)); zs(zs == 0) = 1;
    Z = (Z - mean(Z(refset, :))) ./ zs;
    ir = find(refset);
    if numel(ir) > sig.catchall_max_ref, ir = ir(randperm(numel(ir), sig.catchall_max_ref)); end
    nn = knnsearch(Z(ir, :), Z(left, :), 'K', min(sig.catchall_k, numel(ir)));
    nl = mode(reshape(class(ir(nn)), size(nn)), 2);
    % the neighbours' class is taken only if the event also meets that class's own
    % rules; tiny signals far smaller than Syn, for example, are neither Euks nor Syn
    fz = F(left, :);
    smin_euk = syn_size - sig.euk_min_size_below_syn;
    ok = false(size(nl));
    ok(nl == 1 | nl == 5) = fz.size(nl == 1 | nl == 5) > smin_euk;
    ep = euk_prolike(left); ok((nl == 1 | nl == 5) & ep) = false; %small, dim, no PE: Pro-like, not Euk
    % catch-all Euks carry no Syn-like PE either (rule above)
    hpz = has_pe(left);
    ok(nl == 1) = ok(nl == 1) & ~(hpz(nl == 1) & fz.PEchl(nl == 1) > syn_PEchl - sig.syn_PEchl_tol);
    % no catch-all Euk may have less chl than the Euks the main rules accepted: at
    % any size, chl must reach the lower edge of the main-rule Euks about their
    % chl-size line (catchall_euk_chl_prc percentile of the residuals), and chl per
    % scatter must pass the Euk chl:scatter cut
    if sum(class == 1) >= 50
        e1 = class == 1;
        bl = robustfit(F.size(e1), F.BL3(e1));
        rlo = prctile(F.BL3(e1) - (bl(1) + bl(2) * F.size(e1)), sig.catchall_euk_chl_prc);
        lowe = (nl == 1 | nl == 5) & (fz.BL3 < bl(1) + bl(2) * fz.size + rlo | fz.chlSSC < ccut);
        % Euk with Syn also within euksyn_max_chl_below_euk of the Euk chl-size line
        lowe = lowe | (nl == 5 & fz.BL3 < bl(1) + bl(2) * fz.size - sig.euksyn_max_chl_below_euk);
        ok(lowe) = false;
        info.catchall_euk_floor = [bl' rlo];
    end
    if exist('synpart_cut', 'var'), smax2 = synpart_cut; else, smax2 = syn_size + sig.syn_size_tol; end
    ok(nl == 2) = fz.size(nl == 2) > syn_size - min(sig.syn_size_max_below, max(sig.syn_size_tol_below, sig.syn_size_sd * syn_size_sdv)) & fz.size(nl == 2) < smax2 & ...
        fz.GL2(nl == 2) >= syn_PE_floor & fz.PEchl(nl == 2) > syn_PEchl - sig.syn_PEchl_tol & fz.PEchl(nl == 2) < syn_PEchl + max(sig.syn_PEchl_tol, sig.syn_PEchl_sd * syn_PEchl_sdv);
    ok(nl == 3 | nl == 4) = fz.size(nl == 3 | nl == 4) > syn_size + sig.peeuk_size_above_syn;
    % ... with chl like a cell, as the main PE-euk rule (Nov 2026): faint PE on large detritus
    % (chl:size about -1) took PE-euk neighbours in EN661 06Feb(5) Sample(5)
    if sig.catchall_peeuk_det, ok(nl == 3 | nl == 4) = ok(nl == 3 | nl == 4) & fz.chlSSC(nl == 3 | nl == 4) > sig.detritus_max_chlSSC; end
    % the catch-all adds no Pro (Nov 2026): Pro come only from components named Pro.
    % Catch-all Pro were mostly faint events just above the chl trigger, smaller than
    % the file's Pro (8,153 in 77 underway Pro-setting files, 4% of them old-class Pro,
    % against 91% for main-rule Pro), and most of v2's false Pro in preserved EN655
    % C15N09 (1,337 of 1,603, from the trigger-floor pile). Unassigned events whose
    % neighbours are Pro stay unassigned
    ok(nl == 7) = false;
    if exist('synpart_cut', 'var'), c8 = synpart_cut; else, c8 = syn_size + sig.synpart_min_size_above; end
    ok(nl == 8) = fz.size(nl == 8) > c8;
    % ... and the rest of the Syn-on-particle rules (Syn-like PE:chl and spectrum, about
    % one Syn's chl, below the Syn chl trend); Euk with Syn: the PE-euk size and PE:chl rules
    if exist('synpart_rule', 'var'), sr = synpart_rule(left); ok(nl == 8) = ok(nl == 8) & sr(nl == 8); end
    ok(nl == 5) = ok(nl == 5) & fz.size(nl == 5) > syn_size + sig.peeuk_size_above_syn & fz.PEchl(nl == 5) < syn_PEchl - sig.peeuk_PEchl_below_syn;
    nl(~ok) = 0;
    class(left) = nl;
    left(left) = ok;
    info.n_catchall = sum(left);
    info.catchall_classes = arrayfun(@(x) sum(left & class == x), [1 2 3 4 5 7 8]);
end
info.trace.catchall = left;
% Euk with Syn carries a whole Syn, so its PE (GL2) is at least one Syn's: below the
% euksyn_min_syn_PE_prc percentile of this file's Syn PE it is a Euk with faint PE,
% not Euk with Syn (Nov 2026: AR92, half the Euk with Syn sat below 99% of the file's
% Syn in PE). Applied after every route into class 5 (Syn-check rejects, PE euks,
% catch-all) and before the Euk chl:size cut, which then also checks these Euks
info.euksyn_min_GL2 = NaN; info.n_euksyn_weak = 0;
if sig.euksyn_min_syn_PE_prc > 0 && sum(class == 2) >= 50
    info.euksyn_min_GL2 = prctile(F.GL2(class == 2), sig.euksyn_min_syn_PE_prc);
    weak = class == 5 & F.GL2 < info.euksyn_min_GL2;
    class(weak) = 1; info.n_euksyn_weak = sum(weak);
end
% Pro in runs without Pro settings (Nov 2026): their chl trigger cuts the Pro cluster,
% and the part above it was left unassigned, by the Euk chl gap (flag 2: AE2426
% 08Nov alternating Euk runs, ~30% of the neighbouring Pro runs' Pro lie above the Euk
% run's edge) or as Pro-like (flag 4). Candidates are those unassigned events (Euks
% are not taken: in summer files the small Euks are one cluster reaching Syn size,
% which fixed limits would only slice): no PE, at or above the chl edge, size from
% Syn - pro_eukrun_size_below to Syn + pro_eukrun_size_above, chl:size between the
% detritus cut and pro_eukrun_max_chlSSC, pulse no wider than the Pro limit, not
% the capped trigger pile. As a population they are either the trigger-noise pile
% (edge_pop_is_pile: short chl pulses, or chl collapsed on the edge; flag 3) or
% cells: class 7, with info.pro_lower_bound set since the trigger cuts the cluster.
% Single events with very short pulses (pile_AH_below_ref_event below the reference)
% stay unassigned (flag 3). AE2426 8-10 Nov, 460 Euk runs: Euk-run Pro / the
% neighbouring Pro runs' Pro above the Euk-run edge, median 0.99; Euk unchanged
info.pro_eukrun = struct('n_cand', 0, 'pile', false, 'n_pro', 0, 'stats', {[]});
info.pro_lower_bound = false;
pro_eukrun_pile = false(n, 1);
if ~pro_measured && sig.pro_in_euk_runs && isfield(edge, 'BL3') && any(class == 2) && ~info.edge_pile_capped
    pe_c = class == 0 & ~has_pe & F.BL3 >= edge.BL3 & F.size <= syn_size + sig.pro_eukrun_size_above & ...
        F.size >= syn_size - sig.pro_eukrun_size_below & F.chlSSC > sig.detritus_max_chlSSC & F.chlSSC < sig.pro_eukrun_max_chlSSC & ...
        ~info.trace.wide_pulse & ~pile_ev;
    if isfield(info.trace, 'euksyn_lowchl'), pe_c = pe_c & ~info.trace.euksyn_lowchl; end
    info.pro_eukrun.n_cand = sum(pe_c);
    if sum(pe_c) >= sig.pro_eukrun_min_n
        [isp, st] = edge_pop_is_pile(F, pe_c, edge.BL3, info.pulse_ref_AH, sig);
        info.pro_eukrun.pile = isp; info.pro_eukrun.stats = st;
        if isp
            pro_eukrun_pile = pe_c;
        else
            narrow = pe_c & F.BL3_AH < info.pulse_ref_AH - sig.pile_AH_below_ref_event; %false where the reference is NaN
            class(pe_c & ~narrow) = 7; pro_eukrun_pile = narrow;
            info.pro_eukrun.n_pro = sum(pe_c & ~narrow);
            info.pro_lower_bound = true;
        end
    end
end
% class 10, small Euk-like cells with low chl for their size (Nov 2026; identity
% uncertain). AR77 11-13 Oct: a compact no-PE branch 0.35-0.6 above Syn in size at
% chl:size ~0.6 (bright Euks ~1.4), narrow chl pulses like Syn, 40-80% of the Euk-like
% events while present, almost absent from the matching preserved sample. The Euk
% chl:size cut's pile gate called it Euk in one file and unassigned (flag 5) in the
% next. Fixed limits relative to Syn, so the count follows the cells, not what else
% is in the file: no PE, size Syn + dim_euk_size_above(1) to (2), chl:size between the
% detritus cut and dim_euk_max_chlSSC; at least dim_euk_min_n of them, and cells, not
% the trigger pile. (A limit of 1.0 also took 10-30% of the Euks in summer cruises,
% whose small Euks are one cluster peaking at chl:size 1.0-1.4; 0.8 takes 2-8%. A
% test for a separate low mode instead made the AR77 count flip with the branch's
% shape: 33 jumps over 3x in 65 files, against 6 counting the same events by box.)
info.dim_euk = struct('n', 0, 'pile', false, 'stats', {[]});
if sig.dim_euk_class && isfield(edge, 'BL3') && any(class == 2)
    dc = class == 1 & ~has_pe & F.size >= syn_size + sig.dim_euk_size_above(1) & F.size <= syn_size + sig.dim_euk_size_above(2) & ...
        F.chlSSC > sig.detritus_max_chlSSC & F.chlSSC < sig.dim_euk_max_chlSSC;
    if sum(dc) >= sig.dim_euk_min_n
        [isp, st] = edge_pop_is_pile(F, dc, edge.BL3, info.pulse_ref_AH, sig);
        info.dim_euk.pile = isp; info.dim_euk.stats = st;
        if ~isp, class(dc) = 10; info.dim_euk.n = sum(dc); end
    end
end
% per-file Euk chl:size cut (Nov 2026; from assign_class_discrete, where preserved runs
% get it after their own Pro steps, so discrete_config turns this copy off). Faint
% Syn-sized particles near the chl trigger, with low chl for their size and PE:chl
% between Euk and Syn, can form their own small Euk components (live EN644
% 24Aug2019d Sample(41): 336 such Euks, chl:size 0.42 vs the file's Euks 1.46,
% about 0.3 above the chl edge, no chl gap below them). Reference: the file's clear
% Euks, at least euk_chl.core_above_edge above the BL3 trigger edge (>= core_min_n).
% Euks are unassigned (flag 5) where
%  - chl:size is below the clear Euks' median less below_core (cut no higher than
%    cut_max) and GL2:GL3 below max_GL2GL3, or
%  - BL3 is below the line BL3 = a + ridge_slope*size through the clear Euks by
%    ridge_k robust SDs, or
%  - within pe_edge_win of the chl edge with PE:chl pe_above_core above the clear Euks'
% Applied after the catch-all, so no unassigned event is re-added as Euk.
euk_chl_low = false(n, 1); info.euk_chl_cut = NaN; info.euk_ridge = NaN(1, 3); info.euk_chl_ceiling = NaN; info.euk_chl_pile_peak = NaN;
info.n_euk_lowchl = 0; info.n_euk_dim = 0; info.n_euk_highpe = 0;
if sig.euk_chl_in_v2 && sig.euk_chl.on && isfield(edge, 'BL3')
    ec = sig.euk_chl;
    ec_core = class == 1 & F.BL3 > edge.BL3 + ec.core_above_edge & isfinite(F.chlSSC);
    if sum(ec_core) >= ec.core_min_n
        info.euk_chl_cut = min(ec.cut_max, median(F.chlSSC(ec_core)) - ec.below_core);
        ec_low = class == 1 & F.chlSSC < info.euk_chl_cut & F.GL2GL3 < ec.max_GL2GL3;
        info.n_euk_lowchl = sum(ec_low);
        ec_b = ec.ridge_slope; ec_a = median(F.BL3(ec_core) - ec_b * F.size(ec_core));
        ec_rsd = max(ec.ridge_min_sd, 1.4826 * mad(F.BL3(ec_core) - (ec_a + ec_b * F.size(ec_core)), 1));
        info.euk_ridge = [ec_a ec_b ec_rsd];
        ec_dim = class == 1 & F.BL3 < ec_a + ec_b * F.size - ec.ridge_k * ec_rsd;
        info.n_euk_dim = sum(ec_dim & ~ec_low); ec_low = ec_low | ec_dim;
        ec_pc = median(F.PEchl(ec_core));
        ec_pe = class == 1 & F.BL3 < edge.BL3 + ec.pe_edge_win & F.PEchl > ec_pc + ec.pe_above_core;
        info.n_euk_highpe = sum(ec_pe & ~ec_low); ec_low = ec_low | ec_pe;
        % only where non-phyto particles pile up at the chl trigger (Nov 2026, live): the
        % low chl:size, no-PE particles (Euk or unassigned) must peak in chl within
        % pile_win of the BL3 edge, and only those below the first valley above that peak
        % are cut. Without a pile (HRS2303 05May2023: a separate low-chl cluster 0.6 above
        % the edge, old-class Euks) a per-file cut only splits a wide Euk cluster; above a
        % valley (AR82A, EN649: large low-chl:size cells 1-2 above the edge) the cells are
        % not part of the pile
        pl_ceil = -Inf;
        if ec.pile_gate
            % detritus by chl:size (below detritus_max_chlSSC) is already non-phyto and
            % well apart from the Euks; a cloud of it at the trigger is not a pile of
            % Euk-like particles (HRS2303 2May2023(2) Sample(13): large particles at
            % chl:size ~ -1 made the "pile", and a separate small-Euk group was cut)
            pl_lp = ismember(class, [0 1]) & ~has_pe & F.chlSSC < info.euk_chl_cut & F.chlSSC > sig.detritus_max_chlSSC & ...
                F.BL3 > edge.BL3 - 0.05 & F.SSC > finfo.floor.SSC + 0.1;
            pl_lp = pl_lp | ec_low;
            pl_x = F.BL3(pl_lp) - edge.BL3;
            if numel(pl_x) >= ec.pile_min_n
                [pl_f, pl_g] = ksdensity(pl_x, 'NumPoints', 300, 'Bandwidth', ec.pile_bw);
                [pl_fm, pl_ip] = max(pl_f);
                if pl_g(pl_ip) <= ec.pile_win
                    pl_ceil = Inf;
                    for pl_j = pl_ip + 1:numel(pl_f) - 1
                        if pl_f(pl_j) <= pl_f(pl_j - 1) && pl_f(pl_j) < pl_f(pl_j + 1) && max(pl_f(pl_j + 1:end)) >= ec.pile_valley_rise * pl_f(pl_j)
                            pl_ceil = pl_g(pl_j); break
                        end
                    end
                end
                info.euk_chl_pile_peak = pl_g(pl_ip);
            end
            ec_low = ec_low & F.BL3 - edge.BL3 < pl_ceil;
        else
            pl_ceil = Inf;
        end
        info.euk_chl_ceiling = pl_ceil;
        euk_chl_low = ec_low; class(ec_low) = 0;
    end
end
% rule boundaries for plotting (plot_v2_rules); values as used in this file
R = struct();
R.detritus_chlSSC = sig.detritus_max_chlSSC;
R.chl_min = chl_min; R.pe_min = pe_min;
R.euk_min_size = syn_size - sig.euk_min_size_below_syn;
if exist('euk_cut', 'var'), R.euk_chlSSC_cut = euk_cut; end
if exist('be', 'var'), R.euk_line = be'; end
R.euksyn_max_chl_below_euk = sig.euksyn_max_chl_below_euk;
R.euk_synpe_PEchl = syn_PEchl - sig.syn_PEchl_tol; %Euks with PE above this PE:chl are not Euks
R.peeuk_min_size = syn_size + sig.peeuk_size_above_syn;
R.peeuk_PEchl_max = syn_PEchl - sig.peeuk_PEchl_below_syn;
R.peeuk_max_GL2GL3 = sig.peeuk_max_GL2GL3;
R.peeuk_min_GL2 = sig.peeuk_min_GL2;
R.crypto_dB = info.crypto_dB_cut; if isfield(info, 'crypto_dB_low'), R.crypto_dB_low = info.crypto_dB_low; end
R.crypto_min_PE_above_syn = sig.crypto_min_PE_above_syn;
if exist('sy', 'var') && sum(sy) >= 50
    R.syn_BL2GL3 = median(F.BL2(sy) - F.GL3(sy));
    R.syn_BL2GL2 = median(F.BL2(sy) - F.GL2(sy));
    R.syn_GL2 = median(F.GL2(sy));
end
if exist('synpart_cut', 'var')
    R.synpart_cut = synpart_cut; R.syn_chl_trend = bt'; R.synpart_below_trend = sig.synpart_below_trend;
    R.synpart_PEchl = [syn_PEchl - sig.synpart_PEchl_tol, syn_PEchl + sig.synpart_PEchl_tol];
    R.synpart_BL3 = [syn_BL3 - sig.synpart_BL3_below, syn_BL3 + sig.synpart_BL3_above];
end
R.sat_frac = sig.sat_frac;
if isfield(info, 'rules'), R.syn_check = info.rules.syn_check; end
info.rules = R;
info.syn_ref = struct('PEchl', syn_PEchl, 'GL2GL3', syn_GL2GL3, 'BL3', syn_BL3, 'size', syn_size, 'GL2GL3_range', syn_range);
% Syn spectral subtype (Nov 2026): fixed windows on the bead-normalized excitation
% ratio BL2-GL2 (PUB:PEB proxy), the same in every file and cruise. Kept as a
% separate label so class 2 (and its counts downstream) is unchanged:
% 1/2/3/4 = types A/B/C/D in order of BL2-GL2, 5 = type ambiguous (mixture, below),
% 9 = Syn-bearing event whose spectrum is not usable, 0 = other. The windows are the
% fallback for files with few Syn. In dim Syn the ratio is noisy and spreads type B below the
% first edge (summer: up to half of the dimmest Syn), so for Syn whose BL2 or GL2 is
% less than syn_typeA_min_above_floor above its floor the type-A edge moves left by
% syn_typeA_dim_slope per log10 of shortfall
st = zeros(n, 1, 'uint8');
synlike = ismember(class, [2 5 8]);
% the type uses BL2 and GL2 only, so only those must be on scale (GL3 near its floor,
% common for dim summer Syn, does not matter here)
usable = ~satch('BL2') & ~satch('GL2') & ...
    F.BL2 > finfo.floor.BL2 + sig.syn_type_min_above_floor & F.GL2 > finfo.floor.GL2 + sig.syn_type_min_above_floor;
ex = F.BL2 - F.GL2;
st(synlike) = 9;
typed = synlike & usable;
st(typed) = 1 + sum(ex(typed) >= sig.syn_type_edges(:)', 2);
% type A edge by brightness b (the dimmer of BL2 and GL2 above its floor): the first
% edge where b >= syn_typeA_min_above_floor, further left by syn_typeA_dim_slope per
% log10 below that; events between that edge and the first edge stay B
shortfall = max(0, sig.syn_typeA_min_above_floor - min(F.BL2 - finfo.floor.BL2, F.GL2 - finfo.floor.GL2));
dimA = typed & st == 1 & ex >= sig.syn_type_edges(1) - sig.syn_typeA_dim_slope * shortfall;
st(dimA) = 2;
% per-file mixture (syn_type_mixture, Nov 2026) replaces the windows where there are
% enough typed class-2 Syn: types the file does not need are not used, and events
% whose type is not clear (posterior < syn_mix.pmin) get 5 (ambiguous). Fitted to
% class 2; Euk with Syn (5) and Syn on particle (8) take posteriors from that fit.
% info.syn_mix: method, kept types, positions, widths, elimination path, and soft
% counts (summed posteriors) per class 2/5/8 x type A-D
info.syn_mix = struct('method', 'windows', 'n_fit', sum(typed & class == 2), 'keep', true(1, 4), 'mu', nan(1, 4), 's', nan(1, 4), ...
    'path', {cell(0, 2)}, 'soft', nan(3, 4), 'n_ambiguous', 0);
if isfield(sig, 'syn_mix') && sig.syn_mix.on && sum(typed & class == 2) >= sig.syn_mix.min_n
    bmx = (F.BL2 - finfo.floor.BL2 + F.GL2 - finfo.floor.GL2) / 2;
    it = find(typed);
    [stm, Mx] = syn_type_mixture(ex(it), bmx(it), class(it) == 2, sig.syn_mix);
    st(it) = stm;
    soft = zeros(3, 4); scl = [2 5 8];
    for c = 1:3, soft(c, :) = sum(Mx.post(class(it) == scl(c), :), 1); end
    info.syn_mix = struct('method', 'mixture', 'n_fit', Mx.n_fit, 'keep', Mx.keep, 'mu', Mx.mu, 's', Mx.s, ...
        'path', {Mx.path}, 'soft', soft, 'n_ambiguous', sum(stm == 5));
end
info.syn_type = st;
% flag for unassigned events left out for a specific reason, so they can be
% counted later without entering the phytoplankton classes (class stays 0):
% 1 = Syn-pigment particle with low chl (Euk-with-Syn-like but far below the Euk
%     chl-size line: detritus or grazers with partly degraded Syn and/or Euks)
% 2 = faint particle below the chl gap at the bottom of the Euks
% 3 = trigger-noise pile at the chl edge (capped files, pre-2019)
% 4 = Pro-like event in a file run without Pro settings: no PE, chl at or above the
%     chl trigger edge, no larger than the file's Syn (peak size). Such files do not
%     get class 7, and their chl trigger cuts off dim Pro, so these events count
%     only as a lower bound for Pro (Nov 2026)
% 5 = Euk left out by the per-file Euk chl:size cut (low chl for its size)
fl = zeros(n, 1, 'uint8');
if isfield(info.trace, 'euksyn_lowchl'), fl(class == 0 & info.trace.euksyn_lowchl) = 1; end
fl(class == 0 & fl == 0 & euk_below_gap) = 2;
fl(class == 0 & fl == 0 & pile_cut) = 3;
if ~pro_measured
    if isfield(edge, 'BL3'), at_chl = F.BL3 >= edge.BL3; else, at_chl = F.BL3 > sig.chl_min; end
    fl(class == 0 & fl == 0 & ~has_pe & at_chl & F.size <= syn_size & F.SSC > finfo.floor.SSC + 0.1) = 4;
end
fl(class == 0 & fl == 0 & euk_chl_low) = 5;
fl(class == 0 & pro_eukrun_pile) = 3; %trigger pile found by the Pro-in-Euk-run test
% split-Syn parts are not classified (nothing should reach them; kept explicit). Flag 6: green
% part merged into the preceding event (not counted; the merged event is classified as usual);
% 7: unpaired green part with a Syn-type spectrum (GL2:GL3 above syn_min_GL2GL3, PE no brighter
% than the file's Syn by more than noblue_syn_GL2_above), a separate Syn count (Syn total = class
% 2 + flag 7); 8: unpaired blue part, not counted; 9: other unpaired green parts, not counted.
% dropout_frac_syn (file QC): (merged pairs + flag 7) / (Syn + flag 7); where it is high, Syn chl
% and FSC read low (merged blue parts ~23% and ~44% below whole Syn)
class(noblue | nogreen | mergedg) = 0;
if any(class == 2), gl2s = median(F.GL2(class == 2)); else, gl2s = Inf; end
syn_noblue = noblue & F.GL2GL3 > sig.syn_min_GL2GL3 & F.GL2 < gl2s + sig.noblue_syn_GL2_above;
fl(mergedg) = 6; fl(noblue) = 9; fl(syn_noblue) = 7; fl(nogreen) = 8;
info.n_noblue = sum(noblue);
info.n_syn_noblue = sum(syn_noblue);
info.n_syn_nogreen = sum(nogreen);
info.dropout_frac_syn = (info.n_merged + info.n_syn_noblue) / max(1, sum(class == 2) + info.n_syn_noblue);
info.flag = fl;
info.flag_names = {'Syn-pigment particle, low chl (degraded Syn/Euk, grazer, detritus)' 'below Euk chl gap' ...
    'trigger-noise pile at the chl edge (capped files)' 'Pro-like (no PE, no larger than Syn) in a file without Pro settings: Pro lower bound' ...
    'Euk chl:size cut (low chl for its size)' 'green part of a split Syn, merged into the preceding blue part: not counted' ...
    'unpaired green part of a split Syn (Syn-like PE under green, no SSC/chl pulse): separate Syn count' ...
    'unpaired blue part of a split Syn (no green-laser PE): not counted' 'unpaired green part, not Syn-like: not counted'};
info.syn_type_names = sig.syn_type_names; info.syn_type_edges = sig.syn_type_edges;
info.syn_typeA_min_above_floor = sig.syn_typeA_min_above_floor; info.syn_typeA_dim_slope = sig.syn_typeA_dim_slope;
info.features = F;
end

function s = sat_col(fcsdat, j, frac)
% events at or near the file's maximum in column j (none when the channel is not in the file)
if j == 0, s = false(size(fcsdat, 1), 1); return, end
x = fcsdat(:, j);
s = x >= frac * max(x);
end

function r = syn_ratio_range(mu, syn_GL2GL3, sig)
% [low high] GL2:GL3 of this file's Syn components: means of components named Syn
% with a Syn-type spectrum, mean GL3 clear of its floor, and at least
% sig.syn_comp_min_weight of the fitted events; the Syn median when there are none
use = mu.label == 2 & mu.gl3ok & mu.GL2GL3 > sig.syn_min_GL2GL3 & mu.weight >= sig.syn_comp_min_weight;
if any(use)
    r = [min(mu.GL2GL3(use)) max(mu.GL2GL3(use))];
elseif ~isnan(syn_GL2GL3)
    r = [syn_GL2GL3 syn_GL2GL3];
else
    r = [sig.syn_min_GL2GL3 sig.syn_min_GL2GL3];
end
end

function mu = name_components(g, finfo, sig, pe_min, pro_measured, sdB)
% sdB - chl (BL3) spread of each component for the Pro compactness test
% name mixture components from their means using pigment signatures
mu = array2table(g.mu, 'VariableNames', {'SSC' 'BL3' 'PEchl' 'GL2GL3'});
mu.GL2 = mu.PEchl + mu.BL3;
mu.chlSSC = mu.BL3 - mu.SSC;
mu.label = zeros(height(mu), 1);
mu.gl3ok = false(height(mu), 1);
for c = 1:height(mu)
    m = mu(c, :);
    comp_pe = m.PEchl > sig.pe_min_PEchl && m.GL2 > pe_min;
    % GL2:GL3 says nothing when GL3 is near its floor (dim Syn), so call those Syn;
    % PE euks among them are still picked out event by event below
    comp_gl3ok = m.GL2 - m.GL2GL3 > finfo.floor.GL3 + sig.ratio_min_above_floor;
    % component means average out GL3 noise, so a looser margin serves the Syn spectrum range
    mu.gl3ok(c) = m.GL2 - m.GL2GL3 > finfo.floor.GL3 + sig.comp_ratio_min_above_floor;
    if m.chlSSC < sig.detritus_max_chlSSC
        mu.label(c) = 0;                      %scatter without matching chlorophyll
    elseif comp_pe && (m.GL2GL3 > sig.syn_min_GL2GL3 || ~comp_gl3ok)
        mu.label(c) = 2;                      %Syn PE spectrum
    elseif comp_pe
        mu.label(c) = 3;                      %cryptophyte-type PE
    elseif pro_measured && m.BL3 < sig.pro_max_BL3 && sdB(c) <= sig.pro_max_sd_BL3
        % Pro candidate: dim, no PE, and a compact cluster in chl. Its size is
        % judged against this file's Syn after naming (assign_class_v2, Nov 2026);
        % candidates too large for Pro go back to Euk there
        mu.label(c) = 7;
    elseif pro_measured && m.SSC < sig.pro_max_SSC && m.BL3 < sig.pro_max_BL3
        % a diffuse component in the Pro size/chl window (small Euks and noise in
        % cold water) is left unassigned
        mu.label(c) = 0;
    else
        mu.label(c) = 1;
    end
end
mu.weight = g.ComponentProportion(:);
end

function [c, sd] = peak_centre(x)
% centre = peak of a kernel density of x; sd = half width at half maximum of the
% narrower side / 1.1774 (Gaussian equivalent). Robust to tails on either side.
x = x(isfinite(x));
if numel(x) < 20, c = median(x); sd = max(1.4826 * mad(x, 1), 0.05); return, end
[f, g] = ksdensity(x, 'NumPoints', 400);
[fm, i] = max(f); c = g(i);
lo = find(f(1:i) < fm / 2, 1, 'last'); hi = i - 1 + find(f(i:end) < fm / 2, 1, 'first');
w = [];
if ~isempty(lo), w(end + 1) = c - g(lo); end
if ~isempty(hi), w(end + 1) = g(hi) - c; end
if isempty(w), sd = max(1.4826 * mad(x, 1), 0.05); else, sd = max(min(w) / 1.1774, 0.05); end
end

function cut = spill_valley(r, sig)
% residual (PE:chl above the large Euks' trend) at the valley between the dominant
% density mode (chl spillover) and a strong mode above it (peak at least
% peeuk_spill_peak_rel of the dominant one, at least peeuk_spill_min_sep above it,
% density at the valley at most peeuk_spill_max_valley of the lower of the two
% peaks); Inf where there is no such mode
cut = Inf;
[f, g] = ksdensity(r, 'NumPoints', 400);
[pk, loc] = findpeaks([0 f 0]); loc = loc - 1;
[~, i0] = max(f);
up = loc(g(loc) >= g(i0) + sig.peeuk_spill_min_sep & pk >= sig.peeuk_spill_peak_rel * f(i0));
for k = up(:)'
    [fv, iv] = min(f(i0:k));
    if fv <= sig.peeuk_spill_max_valley * min(f(i0), f(k))
        cut = g(i0 + iv - 1); return
    end
end
end

function ok = pile_near_edge(f, g, e, max_above)
% the density peak of the population below a chl gap (density f on grid g) lies
% within max_above of the chl trigger edge e
[~, i] = max(f);
ok = g(i) <= e + max_above;
end

function [isp, st] = edge_pop_is_pile(F, m, e, ref, sig)
% is the population m (events near the chl trigger edge e) the trigger-noise pile?
% The pile (pre-2019 capped files: 90-95% within 0.1 of the edge) has chl pulses
% 0.24-0.33 (log10 A/H) shorter than the pulse reference, chl spread (robust SD)
% 0.02-0.03 and no size-chl relation; cell clusters near the edge (Syn, Pro, small
% Euks, the AR77 branch, flag 4/5 groups after 2019) have pulses like the reference,
% chl spread 0.07-0.4 and chl rising with size. Pile: median log A/H more than
% pile_AH_below_ref below the reference; without pulse width (no BL3-A column),
% chl SD below pile_max_sd_chl, more than pile_min_frac_edge within 0.1 of the edge,
% and |Spearman rho(size, chl)| below pile_max_rho
x = F.BL3(m); st = struct('n', sum(m), 'dAH', NaN, 'sd_chl', 1.4826 * mad(x, 1), 'f01', mean(x - e < 0.1), 'rho', NaN);
if isfinite(ref) && any(isfinite(F.BL3_AH(m))), st.dAH = median(F.BL3_AH(m), 'omitnan') - ref; end
if sum(m) >= 10, st.rho = corr(F.size(m), x, 'Type', 'Spearman'); end
if isfinite(st.dAH)
    isp = st.dAH < -sig.pile_AH_below_ref;
else
    isp = st.sd_chl < sig.pile_max_sd_chl && st.f01 > sig.pile_min_frac_edge && abs(st.rho) < sig.pile_max_rho;
end
end
