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
%            since dim Pro sit at the edge in Pro runs too); n_lowchl, frac_lowchl
%            (unassigned Pro-like events with chl:size below lp_min_chlSSC, as a
%            fraction of them plus the Pro), uncertain (frac_lowchl above
%            cfg.sig.pro_lowchl_max_frac: Pro not separable from low-chl particles)
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
[cfg.sig.pe_min_PEchl, info_pe] = euk_autofluor_cut(fcsdat, fcshdr, ref, cfg, runtype);
[class, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, assign_pro);
info.euk_autofluor = info_pe;
info.pe_syn_trim = struct('applied', false, 'n', 0, 'gap', -Inf);
if strcmp(runtype, 'PE')
    [class, info.pe_syn_trim] = pe_syn_trim(class, info, cfg);
end
% PE-triggered runs: every recorded event carries some PE, so v2's catch-all counts the
% low-chl particle cloud (no chl trigger; BL3 ~ -1.3) as cell-like and gives it the Syn
% of its neighbours, without the detritus chl:size cut (TN368 15N24(1): 499 catch-all
% "Syn", chl:size median -0.23, PE:chl +0.49, none hand-gated). Catch-all Syn with chl:size
% below pe_catchall_syn_min_chlSSC are unassigned (the 0.6 detritus cut also took real dim
% Syn there: hand-gated Syn chl:size 5th percentile ~0.26)
info.n_pe_catchall_syn_removed = 0;
if strcmp(runtype, 'PE') && isfield(info, 'trace') && isfield(info.trace, 'catchall')
    rm = class == 2 & info.trace.catchall & info.features.chlSSC < cfg.sig.pe_catchall_syn_min_chlSSC;
    class(rm) = 0; info.n_pe_catchall_syn_removed = sum(rm);
end
info.runtype = runtype;
info.large_pro = struct('accepted', false, 'n', 0, 'method', '', 'mu', [], 'sep', NaN, 'ncand', 0, 'lda_n', 0, 'lda_median', NaN, 'lda_frac', NaN, 'n_extended', 0); %(same fields as large_pro returns)
info.pro_swap = struct('applied', false, 'n_to_pro', 0, 'n_to_0', 0);
if pro_run
    [class, info.pro_swap] = pro_component_swap(class, info, cfg);
end
if assign_pro && ~info.pro_swap.applied %(the swap has already found the file's Pro cluster)
    [class, info.large_pro] = large_pro(class, info, cfg, pro_run);
end
% v2 components named Euk that are large Pro by their mean: v2 can find a deep sample's
% large Pro as one clean component and name it Euk because it is larger than its Pro size
% limit (EN720 C07N02 CHL, 71 m: size -1.09, BL3 0.12, chl:size 1.21, GL2:GL3 -0.01;
% 589 of 683 hand Pro named Euk). large_pro missed it there (group fit sep 2.23 < 2.5; the
% fixed boundary's Pro side median 1.56 < 2). Such a component's Euk events become Pro
% if its mean is on the Pro side of lp_lda by lp_comp_min_score, small (size < Syn +
% lp_size_above_syn), dim (BL3 < lp_max_BL3) and with chl:size >= lp_comp_min_chlSSC (deep
% Pro comps 1.21-1.24; low-chl debris comps on the Pro side, 0.71-0.81, took 30-370 events
% in surface CHL runs without hand Pro, e.g. EN661 C07N05, EN668 C014N03) and a mixture
% weight >= lp_comp_min_weight (deep Pro comps 0.11, 0.46; small comps of half hand Euks,
% half unassigned, 0.024-0.051, in surface CHL runs: AR92 C01N02, MVCO MV418N01, EN668 C001N01)
info.n_pro_comp_euk = 0; info.pro_comp_sib = [];
if assign_pro && ~info.pro_swap.applied && ~info.large_pro.accepted && isfinite(cfg.sig.lp_comp_min_score) ...
        && isfield(info, 'gmm') && ~isempty(info.gmm) && isfield(info, 'comp') && istable(info.comp) && isfield(info.trace, 'has_chl')
    F = info.features; mu = info.comp;
    if sum(class == 2) >= 20, syn_size = median(F.size(class == 2)); else, syn_size = cfg.sig.syn_size_default; end
    scm = cfg.sig.lp_lda(1) * mu.chlSSC + cfg.sig.lp_lda(2) * mu.GL2GL3 + cfg.sig.lp_lda(3);
    pc = find(mu.label == 1 & scm > cfg.sig.lp_comp_min_score & mu.SSC < syn_size + cfg.sig.lp_size_above_syn ...
        & mu.BL3 < cfg.sig.lp_max_BL3 & mu.chlSSC >= cfg.sig.lp_comp_min_chlSSC & mu.weight >= cfg.sig.lp_comp_min_weight);
    % v2 can split a deep Pro cluster into two components that differ only in GL2:GL3
    % (HB1907 C04N02 CHL, 60 m: GL2:GL3 -0.10 and +0.18, size/BL3/PE:chl within 0.08; the
    % +0.18 one is on the Euk side of lp_lda, 530 hand Pro). Euk components whose mean is
    % within lp_comp_sib_tol of a converted one in size, BL3 and PE:chl, with chl:size >=
    % lp_comp_min_chlSSC, are converted with it
    if ~isempty(pc) && isfinite(cfg.sig.lp_comp_sib_tol)
        t = cfg.sig.lp_comp_sib_tol; sib = false(height(mu), 1);
        for k = pc(:)'
            sib = sib | (abs(mu.SSC - mu.SSC(k)) < t & abs(mu.BL3 - mu.BL3(k)) < t & abs(mu.PEchl - mu.PEchl(k)) < t);
        end
        sib = find(sib & mu.label == 1 & mu.chlSSC >= cfg.sig.lp_comp_min_chlSSC);
        info.pro_comp_sib = setdiff(sib(:)', pc(:)'); pc = union(pc(:)', sib(:)');
    end
    if ~isempty(pc)
        hc = info.trace.has_chl; ci = zeros(size(class));
        ci(hc) = cluster(info.gmm, [F.size(hc) F.BL3(hc) F.PEchl(hc) F.GL2GL3(hc)]);
        m = class == 1;
        % with them their small end, which v2 unassigned as Pro-like Euk (euk_prolike) before
        % this step: HB1907 C04N02 CHL (60 m), ~480 hand Pro smaller than the Syn peak left
        % at 0, or restored to Euk below (Nov 2026)
        if (~isfield(cfg.sig, 'lp_comp_take_prolike') || cfg.sig.lp_comp_take_prolike) && all(isfield(info.trace, {'euk_prolike' 'after_syncheck'}))
            m = m | (class == 0 & info.trace.euk_prolike & info.trace.after_syncheck == 1);
        end
        m = m & ismember(ci, pc);
        class(m) = 7; info.n_pro_comp_euk = sum(m); info.pro_comp_euk = pc(:)';
    end
end
% the same sibling rule after large_pro: large_pro decides event by event on lp_lda, so a
% component split off a deep Pro cluster by higher GL2:GL3 stays Euk (AR39B C12N07(1) CHL,
% 77 m: component GL2:GL3 +0.38 vs +0.04, size/BL3/PE:chl within 0.1 of the main Pro
% component, chl:size 1.20 vs run Euks 2.06; 314 of 5,195 hand Pro named Euk). Euk
% components within lp_comp_sib_tol of a component that is mostly Pro after large_pro,
% small, dim and with chl:size >= lp_comp_min_chlSSC, become Pro (Nov 2026)
info.n_pro_sib_after_lp = 0; info.pro_sib_after_lp = [];
if info.large_pro.accepted && (~isfield(cfg.sig, 'lp_sib_after_large') || cfg.sig.lp_sib_after_large) && isfinite(cfg.sig.lp_comp_sib_tol) ...
        && isfield(info, 'gmm') && ~isempty(info.gmm) && isfield(info, 'comp') && istable(info.comp) && isfield(info.trace, 'has_chl')
    F = info.features; mu = info.comp; t = cfg.sig.lp_comp_sib_tol;
    if sum(class == 2) >= 20, syn_size = median(F.size(class == 2)); else, syn_size = cfg.sig.syn_size_default; end
    hc = info.trace.has_chl; ci = zeros(size(class));
    ci(hc) = cluster(info.gmm, [F.size(hc) F.BL3(hc) F.PEchl(hc) F.GL2GL3(hc)]);
    np = accumarray(ci(hc), class(hc) == 7, [height(mu) 1]); nn = accumarray(ci(hc), 1, [height(mu) 1]);
    pk = find(np >= 100 & np > 0.5 * nn); sib = false(height(mu), 1);
    for k = pk(:)'
        sib = sib | (abs(mu.SSC - mu.SSC(k)) < t & abs(mu.BL3 - mu.BL3(k)) < t & abs(mu.PEchl - mu.PEchl(k)) < t);
    end
    sib(pk) = false;
    sib = find(sib & mu.label == 1 & mu.chlSSC >= cfg.sig.lp_comp_min_chlSSC & mu.SSC < syn_size + cfg.sig.lp_size_above_syn & mu.BL3 < cfg.sig.lp_max_BL3);
    if ~isempty(sib)
        m = class == 1 & ismember(ci, sib);
        class(m) = 7; info.n_pro_sib_after_lp = sum(m); info.pro_sib_after_lp = sib(:)';
    end
end
% small Euks v2 unassigned as Pro-like (euk_prolike: no PE, smaller than the Syn peak,
% chl below pro_max_BL3, which discrete_config raises to 0.3 for CHL-run Pro) go back
% to Euk where their spectrum is Euk-like: on the Euk side of the Pro-vs-Euk boundary
% (lp_lda) in [chl:size, GL2:GL3]. EN655 C10N16 CHL: 1,330 of 5,325 Syn-sized hand Euks
% (GL2:GL3 ~ +0.1; EN655 Pro ~ -0.2 to -0.6) were cut; user: keep that cluster as Euk.
% CHL runs only: there median agreement rose (unseen cruises 0.939 -> 0.947) with Pro
% unchanged; in Pro runs it only added Euk to runs that already overcount them
info.n_euk_prolike_restored = 0;
if strcmp(runtype, 'CHL') && isfield(info, 'trace') && isfield(info.trace, 'euk_prolike') && isfield(info.trace, 'after_syncheck')
    F = info.features;
    pl = class == 0 & info.trace.euk_prolike & info.trace.after_syncheck == 1 & isfinite(F.GL2GL3) & isfinite(F.chlSSC);
    sc = cfg.sig.lp_lda(1) * F.chlSSC + cfg.sig.lp_lda(2) * F.GL2GL3 + cfg.sig.lp_lda(3);
    back = pl & sc < -cfg.sig.euk_restore_margin;
    class(back) = 1; info.n_euk_prolike_restored = sum(back);
end
% Euks need chl for their size relative to the run's own Euks: faint particles near the
% chl trigger with low chl:size were Euk (EN661 C07N05 CHL, 180 m: 498 auto vs 185 hand
% Euk, chl:size ~0.95 vs Euks 1.77). Cut = median chl:size of the run's clear Euks (chl
% at least euk_chl.core_above_edge above the BL3 trigger edge, >= euk_chl.core_min_n of
% them) less euk_chl.below_core, no higher than euk_chl.cut_max. Euks below it AND with
% GL2:GL3 below euk_chl.max_GL2GL3 are unassigned: the particles have low GL2:GL3, while
% dim real Euks at the low-chl edge of a continuous cloud (EN655 C16N14, TN368 C80N23)
% keep their cloud's GL2:GL3
% dim cryptophytes named Syn: v2 names components by PE spectrum only (GL2:GL3 above
% syn_min_GL2GL3 = Syn), and weak crypto PE reads near the Euks' chl spillover spectrum,
% so a dim crypto cluster can be named Syn. The Syn check rejects it, and v2's re-check
% then splits it event by event at PE:chl = Euk + euksyn_min_PEchl_above_euk (EN727
% C20N02 CHL: one component, size 0.78, BL3 2.25, PE:chl -1.40, GL2:GL3 0.60; 177 to
% class 5, 683 to Euk; user: dim crypto). A Syn-named component that the Syn check mostly
% rejected becomes crypto (3) as a whole when it is much larger than Syn, its PE:chl is
% between the clear Euks' and Syn's, and its chl:size is Euk-like (Nov 2026). Where its
% PE spectrum is not crypto-type it is v2's class 6 population instead (see below)
info.n_crypto_from_syn = 0; info.crypto_from_syn = [];
dc = cfg.sig.dimcrypto;
if dc.on && isfield(info, 'gmm') && ~isempty(info.gmm) && isfield(info, 'comp') && istable(info.comp) ...
        && all(isfield(info.trace, {'has_chl' 'after_gmm' 'after_syncheck'})) && isfield(info.trigger_edges, 'BL3') && sum(class == 2) >= 20
    F = info.features; mu = info.comp;
    core = class == 1 & F.BL3 > info.trigger_edges.BL3 + cfg.sig.euk_chl.core_above_edge & isfinite(F.chlSSC);
    if sum(core) >= cfg.sig.euk_chl.core_min_n
        pc = median(F.PEchl(core)); cc = median(F.chlSSC(core));
        syn_size = median(F.size(class == 2)); syn_pe = median(F.PEchl(class == 2));
        hc = info.trace.has_chl; ci = zeros(size(class));
        ci(hc) = cluster(info.gmm, [F.size(hc) F.BL3(hc) F.PEchl(hc) F.GL2GL3(hc)]);
        rej = info.trace.after_gmm == 2 & info.trace.after_syncheck == 0;
        nr = accumarray(ci(hc), rej(hc), [height(mu) 1]); nn = accumarray(ci(hc), 1, [height(mu) 1]);
        % a tight cluster, not a broad catch-all component: spread of the rejected events in
        % size and PE:chl (EN727 C20N02 0.18 / 0.16; catch-alls that also passed the other
        % tests, AR61B pro_C32N10, MV337N04, TN368 C02N12, EN687 C05N08: 0.53-0.79)
        sdS = zeros(height(mu), 1); sdP = sdS; hr = hc & rej;
        if any(hr)
            sdS = accumarray(ci(hr), F.size(hr), [height(mu) 1], @std, Inf); sdP = accumarray(ci(hr), F.PEchl(hr), [height(mu) 1], @std, Inf);
        end
        k = find(mu.label == 2 & nr >= dc.min_n & nr >= dc.min_rejected_frac * nn & mu.SSC > syn_size + dc.size_above_syn ...
            & mu.PEchl > pc + dc.PEchl_above_euk & mu.PEchl < syn_pe - dc.PEchl_below_syn & mu.chlSSC > cc - dc.chlSSC_below_euk ...
            & sdS < dc.max_sd_size & sdP < dc.max_sd_PEchl);
        if ~isempty(k)
            m = ismember(ci, k) & rej & ismember(class, [1 5]);
            % class 6 (PE euk cluster) where the group's PE spectrum is not crypto-type:
            % BL2:GL3 at or above v2's crypto cut relative to Syn and GL2:GL3 well above
            % Syn. EN727 C20N02: the 860 events are the same population as v2's class 6
            % there (size 0.75 / 0.80, GL2:GL3 0.60 / 0.61, dB +0.05 / +0.03; the run's
            % cryptos dB -0.7 to -0.8), split from it only by v2's PE:chl re-check cut.
            % Class 6 is Euks with PE channels raised by preservation (counted as Euk)
            to = 3; pcl = cfg.sig.peeuk_cluster; syn = class == 2;
            if pcl.on && isfield(info, 'crypto_dB_cut') && isfinite(info.crypto_dB_cut)
                dBm = median(F.BL2(m) - F.GL3(m)) - median(F.BL2(syn) - F.GL3(syn));
                dGm = median(F.GL2GL3(m)) - median(F.GL2GL3(syn));
                if dBm >= info.crypto_dB_cut && dGm >= pcl.min_dGL2GL3_syn, to = pcl.class; end
            end
            class(m) = to; info.n_crypto_from_syn = sum(m); info.crypto_from_syn = k(:)'; info.crypto_from_syn_class = to;
        end
    end
end
% crypto in a crypto-named component dropped by v2's event-level PE test: PE must be
% peeuk_spill_sd robust SDs above the large Euks' PE:chl-vs-chl trend, and a second,
% mid-PE:chl Euk group widens that SD (EN649 C21N18 CHL, 3 m: spill SD 0.38; a clean
% crypto component, 61 events, 53 hand PE euk, GL2:GL3 -0.20, PE:chl -0.71, sat ~2.3 SD
% above the trend; 52 went to 0 and the catch-all put 46 back as Euk). Events v2 had
% as crypto before that test keep it when their spectrum is crypto-type, they are
% clearly larger than Syn and their PE:chl is well above the clear Euks' (Nov 2026)
info.n_crypto_kept = 0;
ck = cfg.sig.cryptokeep;
if ck.on && all(isfield(info.trace, {'after_pro_pe'})) && isfield(info.trigger_edges, 'BL3') && sum(class == 2) >= 20
    F = info.features;
    core = class == 1 & F.BL3 > info.trigger_edges.BL3 + cfg.sig.euk_chl.core_above_edge & isfinite(F.chlSSC);
    if sum(core) >= cfg.sig.euk_chl.core_min_n
        syn_size = median(F.size(class == 2));
        % not the Euk-with-Syn events v2 moved to Euk by its class-5 PE ceiling: v2 read
        % their spectrum as Syn-type (BL2:GL3), and their GL2:GL3 is compressed by GL2
        % saturation (MVCO MV418N01: 197 such events would become dim crypto)
        bright5 = false(size(class)); if isfield(info.trace, 'euksyn_bright'), bright5 = info.trace.euksyn_bright(:); end
        m = info.trace.after_pro_pe == 3 & ~ismember(class, [3 4 5 6]) & ~bright5 & F.GL2GL3 < cfg.sig.peeuk_max_GL2GL3 ...
            & F.size > syn_size + cfg.sig.peeuk_size_above_syn & F.PEchl > median(F.PEchl(core)) + ck.PEchl_above_euk ...
            & F.chlSSC > cfg.sig.detritus_max_chlSSC;
        class(m) = 3; info.n_crypto_kept = sum(m);
    end
end
info.euk_chl_cut = NaN; info.n_euk_lowchl = 0; info.n_euk_dim = 0; info.n_euk_highpe = 0; info.euk_ridge = NaN(1, 3);
if isfield(cfg.sig, 'euk_chl') && cfg.sig.euk_chl.on && isfield(info.trigger_edges, 'BL3')
    F = info.features;
    core = class == 1 & F.BL3 > info.trigger_edges.BL3 + cfg.sig.euk_chl.core_above_edge & isfinite(F.chlSSC);
    if sum(core) >= cfg.sig.euk_chl.core_min_n
        info.euk_chl_cut = min(cfg.sig.euk_chl.cut_max, median(F.chlSSC(core)) - cfg.sig.euk_chl.below_core);
        low = class == 1 & F.chlSSC < info.euk_chl_cut & F.GL2GL3 < cfg.sig.euk_chl.max_GL2GL3;
        % small dim particles (EN661 C07N05: size -2 to -1, BL3 near the trigger edge)
        % pass the chl:size cut; small Euks need relatively more chl, so also cut along a
        % shallower line, BL3 = a + ridge_slope * size through the clear Euks, at
        % ridge_k of their robust SDs (at least ridge_min_sd) below it
        b = cfg.sig.euk_chl.ridge_slope; a = median(F.BL3(core) - b * F.size(core));
        rsd = max(cfg.sig.euk_chl.ridge_min_sd, 1.4826 * mad(F.BL3(core) - (a + b * F.size(core)), 1));
        info.euk_ridge = [a b rsd];
        dim = class == 1 & F.BL3 < a + b * F.size - cfg.sig.euk_chl.ridge_k * rsd;
        info.n_euk_dim = sum(dim & ~low); low = low | dim;
        % near the chl trigger, events with PE:chl well above the clear Euks' (between
        % Euk and Syn; EN644 C02N02, 43 m: 178 small events hand left unassigned) are
        % not Euk
        pc = median(F.PEchl(core));
        pe = class == 1 & F.BL3 < info.trigger_edges.BL3 + cfg.sig.euk_chl.pe_edge_win & F.PEchl > pc + cfg.sig.euk_chl.pe_above_core;
        info.n_euk_highpe = sum(pe & ~low); low = low | pe;
        class(low) = 0; info.n_euk_lowchl = sum(low);
    end
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
% A separate cluster on the BL3 trigger floor called Pro beside the real Pro (EN688
% C11N23 pro_hi, 14 m: v2 Pro component at size -1.96, BL3 -1.64 = edge + 0.14, chl:size
% 0.32 vs Pro 0.71, A/H -0.23 vs 0.00; 7,180 events, none hand Pro; BL3 valley at -1.30
% with 38 events/bin against peaks of 1,494 and 764). Where auto Pro BL3 (0.05 bins,
% 3-bin mean) has a peak within pro_floor.max_above_edge of the edge, a second peak at
% least pro_floor.min_sep above it, both >= pro_floor.min_peak, and a valley below
% pro_floor.max_valley_ratio x the smaller peak, Pro below the valley are unassigned
% (fired on no other held-out pro run, 105 runs)
info.pro_floor = struct('valley', NaN, 'n_removed', 0);
if cfg.sig.pro_floor.on && isfield(info.trigger_edges, 'BL3') && sum(class == 7) >= 200
    v = pro_floor_valley(info.features.BL3(class == 7), info.trigger_edges.BL3, cfg.sig.pro_floor);
    if isfinite(v)
        f = class == 7 & info.features.BL3 < v;
        class(f) = 0; info.pro_floor = struct('valley', v, 'n_removed', sum(f));
    end
end
% Syn with little chl for their size and no Syn-like PE: faint particles on the chl
% trigger floor named Syn (AR77 pro_C13N16: 694 events at the BL3 edge, chl:size 0.17 vs
% Syn 1.51, GL2 1.1 below Syn). Pulse shape, spectrum and FSC/SSC do not separate them
% from hand-gated low-chl Syn (AUC 0.5); PE does: low-chl hand Syn keep their PE (PE:chl
% above the run's Syn), the particles have PE:chl at or below Syn or much less PE. Syn
% with chl:size more than floor_syn.dchl below the run's Syn median and (PE:chl below
% the Syn median or GL2 more than floor_syn.dGL2 below it) are unassigned (needs
% floor_syn.min_syn Syn), in floor_syn.runtypes (pro, CHL: in PE runs the net gain was ~0
% and AR43 PE C05N24 lost 822 hand Syn; held-out sets, 28 Nov 2026). EN644 pro_C23N13:
% 969 removed, none hand Syn (the dim end of the Syn diagonal below the hand gate).
% info.floor_syn.mask marks them whether or not the rule is applied
info.floor_syn = struct('on', false, 'n', 0, 'mask', false(size(class)));
if isfield(cfg.sig, 'floor_syn') && sum(class == 2) >= cfg.sig.floor_syn.min_syn
    F = info.features; fs = cfg.sig.floor_syn; sy = class == 2;
    m = sy & F.chlSSC < median(F.chlSSC(sy)) - fs.dchl & (F.PEchl < median(F.PEchl(sy)) | F.GL2 < median(F.GL2(sy)) - fs.dGL2);
    apply = fs.on && (~isfield(fs, 'runtypes') || ismember(runtype, fs.runtypes));
    info.floor_syn = struct('on', apply, 'n', sum(m), 'mask', m);
    if apply, class(m) = 0; end
end
info.pro = pro_completeness(info, class, cfg);
info.discrete_cfg = cfg.sig;
end

function v = pro_floor_valley(b, e, pf)
% BL3 valley between a trigger-floor peak and a higher Pro peak (NaN if none)
v = NaN; b = b(isfinite(b));
ed = e:0.05:max(b) + 0.05; if numel(ed) < 8, return, end
s = movmean(histcounts(b, ed), 3); x = ed(1:end - 1) + 0.025;
lo = find(x <= e + pf.max_above_edge); if isempty(lo), return, end
[pl, il] = max(s(lo)); il = lo(il);
up = find(x >= x(il) + pf.min_sep); if isempty(up), return, end
[ph, ih] = max(s(up)); ih = up(ih);
[pv, iv] = min(s(il:ih)); iv = iv + il - 1;
if pl >= pf.min_peak && ph >= pf.min_peak && pv <= pf.max_valley_ratio * min(pl, ph) && iv > il && iv < ih
    v = x(iv);
end
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
s = struct('accepted', false, 'n', 0, 'method', '', 'mu', [], 'sep', NaN, 'ncand', 0, 'lda_n', 0, 'lda_median', NaN, 'lda_frac', NaN, 'n_extended', 0);
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
    % (only events above that ceiling: all unassigned no-PE events added false Pro from
    % small particles in 7 hand-zero runs, e.g. EN706 C02N10 679)
    % (whole components, not single events: the HB1907 C15N03 Pro-run Pro component sits
    % at mean size -1.14, so a per-event size cut dropped half of it; v2 lists the
    % components it rejected as Pro in info.pro_comp_rejected; floor groups are small
    % and stay out by the size condition)
    if use_lda && isfield(g, 'pro_max_size_abs') && isfield(info, 'pro_comp_rejected') && ~isempty(info.pro_comp_rejected) && isfield(info, 'gmm')
        mu = info.comp; rj = info.pro_comp_rejected(:);
        rj = rj(mu.SSC(rj) > g.pro_max_size_abs - 0.1 & mu.PEchl(rj) < g.pe_min_PEchl & mu.chlSSC(rj) >= g.lp_min_chlSSC);
        if ~isempty(rj)
            hc = info.trace.has_chl; ci = zeros(size(class));
            ci(hc) = cluster(info.gmm, [F.size(hc) F.BL3(hc) F.PEchl(hc) F.GL2GL3(hc)]);
            euk = euk | (class == 0 & ismember(ci, rj) & F.PEchl < g.pe_min_PEchl);
        end
    end
end
if sum(euk) < g.lp_min_n && use_lda, [class, s] = large_pro_lda(class, F, syn_size, g, s, euk); return, end
% Fit on all the file's Euks (the larger Euks are the reference group: in deep
% samples the small dim "Euks" can be almost all Pro); K = 1-3 by BIC
cand = find(euk & isfinite(F.GL2GL3) & isfinite(F.chlSSC));
% events below the Pro chl:size floor (detritus-like) are left out of the fit: with them in,
% the lowest group was either a detritus group (rejected by lp_min_chlSSC, so the fixed
% boundary ran) or detritus merged with the Pro (accepted), depending on the random start
% (EN687 pro_C05N08: 498-625 vs 1280 Pro across seeds)
if ~isfield(g, 'lp_fit_above_floor') || g.lp_fit_above_floor, cand = cand(F.chlSSC(cand) >= g.lp_min_chlSSC); end
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
n_pro_before = sum(class == 7);
class(cand(side)) = 7;
% the fixed boundary is applied event by event and can cut through the Pro cluster: its
% brighter part (GL2:GL3 and chl:size slightly higher) lands on the Euk side and stays
% Euk (AR39B C13N04 Pro run: 104 of 848 hand Pro). Extend the Pro to the remaining
% candidates within the picked Pro's own spread in [size BL3 GL2GL3] (chi2 quantile
% lp_extend_p; chl:size = BL3 - size, so it adds nothing); the real Euks sit well apart (BL3 ~1.7, GL2:GL3 ~ +0.5).
% Only where the picked group is the run's main Pro (v2 named fewer Pro): a small
% Euk-side sliver next to v2's Pro (EN720 C14N14(1): 290 vs 4700) spreads over the
% Euks, and extending it took 1216 of 2337 hand Euks
s.n_extended = 0;
if isfinite(g.lp_extend_p) && s.lda_n > n_pro_before
    X = [F.size F.BL3 F.GL2GL3];
    P = X(cand(side), :);
    rest = cand(~side);
    if size(P, 1) >= 20 && ~isempty(rest)
        dx = X(rest, :) - median(P, 1);
        d2 = sum((dx / cov(P)) .* dx, 2);
        ext = rest(d2 < chi2inv(g.lp_extend_p, size(X, 2)));
        class(ext) = 7; s.n_extended = numel(ext);
    end
end
s.accepted = true; s.n = s.lda_n + s.n_extended; s.method = 'lda';
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
% ... except bottom slices v2 restored as part of an accepted Pro cluster (info.pro_comp_slice)
slice = false(height(mu), 1);
if isfield(info, 'pro_comp_slice') && ~isempty(info.pro_comp_slice), slice(info.pro_comp_slice) = true; end
if isfield(info, 'pro_comp_adjoin') && ~isempty(info.pro_comp_adjoin), slice(info.pro_comp_adjoin) = true; end %(and Syn-sized parts v2 joined to an accepted Pro cluster)
on_floor = on_floor & ~slice;
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
to_0 = class == 7 & ~to_pro & ~ismember(comp, find(slice)) & ~(d2 < chi2inv(cfg.outlier_p, size(X, 2)));
% ... and with chl:size near the new Pro's (the spread test has no chl:size: EN655
% C15N09 kept 337 low-chl events, chl:size median 0.32 vs Pro 1.10)
cs = F.chlSSC(to_pro);
if numel(cs) >= 20
    to_0 = to_0 | (class == 7 & ~to_pro & ~ismember(comp, find(slice)) & F.chlSSC < median(cs) - g.pro_chlSSC_sd * 1.4826 * mad(cs, 1));
end
class(to_0) = 0; class(to_pro) = 7;
s = struct('applied', true, 'n_to_pro', sum(to_pro), 'n_to_0', sum(to_0), 'pro_comp', find(cand)', 'floor_comp', find(on_floor)');
end

function [class, s] = pe_syn_trim(class, info, cfg)
% PE-triggered runs: a separate faint group at the PE (GL2) trigger, below the Syn
% cluster in chl, passes v2's broad Syn checks (EN617 PE runs: 2,200-4,600 per run,
% left out by the hand gates). Unlike the dim Syn tail of the tuning cruises, which runs
% continuously down the Syn diagonal to the trigger, it is set off from the Syn by a
% density valley in chl. If the chl (BL3) density of this run's Syn has a lower peak
% at least pesyn_gap_min_below below the main peak and the valley between them is
% under pesyn_gap_depth of that lower peak, Syn below the valley and within
% pesyn_near_edge of the PE trigger are unassigned. Runs without such a valley are
% left alone (tuning PE runs: no valley found; TN368: the faint group merges with
% the Syn, not caught)
g = cfg.sig; F = info.features;
s = struct('applied', false, 'n', 0, 'gap', -Inf);
if ~isfield(info.trigger_edges, 'GL2') || ~isfinite(g.pesyn_gap_depth), return, end
e = info.trigger_edges.GL2;
syn = class == 2 & isfinite(F.BL3);
x = F.BL3(syn);
if numel(x) < 200, return, end
[fk, gk] = ksdensity(x, 'NumPoints', 400, 'Bandwidth', g.pesyn_gap_bw);
[~, ip] = max(fk);
lm = find(fk(2:end-1) > fk(1:end-2) & fk(2:end-1) >= fk(3:end)) + 1;
lm = lm(gk(lm) < gk(ip) - g.pesyn_gap_min_below);
if isempty(lm), return, end
[pl, jl] = max(fk(lm)); jl = lm(jl);
[vm, jv] = min(fk(jl:ip)); jv = jl - 1 + jv;
if vm >= g.pesyn_gap_depth * pl || sum(x < gk(jv)) < g.pesyn_gap_min_n, return, end
low = syn & F.BL3 < gk(jv) & F.GL2 <= e + g.pesyn_near_edge;
class(low) = 0;
s = struct('applied', true, 'n', sum(low), 'gap', gk(jv));
end

function p = pro_completeness(info, class, cfg)
% where the Pro cluster sits relative to the chl (BL3) trigger edge
% lower_bound is false only where the cluster is shown complete; fewer than
% pro_edge_min_n Pro cannot show it (calibrate_pro_chl: CHL runs with 52-249
% Pro well above the edge had 1-2% of the Pro run's Pro)
p = struct('n', sum(class == 7), 'peak_above_edge', NaN, 'frac_near_edge', NaN, 'lower_bound', true, 'n_lowchl', 0, 'frac_lowchl', NaN, 'uncertain', false);
% Pro-like events left unassigned for low chl:size: no PE, no larger than Syn +
% lp_size_above_syn, in the Pro's own BL3 range, chl:size below lp_min_chlSSC. Many of
% them next to the Pro mean the sample's Pro are not separable from faded or detrital
% particles and the Pro count is uncertain (EN687 pro_C10N10: 578 such events beside 770
% Pro, frac 0.43; hand gate unusable)
F = info.features; g = cfg.sig;
if p.n >= g.lp_min_group_n
    if sum(class == 2) >= 20, syn_size = median(F.size(class == 2)); else, syn_size = g.syn_size_default; end
    b = prctile(F.BL3(class == 7), [5 95]);
    low = class == 0 & F.PEchl < g.pe_min_PEchl & F.size < syn_size + g.lp_size_above_syn ...
        & F.BL3 >= b(1) & F.BL3 <= b(2) & F.chlSSC < g.lp_min_chlSSC;
    p.n_lowchl = sum(low); p.frac_lowchl = p.n_lowchl / (p.n_lowchl + p.n);
    p.uncertain = p.frac_lowchl > g.pro_lowchl_max_frac;
end
if ~isfield(info.trigger_edges, 'BL3') || p.n < cfg.sig.pro_edge_min_n, return, end
x = info.features.BL3(class == 7) - info.trigger_edges.BL3;
g = linspace(0, 3, 301);
d = ksdensity(x, g, 'bandwidth', 0.05);
[~, k] = max(d);
p.peak_above_edge = g(k);
p.frac_near_edge = mean(x < cfg.sig.pro_edge_band);
p.lower_bound = p.peak_above_edge < cfg.sig.pro_lb_min_peak_above_edge | p.frac_near_edge > cfg.sig.pro_lb_max_frac_near_edge;
end

function [cut, s] = euk_autofluor_cut(fcsdat, fcshdr, ref, cfg, runtype)
% PE detected = PE:chl above the file's Euk autofluorescence: the Euks form the
% main low mode of PE:chl among cell-like events (chl:size above the
% detritus cut); the cut is that mode + cfg.sig.euk_af_above_mode. CHL runs: mode
% ~ -2.05, cut ~ -1.3; PE-triggered runs record only the brightest Euks (EN668
% C002N17: mode ~ -1.2, which the fixed cut called Syn). Without a clear Euk
% mode, cfg.sig.pe_min_PEchl.
s = struct('mode', NaN, 'n', 0, 'size_restricted', false); m = NaN; n = 0;
cut = cfg.sig.pe_min_PEchl;
F = attune_features(fcsdat, fcshdr, ref);
pechl = F.GL2 - F.BL3;
cell_like = F.BL3 - F.SSC > cfg.sig.detritus_max_chlSSC & isfinite(pechl) & F.BL3 > prctile(F.BL3, 10);
% Pro and PE runs: mode among events larger than euk_af_min_size, falling back to all
% cell-like events without a clear mode there. In Pro runs Pro can far outnumber the
% Euks, and preserved Pro carry more green autofluorescence (PE:chl ~ -0.84 in TN368
% C38N24 hi vs Euks -1.38), so their peak was taken for the Euk mode and the cut
% (-0.09) fell inside the Syn (Syn 972 of 2307 hand). In PE runs the Syn (SSC ~ -1.3)
% swamped the Euk peak (no mode; AR92 C16N08, C20N18, C21N19, EN687 C15N05 PE: 0 Euks
% of 815-997 hand). Not in CHL runs: small Euk-poor runs lost all Euks (EN668
% C015N02, EN715 C14N02) and AR92 C25N02 overcounted.
if ismember(runtype, {'pro' 'PE'})
    [m, n] = af_mode(pechl(cell_like & F.SSC > cfg.sig.euk_af_min_size), cfg);
    s.size_restricted = isfinite(m);
end
if ~isfinite(m)
    [m, n] = af_mode(pechl(cell_like), cfg);
end
s.n = n;
if ~isfinite(m), return, end
s.mode = m;
cut = s.mode + cfg.sig.euk_af_above_mode;
end

function [m, n] = af_mode(x, cfg)
% Euk PE:chl mode among cell-like events x (euk_autofluor_cut); NaN without a clear one
m = NaN; n = numel(x);
if n < cfg.sig.euk_af_min_n, return, end
g = linspace(-3.5, 1, 226);
d = ksdensity(x, g, 'bandwidth', 0.08);
pk = find(d(2:end-1) > d(1:end-2) & d(2:end-1) >= d(3:end)) + 1;
pk = pk(d(pk) >= cfg.sig.euk_af_min_peak * max(d) & g(pk) < cfg.sig.euk_af_max_mode);
if isempty(pk), return, end
% the Euk mode: the strongest qualifying peak within euk_af_above_mode of the lowest
% one (not simply the lowest: in AT46 C17N02 a minor bump at -2.36 put the cut at
% -1.61, inside the main Euk mode at -1.8, and all Euks were lost; the window keeps a
% PE-euk peak further up, e.g. -0.84 there, from being taken)
w = pk(g(pk) < g(pk(1)) + cfg.sig.euk_af_above_mode);
[~, b] = max(d(w));
m = g(w(b));
end
