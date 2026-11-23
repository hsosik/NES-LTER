function cfg = discrete_config()
% function cfg = discrete_config()
% classifier_config with the changes needed for preserved (discrete) samples.
% Preserved cells on the bead-normalized scale, relative to live cells from
% the same water (preservation_shift, 6 cruises, Nov 2026):
%   Euk: green/yellow autofluorescence; GL2 +1.3 to +1.65, BL2 +0.5 to +0.8,
%        BL1 +0.15 to +0.8, chl ~0, size -0.02 to -0.1, so PE:chl ~ -2.05
%        (5-95%: -2.56 to -1.42) instead of ~ -3.5
%   Syn: PE:chl +0.15 to +0.27 (median ~ +0.10), size -0.04 to -0.07, spectral
%        ratios within 0.09
%   Pro: size -0.06 to -0.21, chl -0.03 to -0.07, GL2 +0.2 to +0.4
% Hand-gated events, CHL and pro runs from 12 cruises (examine_discrete_events):
%   PE euk PE:chl 5-95% -1.15 to -0.20, GL2 0.42 to 1.54; Euk GL2 95% 0.29
%
% Nov 2026, discrete-sample classifier development

cfg = classifier_config();
cfg.preserved = true;
% PE detected: above the preserved Euks' autofluorescence (Euk 95th percentile -1.42,
% PE euk 5th -1.15); live -2.2
cfg.sig.pe_min_PEchl = -1.3;      %fallback; assign_class_discrete sets it per file from the Euk mode:
cfg.sig.euk_af_above_mode = 0.75;  %... cut = Euk PE:chl mode + this (Euk 5-95% spread -0.5 to +0.6 about the mode)
cfg.sig.euk_af_min_peak = 0.1;     %... mode = lowest KDE peak of cell-like PE:chl at least this fraction of the highest ...
cfg.sig.euk_af_max_mode = -0.8;    %... and below this (Syn ~ +0.1, PE euks -1.15 to -0.2 sit above Euks)
cfg.sig.euk_af_min_n = 200;        %... with at least this many cell-like events
cfg.sig.euk_af_min_size = -1.0;    %... pro and PE runs: among events larger than this (SSC; Pro -1.3 to -1.8, Syn ~ -1.3), all cell-like events if no mode there
% PE euks bright in PE: above preserved Euk GL2 (95th percentile 0.29); live -1
cfg.sig.peeuk_min_GL2 = 0.35;
% Autofluorescent non-cell particles: preserved samples have a cloud with PE:chl
% ~ -1 and little chl for their size (chl:size ~ 0.3), which v2 named Syn
% (EN715 C10N02, EN661 C14N07); with that "Syn" as the file's reference the
% real Euks then failed the Euk-vs-Syn size rule. Hand-gated Syn and Euks have
% chl:size above ~1.0 (5th percentile); Pro has its own SD-based exception. Live 0.04
cfg.sig.detritus_max_chlSSC = 0.6;
cfg.sig.detritus_max_chlSSC_pro = 0.04;
% Pro brightness limit: hand-gated preserved Pro BL3 5-95% -1.42 to +0.26 (live rule
% "dimmer than -0.3"); in CHL runs the chl trigger edge is ~ -0.45, so -0.3 left
% almost no room for a visible Pro cluster. Syn are kept out by their PE (named first)
cfg.sig.pro_max_BL3 = 0.3;
% large, bright Pro (assign_class_discrete large_pro): per-file 2-group split of the
% small no-PE Euks; hand Pro vs Euk at size -1.3 to -0.6: chl:size 1.09 vs 1.70,
% GL2:GL3 -0.10 vs 0.34
cfg.sig.lp_size_above_syn = 0.6;  %candidates smaller than Syn + this ...
cfg.sig.lp_max_BL3 = 1.0;         %... and dimmer than this (hand Pro BL3 99th percentile 0.74)
cfg.sig.lp_min_n = 100;           %... at least this many
cfg.sig.lp_min_d_chlSSC = 0.3;    %Pro group lower in chl:size by at least this ...
cfg.sig.lp_min_d_GL2GL3 = 0.2;    %... and in GL2:GL3 by at least this ...
cfg.sig.lp_min_sep = 2.5;         %... groups this many pooled SDs apart (Mahalanobis) ...
cfg.sig.lp_min_group_n = 50;      %... Pro group at least this many events
cfg.sig.lp_min_chlSSC = 0.6;      %... and its mean chl:size at least this (large-Pro groups 0.78-1.18; a detritus-like group at 0.3 in EN655 C15N09)
% Pro runs where the group fit fails (Pro outnumber Euks): fixed discriminant,
% Pro side if score = a*chlSSC + b*GL2GL3 + c > 0 (hand-gated Pro vs Euk, dcfg6 events)
cfg.sig.lp_lda = [-7.080 -4.316 9.189];
cfg.sig.lp_lda_min_median = 2;    %... applied if the Pro-side events' median score is at least this
cfg.sig.lp_lda_dominant = [0.75 1.2]; %... or if at least this fraction of candidates is on the Pro side and their median at least this
cfg.sig.lp_extend_p = 0.99;        %large_pro_lda: Pro extended to candidates within this chi2 quantile of the picked Pro's spread (Inf: off)
cfg.sig.euk_restore_margin = 0.5;  %v2 Pro-like small Euks back to Euk where the Pro-vs-Euk score (lp_lda) is below -this (Euk side, with a margin)
% Pro runs (assign_class_discrete pro_component_swap): Pro as large as Syn named Euk,
% dimmer groups (trigger floor, faint smear) named Pro (EN655 C15N09)
cfg.sig.pc_min_below_cand = 0.3;        %swap only if every component named Pro is at least this much dimmer (mean BL3) than the candidate
cfg.sig.pc_min_above_edge = 0.8;        %replacement component's mean BL3 at least this above the edge ...
cfg.sig.pc_max_size_above_syn = 0.15;   %... no larger than Syn + this ...
cfg.sig.pc_chlSSC_below_syn = [0.2 1.0]; %... chl:size this far below Syn's (min, max)
% PE-triggered runs (assign_class_discrete pe_syn_trim): a separate faint group below the Syn
% in chl at the PE trigger (EN617), set off by a density valley; the continuous dim Syn tail of
% the tuning cruises has no valley and is kept. (A fixed chl-below-trend cut was tried first
% and cut that real tail: EN687 C10N18 1.02 -> 0.66 of hand Syn.)
cfg.sig.pesyn_near_edge = 0.6;     %only Syn within this of the GL2 trigger edge are unassigned ...
cfg.sig.pesyn_gap_bw = 0.05;       %... KDE bandwidth of the Syn chl (BL3) density ...
cfg.sig.pesyn_gap_min_below = 0.3; %... lower peak at least this below the main Syn peak ...
cfg.sig.pesyn_gap_depth = 0.7;     %... valley under this fraction of the lower peak (Inf: off) ...
cfg.sig.pesyn_gap_min_n = 50;      %... and at least this many Syn below the valley
cfg.sig.pe_catchall_syn_min_chlSSC = 0; %PE runs: v2 catch-all Syn with chl:size below this are unassigned (low-chl particle cloud near chl:size 0). A fixed cut beat cuts relative to the run's Syn core (core chl:size varies 0.73-2.27; the cloud does not): pe_catchall_cut_test, 94 PE runs; 0.1 removed ~45% more hand Syn for the same gain
% Pro in CHL runs (assign_class_discrete): assigned, with a lower-bound flag where
% the chl trigger cuts the Pro cluster (thresholds calibrated in calibrate_pro_chl)
cfg.sig.pro_in_chl_runs = true;
cfg.sig.pro_edge_band = 0.1;               %Pro this close to the chl trigger edge (log10) count as at the edge
cfg.sig.pro_edge_min_n = 400;  %(was 100: CHL runs with 127-364 Pro clear of the edge were wrong or surface: AR39B C6N14 227 vs Pro run 9602, EN688 C10N18 364 at 0.41x, AR92 C26N05 2.8 m 169, EN727 C10N08 127)              %fewer Pro than this: no completeness estimate, lower_bound true
cfg.sig.pro_lb_min_peak_above_edge = 0.3;  %lower bound if the Pro chl mode is closer than this to the edge ...
cfg.sig.pro_lb_max_frac_near_edge = 0.1;   %... or more than this fraction of Pro is at the edge %Pro runs (assign_class_discrete): Pro chl:size 0.27-1.32, so keep the live cut
% SYBR bacteria runs (assign_class_hbac); hand-gated SYBR runs, 6 cruises:
% bacteria chlx (BL3-BL1) 0.35 +- 0.14, pex (GL2-BL1) -0.52 +- 0.18, BL1-SSC >= ~0,
% SSC 5-95% -2.6 to -1.5; Syn PE:chl 5th percentile -0.77, Euk -1.95 to -1.2
cfg.hbac.stain_min = 0;          %stained: BL1 - SSC above this (unstained particles are larger, BL1 near floor)
cfg.hbac.core_max_SSC = -1.4;    %bacteria core for the spillover lines: smaller than this
cfg.hbac.bl1_edge_margin = 0.05; %above the BL1 trigger edge by this
cfg.hbac.min_core = 500;         %fewer core events: default offsets below
cfg.hbac.default_chlx = 0.35; cfg.hbac.default_pex = -0.52; cfg.hbac.default_sd = 0.15;
cfg.hbac.min_sd = 0.08;
cfg.hbac.phyto_sd = 4;           %chl or PE excess: this many robust SDs above the bacteria line ...
cfg.hbac.phyto_min = 0.4;        %... and at least this (log10)
cfg.hbac.syn_min_PEchl = -0.9;   %Syn-like PE:chl (Syn 5th percentile -0.77; Euks -1.95 to -1.2)
cfg.hbac.pro_max_SSC = -1.45;    %chl excess without Syn PE: Pro at or below this SSC, Euk above
cfg.hbac.syn_max_SSC = -0.8;     %Syn by chl excess alone: SSC below this (Syn 95th percentile -1.0)
cfg.hbac.bact_max_SSC = -1.2;    %bacteria no larger than this (SSC 95th percentile -1.5)
cfg.hbac.min_chlSSC = 0.6;       %chl excess needs chl:size above this (as detritus_max_chlSSC in phyto runs) ...
cfg.hbac.min_chlSSC_pro = 0.04;  %... or this for Pro-sized events
% HNA / LNA subtype of bacteria (9), info.bact_type (hbac_hna_split): cut on
% y = BL1 - slope * (SSC + 2) at the run's density valley, else a fixed cut.
% 200 + 175 held-out runs, 25 cruises (hbac_hna_test, Nov 2026): valley median
% -1.26 (5-95% -1.35 to -1.18); no clear valley mostly in summer surface samples
cfg.hbac.hna.slope = 0.35;        %HNA cloud rises with SSC (best slope per run 0.3-0.4)
cfg.hbac.hna.win = [-1.45 -1.10]; %valley search window in y (keeps off the trigger-edge band near BL1 -1.85)
cfg.hbac.hna.flank = 0.35;        %dip = valley density / lower of the maxima within this either side
cfg.hbac.hna.bw = 0.04;           %ksdensity bandwidth
cfg.hbac.hna.min_n = 2000;        %fewer bacteria: fixed cut
cfg.hbac.hna.max_dip = 0.9;       %valley used only if the dip is below this
cfg.hbac.hna.default_cut = -1.26; %fixed cut (median valley)
% Syn defaults (rare-Syn rescue, files without Syn): hand-gated preserved Syn
cfg.sig.syn_PEchl_default = 0.1;   %live -0.2
cfg.sig.syn_size_default = -1.29;  %live -1.24
end
