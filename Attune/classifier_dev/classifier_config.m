function cfg = classifier_config()
% function cfg = classifier_config()
% Settings shared by all cruises for the prototype classifier.
% Nothing here should need to change between cruises.
%
% Nov 2026, classifier redesign prototype

cfg.chans = {'FSC' 'SSC' 'GL1' 'BL1' 'BL2' 'BL3' 'GL2' 'GL3' 'GL4'};
% detector response signal ~ hv^k, from bead checks at two settings (2026_settings_evaluation,
% Sep 2026): SSC 220<->260 8.55, GL1 180<->220 8.55, BL2/BL3/GL4 8.4-8.5, FSC 250<->280 2.1
% (FSC barely responds to hv). The 2019/20 SSC/GL1 sweep agrees (8.3-8.7 between 180 and 260 V)
% but the response steepens with hv (local k ~6.5 at 150 V), so k holds only near the hv used at sea.
cfg.k = struct('FSC', 2.1, 'SSC', 8.5, 'GL1', 8.5, 'BL1', 8.5, 'BL2', 8.5, 'BL3', 8.5, 'GL2', 8.5, 'GL3', 8.5, 'GL4', 8.5);
% which FCB bead population (1 = 0.5, 2 = 1, 3 = 6 um) is on scale in each channel
% (SSC: 1 um; 0.5 and 1 um are equally repeatable, run-to-run SD 0.005 log10, and
% sit a constant 0.84 log10 apart across cruises; thresholds below moved from the 0.5 um scale)
cfg.bead_pop = struct('FSC', 2, 'SSC', 2, 'GL1', 2, 'BL1', 2, 'BL2', 1, 'BL3', 2, 'GL2', 1, 'GL3', 3, 'GL4', 1);
cfg.bead_window_days = 1.5;
% size merge (assign_class_v2): GL1 is side scatter through the OD2 filter on these
% cruises (p.OD2setting = 'GL1'); switch to it above these raw GL1-H counts,
% averaging with SSC in between (batch_mergedSSC2volume convention). GL1 is on the
% 1 um bead scale (bead_pop above) so that it is equivalent to SSC
cfg.gl1_is_od2_scatter = true;
cfg.gl1_merge = [500 1000];
cfg.floor_raw = 10; %raw values below this are treated as not detected

% pigment signatures in bead-relative log10 units (from AR92 summer and EN727 winter, Nov 2026)
cfg.sig.chl_min = -2;              %minimum BL3 to consider an event
cfg.sig.trigger_margin = 0;        %log10 above a fluorescence trigger level (0: count down to the edge, as the old method; was 0.05)
cfg.sig.pe_min_PEchl = -2.2;       %PE detected: GL2-BL3 above this ...
cfg.sig.pe_min_GL2 = -1.8;         %... and GL2 clear of its noise floor (~-2.4)
cfg.sig.syn_min_GL2GL3 = 0.15;     %Syn ~ +0.33, cryptophyte PE ~ -0.1
cfg.sig.detritus_max_chlSSC = 0.04;%Syn and Euks sit at +1.2 to +1.7
cfg.sig.peeuk_min_GL2 = -1;        %PE euks are bright in PE ...
cfg.sig.peeuk_size_above_syn = 0.6;%... larger than the file's Syn ...
cfg.sig.peeuk_PEchl_below_syn = 0.5;%... with less PE per chl than Syn ...
cfg.sig.peeuk_max_GL2GL3 = 0.1;    %... or a cryptophyte-type PE spectrum (winter)
% PE euk split: 4 (cryptophyte-type PE) if BL2:GL3 is this far below the file's Syn
% (old-class 4: -0.44 to -0.79; Syn-like low-PE euks ~0) ...
cfg.sig.peeuk4_dBL2GL3 = -0.3;
cfg.sig.peeuk4_dGL2GL3 = -0.15;     %... or GL2:GL3 this far below Syn where BL2/GL3 are unusable (old 4: -0.21 to -0.34)
cfg.sig.crypto_dB_sd = 2;           %crypto-type spectrum also at least this many robust SDs of the file's Syn spectrum below Syn (summer Syn spread is wide)
cfg.sig.peeuk_spec_PEchl_below_syn = 0; %PE euks entering on spectrum alone need PE:chl at least this far below the Syn median (unless GL2 saturates)
cfg.sig.crypto_syn_low_prc = 1;      %crypto-type spectrum (unsaturated PE euks): below this percentile of the file's Syn spectrum (every Syn pigment type) ...
cfg.sig.crypto_syn_low_margin = 0.1; %... minus this
cfg.sig.peeuk_spill_sd = 3;         %PE present on large particles: PE:chl this many robust SDs above the large Euks' PE:chl-vs-chl trend (chl spillover)
cfg.sig.peeuk_spill_min_n = 50;     %... when there are at least this many large Euks (else fixed pe_min_PEchl / peeuk_min_GL2)
cfg.sig.sat_frac = 0.98;            %raw value >= this fraction of the file's channel maximum counts as saturated
% crypto split into dim, red-shifted (3) and bright (4): accept a per-file 2-group fit on
% [PE in Syn equivalents, BL2:GL3 - Syn] only if ...
cfg.sig.crypto_split_min_n = 40;    %... there are this many crypto-type events with a usable spectrum,
cfg.sig.crypto_min_group_n = 15;    %... each group has this many,
cfg.sig.crypto_max_sd_PE = 0.4;     %... both are compact in PE (SD, log10) ...
cfg.sig.crypto_max_sd_dB = 0.2;     %... and in BL2:GL3 (EN727 / AR95 groups ~0.05-0.15),
cfg.sig.crypto_min_sep = 3;         %... means are this many pooled SDs apart (Mahalanobis),
cfg.sig.crypto_min_dB_offset = 0.1; %... and the dimmer group is the more red-shifted (EN727 -1.0 vs -0.6)
cfg.sig.crypto_min_PE_above_syn = 1; %without an accepted per-file split, crypto-type events with PE < this above Syn (log10, ~10 Syn) are dim (3), brighter are bright (4); EN727 dim groups ~0.25-0.75, bright 1.2-2.2
cfg.sig.crypto_outlier_sd = 4;     %single crypto group: events more than this many robust SDs below (above) its median PE are dim (bright)
cfg.sig.euksyn_max_chl_below_euk = 0.5; %Euk-with-Syn (5) more than this below the Euk chl-size line (log10, ~3x) are Syn on a particle (8)
cfg.sig.euksyn_low_chl_class = 0;   %class for Euk-with-Syn events below that: 0 = leave unassigned (they match neither Euk with Syn nor Syn on particle; the main Syn-on-particle rule can still take those that fit it), 8 = Syn on particle
cfg.sig.euksyn_min_PEchl_above_euk = 1; %Syn-check rejects re-checked as Euks carry Syn PE if PE:chl is this far above the file's Euks (log10)
cfg.sig.syn_PEchl_tol = 0.7;       %Syn events within this of the file's Syn PE:chl ...
cfg.sig.syn_size_tol = 1;          %... and no more than this larger than its Syn ...
cfg.sig.syn_size_tol_below = 0.5;  %... or this smaller (Syn size spread is ~0.2) ...
cfg.sig.syn_size_sd = 4;           %... widened to this many robust SDs of the file's Syn size,
cfg.sig.syn_size_max_below = 0.8;   %... but never more than this below the Syn median (tiny PE signals ~1.5 below are not Syn)
cfg.sig.syn_PE_tol_below = 0.5;     %Syn GL2 no lower than the Syn median minus the larger of this ...
cfg.sig.syn_PE_sd = 4;             %... and this many low-side SDs of Syn GL2
cfg.sig.syn_PE_gap_bw = 0.08;        %Syn PE floor raised to a PE gap (found before the Syn check): KDE bandwidth (log10)
cfg.sig.syn_PE_gap_frac = 0.1;       %... gap = local minimum below the Syn peak with density under this fraction of the peak
cfg.sig.syn_PE_gap_rise = 2;         %... and density rising again below it by at least this factor (a separate faint-PE pile)
cfg.sig.syn_PE_gap_min_pile = 10;    %... with at least this many events below the gap
cfg.sig.syn_PE_gap_min_n = 200;      %... only when there are at least this many Syn
cfg.sig.syn_PEchl_sd = 5;          %... and Syn may have up to this many robust SDs more PE per chl than the file's Syn (dim-chl Syn); detritus cut for Syn uses pro_chlSSC_sd SDs
cfg.sig.syn_GL2GL3_tol = 0.15;     %... within this of the file's Syn PE spectra (bright Syn spread ~0.04) ...
cfg.sig.ratio_min_above_floor = 2; %... checked only where GL3 is 100x above its floor (dim summer Syn at 10-30x are too noisy)
cfg.sig.syn_comp_min_weight = 0.01;%... measured from the range of the file's Syn components (pigment types mix), ignoring components under this weight
cfg.sig.comp_ratio_min_above_floor = 1; %... for that range, component mean GL3 need only be 10x above its floor (means are much less noisy than events)
% particle-associated Syn: one or a few Syn on, or inside, a larger particle
% (aggregate, detritus, or a grazer such as a nanoflagellate that has eaten Syn);
% Syn spectrum and about one Syn's pigment, scatter well above the file's Syn
cfg.sig.synpart_PEchl_tol = 0.3;   %PE:chl and ...
cfg.sig.synpart_GL2GL3_tol = 0.15; %... PE spectrum match the file's Syn
cfg.sig.synpart_BL3_below = 0.5;   %chl within 3x below (partly digested) ...
cfg.sig.synpart_BL3_above = 0.6;   %... to 4x above (several Syn) the file's Syn
cfg.sig.synpart_min_size_above = 0.4; %scatter more than ~2.5x the file's Syn ...
cfg.sig.synpart_size_sd = 3.5;     %... and beyond the file's Syn size spread (robust SDs), so the
cfg.sig.synpart_below_trend = 0.3;  %... and chl this far below the core Syn chl-vs-size trend extended to its size (large Syn follow the trend)
                                   %    large end of a broad summer Syn cluster stays Syn; no upper limit
cfg.sig.euk_chlSSC_tol = 0.7;      %Euk events no more than this below the file's Euk chl:size
cfg.sig.euk_min_size_below_syn = 0.5; %Euks no smaller than ~3x less scatter than the file's Syn
cfg.sig.euk_rescue_size_above_syn = 0.5; %events in detritus-named components return to Euk if in the Euk chl:size band, no PE, and this much larger than Syn
cfg.sig.nonphyto_min_n = 50;       %Euk chl:size cut uses the midpoint to the non-phyto particles when there are this many of them
cfg.sig.edge_pro_band = 0.4;        %non-Pro-setting files: no-PE events within this of the chl trigger edge (pile is ~0.4 thick in AR95 19Oct Euk(20); cuts no old-class Euks) ...
cfg.sig.edge_pro_size_above_syn = 0.1; %... and no larger than the file's Syn (+ this) are the Pro tail, left unassigned
cfg.sig.euk_gap_size_above = 1;     %chl gap at the bottom of the Euks: small no-PE particles (smaller than Syn + this) ...
cfg.sig.euk_gap_bw = 0.05;          %... KDE bandwidth of their chl (log10)
cfg.sig.euk_gap_min_above_edge = 0.2; %... Euk peak at least this far above the chl trigger edge
cfg.sig.euk_gap_frac = 0.1;         %... gap = local minimum below the peak with density under this fraction of the peak ...
cfg.sig.euk_gap_rise = 1.5;         %... and density rising again below it (pile at the trigger) by at least this factor ...
cfg.sig.euk_gap_min_pile = 20;      %... with at least this many events below the gap
cfg.sig.euk_gap_min_n = 200;        %... only with at least this many candidate particles
cfg.sig.syn_PEchl_default = -0.2;  %used when a file has no Syn
cfg.sig.syn_size_default = -1.24;
cfg.sig.pro_max_SSC = -1.44;        %smaller than Syn
cfg.sig.pro_max_BL3 = -0.3;        %dimmer than Syn and Euks
cfg.sig.pro_GL2_above_euk = 0.3;   %Pro have no PE: GL2 no more than this above the file's Euk median GL2 (noise level) ...
cfg.sig.pro_GL2_default_above_floor = 1.5; %... or above the GL2 floor when a file has too few Euks (AR92: Euk GL2 ~1.1-1.2 above floor)
cfg.sig.pro_GL2_sd = 4;            %... and more than this many robust SDs above the Pro cluster's median GL2
cfg.sig.pro_chlSSC_sd = 4;         %Pro fail the detritus chl:size cut only if also this many robust SDs below the Pro cluster's median chl:size
cfg.sig.pro_max_sd_BL3 = 0.4;       %Pro component must be compact in chl (SD; real Pro 0.09-0.29 over 14 Pro files, spurious AR95 cluster 0.50)

cfg.gmm_K = 3:9;
cfg.gmm_max_fit = 10000;
cfg.outlier_p = 0.999;
% catch-all: unassigned events with cell-like PE or chl take the majority class of
% this many nearest assigned neighbours (standardized size, chl, PE, BL2, GL3)
cfg.sig.catchall_k = 10;
cfg.sig.catchall_euk_chl_prc = 1;        %catch-all Euks need chl above this percentile of main-rule Euk residuals about their chl-size line
cfg.sig.cell_PE_below_syn = 0.5;   %PE like a cell: GL2 no more than this below the file's median Syn GL2 (~1/3 of a Syn)
cfg.sig.catchall_max_ref = 30000;   %subsample of assigned events used as neighbours


% Syn spectral subtypes (Nov 2026, syn_spectral_types.m): the bead-normalized
% excitation ratio BL2-GL2 (PE excited at 488 vs 532 nm, a PUB:PEB proxy) has three
% modes that recur across 20 cruises, 2020-2026: -0.23 (A, every file), 0.00 (B) and
% +0.45 (C, summer/fall only); troughs at -0.10 and +0.19. Labelled per event in
% info.syn_type for Syn (2), Euk with Syn (5) and Syn on particle (8):
% 1 = A, 2 = B, 3 = C, 9 = spectrum not usable (saturated, or BL2/GL2/GL3 near floor), 0 = other classes
cfg.sig.syn_type_edges = [-0.10 0.19];
cfg.sig.syn_type_names = {'A (low PUB:PEB)' 'B (mid)' 'C (high PUB:PEB)'};
cfg.sig.syn_type_min_above_floor = 1; %BL2, GL2 and GL3 at least this far above their floors
