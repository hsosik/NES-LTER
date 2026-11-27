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
% cruises (p.OD2setting = 'GL1'); switch to it above this GL1 level, averaging with
% SSC in between. GL1 is on the 1 um bead scale (bead_pop above) so that it is
% equivalent to SSC. The window is in 1 um bead units (GL1-H / bead GL1-H), not raw
% counts, so it stays where both channels are linear when voltages change: GL1 reads
% high below ~10^2 raw (OD2 low-signal excess), SSC compresses above ~10^5.2 (220 V)
% to ~10^5.6 (240 V). 0.5-2x bead is linear on both at SSC/GL1 220/220 (all cruises
% to 2026; raw GL1 ~560-2200, close to the earlier 500-1000 counts) and 240/180
% (AR104). Nov 2026, classifier_dev/pt_transfer/README.md
% Callers set this per cruise: strcmpi(attune_config(cruise).OD2, 'GL1') (false before
% Aug 2019, when OD2 sat on SSC or nowhere); assign_class_v2 also skips the merge when GL1
% is not in the file or ref.OD2 (cfg.cruise set) is not 'GL1'. Optional cfg.cruise
% lets bead_reference record the OD2 placement in ref.
cfg.gl1_is_od2_scatter = true;
cfg.gl1_merge_bead = [0.5 2];
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
cfg.sig.peeuk_spill_method = 'valley'; %'mad': floor = trend + peeuk_spill_sd robust SDs over all large Euks; 'valley': no higher than the valley below a separate upper (PE) mode of the residuals, where there is one (PE euks / particles in the Euk set widen the SD)
cfg.sig.peeuk_spill_peak_rel = 0.3; %'valley': upper mode peak at least this fraction of the dominant (spillover) peak ...
cfg.sig.peeuk_spill_min_sep = 0.6;  %... at least this far above it (log10 PE:chl; two Euk modes without PE: 0.4-0.6) ...
cfg.sig.peeuk_spill_max_valley = 0.6; %... with the valley density at most this fraction of the lower peak (PE modes 0.36-0.49; Euk modes 0.68-0.87)
% class 6, PE euk cluster: Euk-with-Syn (5) events forming one compact population with
% their own PE. PE:chl and GL2:GL3 robust SDs of the class-5 events (usable spectra) at
% most max_sd_*; among those within win_sd of the centre (at least min_n) GL2 rising
% with BL3 (robust slope within slope, correlation at least min_r) and median GL2:GL3
% at least min_dGL2GL3_syn above the file's Syn: they take class (6)
cfg.sig.peeuk_cluster = struct('on', true, 'class', 6, 'min_n', 50, 'max_sd_PEchl', 0.18, 'max_sd_GL2GL3', 0.06, ...
    'win_sd', 3, 'slope', [0.7 1.3], 'min_r', 0.7, 'min_dGL2GL3_syn', 0.15);
cfg.sig.peeuk_sat_spectrum = true; %true: saturated events with a usable BL2:GL3 go by spectrum: GL2-saturated ones are crypto only with a crypto-type spectrum, and BL3-saturated ones (PE:chl overstated) with a Syn/spillover-like spectrum are not PE euks
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
cfg.sig.euksyn_min_syn_PE_prc = 1;   %Euk with Syn (5) needs GL2 at least this percentile of the file's Syn GL2 (one Syn's PE); below -> Euk (1). 0 = off
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
cfg.sig.edge_pro_band = 0.4;        %non-Pro-setting files: band above the chl trigger edge where no-PE events pile up (~0.4 thick in AR95 19Oct Euk(20)); defines the trigger-noise pile, and Syn rescue seeds lie above it
% (the Pro-tail cut, edge_pro_size_above_syn, was removed Nov 2026: it sliced real small-Euk clusters in summer files; see assign_class_v2 'Euks carry no PE')
cfg.sig.edge_pile_max_frac = 0.5;   %non-Pro files: where no-PE events in that band are more than this fraction of chl events (pre-2019 trigger-noise pile) ...
cfg.sig.edge_pile_fit_ratio = 1;    %... they enter the mixture fit sample at most this many times the other events
cfg.sig.euk_gap_size_above = 1;     %chl gap at the bottom of the Euks: small no-PE particles (smaller than Syn + this) ...
cfg.sig.euk_gap_bw = 0.05;          %... KDE bandwidth of their chl (log10)
cfg.sig.euk_gap_min_above_edge = 0.2; %... Euk peak at least this far above the chl trigger edge
cfg.sig.euk_gap_frac = 0.1;         %... gap = local minimum below the peak with density under this fraction of the peak ...
cfg.sig.euk_gap_rise = 1.5;         %... and density rising again below it (pile at the trigger) by at least this factor ...
cfg.sig.euk_gap_min_pile = 20;      %... with at least this many events below the gap
cfg.sig.euk_gap_min_n = 200;        %... only with at least this many candidate particles
cfg.sig.euk_chl_in_v2 = true;       %per-file Euk chl:size cut inside assign_class_v2 (from assign_class_discrete; live files). discrete_config keeps it off: preserved runs get it after their Pro steps
cfg.sig.euk_chl.on = true;          %... settings as discrete_config's: clear Euks this far above the BL3 trigger edge ...
cfg.sig.euk_chl.core_above_edge = 0.8;
cfg.sig.euk_chl.core_min_n = 30;    %... at least this many ...
cfg.sig.euk_chl.below_core = 0.6;   %... Euks need chl:size above their median less this ...
cfg.sig.euk_chl.cut_max = 1.2;      %... with the cut no higher than this ...
cfg.sig.euk_chl.max_GL2GL3 = 0;     %... removing only Euks also below this GL2:GL3
cfg.sig.euk_chl.ridge_slope = 0.5;  %also unassign Euks below the line BL3 = a + slope*size through the clear Euks ...
cfg.sig.euk_chl.ridge_k = 3.5;      %... by this many of their robust SDs ...
cfg.sig.euk_chl.ridge_min_sd = 0.1; %... the SD at least this
cfg.sig.euk_chl.pe_edge_win = 1.0;  %also unassign Euks within this of the BL3 trigger edge ...
cfg.sig.euk_chl.pe_above_core = 0.8;%... with PE:chl this far above the clear Euks' median
cfg.sig.euk_chl.pile_gate = true;   %... only where the low chl:size no-PE particles pile up at the chl trigger (assign_class_v2 only; discrete_config's cut has no such field) ...
cfg.sig.euk_chl.pile_win = 0.3;     %... their chl density peaking within this of the BL3 edge ...
cfg.sig.euk_chl.pile_min_n = 30;    %... with at least this many of them ...
cfg.sig.euk_chl.pile_bw = 0.08;     %... (KDE bandwidth, log10) ...
cfg.sig.euk_chl.pile_valley_rise = 1.5; %... and only below the first valley above the pile (density rising again by this factor above it)
cfg.sig.euk_gap_peak_local = true;  %... Euk peak = highest local maximum above the edge window when the window's top density is only the trigger pile's slope (Nov 2026, AR39B)
cfg.sig.euk_gap_pile_max_above_edge = 0.7; %... and only where the population below the gap peaks within this of the chl edge (trigger noise; underway piles 0.37-0.60, a real small-Euk cluster in MVCO ~1.0)
cfg.sig.pro_in_euk_runs = true;        %runs without Pro settings: unassigned no-PE events at the chl edge that are cells (not the trigger pile) become Pro, a lower bound (assign_class_v2) ...
cfg.sig.pro_eukrun_size_above = 0.25;  %... no larger than Syn + this ...
cfg.sig.pro_eukrun_size_below = 1.1;   %... no smaller than Syn - this ...
cfg.sig.pro_eukrun_max_chlSSC = 1.4;   %... chl:size below this ...
cfg.sig.pro_eukrun_min_n = 50;         %... at least this many
cfg.sig.pile_AH_below_ref = 0.15;      %trigger-noise pile: median chl pulse log10(A/H) more than this below the pulse reference (pile 0.24-0.33 below; cells within 0.1) ...
cfg.sig.pile_AH_below_ref_event = 0.25;%... single events this far below stay unassigned when the rest are cells
cfg.sig.pile_max_sd_chl = 0.06;        %... without pulse width: chl robust SD below this (pile 0.02-0.03, cells 0.07-0.4) ...
cfg.sig.pile_min_frac_edge = 0.7;      %... more than this fraction within 0.1 of the chl edge ...
cfg.sig.pile_max_rho = 0.2;            %... and no size-chl relation (|Spearman rho| below this; cells 0.2-0.9)
cfg.sig.dim_euk_class = true;          %class 10, small Euk-like cells with low chl for their size (AR77 branch): no PE ...
cfg.sig.dim_euk_size_above = [0.25 1.0]; %... size from Syn + (1) to Syn + (2) ...
cfg.sig.dim_euk_max_chlSSC = 0.8;      %... chl:size below this (AR77 branch ~0.6; summer small Euks peak at 1.0-1.4) ...
cfg.sig.dim_euk_min_n = 200;           %... at least this many, and cells (not the trigger pile)
cfg.sig.syn_ref_max_size = -0.5;    %file's Syn reference from Syn-named events no larger than this (real Syn peak -1.5 to -0.75; broad large PE components named Syn otherwise move it) ...
cfg.sig.syn_ref_min_n = 20;         %... when there are at least this many (else all Syn-named events)
cfg.sig.syn_ref_above_edge = true;  %... from those clearly above the chl trigger edge (edge_pro_band) when there are that many (PE+chl trigger files)
cfg.sig.catchall_peeuk_det = true;  %catch-all PE euks need chl:size above detritus_max_chlSSC, as main-rule PE euks
cfg.sig.noblue_on = true;           %split Syn (laser dropout): pairs merged (green part flag 6); unpaired green parts flag 7 (Syn-like, separate Syn count) / 9; unpaired blue parts flag 8 ...
cfg.sig.noblue_min_GL2BL2 = 1.5;    %... raw log10 GL2-A minus BL2-A above this (Syn 0.4-0.9; EN661/EN649 dropouts 1.5-3.5; other held-out files: 0 at the 95th percentile, at most 282)
cfg.sig.noblue_max_SSCW = 30;       %... and SSC-W below this (raw; Syn ~200) ...
cfg.sig.noblue_max_BL3W = 50;       %... and BL3-W below this (Syn ~270; in the EN661 episode dropout widths spread 0-50)
cfg.sig.noblue_syn_GL2_above = 1;  %unpaired Syn-like green parts (flag 7, separate Syn count): Syn-type GL2:GL3 and GL2 no more than this above the file's Syn median; other dropouts flag 7
cfg.sig.nogreen_min_BL2 = 3.5;      %blue-laser part of a split Syn (unpaired: flag 8, not counted): raw log10 BL2-A above this ...
cfg.sig.nogreen_min_BL2BL3 = -0.8;  %... BL2-A minus BL3-A above this (Syn ~+0.7; Euks ~-2) ...
cfg.sig.nogreen_max_GL2W = 50;      %... GL2-W below this (Syn ~260), and GL2-A below BL2-A by more than noblue_min_GL2BL2
cfg.sig.dropout_min_split = 20;     %split Syn are merged (and leftover parts flagged) only where at least this many pairs are found: a blue part ...
cfg.sig.dropout_pair_dt = 1;        %... followed within this many Time units (1 ms; both parts carry the same Time stamp) ...
cfg.sig.merge_min_GL2BL2 = 1.0;     %... by an event with raw GL2-A above BL2-A by more than this ...
cfg.sig.split_min_GL2 = 3;          %... and raw log10 GL2-A above this (also for unpaired green parts)
cfg.sig.syn_ref_spec_min_n = 20;    %Syn reference PE:chl from those with a measured Syn-type spectrum (GL3 clear of floor, GL2:GL3 > syn_min_GL2GL3) when at least this many
cfg.sig.syn_PEchl_default = -0.2;  %used when a file has no Syn
cfg.sig.syn_size_default = -1.24;
cfg.sig.syn_rescue = true;         %no component named Syn: events near these defaults with a Syn PE spectrum start as Syn (rare Syn; assign_class_v2) ...
cfg.sig.syn_rescue_min_n = 20;     %... if there are at least this many
cfg.sig.pro_max_SSC = -1.44;        %smaller than Syn (fixed limit, used only where a file has no Syn component; Pro components are otherwise judged against Syn, below)
cfg.sig.pro_min_size_below_syn = 0.1; %Pro component size at least this far below the Syn components' (real Pro 0.12-1.0 below in 108 Pro-setting files; fall Pro 0.15-0.5 below)
cfg.sig.pro_max_size_abs = -1.15;   %... and Pro component size no more than this (largest real Pro cluster in the 108-file survey: -1.25)
cfg.sig.pro_max_size_below_syn = 1.1; %... and no more than this far below (tiny high-chl artifacts 1.2-1.5 below: AR104 23Oct wPro(9), AR79)
cfg.sig.pro_max_AH_above_ref = 0.2; %Pro-setting files: chl pulse shape log10(BL3 A/H) no more than this above the file's reference (median of small PE-bearing events, mostly Syn: their 99th percentile is 0.06-0.15 above); wider-pulse faint signals (0.3-0.4) are not Pro and stay out of the mixture fit
cfg.sig.pro_AH_sd = 4;             %... and Pro A/H no more than this many SDs (narrow-side spread) above the Pro median (summer Pro: broad A/H, spread ~0.1)
cfg.sig.pulse_ref_size_above = 0.5; %... reference events no larger than syn_size_default + this
cfg.sig.pro_max_chlSSC_above_syn = 0; %Pro component chl:scatter below the Syn components' (real Pro 0.1-1.0 below; the tiny artifacts at or above Syn)
cfg.sig.pro_max_BL3 = -0.3;        %dimmer than Syn and Euks
cfg.sig.pro_GL2_above_euk = 0.3;   %Pro have no PE: GL2 no more than this above the file's Euk median GL2 (noise level) ...
cfg.sig.pro_GL2_default_above_floor = 1.5; %... or above the GL2 floor when a file has too few Euks (AR92: Euk GL2 ~1.1-1.2 above floor)
cfg.sig.pro_GL2_sd = 4;            %... and more than this many robust SDs above the Pro cluster's median GL2
cfg.sig.pro_chlSSC_sd = 4;         %Pro fail the detritus chl:size cut only if also this many robust SDs below the Pro cluster's median chl:size
cfg.sig.pro_max_sd_BL3 = 0.4;       %Pro component must be compact in chl (SD; real Pro 0.09-0.29 over 14 Pro files, spurious AR95 cluster 0.50)
cfg.sig.pro_max_chlSSC_below_syn = 1; %Pro component chl:scatter no more than this below the Syn components' (real Pro 0.1-1.0 below over 108 Pro files; AR95 19Oct AltPro(30) false group 1.1-1.2)
cfg.sig.pro_chlSSC_edge_start = 0.2; %... that limit holds for components within this of the chl trigger edge (BL3 mean); further above it grows ...
cfg.sig.pro_chlSSC_edge_slope = 0.5; %... by this per log10 unit of height above (edge + start) ...
cfg.sig.pro_max_chlSSC_below_syn_far = 1.3; %... up to this (AR79 10Nov Pro(14): real Pro 1.03 below Syn, 0.55 above the edge; AR95 19Oct AltPro(30) false group 1.13 below, 0.20 above)
cfg.sig.pro_max_sd_size = 0.45;     %Pro component size spread (robust SD of its narrow-pulse events) no more than this (real Pro 0.10-0.38; EN727 C19N13 false column ~0.5)
cfg.sig.pro_slice_size_tol = 0.25;  %a Pro candidate failing only the chl:scatter limit is still Pro if an accepted Pro component is within this in size and brighter in chl (the cluster's bottom slice at the chl trigger); also a Syn-sized candidate (within pro_min_size_below_syn of Syn, no larger) next to an accepted Pro component

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
% excitation ratio BL2-GL2 (PE excited at 488 vs 532 nm, a PUB:PEB proxy) has
% modes that recur across 20 cruises, 2020-2026: -0.23 (every file), 0.00 and
% +0.45 (summer/fall only); troughs at -0.10 and +0.19. Labelled per event in
% info.syn_type for Syn (2), Euk with Syn (5) and Syn on particle (8).
% Nov 2026 (output/syn_type_d): a fourth mode at -0.48 below the old A, with PE:chl
% ~0.3 above it and the same emission ratios, on the fresher northern shelf
% (AE2426: 4.5% of Syn in 41% of files; AR77; also EN657, EN712, EN715, EN727, AR95).
% Types are now lettered in order of BL2-GL2 (old A/B/C became B/C/D):
% 1 = A (-0.48), 2 = B (-0.23), 3 = C (0.00), 4 = D (+0.45),
% 9 = spectrum not usable (BL2 or GL2 saturated or near floor), 0 = other classes.
% Type-A edge by brightness b = min(BL2, GL2) above floor. In dim Syn the ratio is
% noisy: summer Syn (mostly 1.5-2.25 above floor) spread type B below -0.40 (39-64%
% of events at b 1-1.25, ~0 above 2.5), 2-7% of all Syn, while real type A keeps a
% flat share when bright. So the first edge holds for b >= syn_typeA_min_above_floor
% and moves left by syn_typeA_dim_slope per log10 below it (-0.50 at b 1.75, -0.60
% at 1.5). Files grouped by bright type A (10 cruises, ~2,000 files): type A kept
% 0.90, false A 0.64% of Syn; a fixed limit (dim events never A) kept 0.77 at 0.28%
% (summer type A is dim: EN706, EN720, AR92 inshore). A PE:chl check was tried first
% and dropped: type A PE:chl is +0.3 above B in AE2426/AR77 but only 0-0.07 above it
% in summer type A (AR92, EN687) and in preserved samples
cfg.sig.syn_type_edges = [-0.40 -0.10 0.19];
cfg.sig.syn_type_names = {'A (lowest PUB:PEB)' 'B (low PUB:PEB)' 'C (mid)' 'D (high PUB:PEB)'};
cfg.sig.syn_typeA_min_above_floor = 2; %type A: first edge where min(BL2, GL2) is at least this far above its floor (log10) ...
cfg.sig.syn_typeA_dim_slope = 0.4;     %... edge this much further left per log10 below that
cfg.sig.syn_type_min_above_floor = 1; %BL2 and GL2 at least this far above their floors
% Per-file mixture typing (Nov 2026; syn_type_mixture, prototype output/syn_type_d/
% mix_proto4.m). The windows above count the noisy tails of dim Syn as neighbouring
% types (summer files: 5-7% "type A", all B spread) and call any shoulder a type. The
% mixture fits types A-D along BL2-GL2 per file, with the noise width set by
% brightness (mean of BL2 and GL2 above floor; band widths fitted jointly over 1,858
% files), and keeps only the types the file needs. Summer type A share 0.1-0.6% (was
% 5-7%), AE2426 2.7% (unchanged), type A kept in AE2426 194/301 files and AR92 inshore.
% Events whose most probable type has posterior below pmin get code 5 (ambiguous):
% 3-6% of Syn in bright-Syn cruises, 22-34% in dim summer cruises. Files with fewer
% than min_n typed Syn (class 2) use the windows and sloped edge above
cfg.sig.syn_mix.on = true;
cfg.sig.syn_mix.min_n = 500;          %typed class-2 Syn needed for the mixture (else windows)
cfg.sig.syn_mix.mu0 = [-0.48 -0.23 0.00 0.45]; %type positions A-D ...
cfg.sig.syn_mix.win = [0.12 0.10 0.10 0.12];   %... each may move this far per file
cfg.sig.syn_mix.sgrid = 0.03:0.01:0.08;        %intrinsic sd of each type per file (wider ones swallowed neighbours)
cfg.sig.syn_mix.band_edges = [1 1.25 1.5 1.75 2 2.25 2.5 2.75 3]; %brightness bands (mean of BL2, GL2 above floor) ...
cfg.sig.syn_mix.band_sd = [0.377 0.21 0.21 0.175 0.13 0.10 0.085 0.07 0.065]; %... total sd of a type in each (last: >= 3), joint fit
cfg.sig.syn_mix.s0ref = 0.06;         %typical intrinsic sd inside band_sd (noise = band_sd^2 - s0ref^2)
cfg.sig.syn_mix.alpha = 200;          %band weights pulled toward the file's weights by this many pseudo-events
cfg.sig.syn_mix.pen = 30;             %drop a type while that costs less than pen*log(n) in log-likelihood
cfg.sig.syn_mix.wmin = 0.01;          %... and always below this share ...
cfg.sig.syn_mix.shoulder_max = 0.10;  %... or below this share with its position pinned at the limit next to a kept type
cfg.sig.syn_mix.bg_mu = 0.1; cfg.sig.syn_mix.bg_sd = 0.5; cfg.sig.syn_mix.bg_max = 0.05; %broad outlier component, at most bg_max per band
cfg.sig.syn_mix.pmin = 0.8;           %type assigned where its posterior is at least this; else code 5 (ambiguous)
cfg.sig.syn_mix.xbins = -1:0.02:1.2;  %ratio bins of the fit
