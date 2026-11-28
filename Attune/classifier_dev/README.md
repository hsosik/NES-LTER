# Classifier redesign: exploration

Goal: replace the per-cruise `assign_class_<cruise>.m` functions with one classifier that works across cruises, seasons and instrument settings.

## Scripts

| script | what it does |
|---|---|
| `collect_bead_runs.m` | all QC'd FCB bead runs from every cruise (voltage + log10 peak position per channel) → `output/bead_runs.mat` |
| `extract_cluster_centers.m` | for a time-spread sample of files per cruise: median position of existing Syn/Euk/Pro classes, detector voltages, nearest bead run, lat/lon/temperature → `output/cluster_centers.mat` |
| `fit_gain_curves.m` | fits signal ∝ voltage^k per detector from bead runs |
| `analyze_normalization.m` | compares raw, voltage-corrected and bead-relative cluster positions; writes `output/normalization_spread.csv` and figures |

Run from this folder, e.g. `extract_cluster_centers(16, {}, 'output/cluster_centers.mat')` then `analyze_normalization`.

## Test 1 (Nov 2026): can normalization put all cruises on one scale?

Data: 648 files from 39 cruises (2019–2026), 332 bead runs. Labels are the existing class files.

1. **Detector response is a power law, signal ∝ hv^8.5.** Beads give k = 8.5 (BL2, GL2). Syn in Pro/standard alternating samples gives k = 8.4 ± 0.6 (BL3), measured on the same water at different voltages. A 10% voltage change shifts a cluster ~2.2×.
   - Bead checks at two settings (`2026_settings_evaluation`, Sep 2026) confirm k ≈ 8.5 for SSC, GL1, BL2, BL3 and GL4. FSC is the exception: k ≈ 2.1.
   - An SSC/GL1 bead sweep from 150 to 270 V (2019/20) shows the response is not a pure power law. The local k rises from about 6.5 at 150 V to 8.5–9 above 200 V, so a single k holds only near the voltages used at sea.
   - An earlier SSC k of 7.1 came from Syn in only 4 off-voltage files and has been replaced with 8.5.
2. **Voltage correction removes the settings effect.** The Pro-setting vs standard offset in Syn BL3 drops from 0.56 to 0.09 log10; for Euk it drops from 0.71 to 0.06.
3. **Scattering becomes cruise-independent once bead-relative.** Between-cruise SD of Syn position:
   - SSC: 0.39 → 0.07 log10
   - GL1 (SSC through the OD2 filter): 0.74 → 0.06 log10
4. **Fluorescence brightness does not become cruise-independent.** Between-cruise SD stays at 0.4–0.5 log10, a ~1.5 decade range. That is biology, not instrument: Syn GL2 (bead-normalized) vs temperature gives r = −0.84. Cold-water Syn are brighter and larger.
5. **Spectral ratios are nearly invariant.** Syn falls on one line in chl–PE space across all cruises and seasons (`output/syn_euk_chl_pe.png`). Between-cruise SD:
   - PE/chl (GL2−BL3): 0.12
   - GL2/GL3: 0.03
   - pigment per scatter (BL3−SSC, GL2−SSC): 0.36–0.46, varying with temperature

### Design implications

- **Population identity** comes from spectral shape (log ratios GL2/BL3, GL2/GL3, BL1, …). These ratios are universal across cruises and seasons, so they can be fixed in a shared config.
- **Size** comes from bead-relative SSC/GL1. It is instrument-invariant, so remaining variation is real cell size.
- **Brightness and size position** change with season and water mass. Clusters must be located per file (mixture model or density peaks seeded from neighbouring files), not with fixed thresholds.
- Correct to a reference voltage first using the header `hv`, so Pro/standard alternating settings and future setting changes need no special cases.

## Prototype classifier (Nov 2026)

`assign_class_v2(fcsdat, fcshdr, ref, cfg, pro_measured)` with `ref = bead_reference(beadstat, time, cfg)` and `cfg = classifier_config()`. There are no cruise- or date-specific settings. Steps:

1. **Common scale:** all channels go onto the bead-relative scale (`attune_features`), using voltages and trigger levels from the FCS header.
2. **Pigment check:** keep events whose chlorophyll is above the BL3 trigger edge, or whose PE is clearly detected (dim summer Syn trigger on GL2 only).
3. **Size:** SSC, switching to GL1 (scatter through the OD2 filter) above 1000 raw GL1-H counts and averaging the two between 500 and 1000, as in `batch_mergedSSC2volume`. SSC and GL1 are each relative to their own 1 µm bead peak and taken as equivalent, with no fitted offset.
4. **Per-file clusters:** a Gaussian mixture in [size, BL3, PE:chl, GL2:GL3], with the number of components chosen by BIC. Components are named by pigment signature rules in `cfg.sig`.
5. **Event-level checks, relative to each file's own clusters:**
   - too little chlorophyll for the particle size → 0
   - Syn must sit near the file's Syn PE:chl and size
   - Euks must sit near the file's Euk chlorophyll per size
   - PE-containing euks are picked cell by cell: PE present, larger than Syn, and less PE per chl than Syn or a cryptophyte-type spectrum

`compare_v2` and `plot_compare_counts` compare against the existing class files. `plot_v2_vs_old` draws side-by-side cytograms.

### Agreement with existing classes, 40 files each

| | Syn | Euk (no PE) | PE euk | Pro |
|---|---|---|---|---|
| AR92 (Aug, warm): files within 10% | 39/40 | 39/40 | 4/35 | 10/11 |
| AR92 median v2/old | 0.97 | 0.99 | 0.66 | 0.97 |
| EN727 (Jan, cold): files within 10% | 40/40 | 40/40 | 8/40 | – |
| EN727 median v2/old | 1.00 | 1.02 | 0.86 | – |

Syn here excludes the new particle-associated Syn class (code 8, see the winter Syn smear section): 2% of v2 Syn-like events in AR92, 8% in EN727. Before that class and the tighter Syn check, EN727 Syn was 26/40 within 10% (median 1.07).

Euk minimum size rule: no known Euks are much smaller than Syn (*Ostreococcus* only slightly), so Euks must be within 0.5 log10 of the file's Syn size (`euk_min_size_below_syn`).
- This removed the tiny winter particles at the trigger edge (scatter ~10× below Syn).
- Pro is not an alternative for them in cold water.
- Winter Euks went from 16/40 to 40/40 files within 10%.

Open questions:
- ~~**Winter Syn excess (+7%)**~~: resolved by the particle-associated Syn class and the tighter Syn check (see the winter Syn smear section).
- **Unassigned class:** class 0 mixes electronic noise, non-phytoplankton particles (detritus, bacteria) and below-detection events. Splitting it would make the exclusions inspectable.
- **PE euks:** they are rare, and the old and new definitions differ near the Syn line.
- **Pro near the chlorophyll trigger edge:** ~3% of old Pro events are unassigned.

### Noise vs small Euks using all parameters (EN727, 25 files)

`examine_disputed`, `discriminate_noise` and `check_disputed_ratios` analyse events in the small/dim region where the classifiers overlap.

- **Agreed Euks and agreed noise separate almost perfectly using parameters neither classifier uses** (bagged trees, cross-validated by file):
  - AUC 0.998, from BL1, BL2, GL4, GL1, pulse area/height ratios and FSC-W.
  - The most informative are GL4, BL2, GL1, BL3 A/H, SSC A/H and FSC-W.
  - Noise has broadband fluorescence (BL2−BL3 −1.3 vs −2.6 for Euks) and flatter pulses (BL3 A/H 0.009 vs 0.087).
  - **These are good candidates for an explicit noise-rejection step.**
- **Disputed events** (old 0, v2 Euk; ~12% of v2 Euks in winter):
  - They are tiny: scatter ~10× below Syn.
  - They sit on the BL3 trigger edge, with PE and BL1 at the detection floor.
  - Their FSC pulses are the shortest of any group.
  - Their concentration tracks Syn and Euks (Spearman r = 0.88 and 0.83) and moves opposite to agreed noise (r = −0.50). Trigger settings were constant for the whole cruise.
  - A model trained on agreed classes scores 81–96% of them as noise, but they lie outside both training distributions, so that is extrapolation.
  - The GL4/BL3 chlorophyll ratio does not discriminate near the detection floor.
  - Conclusion: probably the bright tail of a population mostly below detection. They can't be identified reliably, and their counts depend on distance to the trigger.

### Winter Syn smear (`examine_syn_smear`, EN727, 25 files)

Events v2 calls Syn but the old classifier left as 0 (8% of v2 Syn, 1–24% per file). Offsets are relative to each file's agreed Syn (`output/syn_smear_EN727.png`). The smear has two parts:

1. **Large-scatter tail (~70% of the smear): Syn pigment on a bigger particle.**
   - Spectrum identical to core Syn: PE:chl +0.03, GL2/GL3 −0.01, GL4/BL3 +0.08, BL2/BL3 0.00, BL1/BL3 +0.04.
   - Pulse shapes (A/H) are identical.
   - Pigment amount only +0.1–0.2 log10 (1.3–1.5×), while SSC is +0.5–1.1 (3–12×). FSC and GL1 rise too, so these are really larger particles, not an SSC artefact.
   - Not doublets, which would show +0.3 in both pigment and scatter. FSC-W is only slightly longer (+0.1–0.2).
   - Most consistent with one Syn cell on, or alongside, a non-fluorescent particle (aggregate or detritus).
   - **Not a separate cluster:** Syn-spectrum events form a continuous tail out to +2.2 in scatter (~1000 per 0.25 log bin, slowly declining). The v2 cut at `syn_size_tol` = 1 is arbitrary; about 3000 more such events lie beyond it.
   - The smear count barely tracks Syn concentration (Spearman r = 0.02) and only weakly tracks non-phyto particles of the same size (r = 0.25).
2. **Small dim events (~30% of the smear): not Syn.**
   - Scatter 10–100× below Syn, chl −1.1, and a different PE spectrum (GL2/GL3 −0.23, PE:chl +0.34, BL1/BL3 +0.74).
   - They enter through the PE-detected path and pass the Syn checks because those have no lower size bound or GL2/GL3 tolerance.

**Implemented:**
- **Class 8, particle-associated Syn.** Picked event by event after clustering:
  - PE:chl within 0.3 and GL2/GL3 within 0.15 of the file's Syn.
  - Chl between 3× below and 4× above the file's Syn, allowing for partly digested Syn and several Syn per particle.
  - Scatter above the file's Syn by more than 0.4 log10 and more than 3.5 robust SDs of the file's Syn size, with no upper limit and no chl:size detritus rule.
  - The particle may be an aggregate, detritus, or a grazer such as a nanoflagellate that has eaten Syn. These can't be told apart with these measurements.
  - Each event holds roughly one Syn, so it could count toward Syn abundance, but its scatter doesn't describe Syn size or carbon. Not yet used by `generate_attune_table`.
- **Tighter Syn check:** no more than 0.5 log10 smaller than the file's Syn (`syn_size_tol_below`), and GL2/GL3 within 0.15 of the file's Syn (`syn_GL2GL3_tol`).
  - The GL2/GL3 check applies only where GL3 is more than 2 decades above its floor (`ratio_min_above_floor`).
  - Dim summer Syn sit only 1–1.5 decades above the floor, where the ratio is noisy. Applying the check there removed 22% of AR92 Syn.
- **Result (40 files each):**
  - EN727 Syn: 40/40 files within 10% of old (median 1.00), up from 26/40.
  - Class 8 holds 9,400 events (8% of Syn-like events); 87% of them were unassigned by the old classifier.
  - AR92: 1.9% of old Syn moved to class 8 and 1.6% became unassigned; Syn median v2/old is 0.97.

**Summer check (`examine_synpart`, AR92 and EN727, 25 files each, `output/synpart_*.png`).** Is class 8 a separate group or just the large end of the Syn cluster? Two tests:
1. Within Syn, chl rises with scatter (slope 0.84 in AR92, 0.62 in EN727). If class 8 were big Syn, its chl would continue that trend.
2. The Gaussian shape of the file's Syn scatter predicts how many Syn should fall beyond the class-8 cut.

With the original fixed cut (+0.4 above the file's Syn):
- **AR92:** 82% of class 8 sat at +0.3–0.7, with chl near the Syn trend, FSC barely raised, and counts matching the Gaussian Syn tail (observed/expected 1.1). These were big Syn, not particle-associated. Summer Syn scatter is broader (robust SD 0.20 vs 0.14 in winter).
- **EN727:** 67× more events than the Syn tail predicts, with flat chl (+0.15) and FSC rising with SSC from the start. This is a genuinely separate group.

**Fix:** the class-8 cut is now the file's Syn size + max(0.4, 3.5 robust SD of Syn size) (`synpart_size_sd`). After the fix:

| | AR92 | EN727 |
|---|---|---|
| class 8 / Syn (median per file) | 0.4% | 7% |
| class 8 vs Gaussian Syn tail | 17× | 330× |
| class 8 chl offset (vs Syn trend) | +0.21 (trend +0.86) | +0.18 (trend +0.61) |
| Syn files within 10% of old (median v2/old) | 39/40 (0.99) | 37/40 (1.00) |

- In both seasons class 8 now has about one Syn's chl, while FSC and GL1 rise with SSC (real larger particles).
- In AR92 a subset has FSC-W 30–100× longer than Syn, consistent with elongated particles such as aggregates or chains.
- Highest class-8 fraction: one low-Syn EN727 file (345 class 8 vs 744 Syn).

### Pro robustness (`check_pro`, 45 Pro-setting files from AR92, AR95, HRS2609)

- **Warm water (>22 °C):** v2 is within ~5% of the old counts.
  - Pro is separable from broadband-fluorescing particles: those form their own BL2/GL4 diagonal.
  - **But 20–43% of Pro events sit within 0.2 log10 of the BL3 trigger edge** (median only ~0.25 above it), so the population is truncated by the detection limit.
  - At the edge every other channel is at its floor, so dim Pro cannot be told from edge pile-up by any parameter.
- **AR95 (October, 15–17 °C):** in 7 files v2 counts 400–2400 Pro where the old classifier found none.
  - These "Pro" are larger and brighter than summer Pro, overlapping Syn in size and reaching Syn-level chlorophyll.
  - A subset has an unusual BL3 pulse shape (A/H ~0.35–0.4). Probably not Pro.
  - The Pro rule uses fixed thresholds, not limits relative to the file's Syn.
- **HRS2609 "less sensitive" setting (21 °C):** 361 v2 Pro where the old classifier found 0, with 76% within 0.2 of the edge, so likely edge pile-up.

Candidate fixes:
1. Pro must be smaller and dimmer in chlorophyll than the file's own Syn (event level, like the Euk rule).
2. Report a per-file Pro detection-quality metric (fraction of Pro near the trigger edge), and flag files where Pro isn't resolved from the edge instead of counting the pile-up.
3. Optionally, a QC flag for Pro reported in cold water using underway temperature.


### SSC bead reference: 1 µm instead of 0.5 µm (Nov 2026)

`cfg.bead_pop.SSC` is now the 1 µm FCB bead.
- In all QC'd bead runs (213 with SSC medians), the 0.5 and 1 µm peaks are equally repeatable in SSC: within-cruise run-to-run SD 0.005 log10 for both.
- The two peaks sit a constant 0.84 ± 0.01 log10 apart across 20 cruises.
- So the switch is a constant shift. The four absolute thresholds on the SSC scale moved by 0.84: `detritus_max_chlSSC`, `syn_size_default`, `pro_max_SSC`, and the GL1-merge SSC threshold in `assign_class_v2`. The Syn-relative rules don't change.
- Rerun of the 40-file comparisons: all counts within 0.5% of the 0.5 µm run (AR92 Syn 39/40 within 10% of old, median 0.99; EN727 37/40, median 1.00). The previous outputs are kept as `output/compare_v2_*_bead05.mat`.


### Per-sample cytograms and agreement range (Nov 2026)

- `plot_class_cytograms`: three rows (existing classes, v2 classes, disagreements only) × six views. The views are size–chl, chl–PE, GL3–PE, size–PE:chl, chl–BL2 (broadband fluorescence) and size–FSC-W, with trigger edges dashed.
- `score_agreement` measures per-event agreement: events with the same class / events either method calls phytoplankton.
- `show_agreement_examples` draws samples evenly spaced from highest to lowest agreement into `output/examples/`.
- **Agreement over 60 files** (15 each from AR92, AR95, EN727, HRS2609): median 0.955, range 0.81–0.99.

**Bug found from the plots:** the lowest-agreement AR92 file (0.78) lost 2,600 dim Syn just above the GL2 trigger edge. The mixture model put them in their own component, whose GL2/GL3 (0.08) read as cryptophyte PE because GL3 sits near its floor. The component was named PE euk, then rejected cell by cell as too small.
- **Fix:** component naming ignores GL2/GL3 when the component's GL3 is less than `ratio_min_above_floor` above its floor. The file is now at 0.954.
- The pre-fix figure is `output/examples/before_fix_lowest_AR92.png`.

**Where the remaining disagreements sit:**
1. **By definition:** class 8 (particle-associated Syn) and the Syn smear. These dominate the lowest-agreement winter files. Hand labels can't settle them; they need a decision about what the class should be.
2. **Boundaries through continuous clouds at a detection edge:** Pro vs noise at the SSC/BL3 trigger (e.g. 1,100 Pro → none and 560 none → Pro in an AR92 Pro-setting file), dim Syn at the GL2 edge, small dim Euks at the BL3 edge. There's no gap, so a human would draw a line as arbitrary as either rule.
3. **Rare groups that are separated in some view:** PE euks along the cryptophyte line, Euks vs detritus at large size. A labeller with all six views could do well here, and the methods disagree most here in relative terms.
4. **Scatter ceiling (fixed):** large events piled up at size ≈ +1 in EN727 and AR95. See *SSC roll-off and the GL1 merge* below.

### Gaps found

- Bead peak medians are missing for 19 older cruises: those with only `beadstat_2021`, plus AR34B, AR63 and RB1904, which have no usable bead table. Re-running bead processing (wrapper step 4) on them would extend bead normalization to the full archive. This matters most for early cruises with OD2 on SSC, where voltage correction alone can't fix the scatter scale.
- Euk median GL2 sits near the noise floor (~100), so Euk PE ratios are noise-dominated. Noise handling must account for the trigger/detection floor in each channel.

## SSC roll-off and the GL1 merge

- The pile-up at size ≈ +1 was not SSC clipping. Raw SSC-H peaks at 5–9×10⁵ and only one event per file is at the top.
- SSC-H rolls off gradually from about 3–10× the 1 µm bead, while GL1 (OD2 scatter) stays linear to +2.7. Overlap slope of SSC on GL1 is 0.67–0.87 (`examine_ceiling.m`, `output/ceiling_*.png`).
- The old prototype merge only replaced events within 0.05 of the SSC maximum, so it almost never switched. Its fitted offset was also biased by the compressed part of the overlap.
- Now the merge follows `batch_mergedSSC2volume`: GL1 on its own 1 µm bead scale (`cfg.bead_pop.GL1 = 2`), used above 1000 raw counts and averaged with SSC between 500 and 1000. In the linear range, SSC minus GL1 (both 1 µm-normalized) is −0.02 to −0.11, so the bead equivalence holds.
- The GL1 switch assumes OD2 is on GL1 (`cfg.gl1_is_od2_scatter`). It needs a per-cruise setting for cruises with OD2 on SSC or no OD2.
- Rescored 60 files: agreement median 0.958 (was 0.955). v2/old Euk counts are within 10% in all 60 files (median 0.99–1.01 by cruise). The detritus chl:size rule didn't remove large Euks.

## Syn pigment types and the GL2:GL3 check

- In example 1 (AR92 15Aug(3) Sample(18)), about 50 bright Syn were dropped by the Syn spectrum check (GL2:GL3 more than 0.15 from the file's Syn median, 0.30). Their PE:chl and scatter were Syn-like. Pulse area/height was about 15% longer than equally bright kept Syn, consistent with elongated, dividing or paired cells (`debug_example1.m`).
- GL2:GL3 rises with chl across this file's Syn (0.25 to 0.40). The mixture components show why: two groups, about 0.26 and 0.34–0.43, i.e. a mixture of Syn pigment types, not one population with a trend.
- A brightness-trend line was tried first. It recovered about half and changed Syn counts by at most 0.5%, but it still assumes one population.
- Now the check accepts GL2:GL3 within ±0.15 of the range spanned by the file's Syn components (Syn-type spectrum, weight ≥ 0.01, mean GL3 ≥ 1 decade above its floor; `syn_ratio_range` in `assign_class_v2`). It falls back to the Syn median when no component qualifies. Class 8 uses the same range.
- The component floor margin is looser than the event-level one (2 decades) because component means are much less noisy. In this summer file Syn GL3 sits only ~1.4 decades above its floor, and with the 2-decade margin no component qualified.
- Example 1: dropped bright Syn 50 → 4; range 0.26–0.40.
- 60 files: agreement unchanged or up by ≤0.002; Syn counts up 0–0.4% per file (one EN727 file +2.4%); PE euk counts unchanged, so cryptophyte-type PE is still excluded.
- `assign_class_v2` now also returns `info.trace` (class after each step) and `info.floor` for debugging.

## Trigger margin set to 0

- v2 used to require BL3 (and GL2, for PE triggering) more than 0.05 log10 above the trigger edge before considering an event. In example 4 (AR92 15Aug(23) Pro(5), trigger SSC AND BL3) that margin alone removed 1,716 of the 1,740 Pro dropped by v2. Those were the dim end of the Pro population (BL3 −1.43 to −1.39, edge −1.44), with SSC like kept Pro (`debug_example4.m`).
- `cfg.sig.trigger_margin` is now 0, so events count down to the trigger edge as in the old method. In that file, Pro dropped by v2 fell from 1,740 to 34; v2 Pro is now about 1.5% above the old count, mainly from low-PE:chl Pro the old gate excluded.
- 60 files: median agreement AR92 0.971 → 0.980 (min 0.905 → 0.956), AR95 0.951 → 0.952, EN727 0.922 → 0.920, HRS2609 0.971 → 0.972.
- EN727 24Jan(1) Sample(1) Syn is now 936 vs 785 in the old classes (+19%; it was +10% with the margin). The extra events are dim Syn near the chl trigger edge. EN727 27Jan(7) Sample(11) is +11%. Not yet examined.
- Counts near the trigger edge are detection-limited either way. A per-file flag for the fraction of each class within ~0.1 log10 of the edge would make that visible in the data product.

## Pro rules and fit reproducibility

- **Reproducible fits.** The mixture fit depends on its random start. `rng(1)` gave a different generator on parallel workers than in serial runs, so scores (parallel) and plots (serial) disagreed. `assign_class_v2` now uses `rng(seed, 'twister')` (`cfg.seed`, default 1).
- **Pro have no PE.** Pro events must have GL2 no more than 0.3 above the file's median Euk GL2 (`pro_GL2_above_euk`), or 1.5 above the GL2 floor when there are fewer than 20 Euks. In AR95 18Oct AltPro(6) Pro(5) this removed half of 257 spurious Pro; in AR92 Pro files it removes <0.6%.
- **Pro need a compact cluster.** A component in the Pro size/chl window is named Pro only if its chl SD ≤ 0.4 (`pro_max_sd_BL3`); otherwise it is left unassigned. Real Pro components in 14 Pro-setting files: chl SD 0.09–0.29; the spurious AR95 18Oct component: 0.50, and its Pro count swung 0–172 across starts. Now 0 Pro for all starts there, real Pro counts unchanged (`examine_pro_components.m`). Calibrated on one spurious case.
  - Refitting with the chosen K from other starts did not separate it (CV 0.14 vs 0.13 for a real file); the instability came from BIC picking different K, and redoing that triples run time.
- **AR95 21Oct AltPro(24) Pro(10):** v2 finds a stable, compact Pro cluster (~1,000 events, chl 0.9 below Syn, smaller, no PE) that the old classes do not have.
- **Pro tail in Euk-setting files.** In files run without Pro settings, the bright tail of Pro piles up on the chl trigger edge. No-PE events within 0.4 of the BL3 edge and no larger than the file's Syn (+0.1) are left unassigned (`edge_pro_band`, `edge_pro_size_above_syn`). AR95 19Oct AltPro(19) Euk(20): v2 Euk 6,147 → 5,390 (old 5,422); HRS2609 17Aug(6) Euk(8): 2,347 → 1,839 (old 1,815). The band of 0.4 removes none of the old-class Euks in the AR95 file.
  - *Removed later (Nov 2026): it sliced the bottom off real small-Euk clusters in summer files. See "Early cruises (before Aug 2019)" at the end.*
- 60 files after these: median agreement AR92 0.979, AR95 0.952, EN727 0.924, HRS2609 0.971.

## Low-PE (3) and high-PE (4) euks

- PE euks form two groups in all four cruises (`examine_peeuk.m`): a diffuse low-PE group with a Syn-like PE spectrum (GL2:GL3 ~ Syn, PE:chl −2 to −1) and compact cryptophyte-type groups (GL2:GL3 0.15–0.35 below Syn, PE:chl −0.8 to −0.4).
- GL2 saturates for 51–65% of PE euks in EN727/AR95 (5–14% in summer). BL2 (blue-laser PE) stays on scale for 84–94% of those, and BL2:GL3 separates cryptophyte-type PE (0.45–1.0 below Syn) without using GL2 (`examine_bl2.m`, `examine_peeuk_split.m`).
- `assign_class_v2` now splits PE euks: 4 if BL2:GL3 is more than 0.3 below the file's Syn (GL2:GL3 more than 0.15 below where BL2/GL3 are unusable; GL2-saturated with no usable BL2/GL3 → 4), else 3. The old scripts' 3 (`in_nonsyn_lowPE`) / 4 (`in_nonsyn_hiPE`) were split mainly by PE amount.
- Against the old 3/4 (PE euks both methods call): old 4 → new 4 for 96–99%. Old 3 → new 4: AR92 83 of 418, AR95 390 of 682, EN727 863 of 962, HRS2609 280 of 669.
  - EN727/AR95: a tight cluster on the cryptophyte line (BL2:GL3 ~1.0 below Syn) with moderate PE (PE:chl ~ −0.8); old method called it low PE by amount.
  - Summer: dim events (PE:chl ~ −1.4) whose BL2:GL3 is ~0.5 below Syn; low PE, so the spectral ratio may be noisy there.
- Agreement with 3/4 merged is unchanged (AR92 0.979, AR95 0.952, EN727 0.924, HRS2609 0.971); strict 3/4 lowers EN727 to 0.908.
- Decision (Nov 2026): keep the split spectral only, with no PE:chl or PE-amount requirement. Class 4 means cryptophyte-type PE, including the moderate-PE winter cluster and dim summer events; class 3 means Syn-like PE spectrum.

## Large Euks absorbed by a detritus component

- In HRS2609 16Aug(13) Euk(9), 54 large bright Euks (size ~1, chl ~2) were in a broad mixture component whose mean has little chl for its scatter (chl:size −0.06), so the component was named detritus (`debug_euk.m`).
- Events in detritus-named components (not outliers) now return to Euk if they sit in the file's Euk chl:size band, have no PE, and are ≥0.5 larger than the file's Syn (`euk_rescue_size_above_syn`). That file: 50 of 54 restored; the other 4 are below its Euk chl:size cut.
- 60 files: Euk counts change 0–0.3% except HRS2609 (up to +6.5%); HRS2609 Euk v2/old median 0.994 → 1.005; agreement unchanged or up.
- Plot colours: low-PE euk (3) cyan, high-PE euk (4) magenta.

## Pro PE cutoff relative to the Pro cluster

- The Pro PE rule (GL2 ≤ Euk median + 0.3) drew a horizontal line through the top of the Pro cloud in AR92 18Aug ProAlt(2) Pro(13), dropping 123 Pro whose GL2 just continued the cluster's distribution.
- Now a Pro event is dropped only if GL2 is also more than 4 robust SDs above the Pro cluster's median GL2 (`pro_GL2_sd`). That file: old Pro dropped 237 → 115; the rest are larger events the fit put in a detritus component.
- Low-PE euks (3): in the 60 scored files v2 calls at most ~100 per sample, as a diffuse group, not a tight cluster (`output/examples/lowPE*.png`).

## Dim crypto (3), bright crypto (4), Euk with Syn (5)

Replaces the low-/high-PE split above. The old scripts intended 3 = dim cryptos (possibly with Syn grazers) and 4 = bright cryptos.
- **5 = Euk with Syn:** Syn-like PE spectrum, PE ~1–1.6 Syn equivalents, chl normal for Euks of that size (`examine_lowpe.m`). Not chance coincidences (low-PE per Euk does not scale with Syn abundance across files). Possible grazers with ingested Syn, or Syn attached; cytometry cannot tell which. Not included in Euk totals: classes are mutually exclusive, and any sums are decided case by case downstream (decision Nov 2026).
  - Events more than 0.5 log10 below the file's Euk chl-size line go to 8 (Syn on a particle): two groups in EN727 (1.2 and 0.24 below) and AR92 (0.8 below, 22%) (`examine_euksyn_en727.m`).
- **3 = dim crypto, 4 = bright crypto:** cryptophyte-type spectrum (BL2:GL3 > 0.3 below Syn). In EN727/AR95 two compact groups: dim (PE ~0.25–0.75 above Syn, BL2:GL3 ~0.85–1.0 below Syn, side scatter ~0) and bright (PE 1.2–2.2, BL2:GL3 ~0.6 below, side scatter 0.2–1.2) (`output/crypto_groups.png`).
  - Per-file 2-group fit on [PE in Syn equivalents, BL2:GL3 − Syn], accepted only if each group ≥15 events, compact (PE SD ≤ 0.4, BL2:GL3 SD ≤ 0.2), ≥3 pooled SDs apart, and the dimmer group more red-shifted. Accepted in all 3 crypto-rich EN727 test files and AR95 22Oct(10); not in summer.
  - Without an accepted split: crypto-type events with PE > 0.5 log10 (~3 Syn) above Syn are 4; dimmer ones are Euks (1).
- 60 files (old 3/4 → v2 3/4/5): AR92 785/91 → 0/142/353; AR95 791/2712 → 208/2291/274; EN727 1054/2334 → 848/1637/110; HRS2609 1111/193 → 0/264/396. Strict agreement EN727 0.907 → 0.918; others unchanged.
- Plot colours: dim crypto cyan, bright crypto magenta, Euk with Syn olive.

## Pro and the detritus chl:size cut

- The fixed detritus cut (chl:size < 0.04) sliced the lower-right edge of the Pro cloud, where the largest Pro sit at the chl trigger edge: in AR92 15Aug(18) Pro(7) it dropped 466 of 477 lost Pro (size ~−1.31 vs −1.53 for kept Pro).
- Pro now fail it only if also more than 4 robust SDs below the Pro cluster's median chl:size (`pro_chlSSC_sd`). That file: old Pro dropped 477 → 12; v2 Pro 27,621 vs old 26,944.

## Large Euks, bright cryptos near GL2 saturation, and the dim/bright PE floor

- **Large Euks in a PE-euk component** (AR95 22Oct(5) Sample(24)): a broad large-particle component held bright cryptos and large Euks and was named PE euk; the Euks (no PE) failed the PE-euk test and became 0. Events in a PE-euk component without PE are now Euks. Lost large Euks 63 → 3.
- **Euk chl:size cut** is now the looser of (Euk median − 0.7) and the midpoint between the Euk median and the median of the file's large non-phyto particles (`nonphyto_min_n` ≥ 50), so variable Euks only need to be closer to the Euks than to the non-phyto cloud.
- **Bright cryptos near GL2 saturation** (AR95 18Oct AltPro(16) Euk(12)): GL2 ~1.4–1.5 gives a Syn-like PE:chl and GL2:GL3 just above 0.1, so 101 old bright cryptos failed the PE-euk test. The test now also accepts BL2:GL3 > 0.3 below Syn. v2 bright crypto 287 → 382 of 393 old.
- **Dim/bright fallback**: when the per-file 2-group split is not accepted, crypto-type events with PE < 1 log10 (~10 Syn) above Syn are dim (3), brighter are bright (4) (`crypto_min_PE_above_syn`). (An earlier version wrongly sent these dim events to Euk.)
- 60 files, old dim/bright → v2 dim/bright/Euk+Syn: AR92 785/91 → 268/127/353; AR95 791/2712 → 1039/2910/274; EN727 1054/2334 → 1119/1703/113; HRS2609 1111/193 → 821/229/396.
- Median agreement AR92 0.981, AR95 0.953, EN727 0.908 (was 0.918; lowest file 0.783, more v2 adds), HRS2609 0.970. Euk v2/old median 1.007–1.029.

## Dim-chl Syn: Syn checks relative to the Syn cluster

- Dim-chl Syn near the chl trigger edge, or detected by PE alone in OR_GL2 files, were dropped by fixed-width checks: HRS2609 SummerSyn(5) Pro(21) lost 90 old Syn (48 detritus cut, 16 PE:chl, 13 too small); HRS2609 16Aug(7) Euk(19) lost 282 (118 PE:chl — chl ~5x below the main Syn with normal PE, PE:chl +0.14 vs −0.56 — 56 too small, 52 detritus cut).
- Now: Syn fail the detritus chl:size cut only if also >4 robust SDs below the Syn cluster's median chl:size; Syn may have up to max(0.35, 5 robust SDs) more PE:chl than the file's Syn (the limit below, which keeps PE euks out, stays fixed at 0.35); the lower size limit is max(0.5, 4 robust SDs) below the Syn median (`syn_PEchl_sd`, `syn_size_sd`).
- Old Syn lost: 90 → 19 and 282 → 63 in those files.
- 60 files after the Syn change: median agreement AR92 0.985, AR95 0.952, EN727 0.908, HRS2609 0.974; median events dropped per file AR92 215 → 83, HRS2609 275 → 93.

## Single crypto group

- AR95 18Oct AltPro(16) Euk(12): the crypto-type events form one compact group (2-group split rejected, groups 1.0 SD apart), and the per-event PE floor cut it into "dim" and "bright" halves.
- Now, without an accepted split: if the crypto-type events form one compact group (robust SD of PE ≤ 0.4 and of BL2:GL3 ≤ 0.2), the whole group is labelled by its median PE (≥ 1 log10 above Syn → bright); only diffuse groups (summer) are split by the PE floor event by event. That file: all bright, 382 of 393 old bright cryptos.
- Single-group outliers (Nov 2026): within a single compact crypto group, events more than 4 robust SDs below (above) the group's median PE are dim (bright) (`crypto_outlier_sd`). AR95 21Oct AltPro(5) Pro(2): 5 low-PE events (PE ~0 vs bright median ~1.35) that had been called bright are now dim; all 80 old bright cryptos stay bright.
- Decision (Nov 2026): option (a). Dim events with a Syn-like PE spectrum are Euk with Syn (5) or, with much less chl than Euks, Syn on particle (8); dim crypto (3) requires a cryptophyte spectrum. Many v1 dim cryptos therefore become 5 (e.g. 14 of 26 in AR95 21Oct AltPro(5)).
- Two-sided single-group fix, 60 files: AR95 dim 437 → 505, bright 3512 → 3444; other cruises within ±10; agreement unchanged.

## Syn-check rejects re-checked as Euks

- AR95 21Oct AltPro(5) Pro(2): 4 v1 dim cryptos were unassigned. The mixture put them in a Syn component, the Syn check rejected them (PE:chl −2.2 to −1.1 vs Syn −0.29), and nothing re-examined them; three were just under the PE-euk size minimum (Syn + 0.6) and one just under the has-PE PE:chl cutoff. They have Euk-level chl for their size, ~0.1–1 Syn's PE and a Syn-like GL2:GL3, like the file's Euk-with-Syn events.
- Events rejected by the Syn check are now re-checked as Euks: chl within 0.5 of the Euk chl-size line, not smaller than Syn − 0.5, above the detritus cut → Euk with Syn (5) if PE is detected, PE:chl > Euk median + 1 (`euksyn_min_PEchl_above_euk`) and the spectrum is Syn-like (or GL3 too low to tell); otherwise Euk (1). That file: all 4 → 5; 42 rejects → 5 and 8 → 1 overall.

## Catch-all for cell-like events

- Principle (Nov 2026): events with PE like a cell and/or chl like a cell for their size should not end up unassigned.
- Final step in `assign_class_v2`: unassigned events with cell-like PE (GL2 no more than 0.5 below the file's median Syn GL2, ~1/3 of a Syn; `cell_PE_below_syn`) or cell-like chl (above the detritus cut and, for particles well larger than Syn, above the file's Euk chl:size cut) take the majority class of their 10 nearest assigned neighbours (standardized size, chl, GL2, BL2, GL3), only if they meet that class's own rules: Euk/Euk-with-Syn no more than 0.5 smaller than Syn; Syn inside the Syn size window; cryptos larger than Syn + 0.6; Pro only in Pro files; Syn-on-particle above its size cut. Otherwise they stay unassigned. The Pro tail at the chl trigger in non-Pro files stays unassigned. *(Later: the Pro-tail cut is gone; the catch-all now keeps out the capped trigger-noise pile instead, and catch-all Euks may not carry Syn-like PE. See "Early cruises (before Aug 2019)".)*
- First version counted any detected PE as cell-like and ignored the class rules: in EN727 it added ~500 faint-PE non-phyto particles per file and ~2,300 tiny signals (size −3 to −2) as Euk/Syn; EN727 median agreement fell 0.908 → 0.822.
- Syn size window now capped at 0.8 below the Syn median (`syn_size_max_below`): the mixture had put tiny PE signals in a Syn component, widening the robust-SD window. EN727 26Jan(13) Syn: 1,811 → 1,135 (old 1,093).
- Open: EN727 24Jan(1) Sample(1) Syn 1,172 vs old 785 (winter Syn smear file).
- 60 files with the constrained catch-all: median agreement AR92 0.982, AR95 0.951, EN727 0.902, HRS2609 0.977 (before the catch-all 0.985 / 0.952 / 0.908 / 0.974). Median old-phyto events dropped per file: AR92 83 → 22, AR95 96 → 8, EN727 74 → 19, HRS2609 93 → 17. Euk v2/old 1.020–1.060; Syn v2/old 1.000–1.017.

## Syn far from the cluster

- AR95 18Oct AltPro(6) Pro(5): the catch-all called 71 events Syn with chl at the trigger edge and PE at noise level (the former spurious Pro); its Syn test only checked size. EN727 24Jan(1): the main rules kept ~450 winter-smear events as Syn. Feedback: upward-widened Syn checks kept smear in Syn, which inflated the Syn size spread (robust SD 0.37 vs 0.13) and moved the Syn-on-particle cut up to +0.35.
- Now: Syn must have GL2 ≥ Syn median − max(0.5, 4 low-side SDs) (`syn_PE_tol_below`, `syn_PE_sd`); the Syn-on-particle cut uses the core Syn centre and spread from the low side of the Syn size distribution (p15.9 + SD, SD = p15.9 − p2.3); Syn larger than that cut that are not Syn on a particle become unassigned; catch-all Syn also need Syn-level PE and a Syn-like PE:chl.
- AR95 18Oct(6): Syn added 91 → 0. EN727 24Jan(1): Syn 1,172 → 702 (old 785), Syn on particle 604.

## Syn on particle vs large Syn, and the Syn centre

- AR95 21Oct AltPro(29) Euk(15): Syn on particle extended up the main chl-vs-size diagonal; 164 large Syn on the upper end of the Syn cluster, whose chl rises with scatter at the Syn slope, were relabelled Syn → Syn on particle.
- Syn on particle (8) now also requires chl more than 0.3 below the core Syn chl-vs-size trend extended to the event's size (`synpart_below_trend`; robust fit on core Syn within 2 SD of the size peak). Large Syn on the trend stay Syn; only large events below the trend that fail the class-8 test become unassigned. That file: Syn → Syn on particle 164 → 45, agreement 0.957 → 0.966.
- The Syn centre (and spread) for the Syn size checks and the class-8 cut is now the peak of a kernel density of Syn size (spread = narrower half-width at half maximum / 1.1774; `peak_centre`). The median is pulled up by the winter smear (EN727 24Jan(1): small core Syn rejected as "too small" and re-checked into Euk with Syn); a low-side estimate is pulled down by tiny PE signals (EN727 26Jan(13): Syn 531 vs old 1,093, Euk 3,610 vs 2,017).
- Test files (v2 Syn vs old): EN727 24Jan(1) 848/785, EN727 26Jan(13) 1,120/1,093, AR95 21Oct(5) 8,629/8,693, HRS2609 16Aug(7) 16,396/16,432, AR92 15Aug(3) 37,127/37,133.
- 60 files after these two changes: median agreement AR92 0.983, AR95 0.951, EN727 0.902, HRS2609 0.977; Syn within 10% of old in all 60 files; Syn v2/old median 0.993–1.025; Euk v2/old 1.020–1.062; Syn on particle AR92 1,655, AR95 1,905, EN727 3,096, HRS2609 2,603; median old-phyto events dropped per file 20–60.

### Catch-all Euks held to the main-rule chl floor (Nov 2026)

The catch-all was adding small, low-chl events to Euk. In AR95 18Oct AltPro(6) Pro(5) it added 118 Euks, about 1.3 log units below the Euk chl-size line, near the chl trigger. These events passed because for particles near Syn size only the detritus chl:scatter cut was applied. Now no catch-all Euk (or Euk with Syn) may have less chl than the Euks the main rules accepted. At any size its chl must be above the `catchall_euk_chl_prc` (1st) percentile of the main-rule Euks' residuals about their chl-size line, and it must also pass the Euk chl:scatter cut. In that file, v2 Euk additions fell from 132 to 34.

Effect on the 60 files:

| Cruise | Median agreement | Euk v2/old median |
|---|---|---|
| AR92 | 0.983 → 0.984 | 1.020 → 1.012 |
| AR95 | 0.951 → 0.953 | 1.041 → 1.033 |
| EN727 | 0.902 → 0.907 | 1.062 → 1.055 |
| HRS2609 | 0.977 → 0.977 | 1.024 → 1.022 |

Syn counts are unchanged. The largest Euk overshoot (v2/old maximum) in AR95 is now 1.16 (was 1.27).

### Syn-level PE events no longer re-checked into plain Euk (Nov 2026)

Some Syn-sized events had Syn-like PE per chl but a PE spectrum (GL2:GL3) outside the file's Syn range by more than `syn_GL2GL3_tol`. Some sat around 0.5 (Syn are about 0.3–0.4) and some around 0.0–0.15. The Syn check rejected them on spectrum alone. Because their spectrum also failed the Euk-with-Syn test, the re-check of rejected Syn made them plain Euk. Euks have about 1000x less PE per chl, so these were plain Euks with Syn-level PE: 286 across the 60 files, almost all called Syn in the old classes.

Now, re-checked events with PE per chl within `syn_PEchl_tol` of Syn are not made plain Euk. They stay unassigned for the following steps. In the 60 files, 482 events were held back this way:

| Final class | Events |
|---|---|
| Syn (catch-all) | 218 |
| Crypto (PE-euk step; larger events with crypto spectrum, called PE euk in the old classes too) | 238 |
| Unassigned | 25 |
| Euk with Syn | 1 |

About 35 plain Euks with Syn-level PE remain. They come straight from GMM Euk components, not from this re-check.

Agreement and Euk and Syn counts change by less than 0.003 in every cruise. EN727 median Syn v2/old goes from 1.025 to 1.028.

### Syn PE floor at a PE gap (Nov 2026)

In OR-trigger files (BL3 or GL2), faint events that clear only the GL2 trigger pile up at its edge. In AR92 14Aug(1) Sample(1) the initial clustering put this pile in Syn. The Syn check removed most of it as too small, but its Syn-sized part (about 170 events at GL2 −1.8 to −1.3, chl below the chl trigger) stayed in Syn. Once the pile has been cut down, those leftovers look like a flat tail with no trough, and the statistical floor (4 low-side SDs, at −2.33 here) removed none of them.

Two alternatives were tried and rejected:

- **Euk no-PE noise floor:** a percentile of Euk GL2. Using all Euks, it wiped out Syn in summer files: at that step the Euks still hold cryptos, and large Euks show PE-channel signal from chl spillover. Using only small Euks, the noise level (about −2.1) sits below the pile, so it removed none of it, and in 3–4 files the percentile was inflated enough to cut real Syn.
- **Gap found after the Syn check:** never fired.

The rule now finds the gap in the Syn PE distribution before the Syn check, while the pile is intact. It looks for a local minimum below the Syn peak with density under `syn_PE_gap_frac` (0.1) of the peak, rising again below it by `syn_PE_gap_rise` (2) into a pile of at least `syn_PE_gap_min_pile` (10) events. Bandwidth is 0.08. The Syn PE floor is the higher of the statistical floor and the gap. Euks are not involved; the gap is found from events the initial clustering called Syn.

In AR92 14Aug(1) the gap sits at −1.26. Syn goes from 7,333 to 7,151, and agreement from 0.960 to 0.963. The old classes had 71 Syn below −1.26 (their own faint tail).

Across the 60 files the gap is found in 24. In most of these it lies below all Syn and has no effect. Elsewhere it removes 1–35 Syn, and 104 in HRS2609 less_sensitive(1). Cruise medians of agreement, Euk and Syn are unchanged; AR92 minimum agreement goes from 0.960 to 0.963. In several AR95 and EN727 files the gap lies only 0.6–1 log units below the Syn median, but it removes at most 9 events there.

### Crypto spectrum relative to the Syn spread; no PE-euk entry on spectrum alone (Nov 2026)

- **Problem:** in your old classes, dim crypto (3) look the same in every file: PE:chl 1.3–1.5 log below Syn, PE about one Syn's worth, size 1.1–1.6 above Syn. Their spectrum is crypto-shifted only in EN727 (BL2:GL3 median 0.88 below Syn). In AR92, AR95 and HRS2609 it sits only 0.05–0.27 below Syn.
- **Why the fixed cut failed:** the v2 crypto cut (BL2:GL3 0.3 below Syn) is 3.7 Syn SDs in EN727 but only 1–1.8 SDs in the AR92 and HRS2609 summer files, where Syn pigment types mix. There v2 sent the old dim crypto to Euk with Syn or Euk. It also pulled Syn-like events (PE:chl at Syn level, about 1 log larger than Syn) into dim crypto: in HRS2609 17Aug(11), 182 of v2's 233 dim crypto were old Syn or unassigned.
- **Changes:**
  - **Relative cut:** the crypto-type spectrum must be below Syn by at least the larger of 0.3 and `crypto_dB_sd` (2) robust SDs of the file's Syn BL2:GL3. The same applies to the GL2:GL3 fallback. At 3 SDs, half of AR92's bright cryptos would have failed.
  - **No entry on spectrum alone:** events entering the PE-euk step on spectrum alone need PE:chl below the Syn median (`peeuk_spec_PEchl_below_syn` = 0), unless GL2 saturates (bright cryptos).
  - **Your decision:** summer and fall low-PE:chl euks with a Syn-like spectrum are Euk with Syn (5).
- **Effect over the 60 files:**

  | Cruise | v2 dim crypto | v2 bright crypto | Euk with Syn |
  |---|---|---|---|
  | AR92 | 287 → 32 | 128 → 129 | 609 → 672 |
  | AR95 | 510 → 466 | 3,637 → 3,635 | 656 → 665 |
  | EN727 | 1,094 → 932 | 1,766 → 1,880 | 244 → 247 |
  | HRS2609 | 995 → 247 | 229 → 211 | 759 → 881 |

  Agreement medians are unchanged (HRS2609 0.977 → 0.978).
- **Still open:** many old dim crypto in summer files are plain Euk in v2. They fail the absolute PE thresholds `pe_min_PEchl` (−2.2) and `peeuk_min_GL2` (−1.0), sitting at PE:chl about −2.4 and GL2 −0.8 to −1.2. That is about 10x more PE per chl than plain Euks (about −3.5).

### PE euks: chl-spillover floor, saturated GL2 = bright crypto, spectrum vs every Syn type (Nov 2026)

The old dim crypto class is known to be a mixture, so agreement with it is not a target.

**Changes:**

- **PE present on large particles.** PE:chl must be at least `peeuk_spill_sd` (3) robust SDs above the large Euks' PE:chl-vs-chl trend. Large Euks show PE-channel signal only from chl spillover, and their PE:chl falls with chl (about −3 at chl 0.5, −4 at chl 2.5). This replaces the fixed `pe_min_PEchl` (−2.2) and `peeuk_min_GL2` (−1.0) cuts, which ran through a separate PE band 1–1.5 log above the Euk cloud in summer files (`output/examples/pe_spillover.png`). The fixed cuts remain the fallback when a file has fewer than 50 large Euks.
- **Saturated GL2 means bright crypto.** PE euks at or near GL2 saturation are bright crypto regardless of spectrum. No other class carries that much PE, and saturation compresses PE:chl and GL2:GL3 toward Syn-like values. The dim/bright two-group fit still uses only events with a crypto-type spectrum, so the EN727 split is still accepted.
- **Spectrum vs every Syn type.** Syn can be a mix of pigment types. HRS2609 17Aug(11) has a second type about 0.7 below the main one in BL2:GL3, and v2's dim crypto there had exactly that spectrum (`output/examples/dim_crypto_spectra.png`). An unsaturated PE euk is crypto-type only if its spectrum is below the 1st percentile of the file's Syn spectrum minus 0.1 (`crypto_syn_low_prc`, `crypto_syn_low_margin`). Otherwise it is Euk with Syn. The exception is events with more than 10 Syn's worth of PE (`crypto_min_PE_above_syn`): they keep the median-relative cut, because unsaturated bright cryptos in HRS2609 sit within the Syn spectral spread (without this, HRS2609 bright crypto fell from 211 to 58).

**Results, 60 files:**

| Cruise | Dim crypto (v2) | Bright crypto (v2) | Euk with Syn | Euk v2/old (median) | Agreement (median) |
|---|---|---|---|---|---|
| AR92 | 32 → 5 | 129 → 140 | 672 → 1,181 | 1.009 → 0.994 | 0.984 → 0.983 |
| AR95 | 466 → 425 | 3,635 → 3,676 | 665 → 1,065 | 1.033 → 1.021 | 0.953 → 0.945 |
| EN727 | 932 → 947 | 1,880 → 1,879 | 247 → 362 | 1.053 → 1.048 | 0.907 → 0.904 |
| HRS2609 | 247 → 26 | 211 → 217 | 881 → 1,711 | 1.022 → 1.001 | 0.978 → 0.974 |

- In summer files, the low-PE:chl euks with a Syn-type spectrum are now Euk with Syn, as decided.
- The small agreement drops come from the old dim crypto class, which mixes groups v2 now separates.

### PE euks: spill floor at the residual valley, saturated events by spectrum (Nov 2026)

**Spill floor (`peeuk_spill_method = 'valley'`).** The large-Euk fit set can hold PE euks and particles. In AT46 C17N02 (preserved CHL) it held 736 hand Euks, 290 hand PE euks and 434 particles. The robust SD came out at 0.47 instead of the Euks' 0.20, so the floor sat above every crypto (12 auto vs 383 hand PE euks).
- **Rule.** The floor is still trend + 3 robust SDs, but no higher than the density valley below a separate upper mode of the residuals. The upper mode's peak must be at least 0.3 of the dominant peak, at least 0.6 above it, and the valley at most 0.6× the lower peak (`peeuk_spill_peak_rel`, `peeuk_spill_min_sep`, `peeuk_spill_max_valley`).
- **Not adopted.** Trim-and-refit fails at ~50% contamination (SD stays 0.40–0.47). Fitting on the lowest residual mode put the floor inside a second Euk mode with no PE in summer runs: EN655 C10N13 CHL went from 2 to 969 PE euks (hand 3). Two Euk modes sit 0.4–0.6 apart with a shallow valley (0.68–0.87 of the lower peak). PE modes sit 0.8–1.1 above with a deep one (0.36–0.49).
- **Held-out (607 CHL/PE/pro runs).** The rule fires in 2 runs: AT46 C17N02 agreement 0.490 → 0.600 (PE euks 12 → 408, hand 383) and EN661 C12N03 0.703 → 0.757 (0 → 61, hand 60). No other run changes. Live (7 cruises × 40 files): no file changes.

**Saturated events (`peeuk_sat_spectrum = true`).** Chl spillover and Syn PE share the same BL2:GL3 spectrum.
- **GL2 saturated, BL2/GL3 on scale.** Such events are crypto only with a crypto-type spectrum.
- **BL3 saturated.** Chl is capped, so PE:chl is overstated. An event with a usable non-crypto spectrum is not a PE euk.
- **Example.** In EN627 C45N05 about 10 BL3-saturated bright Euks had been forced to bright crypto. About 20 more hand-gated "Euks" there carry crypto-type PE far above spillover and stay crypto.
- **Held-out (on top of the valley floor).** 151 runs change: 138 up, 13 down (worst −0.0008). PE euks within 20% of hand: tuning CHL 0.41 → 0.44, 2a 0.50 → 0.53.
- **Live.** Cruise means are flat or up (AR95 0.9338 → 0.9344, EN727 0.9281 → 0.9280).

**Set means, both changes (held-out, mean agreement):** tuning 0.8884 → 0.8889, 2a 0.8911 → 0.8920, 2a′ 0.9228 → 0.9230, 2b 0.9059 → 0.9059.

**Class 6, PE euk cluster (`peeuk_cluster`, on); counted as Euk since 27 Nov 2026.** A compact population of class-5 events with apparent PE of their own. In AT46 C17N02, 270 class-5 events form one tight cluster.
- **Spectrum.** After removing Euk spillover, GL2:GL3 is 0.62 (Syn 0.29) and BL2−GL2 is −0.42.
- **PE tracks chl.** PE rises with chl (slope 0.94), with PE:chl robust SD 0.13. Poisson ingestion of the file's Syn would give 0.23–0.41.
- **Rule.** Class-5 events with usable spectra become class 6 when their PE:chl and GL2:GL3 robust SDs are at most 0.18 and 0.06, the events near the centre have GL2 rising with BL3 (slope 0.7–1.3, r ≥ 0.7), and their median GL2:GL3 is at least 0.15 above the file's Syn. Tightness and slope alone are not specific: in many preserved runs class 5 is tight with PE tracking chl; the spectrum is what is shared.
- **Where it fires.** 5 of 607 held-out runs: AT46 C17N02, EN727 C20N02, AR43 C02N24 PE, MVCO MV417N04, MV432N04 PE. All five share GL2:GL3 0.55–0.61 and BL2−GL2 −0.34 to −0.42. It never fires in live files (7 cruises × 40 files), where class 5 has PE:chl SD 0.4–1.6 and PE independent of chl. Not Crocosphaera-like: PE:chl ~13x below Syn, cells larger than Euks, 4.7 °C water (AT46).
- **What it is (27 Nov 2026): Euks whose PE channels are raised by preservation, not a PE population.** Each class-6 sample was compared with the live file nearest in time (AT46 C17N02 75 min, AR43 C02N24 1 min, EN727 C02N02 53 min, EN727 C20N12 1 min). The live sample has the same Euks (same size and chl) with GL2 1.5–2.3 log10 lower and GL2:GL3 −0.3 to −0.5 (preserved 0.53–0.61), and no events where class 6 sits. The whole preserved Euk cloud shifts up on GL2 (+1.3 to +2.5) while Syn moves ~0.3. Class 6 is the top of that shifted cloud, where the extra signal is tight and tracks chl. Hand gates in AT46 C17N02 follow the same shift (383 hand PE euks).
- **Downstream.** Counted as Euk: in the Euk columns (`Euk_celltypes = [1 6]`) with its own `Euk_PE_shift` columns (`generate_attune_table`); in `euk` and `euk_leqXum` with its own `euk_pe_shift_*` subset in `preserved_table` and the `merge_discrete_tables` extra of the same name; as Euk (1) in the held-out scoring, `live_holdout_*` and `score_agreement`'s merged agreement. Kept as its own code in the class files, `n6_peeuk_cluster` in the class log, and the review (brown, "Euk, PE shift (preserved)").
- **EN727 C20N02 and `dimcrypto`.** The 860 events the discrete `dimcrypto` rule moved to dim crypto are the same population as v2's class 6 there: size 0.75 vs 0.80, GL2:GL3 0.60 vs 0.61, BL2−GL2 −0.32 vs −0.34, dB +0.05 vs +0.03 (the run's cryptos: dB −0.7 to −0.8). v2's re-check of Syn-check rejects split the cloud at a PE:chl cut (below: Euk, above: class 5, then 6). `dimcrypto` now sends its group to class 6 when the group's dB is at or above the crypto cut and its GL2:GL3 is at least 0.15 above Syn; crypto-type groups still go to 3. EN727 C20N02: class 6 192 → 1,052 events, dim crypto 901 → 41.

**Dim/bright crypto split: robust compactness (`crypto_robust_compact`, on; 27 Nov 2026).** The two-group fit's compactness check now uses each group's robust SD (1.4826 MAD of its hard-assigned members) instead of the Gaussian SD. In EN727 27Jan2025(4) Sample(41) a few dim stragglers gave the dim group a dB SD of 0.217 (limit 0.2; robust ~0.1). That rejected a clear split (separation 4.3). The single-group fallback's median PE then fell at 0.999, just under the 10-Syn floor, and 39 unsaturated bright cryptos were called dim. With the change: 31 dim / 100 bright (was 68 / 63).
- **Held-out:** hand low-/high-PE euk vs auto dim/bright match: tuning 0.908 → 0.932, 2a 0.945 → 0.947, 2a′ 0.798 (same), 2b 0.825 → 0.824. Split accepted in 23 → 29, 25 → 28, 0 → 0 and 8 → 9 runs. Overall agreement unchanged in every set.
- **Live (7 cruises × 40 files):** split accepted in more files (EN727 36 → 40, AR77 1 → 7, AR95 4 → 5); no file's agreement moves by more than 0.005.

**Class 6 across all preserved CHL runs (27 Nov 2026).** Class 6 fires in 31 of 2,086 runs.
- **Season.** On the NES transect it fires in 3.3% of winter runs and 0.4% of summer runs.
- **Place.** It concentrates on the winter inner shelf (> 40.8°N, Dec–Mar: 12 of 120 runs, 10%). Every cluster of about 250 or more events is there: EN727 cast 20, EN627 C12/C07, AT46 C17 and AR99 C03.
- **Storage time is not the cause.** Storage runs from 30 to 2,400 days, and the large-Euk shift does not track it (Spearman −0.04).
- **The shift itself is everywhere.** Large preserved Euks sit about 0.1 above Syn in GL2:GL3 in nearly all runs (winter median 0.12, summer 0.09). Class 6 appears where one abundant Euk type shifts strongly and uniformly.

## Cryptos: spectral dim/bright split, catch-all spectrum test, Syn clumps, Euk-with-Syn PE ceiling (27 Nov 2026)

**Dim (3) / bright (4) by spectrum in live files** (`crypto_split_method = 'spectrum'`, since replaced by `'gmm2d'` for live and preserved: see the two-step rule below).
- **What changes.** Live crypto-type events with a usable BL2:GL3 are dim when BL2:GL3 is more than `crypto_red_dB` (−0.77) below the file's Syn, else bright. GL2-saturated cryptos are always bright. Events without a usable spectrum are split by the 10-Syn PE floor.
- **Why.** In 60 live files each from EN727, AR77 and AR95, the crypto-type events form a red-shifted group (BL2:GL3 0.85–1.0 below Syn) and a main group (0.52–0.67 below). Each is compact in spectrum (robust SD 0.04–0.08) and flat across its PE range (slope with PE −0.08 to +0.04).
  - In AR77 the two groups separate in 44 of 57 files.
  - PE amount does not separate them. AR77's red group holds 2–6 Syn of PE, and its main group runs from 0.1 to 100 Syn without a gap. AR95 has a red group at about 5 Syn. The old PE-based split cut AR77's main group in two: two-thirds of its "dim" were the faint end of the main group.
- **Saturated events.** GL2-saturated cryptos read low in BL2:GL3 as they brighten. In AR95 they spilled below −0.77, while the unsaturated main peak sat at −0.65, clear of the cut. They are bright regardless.
- **Preserved samples (first reading, corrected below).** The spectral gap seemed to close there, and hand low-/high-PE gates in 66 held-out runs follow PE (AUC 0.99) more than spectrum (0.85). The closed gap came from GL2-saturated events.
- **Live scoring (7 cruises × 40 files).** Only dim/bright change; other classes change by at most 90 events in 40 files, through the catch-all.
  - Dim (3): EN727 2,605 → 2,607, AR77 1,514 → 965, EN706 469 → 87, HRS2609 283 → 87.
  - Agreement with the old (PE-split) classes falls by at most 0.0016.
- **Per-file or pooled cut (code kept, off: `crypto_red_valley`).**
  - A per-file cut from two Gaussian components on BL2:GL3 tracks EN727's drift (groups move by ~0.1 over the cruise). But it jumps between −0.70 and −0.90 file to file in AR77, and is rejected in 52 of 60 AR95 files, where the red group is small and often absent.
  - Pooling over ±12 h smooths it (median jump 0.005). With saturated events kept bright, the fixed −0.77 falls in the gap in all three cruises, so no pooling pass was added to `process_cruise`.

**Two-step rule, live and preserved (`crypto_split_method = 'gmm2d'`, 28 Nov 2026).** One rule for both; `discrete_config` no longer keeps `'gmm'`.
1. GL2-saturated cryptos are bright.
2. A red-shifted spectrum is dim: BL2:GL3 more than `crypto_red_dB` (−0.77) below the file's Syn.
3. Among the unsaturated main-spectrum events, a two-component fit on PE is accepted when each group has at least `crypto_min_group_n` events, both are compact (robust SD ≤ `crypto_max_sd_PE`) and their means are at least `crypto_min_sep` pooled SDs apart. The lower-PE group is then dim. Events without a usable spectrum take the nearer PE centre.
4. Catch-all cryptos follow the same steps (`info.crypto_pe_gmm`). The split is recorded in `info.crypto_split` (`pe_split`, `pe_sep`, `pe_centres`, `pe_sd`, `pe_n`, `n_pe_dim`).
- **Why preserved changed.** Without GL2-saturated events the preserved spectra show the same red (−0.91) and main (−0.61 to −0.64) groups as live. The per-file `'gmm'` fit split AR77's single main group by PE. Spectrum alone missed preserved runs with a clear low-PE group in the main spectrum (EN668 C001N01, MV428N03, AR92 C25N02), hence step 3. A one-step 2-D fit on [PE, spectrum] mis-split AR77 Sample(43) live (15 main-group events dim), hence the order.
- **Live (7 cruises × 40 files).** Step 3 fires in 3 AR95 files (dim 721 → 872) and 1 EN706 file (87 → 108); all else matches the spectral rule. Agreement with the old classes moves by at most 0.0001.
- **Preserved held-out, hand low-/high-PE match** (old per-file fit / spectral only / two-step): tuning 0.932 / 0.816 / 0.910; 2a 0.948 / 0.780 / 0.911; 2a′ 0.798 / 0.963 / 0.963; 2b 0.824 / 0.760 / 0.855; all 0.913 / 0.794 / 0.900. Merged crypto agreement is identical for all three. The two-step rule sits slightly below the old fit because it does not split a single gap-free spectral group by PE (mainly AR77); the hand gates do.
- **Examples.** EN668 C001N01: 426 of 483 hand low-PE are dim. MV428N03: 269 of 314. AR92 C25N02: 146 of 227.

**Crypto spectrum against one Syn type (`crypto_syn_ref_window = [-0.40 -0.10]`, 28 Nov 2026).** Step 2 measures BL2:GL3 against the median of Syn with BL2−GL2 in that window (type B since type D was added), not all Syn. The cut is then `crypto_red_dB_ref` (−0.70). Files with fewer than `crypto_syn_ref_min_n` (50) such Syn keep the all-Syn reference and −0.77. The offset is in `info.crypto_ref_offset`, the count in `info.crypto_ref_n`; catch-all cryptos use the same offset.
- **Why.** The all-Syn median BL2:GL3 moves with the file's Syn type mix (0.09–0.75 over 598 live files, 20 cruises). Against all Syn, crypto dB drops as the file's share of high-PUB Syn (BL2−GL2 > 0.19) rises; against the window it does not. Median crypto dB (all Syn / window) by high-PUB share: below 5%, −0.86 / −0.84; 5–30%, −0.69 / −0.49; above 30%, −0.93 / −0.51. AR70B 17Nov2022 Sample(50) has 4,670 high-PUB Syn of 8,148, so the all-Syn reference is 0.64 too high. Its cryptos (absolute BL2−GL3 −0.20) are main-type but read as red.
- **Crypto types.** Two spectral groups, not three: per-file peaks at absolute BL2−GL3 −0.77 and −0.34 (SD 0.042 / 0.037; summer main −0.47). They differ in emission (GL2−GL3 −0.15 vs +0.19), not excitation (BL2−GL2 one band, −0.60 to −0.52), and merge below about Syn PE. Saturated events (GL2, BL2, GL3) were left out.
- **Window Syn.** At least 100 events in 99.3% of the 598 files. Their BL2−GL3 robust SD across files is 0.047 (all Syn 0.119). Against them the per-file red crypto peaks spread 0.022 (all Syn 0.059), and the event valley at −0.70 is 0.12 of the smaller peak deep (all Syn 0.31 at −0.77).
- **Live (8 cruises × 40 files), dim with all-Syn / window reference.** AR95 872 → 653, HRS2609 86 → 23, AR70B 819 → 638 (Samples 50 and 28 now all bright); others within 20. Agreement with the old classes moves by at most 0.0002.
- **Preserved held-out, hand low-/high-PE match** (all Syn / window): tuning 0.910 / 0.935; 2a 0.911 / 0.931; 2a′ 0.963 / 0.963; 2b 0.855 / 0.831; all 0.900 / 0.913. Hand high-PE called dim: 267 → 96. Merged agreement unchanged.

**Saturated cryptos following the red group (`crypto_sat_follow`, off).** GL2-saturated cryptos with a red BL2:GL3 become dim when the file's unsaturated red group reaches GL2 saturation (90th-percentile PE within `crypto_sat_reach` of it, at least `crypto_sat_min_n` events), the main group does not, and the saturated events' median is below the cut. It was written for AR70B Samples 50 and 28, which looked like a single red group running into saturation. That was the all-Syn reference artifact above. With the window reference the rule fires in one preserved 2b run (39 events) and no live file, so it is off; with the all-Syn reference it lowered the preserved match from 0.900 to 0.880.

**Catch-all cryptos.** Unassigned events whose neighbours are cryptos now need a crypto-type spectrum (`catchall_crypto_spectrum`): BL2:GL3 below the crypto cut, GL2:GL3 below `peeuk_max_GL2GL3` only where BL2:GL3 is unusable, or GL2 saturated. They take dim/bright by the spectral rule, not their neighbours' label.
- AR77 11Oct2023(12) Sample(43): 223 catch-all cryptos had a Syn spectrum (dB −0.13) and Syn PE:chl; 52 remain, with dB −0.61.
- A BL3-saturated event with a Syn-like BL2:GL3 had passed on GL2:GL3 alone.

**Class 11, Syn clump (`catchall_syn_clump_class = 11`).** Those catch-all events with a Syn-type spectrum and Syn PE:chl: PE and chl about 4–5 Syn's worth on the Syn line, about 1.5 log10 larger than Syn. They are not counted as Syn.
- **Where it is carried:**
  - `n11_syn_clump` in the class log;
  - `Syn_clump` columns in `generate_attune_table`;
  - `syn_clump` in `preserved_table` and the `merge_discrete_tables` extra;
  - the review config, QC, viewer (grey) and guide.
- **Effect:**
  - Live: crypto counts fall 1–10% (AR77 9,341 → 8,694; AR95 12,237 → 11,573), with 31–735 class 11 per 40 files; most were unassigned in the old classes.
  - Held-out: means change by at most 0.0001 (tuning 0.8889 → 0.8888, 2a 0.8920 → 0.8919, 2a′ and 2b unchanged); 16 of about 2,250 class-11 events were hand Syn.

**Euk with Syn PE ceiling (`euksyn_max_syn_PE_prc`, not set: off).** Class 5 above the file's Syn GL2 99th percentile (BL2 against the Syn BL2 percentile where GL2 saturates) becomes Euk (1). It applies after class 6, Syn on particle and the catch-all.
- **The `cryptokeep` fix.** Discrete `cryptokeep` now skips events moved by the ceiling. Their GL2:GL3 is compressed by GL2 saturation: MVCO MV418N01 would otherwise have turned 197 of them into dim crypto.
- **Live held-out.** PE euk totals move onto the old classes: T 9,314 → 9,042 (old 9,071), L1 45,023 → 43,774 (43,818), L2 21,970 → 21,269 (21,228). Merged agreement changes by at most 0.0005; 1–2 files per set flip pass/fail near the ±20% limit.
- **Preserved held-out.** Set means: tuning 0.8887 (unchanged), 2a 0.8918 → 0.8920, 2a′ 0.9230 (unchanged), 2b 0.9058 → 0.9060. The 3,319 moved events are hand Euk 1,099, PE euk 339, Syn 108 and none 1,773.
- **Examples.**
  - AR77 pro_C13N16: only 3 of 39 move. This pro run's Syn GL2 spread is wide (p99 0.9 above the median), and most of its class 5 sits at or below the Syn median.
  - MVCO MV418N01: 201 of 427 move, 199 of them GL2-saturated and 197 hand-unlabelled.

## Syn spectral types across cruises (`syn_spectral_types.m`, Nov 2026)

- **Method.** 148 files from 20 cruises (all those with `beadstat_2026`, up to 8 standard files each) were classified with v2. Syn with usable BL2, GL2 and GL3 (not saturated, at least 1 log above floor) were described by bead-normalized ratios:
  - BL2 − GL2: PE excited at 488 vs 532 nm, a proxy for PUB:PEB excitation
  - GL2 − GL3: PE emission ratio
  - GL4 − GL2
  - PE:chl
- **Outputs.** Per-file 1–3 component mixtures are in `output/syn_spectral_types.mat` (C = groups, P = 1000 events per file); the figure is `output/syn_spectral_types.png`.
- **Result.** BL2 − GL2 has three sharp modes that recur across cruises and years, at about **−0.23, 0.00 and +0.45**. Within-group SD is about 0.07, and the emission ratio GL2 − GL3 is about 0.35 for all three.
- **Season.** Winter and spring cruises (SST about 8–12 °C) have only the −0.23 and 0.00 modes. Summer and fall cruises (SST over 20 °C) add the +0.45 mode, a higher PUB:PEB type consistent with 3c-like or offshore Syn. The −0.23 mode also has higher GL4:GL2.
- **Lower-emission group.** A smaller group with lower emission ratio (GL2 − GL3 about 0.1–0.2) appears in some warm-water files.
- **Caveats:**
  - The per-file mixture nearly always picks K = 3, so it over-splits; the pooled event histogram shows the same modes independently.
  - Naming the modes as pigment types (3a, 3b, 3c, 3d) needs reference cultures run with beads.
  - Bead-lot differences in the beads' own spectra could shift a cruise slightly, but the modes line up across 2020–2026.

### Syn subtype label in v2 (Nov 2026)

*Superseded 26 Nov 2026 by four types A/B/C/D: the A/B/C below are now B/C/D (see "Fourth Syn type" at the end).*

- **Labels.** `assign_class_v2` returns `info.syn_type`, one label per event, for Syn (2), Euk with Syn (5) and Syn on particle (8):
  - 1 = type A, 2 = B, 3 = C
  - 9 = spectrum not usable (BL2, GL2 or GL3 saturated or within `syn_type_min_above_floor` = 1 of the floor)
  - 0 = other classes
- **Boundaries.** The types are fixed windows on the bead-normalized BL2 − GL2: A < −0.10 ≤ B < +0.19 ≤ C (`syn_type_edges`). These are the troughs between the pooled modes at −0.23, 0.00 and +0.45. The class codes are unchanged; the subtype is a separate label, so class 2 counts and downstream code are not affected.
- **Outputs.** `score_agreement` adds `new_synA`, `new_synB`, `new_synC` and `new_syn_untyped`, counted for class 2. `plot_class_cytograms` adds a GL2 vs BL2 panel, with v2 Syn coloured by type and dashed lines at the boundaries.
- **60 files.** Classes are identical to the previous run.

  | Cruise | A (median share per file) | B | C | Untyped | Files with more than 5% C |
  |---|---|---|---|---|---|
  | AR92 (Aug) | 0.40 | 0.06 | 0.38 | 0.08 | 14 of 15 |
  | AR95 (Oct) | 0.16 | 0.69 | 0.04 | 0 | 6 of 15 |
  | EN727 (Jan) | 0.55 | 0.44 | 0.002 | 0 | 0 of 15 |
  | HRS2609 (Aug) | 0.27 | 0.05 | 0.47 | 0.08 | 13 of 15 |

  Untyped Syn (about 8% in the summer files) are mostly dim Syn near the GL3 or BL2 floor.

### Chl gap at the bottom of the Euks (Nov 2026)

- **Problem.** In AR92 14Aug(1) Sample(1), v2 added 226 Euks that sat on the pile of faint events at the chl trigger: median 0.24 above the edge, below every agreed Euk in chl:scatter. The chl histogram of small no-PE particles there shows a pile at the edge, a trough about 0.35 above it, then the Euk cluster.
- **Rule.** Same approach as the Syn PE gap. Take small no-PE particles (Euks and unassigned, smaller than Syn + 1) and find a deep minimum in their chl (BL3) between the Euk peak (at least 0.2 above the edge) and a pile below it:
  - minimum density under 0.1 of the peak (`euk_gap_frac`)
  - density rising again below it by at least 1.5x (`euk_gap_rise`)
  - at least 20 events below the minimum (`euk_gap_min_pile`)
  - KDE bandwidth 0.05

  Where the gap is found, small no-PE Euks below it are left unassigned, and the catch-all may not re-add them. Files without a gap keep the existing rules.
- **AR92 14Aug(1).** The gap sits 0.33 above the edge. Euks added relative to the old classes drop from 226 to 91; the remaining ones are larger (about Syn + 0.9) with low chl:scatter.
- **60 files.** The gap is found in 43 of 60 files.

  | Cruise | Agreement median | Agreement min | Euk v2/old median |
  |---|---|---|---|
  | AR92 | 0.983 → 0.983 | 0.957 → 0.960 | 0.994 → 0.991 |
  | AR95 | 0.945 → 0.951 | 0.818 → 0.848 | 1.021 → 0.989 |
  | EN727 | 0.904 → 0.920 | 0.780 → 0.787 | 1.048 → 1.023 |
  | HRS2609 | 0.974 → 0.976 | 0.958 → 0.958 | 1.001 → 0.996 |

  Euk v2/old across all files now ranges 0.96–1.19 (was up to 1.25).
- **To check.** In Pro-setting files (AR92 and AR95 `Pro` files) the gap is 1.1–1.7 above the chl edge, because the pile there is the faint Pro tail. Euk counts in those files stay within a few percent of the old counts, for example AR95 19Oct AltPro(30) Pro(10) at 3,046 vs 3,123.

### Euk with Syn / Syn on particle side paths fixed; low-chl Euk with Syn unassigned and flagged (Nov 2026)

- **Bug: side paths bypassed the class rules.** Two later steps could assign Euk with Syn or Syn on particle without checking the class's own rules. In HRS2609 17Aug(11), 61 of 286 Euk with Syn and 159 of 553 Syn on particle came in this way.
  - **Syn-reject re-check.** It assigned Euk with Syn using only the Euk minimum size (Syn − 0.5) and a Syn-type spectrum. It now also requires the PE-euk rules: size larger than Syn + 0.6 and PE:chl at least 0.5 below Syn. Events that fail stay unassigned for later steps.
  - **Catch-all, Euk with Syn.** It now requires the PE-euk size and PE:chl rules and chl within `euksyn_max_chl_below_euk` of the Euk line.
  - **Catch-all, Syn on particle.** It checked only the size cut, and every event it added failed at least one Syn-on-particle rule. It now requires the full main rule (`synpart_rule`).
- **Low-chl Euk with Syn.** These sit more than 0.5 below the Euk chl-size line.
  - Their PE:chl is 0.5–2 below Syn and their chl is 1–2 log above one Syn's, so they fit neither Euk with Syn nor Syn on particle. They are likely detritus or grazers with partly degraded Syn and/or Euks.
  - They are now left unassigned (`euksyn_low_chl_class` = 0; it was 8). The main Syn-on-particle rule can still take those that fit it, and the catch-all may not re-add them.
  - Plot: `output/examples/euksyn_lowchl_unassigned.png`.
- **Effect over the 60 files.** Euk, Syn and crypto counts are essentially unchanged.

  | Cruise | Euk with Syn | Syn on particle | Agreement median | Agreement min |
  |---|---|---|---|---|
  | AR92 | 1,167 → 947 | 1,948 → 1,097 | 0.983 → 0.985 | 0.960 → 0.966 |
  | AR95 | 1,030 → 758 | 2,378 → 1,490 | 0.951 → 0.955 | 0.848 → 0.849 |
  | EN727 | 363 → 249 | 3,457 → 2,506 | 0.920 → 0.929 | 0.787 → 0.807 |
  | HRS2609 | 1,679 → 1,372 | 3,222 → 2,025 | 0.976 → 0.978 | 0.958 → 0.963 |

- **New flag variable.** `info.flag` marks why some unassigned events were left out; class stays 0:

  | Flag | Meaning | 60-file total | Median per file |
  |---|---|---|---|
  | 1 | Syn-pigment particle with low chl | 1,660 | 8–29 |
  | 2 | Faint particle below the Euk chl gap | 20,102 | 0–402 |
  | 3 | Pro tail at the chl trigger (non-Pro files) | 376 | 0–4 |

  `score_agreement` reports per-file counts as `flag_synpig`, `flag_below_gap` and `flag_pro_tail`.
  - *Flag 3 now means the trigger-noise pile at the chl trigger in capped files (before Aug 2019), since the Pro-tail cut was removed. The `score_agreement` column is now `flag_edge_pile`; agreement files saved before Nov 2026 still have `flag_pro_tail`.*
- **Planned for integration.** The per-sample class files would carry `syn_type` and `flag` as extra per-event variables next to `class`, with `notes` describing the codes. Code that reads only `class` is unaffected.

### New hv exponents and bead stats; Pro chl:scatter rule (Nov 2026)

- **New hv exponents and bead stats.** Two changes from sibling sessions:
  - commit c5c6980: hv exponents SSC 8.5, FSC 2.1
  - commit c9f722b: new `beadstat.mat` read by `bead_reference`

  Rescoring the 60 files changed almost nothing: agreement medians within 0.002, Syn and Euk ratios within 0.003. The exception is AR95 19Oct AltPro(30) Pro(10), where agreement fell from 0.961 to 0.905.
- **Why AR95 19Oct AltPro(30) changed.** Neither change was the real cause:
  - With the old SSC exponent the result was the same.
  - The bead positions moved 0.011 log or less.

  The small shift tipped the initial mixture fit into forming a compact, low-chl group (4.5% of events, chl SD 0.15) and naming it Pro. That gave 579 events the old classes leave unassigned. They are the bottom of the faint smear in an SSC-triggered file: possibly some real Pro mixed with detritus and small particles, as seen near shore in well-mixed seasons. The decision was not to count them as Pro.
- **Fix: Pro chl:scatter rule.** A Pro component's chl:scatter must be within `pro_max_chlSSC_below_syn` (1.0) of the Syn components' (weighted mean).
  - Real Pro clusters in 14 Pro-setting files sit 0.28–0.82 below Syn; the false group sits 1.2 below.
  - AltPro(30) agreement returns to 0.958, with no Pro.
  - Every other file is unchanged.
- **Open case.** AR95 21Oct AltPro(24) Pro(10) has 1,296 v2 Pro against none in the old classes, and it is AR95's lowest agreement (0.849). It is a separate compact cluster, 0.27 smaller than Syn, chl:scatter 0.64 below Syn, no PE. It looks like real Pro the old classifier missed. Plot: `output/examples/prochk_1_AR95.png`.
- **Pro catch-all fix.** You judged that AR95 21Oct AltPro(24) has real Pro, but not the separate group with higher PE (chl −1 to −0.5, GL2 about −1.6). That group came from the catch-all, which accepted Pro whenever the file was Pro-setting.
  - The catch-all now requires the main Pro rules: GL2 at or below the file's Pro PE ceiling, SSC < `pro_max_SSC`, BL3 < `pro_max_BL3`, and chl:scatter within `pro_chlSSC_sd` robust SDs of the file's Pro.
  - In AltPro(24), Pro goes 1,296 → 1,110 and agreement 0.849 → 0.865.
  - The other Pro-setting files lose 0–150 catch-all Pro each (at most about 3%), all moving closer to the old counts.
  - AR95 median agreement 0.955 → 0.958. Nothing else changed.
  - A PE-gap ceiling for Pro was tried and removed, because these events never passed the main rule.

### Syn typing needs only BL2 and GL2 (Nov 2026)

- **Change.** The Syn type uses BL2 − GL2 only. It previously also required GL3 to be at least 1 log above its floor. Dim summer Syn often miss that, so 23–30% of Syn were untyped in AR92 15Aug(13) and HRS2609 15Aug(4), almost all only because of GL3. Now only BL2 and GL2 must be on scale (not saturated, at least `syn_type_min_above_floor` above floor).
- **Effect.** The median untyped share is now 0 in every cruise, and at most 0.3% in any file (was up to 61%).

  | Cruise | Syn A | Syn B | Syn C |
  |---|---|---|---|
  | AR92 | 105,619 → 114,329 | 17,542 → 22,574 | 92,922 → 110,657 |
  | HRS2609 | 83,842 → 88,728 | 15,222 → 18,563 | 112,673 → 129,516 |

- **Other changes in this rescore.** Small class changes also appeared, for example AR95 18Oct AltPro(11) Pro(8) 0.938 → 0.922 and HRS2609 below-gap flags 2,819 → 2,246. They come from concurrent edits to `assign_class_v2.m` by the pre-2019 channel-mapping session, not from the typing change.

## Early cruises (before Aug 2019)

### Channel map from `attune_config` (Nov 2026)

- **Problem.** Before the 16 Aug 2019 filter change (`attune_config(...).layout == 'v2019'`: EN608, AR28B, AR29, EN617, AR31A, EN627, AR34B, RB1904, TN368) the green-laser parameters are shifted, and there is no green side scatter:

  | Physical channel | File parameter before Aug 2019 |
  |---|---|
  | GL1 (530/30 + OD2 green scatter) | none |
  | GL2 (575/36, PE) | `GL1-H` |
  | GL3 (620/15) | `GL2-H` |
  | GL4 (680/40) | `GL3-H` |

  v2 read raw `<ch>-H` names, so on these cruises it read the wrong channels.
- **Fix.**
  - `bead_reference` stores `ref.chan` and `ref.layout` from `attune_config` for the sample time, plus `ref.OD2` when `cfg.cruise` is set.
  - `attune_features` reads each physical channel through `ref.chan`. A missing channel (physical GL1 before 2019) becomes a NaN column with no floor. Trigger names are translated too (file `OR_GL1` is physical GL2 before 2019; raw names in `info.trigger_file`). New `info.col` gives the file column for each physical channel; the saturation checks use it.
  - The GL1 size merge runs only if `cfg.gl1_is_od2_scatter` is true, GL1 has values, and `ref.OD2` (when known) is `'GL1'`; `info.size_from_gl1` records it. Callers set `cfg.cruise = cruise` and `cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1')`.
  - `beadstat` needed no change: the FCB bead pipeline and the PT transfer already write physical channel names (GL1 = NaN before 2019).
- **Effect.** After 2019 every class, Syn type and flag is identical. Early-cruise agreement with the old classes, median over 12 files: TN368 0.69 → 0.89, AR34B 0.77 → 0.82, RB1904 0.76 → 0.93, EN627 0.54 → 0.90 (raw names → channel map).

### Trigger-noise pile and rare Syn (Nov 2026)

- **Problem.** In AR34B, RB1904 and EN627 (no OD2, low trigger thresholds), no-PE events within 0.4 of the chl trigger are 60–98% of the chl-bearing events (non-Pro files after 2019: under 25%; `measure_edge_pile.m`). The mixture spent its components on them: real Syn were absorbed, and a few bright cryptos were named Syn. Some files fell to 0.03–0.25.
- **Changes in `assign_class_v2`:**
  1. **Capped fit.** In non-Pro files where the pile is over `edge_pile_max_frac` (0.5) of chl events, it enters the fit sample at most `edge_pile_fit_ratio` (1) times the other events, with enough events left to fit the largest mixture (`info.edge_pile_frac`, `info.edge_pile_capped`).
  2. **Pile labels.** In those files pile events may only be Euk or unassigned and stay out of the catch-all. Euk-labelled pile events are unassigned at any size, as the old per-cruise classifiers did (raw BL3 < 200 with a trigger of 100). Flag 3.
  3. **Rare Syn.** If no component is named Syn, events with PE, a Syn-type spectrum, and PE:chl and size near the `cfg` Syn defaults start as Syn (at least `syn_rescue_min_n`, 20), then pass the usual Syn checks (`info.syn_rescued`). They must have chl above the chl trigger + 0.4: in RB1904 12May Sample(1) the PE-triggered pile below the chl trigger otherwise became the Syn reference (2,187 false Syn on particle; saturated bright cryptos called dim).
  4. **Syn PE gap.** The gap search uses only Syn-named events no larger than Syn + `syn_size_tol`. A bright-crypto component named Syn had put the gap above the real Syn and removed every one of them (AR34B 16Apr c, 22Apr d).
  5. **Guards.** If no mixture can be fitted, v2 returns all 0 instead of crashing (AR28B 03Apr2018A Sample(10) is 100% pile). Both "nothing to classify" exits now set `info.syn_type`, `info.flag` and `info.features`.
- **Effect.** Median (lowest) agreement over 12 files: AR34B 0.92 (0.78), EN627 0.92 (0.83), RB1904 0.93 (0.86), TN368 0.89 (0.56). Changes 1–3 never trigger after 2019; change 4 moved 1–3 events in 3 of 48 files there.

### Pro-tail cut removed; Euks carry no Syn-like PE (Nov 2026)

- **Why.** In summer files the small no-PE events near the chl trigger are a real Euk cluster that peaks 0.4–0.9 above the trigger, and the Pro-tail cut sliced off its lower part: about 3,500 old-class Euks in 12 TN368 files, and some in AR92 and HRS2609 (`measure_pro_tail.m`, `assess_pro_tail_removal.m`). In winter files the Euk-gap rule already handles most of the pile.
- **What the cut was also doing.** Outside EN727 most of the events it removed were faint Syn-like events detected by the PE trigger alone (chl below the chl trigger, PE:chl like Syn) that a Euk component had absorbed (AR95, HRS2609, AR34B, RB1904: 40–90 per file; `show_pro_tail_examples.m`).
- **Changes.**
  - The Pro-tail cut is removed (`edge_pro_size_above_syn` gone; `edge_pro_band` stays for the pile and the Syn rescue).
  - New rule: an event named Euk that has PE with PE:chl above the file's Syn PE:chl minus `syn_PEchl_tol` is unassigned, in the main rules and the catch-all (`info.n_euk_synpe`, `info.rules.euk_synpe_PEchl`).
  - Flag 3 now marks only the capped trigger-noise pile; `score_agreement` reports it as `flag_edge_pile`.
  - `score_agreement` now sets `cfg.cruise` and `cfg.gl1_is_od2_scatter` from `attune_config`, as the processing wrapper will.
- **Effect.** Median agreement on 96 non-Pro files (`rescore_euk_pe_rule.m`), before → now: TN368 0.886 → 0.909 (lowest 0.56 → 0.62), AR92 0.983 → 0.987, others within 0.002 (AR34B 0.922 → 0.921, EN727 0.930 → 0.928, AR95 0.971 → 0.970). No file is worse by more than 0.005. Pro-setting files: 19 of 24 identical, the rest within 4–33 events.
- **Left as Euk.** No-PE events at the chl trigger that the old classes left out (`scan_edge_euks.m`, `show_edge_euk_examples.m`): 57 of 6,781 Euks in EN727 25Jan(9), where the pile has no gap deep enough for the gap rule; 0–9 in every other file.

### Syn on particle needs chl above the chl trigger (Nov 2026)

- **Problem.** In EN617 20Jul2018A Sample(1) the chl trigger cuts the Syn cluster in half (53% of Syn have chl below it). 2,201 of 2,269 Syn on particle there had chl below the chl trigger: a band of faint large events detected by the PE trigger alone, two thirds of them also within 0.3 of the PE trigger. Their chl is not measured, so "less chl than one Syn plus a Euk" cannot be judged.
- **Change.** `synpart_rule` also requires chl above the chl trigger edge (`F.BL3 > edge.BL3`, where a chl trigger exists).
- **Effect** (`compare_synpart_edge.m`, `output/early_layout_robust/synpart_edge.mat`; 96 non-Pro files, 24 Pro-setting files, 6 EN617 files):
  - 33 of 126 files change; none gets worse; the largest gain is +0.0095 (HRS2609 SummerSyn(12) Euk(9)). Median agreement: TN368 0.909 → 0.910 (lowest 0.617 → 0.626), AR92 0.987 → 0.987; others unchanged.
  - Syn on particle, total before → after: TN368 940 → 496, HRS2609 5,118 → 4,559, AR92 846 → 738, AR34B 326 → 297, RB1904 2,310 → 2,242, EN617 4,420 → 2,133. Pro-setting files, EN627, EN727 and AR95 are unchanged.
  - EN617 20Jul A: 2,269 → 65; the band is now unassigned. The 65 left are its upper fringe just above the chl trigger (scatter 1–2 log10). EN617 25Jul A Sample(99) keeps all 1,825, which are well clear of both triggers, with Syn PE:chl and Euk-like FSC-W; they look genuine.

### Validation and review scripts

| script | what it does |
|---|---|
| `validate_early_layout.m` | agreement with the old classes for one early cruise, with the `attune_config` cruise settings; optional second classifier for before/after |
| `run_early_layout_validation.m` | runs the above for all early cruises; `'regress'` compares post-2019 files event by event |
| `summarize_early_layout.m` | one row per cruise from an output folder |
| `review_early_file.m` | review cytograms (`plot_class_cytograms`) and rule panels (`plot_v2_rules`) for one early-cruise file |
| `measure_edge_pile.m` | trigger-noise pile fraction per file |

Results: `output/early_layout/` (channel map only), `output/early_layout_robust/` (review figures in `review/`). There, `agreement_<CRUISE>_layout.mat` for TN368, AR34B, RB1904 and EN627 predate the Pro-tail removal; their scores with the final code are in `euk_pe_rule.mat`. No old class files exist for EN608, AR28B, AR29, EN617 and AR31A, so those have counts only, re-run with the final code (6 files each, no errors). EN617 Syn on particle is reviewed in the section above (`review/EN617_*`).

### Pro-like small Euks; event-level Pro pick tried and switched off (Nov 2026)

- **Pro-like events are no longer Euk.** No-PE events that are smaller than the Syn peak and dimmer than `pro_max_BL3` are now Pro-like, not Euk (`info.trace.euk_prolike`); the catch-all may not re-add them as Euk.
  - In files run without Pro settings they become class 0, and flag 4 (Pro-like) when at or above the chl edge. That feeds the Pro lower bound.
  - Small Euks brighter than Pro can be are kept, which protects Ostreococcus-like cells.
  - In AR104 wPro runs, 17–71 Euks per file move out.
  - 60-file Euk v2/old ranges: AR92 0.96–1.06 → 0.96–1.00, AR95 0.98–1.05 → 0.98–1.01.
- **Event-level Pro pick: tried, now OFF** (`pro_pick_min_n` = Inf). It applies only to Pro-setting files where the mixture names no Pro. It picks small no-PE, low-chl events with chl:scatter within 1.0 of Syn, if their size distribution has its own peak.
  - In AR104 wPro runs (run as Pro settings) there is no such peak. About 1,300 small faint events per file spread evenly into the trigger pile, mostly 1.1–1.5 below Syn in chl:scatter, so the pick correctly found nothing. Plot: `output/examples/ar104_wpro3_rules.png`.
  - In AR95 it invented Pro where the old classes have none: 677 in 19Oct AltPro(30), the faint-detritus case judged not Pro, and 493 in 18Oct AltPro(11).
- **AR104 Pro flag.** `FCSfileinfo` on the NAS has pro_measured = 0 for 2,848 of 2,856 files. It predates the voltage-based `attune_pro_measured` rule, under which the wPro runs (AND_SSC 500 at 240 V / AND_BL3 300 at 410 V) count as Pro-capable.
- **Concurrent change from another session.** Syn on particle now requires chl above the chl trigger edge (uncommitted edit by the pre-2019 / wrapper sessions, for EN617). That is why Syn on particle fell in OR_BL3-triggered summer files: AR92 1,094 → 950, HRS2609 2,025 → 1,817.

### Pro components judged against Syn: size and chl:scatter (Nov 2026)

- **Survey.** `survey_pro.m` classified 12 Pro-capable Pro-named files per cruise (voltage rule, `attune_pro_measured`) in AR79, AR82A, AE2426, AR92, AR95, HRS2609 and AR104, and listed every mixture component that is dim and has no PE (`output/pro_survey_base.mat` before, `output/pro_survey.mat` after). `comp_detail.m` lists all components of chosen files.
  - Real Pro components (holding the old-class Pro) sit **0.12–1.0 smaller than the Syn components** and **0.1–1.0 below Syn in chl:scatter**.
  - Summer Pro (AR92, HRS2609) are 0.25–0.7 smaller than Syn, with size −1.5 to −2.0. Fall and spring Pro (AR79, AR82A, AE2426, AR104 offshore) are often only 0.15–0.5 smaller, with size about −1.40.
- **The AR104 problem was the fixed size cut, not chl:scatter.** A component could be Pro only if its size was below `pro_max_SSC` (−1.44). Fall/spring Pro clusters at −1.38 to −1.43 were named Euk, then removed as Pro-like Euks, and the catch-all gave back only part of them.
  - AR104 23Oct wPro(1): the Pro cluster (half the chl events, 0.8–0.9 below Syn in chl:scatter, well within the 1.0 rule) gave 0–4,900 Pro per file.
  - AR82A 06Apr alternating_shelf_pro: 0 v2 Pro against 19,375 old.
- **Change** (`assign_class_v2`, after `name_components`):
  - `name_components` names any compact, dim, no-PE component a Pro candidate, whatever its size.
  - Candidates are then judged against the file's Syn components (weighted means):
    - size at least `pro_min_size_below_syn` (0.1) below Syn, else Euk as before;
    - size no more than `pro_max_size_below_syn` (1.1) below Syn;
    - chl:scatter below Syn's (`pro_max_chlSSC_above_syn` = 0);
    - chl:scatter within `pro_max_chlSSC_below_syn` (1.0) of Syn, unchanged.
  - Rejected candidates: `info.pro_comp_rejected`. Files without a Syn component keep the fixed `pro_max_SSC` window.
  - Catch-all Pro must be smaller than the Syn peak (was SSC < `pro_max_SSC`).
  - The two new upper limits remove tiny high-chl artifacts: 1.2–1.5 smaller than Syn, at or above Syn in chl:scatter. They were named Pro in AR104 23Oct wPro(9) Sample(1) (9,154 Pro) and gave 300–430 Pro per file in AR79.
- **Effect.**
  - **60-file check** (AR92, AR95, EN727, HRS2609): median and lowest agreement per cruise unchanged. 5 files change by at most 0.0005 and 23 Pro.
  - **Survey files.** Pro-setting files whose Pro the size cut had lost now match the old classes:

    | file | old | before | after |
    |---|---|---|---|
    | AR82A 06Apr alternating_shelf_pro | 19,375 | 0 | 20,035 |
    | AR82A 07Apr Pro_only(2) Sample(34) | 19,204 | 7,540 | 18,378 |
    | AR82A 09Apr alternating_shelf_pro(3) Pro(5) | 14,956 | 9,431 | 16,124 |
    | AR82A 08Apr alternating_shelf_pro(9) Pro(16) | 4,202 | 2,636 | 4,655 |
    | AE2426 08Nov alternating_shelf_pro(7) | 3,729 | 1,886 | 4,270 |

    AR79 10Nov Pro(14) goes 11,905 → 15,984, with old 0, as in the first file of the run. Summer files are unchanged.
  - **AR104 offshore** (3 files per Pro-run group, `output/ar104_pro_scan.mat`):
    - 23 Oct 05:00 to 24 Oct 00:00 is now a continuous 15,000–33,000 Pro per file. Before, about 10 files there fell to 0–7,600 (23Oct wPro(1) Sample(50): 0 → 19,106; 22Oct wPro(13) Sample(14): 5,965 → 16,256).
    - 20 Oct 01:00–06:00 is unchanged (about 500–5,800).
    - The shelf (22 Oct) and 24 Oct after 01:00 stay at 0.
  - **Cases not to count:** AR95 19Oct AltPro(30) Pro(10) stays at 0 Pro (the false group is 1.13 below Syn in chl:scatter, 0.20 above the chl edge). AR95 21Oct AltPro(24) keeps its real cluster (about 1,100).
- **Chl:scatter limit not loosened.** With the size fix, no AR104 offshore Pro component falls outside the 1.0 limit. The components that do are near-edge groups (0.1–0.3 above the chl trigger, 1.0–1.2 below Syn) with no old-class Pro, in AE2426, AR95 and HRS2609. A cruise-level Pro reference is not needed for these files.
- **Open (unchanged by this change; old classes have 0 Pro):**
  - AR95 18Oct AltPro(11): 350–550 Pro in about a third of its files, from a small compact group 0.6–0.8 smaller than Syn, 0.6–0.8 below in chl:scatter, 0.2–0.7 above the edge.
  - Similar small groups of 360–1,600 Pro in AR95 18Oct AltPro(7) and 19Oct AltPro(29), AE2426 08Nov (5)/(9) and 09Nov (3)/(12), AR92 18Aug ProAlt(15) Pro(6), and HRS2609 less_sensitive (6).
  - AR79 09Nov Pro(11): 553 old Pro, 0 in v2. They sit in a broad mixed component (size SD 1.0).

### Pro: chl pulse shape, robust compactness, chl:scatter limit by height above the trigger, size ceiling (Nov 2026)

- **Wide-pulse faint signals among Pro.** In sparse-Pro files the Pro region holds a second group that is not Pro:
  - AR95 21Oct AltPro(24): 24% of v2 Pro; the diffuse small signals you pointed out.
  - It is 150–550 events per file, also under rich Pro clusters (AR104, AR95 AltPro(13)).
- **What separates them.** GL4 (chl excited by the green laser), BL1, BL2 and FSC don't separate it from the compact Pro core: all sit at noise for these events.
  - The **BL3 pulse shape**, log10(BL3-A / BL3-H), does. The compact core sits at about 0, like Syn; the second group forms a separate band at +0.3 to +0.45.
  - Its SSC pulse shape is normal, so it is not two particles passing together.
  - It sits at the smallest sizes, near the SSC trigger.
  - Plots: `output/examples/pro_pulse_shape_channels.png`, `pro_pulse_shape_en727_ar95.png`.
  - Every cruise from EN608 (2018) on records -A, -H and -W for all channels.
- **Changes:**
  - `attune_features` adds `F.BL3_AH` = log10(BL3-A / BL3-H), read through the same channel map. It needs no bead correction: area and height scale alike with voltage.
  - Reference: `info.pulse_ref_AH`, the median `BL3_AH` of small PE-bearing events (mostly Syn; their 99th percentile is 0.06–0.15 above it).
  - Pro-setting files only:
    - Small, dim, no-PE events more than `pro_max_AH_above_ref` (0.2) above the reference are left out of the mixture fit, though they are still assigned with it (`info.trace.wide_pulse`, `info.n_fit_wide_excluded`).
    - Class-7 events, and catch-all Pro, are removed above an **adaptive cut** `info.pro_AH_max` = max(reference + 0.2, Pro median + `pro_AH_sd` (4) × the Pro's narrow-side spread) (`info.trace.pro_wide_pulse`, `info.n_pro_wide`). The narrow-side spread (median − 15.87th percentile) can't be stretched by a wide group above.
  - A fixed 0.2 cut sliced the tail of the broad, single-peaked A/H distribution of summer Pro near the chl trigger (AR92, HRS2609: 280–430 old-class Pro per file). With the adaptive cut that is 7–24. Tight fall Pro (spread about 0.05) keep a cut near 0.25, below their wide group. Plot: `output/examples/pro_pulse_shape_summer.png`.
- **Robust compactness.** The Pro compactness test (`pro_max_sd_BL3`, 0.4) now uses the robust SD (MAD) of chl among the narrow-pulse events assigned to each component (`info.comp.sd_BL3`). It falls back to the fitted Gaussian SD when a component has fewer than 20 such events. The fitted SD is inflated by a tail toward the chl trigger: AR95 19Oct AltPro(29), a compact cluster, had a fitted SD of 0.43.
- **Chl:scatter limit grows with height above the chl trigger edge.** It is `pro_max_chlSSC_below_syn` (1.0) for components within `pro_chlSSC_edge_start` (0.2) of the edge. Above that it rises by `pro_chlSSC_edge_slope` (0.5) per log unit, up to `pro_max_chlSSC_below_syn_far` (1.3) (`info.comp.pro_chl_lim`).
  - AR79 10Nov Pro(14): its main cluster is 1.03 below Syn but 0.55 above the edge. It was rejected at 1.0, giving 0 Pro; now 19,584.
  - Near-edge groups stay out: AR95 19Oct AltPro(30) (1.13 below, 0.20 above), and AE2426, AR95, HRS2609 and AR104 groups (1.03–1.16 below, 0.16–0.27 above).
- **Absolute Pro size ceiling** `pro_max_size_abs` (−1.15, bead-relative size, the same scale in every cruise). The largest real Pro cluster in the survey is at −1.25. On the AR104 24 Oct shelf, Syn sit at −0.41, so a broad dim group at −0.88 was "smaller than Syn" and got 1,266 Pro; now 0.
- **Effect** (`output/pro_survey.mat`, `output/ar104_pro_scan.mat`):
  - **60-file check**, before → after this section:

    | Cruise | Median agreement | Lowest file |
    |---|---|---|
    | AR92 | 0.9855 → 0.9852 | 0.969 → 0.969 |
    | AR95 | 0.9582 → 0.9552 | 0.864 → 0.881 |
    | EN727 | unchanged | unchanged |
    | HRS2609 | 0.9783 → 0.9777 | 0.963 → 0.963 |

    The AR95 median drop comes from wide-pulse events that the old classes counted as Pro (19Oct AltPro(14): 186 of them; 19Oct AltPro(9); 21Oct AltPro(5)). AltPro(24) rises because its wide group is gone (1,109 → 901 Pro).
  - **Survey of 84 Pro-setting files:**
    - The diffuse groups judged not to be clusters go to 0: AR95 18Oct AltPro(7) and AltPro(11), AE2426 09Nov (3) and (12).
    - AR95 19Oct AltPro(29) and 21Oct AltPro(24) keep their compact clusters (657 and 847).
    - Rich Pro files mostly lose 0–3%. Exceptions: AR82A 09Apr shelf_pro(3) 16,124 → 15,047 (old 14,956); AR82A 09Apr(8) 1,431 → 709 (old 1,162).
    - AR79 09Nov Pro(9) (old 1,186) and AR82A shelf_pro(8) (old 264) go to 0. Their Pro-like events are in a broad component at Syn size, or in no component.
  - **AR104:** the offshore run (23 Oct 05:00 – 24 Oct 00:00) loses 2–10% per file (wide-pulse events) and stays continuous at about 14,000–32,000 per file. The shelf stays at 0, except 24 Oct 01:09 wPro(4) Sample(1), with 1,176 Pro one hour after the run's last strong files.
  - **Preserved runs** (discrete wrapper, 12 cruises): median agreement Pro runs 0.903 → 0.908, CHL runs 0.914 → 0.923.
- **Event-level Pro pick still OFF.** With `pro_pick_min_n` = 100 it finds the small compact clump in AR95 19Oct AltPro(30) (433 Pro). But it also takes a sparse trigger-edge group in 18Oct AltPro(7) (222) and tiny artifacts in AR79 Pro(11) (71): it lacks the component rules (size limits, chl:scatter, height above the edge).
- **Open:**
  - EN727 preserved C19N13 (no real Pro, by your ruling, like C19N02, C09N01 and C06N06): v2 still names a Pro component on the chl trigger floor at the smallest sizes, about 395 events.
  - Pro-like flag 4 (the Pro lower bound in files without Pro settings) has no pulse-shape filter yet.

### Catch-all adds no Pro; bottom slice of a Pro cluster (Nov 2026)

- **Catch-all Pro removed.** Class 7 now comes only from components named Pro; unassigned events whose neighbours are Pro stay unassigned.
  - In 84 surveyed Pro-setting files the catch-all had added 8,153 Pro in the 77 with old classes (only 4% of them old-class Pro, against 91% for main-rule Pro) and 1,942 in AR104. Most were faint events 0.1–0.3 above the chl trigger, 0.2–0.6 smaller than the file's Pro.
  - In preserved EN655 C15N09, run through v2 alone, they were 1,337 of v2's 1,603 Pro, from the trigger-floor pile.
  - Losses: AR79 10Nov Pro(14) 19,584 → 13,260 (a Syn-sized cluster named Euk, which the catch-all had partly given back as Pro), AR92 18Aug ProAlt(15) 1,229 → 723 (old 0), AR104 total −2%.
- **Bottom slice of a Pro cluster.** Where the chl trigger cuts a Pro cluster, the mixture can give its lowest part a component of its own: the same size, lower chl, so lower chl:scatter. Near the edge such a slice fails the 1.0 limit.
  - HRS2609 18Aug(1) Pro(19): slice 1.04 below Syn, 0.15 above the edge, 743 of its 842 events old-class Pro. Only the catch-all had rescued it.
  - A candidate that fails only the chl:scatter limit, and is no more than `pro_max_chlSSC_below_syn_far` below Syn, is now Pro if an accepted Pro component is within `pro_slice_size_tol` (0.25) in size and brighter in chl (`info.pro_comp_slice`).
  - The near-edge false groups have no Pro cluster above them and stay out: AR95 19Oct AltPro(30), AE2426 09Nov(10).
- **Effect of both, 60-file check** (start of the Nov 2026 Pro work → now):

  | Cruise | Median agreement | Lowest file |
  |---|---|---|
  | AR92 | 0.9855 → 0.9852 | 0.969 → 0.969 |
  | AR95 | 0.9582 → 0.9573 | 0.864 → 0.882 |
  | EN727 | unchanged | unchanged |
  | HRS2609 | 0.9783 → 0.9779 | 0.963 → 0.963 |

  HRS2609 18Aug(1) Pro(19) goes 0.973 → 0.978. Slices also return in AR95 21Oct AltPro(3) (5,351 → 5,580, old 5,591) and HRS2609 17Aug(5). AR104 is unchanged by the slice rule. Tables: `output/pro_survey.mat`, `output/ar104_pro_scan.mat`.
- **Open: Pro clusters at Syn size are named Euk** (candidate size ≥ Syn − `pro_min_size_below_syn`, 0.1). Examples:
  - AR79 10Nov Pro(14): a second cluster 0.09 smaller than Syn, about 11,000 events.
  - Some AR104 23Oct wPro(1) groups, 0.04–0.17 smaller than Syn.
  - Preserved EN655 C15N09: the true Pro cluster, 0.02 smaller than Syn. The discrete wrapper swaps it back.

### Pro size-spread cap (Nov 2026)

- **Problem.** EN727 preserved Pro run C19N13 has no real Pro (user ruling, as for C19N02, C09N01 and C06N06), but v2 named two components Pro. Both are pressed against the SSC trigger wall (AND_SSC 900 in these preserved runs):
  - a compact corner group: size −2.12, 0.20 above the chl edge, 0.22 above the SSC edge;
  - a broad column above it: size −1.70, size SD 0.48.
- **Height above the SSC trigger doesn't separate them** (`ssc_wall` survey of the 84 underway Pro-setting files):
  - Real Pro components sit 0.24–0.90 above the SSC trigger edge, with up to 23% of their events within 0.15 of it.
  - The corner group sits 0.22 above it, with 34% of its events there.
- **Size spread does separate the column.** Real underway Pro components have a size SD of 0.10–0.38. The only two above 0.45 are AR92 18Aug ProAlt(15) (old classes 0 Pro) and a 383-event side component in AR79 10Nov Pro(6).
- **Change.** A Pro candidate whose size spread exceeds `pro_max_sd_size` (0.45) is rejected (label 0) and can't be restored as a bottom slice. The spread is the robust SD of the narrow-pulse events assigned to it, like `sd_BL3` (`info.comp.sd_size`).
- **Effect:**
  - 60-file check: unchanged.
  - Survey: AR92 18Aug ProAlt(15) 723 → 0; AR79 10Nov Pro(6) 35,013 → 34,630 (old 33,566).
  - AR104: two offshore files lose 4–6% (22Oct wPro(8) Sample(50), wPro(9) Sample(1)); total −0.4%.
  - EN727 C19N13, v2 alone: 468 → 198 Pro. The compact corner group remains.
  - EN727 C19N13, through the preserved wrapper (`assign_class_discrete`, 'pro'): 0 Pro, as for C19N02, C09N01 and C06N06. EN655 C15N09: 946 (hand 799).
- **Cost in AR92 preserved C21N05** (CHL run, 85.8 m; redesign-wrapper check):
  - Effect: Pro 3,962 → 3,525/ml (hand 4,131). AR92 per-Niskin Pro median 1.067 → 1.050 of hand, lowest tenth 0.943 → 0.900; other runs unchanged.
  - The rejected component is not the deep Pro cluster. It is a 133-event smear from size −1.80 to −0.43, with only 32% within ±0.25 of its peak; the wrapper's large-Pro step still finds the cluster.
  - No spread measure separates it from the C19N13 false column, which is just as lopsided. Compared: robust SD, half-widths either side of the peak, fraction near the peak.
  - Kept as is (user decision, Nov 2026).

### Syn-sized part of a Pro cluster (Nov 2026)

- **Problem.** AR79 10Nov Pro(14): the mixture split the Pro cluster so that its centre (about 11,000 events) sat in a component only 0.09 smaller than Syn. The rule "candidate size ≥ Syn − `pro_min_size_below_syn` → Euk" named it Euk. The Pro-like step then left 9,010 of its events unassigned, so only the cluster's outer ring was Pro.
  - The old classes have 0 Pro in this file. Their Pro gate failed here: the same cluster gives 18,131 old-class Pro in the run's previous file, 10Nov Pro(12).
  - Plot: `output/examples/ar79_10nov_pro14.png`.
- **Change.** A too-large Pro candidate is still Pro when:
  - it is no larger than the Syn components;
  - it passes the other Pro rules (`pro_max_size_abs`, the chl:scatter limits, `pro_max_sd_size`);
  - an accepted Pro component is within `pro_slice_size_tol` (0.25) in size.

  These are listed in `info.pro_comp_adjoin`. A lone Syn-sized group stays Euk.
- **Effect:**
  - AR79 10Nov Pro(14): Pro 13,260 → 24,428; Euk 3,157 → 673 (old 675).
  - 60-file check and every other survey file: unchanged.
  - AR104: two offshore files gain the Syn-sized parts of their clusters (19Oct testPro(10) Sample(2) 4,119 → 4,720; 22Oct wPro(12) Sample(50) 16,343 → 19,130).
  - Preserved runs through the wrapper: unchanged (EN727 no-Pro runs 0; EN655 C15N09 946). The wrapper's Pro swap leaves `pro_comp_adjoin` components alone.

### Event-level Pro pick removed (Nov 2026)

- The switched-off pick code and its four `pro_pick_*` settings are removed from `assign_class_v2.m` and `classifier_config.m`. With the pick off, the classes don't change: the six example agreements are identical.
- `info.pro_pick` is no longer returned. `discrete/assign_class_discrete.m` still sets `pro_pick_min_n` and checks `isfield(info, 'pro_pick')`. Both are harmless now, but its `pro_pick_trim` step can no longer run and could be removed by the discrete-sample session.

### Pulse shape for small Syn / Euk and noise: checked, no rule added (Nov 2026)

- **Wide chl pulses** (BL3 A/H more than 0.2 above the Syn reference) don't mark noise among small Syn or Euks in any cruise checked: AR92, HRS2609, EN727, AR95, AT46, AR99, EN627, AR34B, TN368. Only large Euks are wider, as expected. The wide-pulse criterion stays a Pro-only rule.
- **The side-scatter pulse (SSC A/H)** depends on instrument configuration: 20–35% of normal Syn look "wide" in pre-2019 cruises. It is not usable as a criterion.
- **Narrow chl pulses** (more than 0.2 below the reference) mark electronic trigger noise in pre-2019 files: 63–71% of the flagged trigger-noise pile in EN627 and AR34B, median −0.3. A few reach the dim end of Syn (AR34B: up to 8% of the dimmest tenth).
- **A narrow-pulse exclusion rule was tried and removed.** It cut 4–5% of real summer Pro (AR92 Pro(5) 32,305 → 30,964; HRS2609 Pro(21) 19,730 → 18,853), and dropping the pre-2019 noise from the fit destabilized the mixture (EN627 4Feb2019b Syn 5,636 → 2,050; 1Feb2019a Euk 429 → 1,463).
- **Euk loss in the held-out winter CHL runs** (AT46, AR99) is unrelated to pulse shape: the small Euks there have normal pulses.

### Syn reference from real Syn; Euk chl gap only at the trigger edge (Nov 2026)

Found in the held-out check of preserved winter CHL runs, where most hand-gated Euks were lost because the file's Syn reference was wrong. Three changes, all in `assign_class_v2.m` and set in `classifier_config.m`:

- **Syn reference from Syn-named events no larger than size −0.5** (`syn_ref_max_size`, with at least `syn_ref_min_n` = 20 of them).
  - A broad PE-bearing component of large particles can be named Syn. Its events outnumbered the real Syn, the reference moved to size +0.2 to +0.5, and "no Euk much smaller than Syn" removed the small Euks.
  - Real Syn peak at −1.5 to −0.75 in 328 underway files (winter EN627 the largest). Larger Syn-named components are mostly the large-Syn tail; naming them Euk instead changed 89 underway files and dropped Syn events, so the components keep their names and only the reference is restricted.
- **Syn reference PE:chl from Syn with a measured Syn-type spectrum** (GL3 clear of its floor, GL2:GL3 above `syn_min_GL2GL3`), when there are at least `syn_ref_spec_min_n` = 20.
  - In MVCO MV337N04 a cloud of dim PE particles (PE:chl about −0.7, GL3 at its floor so no spectrum test) was named Syn and outnumbered the 42 real Syn (PE:chl +0.35) ten to one.
  - GL2:GL3 alone doesn't separate the two (the cloud's GL2:GL3 is noise), nor does BL2:GL2 (the cloud's PE spectrum is Syn-like); PE amount and PE:chl do.
- **Euk chl gap only where the pile below it sits at the trigger edge**: its density peak within `euk_gap_pile_max_above_edge` = 0.7 of the chl edge.
  - In MVCO the gap window took in the large Euks as the "Euk peak", and the whole small-Euk cluster (peak about 1.0 above the edge) fell below the minimum between the two.
  - Underway trigger-noise piles peak 0.37–0.60 above the edge, so their gaps are kept.

**Effect** (hand Euk recovered, preserved runs through the discrete wrapper):

| run | before | after |
|---|---|---|
| EN627 C45N05 | 83 of 1,646 | 1,586 (Syn ref +0.28 → −1.21) |
| AR99 C02N16 | 321 of 1,158 | 1,056 (Syn ref +0.19 → −1.19) |
| MVCO MV337N04 | 989 of 1,843 | 1,791 (Syn 488 → 56, all 42 hand Syn) |

- AT46 C17N02 is unchanged; it has the opposite problem (too many Euks).
- **60-file check:** medians AR92 0.9852, AR95 0.9573 → 0.9575, EN727 0.9293, HRS2609 0.9779 → 0.9782; no file moves by more than 0.002.
- **Underway survey (328 files):** 10 change. AR88 26Apr(3) Sample(30) Euk 535 → 3,613 (old 3,577) and EN627 1Feb2019a Sample(1) Euk 429 → 1,444 (both v2 and the old classes had the reference at +0.47 there). Seven files lose 20–120 of 1,700–9,300 Syn; EN688 EXPORTS C01N24 loses 472 of 1,374 Syn, but it is a preserved PE run filed with the underway data and only changes when run through v2 alone. Through `assign_class_discrete`, the right path for it, Syn stay at 781 with the same Syn reference.

### Fourth Syn type; types relettered in order of BL2 − GL2 (26 Nov 2026)

Investigation in `output/syn_type_d` (scans of 1,000 underway files each from AR77 and AE2426, plus the 20-cruise reference set).

- **New mode.** A separate Syn mode sits at BL2 − GL2 ≈ −0.48, below the former type A (−0.23).
  - Its trough with the former A is at −0.39 to −0.41, and its spread (SD ~0.05) matches the other types.
  - Compared with the former A, it has PE:chl about 0.3 higher (absolute ~0.32), the same chl, the same emission ratios (GL2−GL3, GL4−GL2) and is +0.13 larger.
  - It is not saturation, the trigger edge, doublets, BL2 compression or an HV or laser-delay effect: it sits in the same files as the other types, at constant HV.
- **Where.** It is on the fresher northern shelf.
  - AE2426: a separate mode in 407 of 1,000 files, a median 4.5% of Syn in those files (max 8%).
  - AR77: 177 of 1,000 files, about 1.7%.
  - Also in reference files from EN657, EN712, EN715, EN727 and AR95.
- **Relettering.** Types are now lettered in order of BL2 − GL2: 1 = A (new, −0.48), 2 = B (former A, −0.23), 3 = C (former B, 0.00), 4 = D (former C, +0.45). 9 and 0 are unchanged. Edges: `syn_type_edges = [-0.40 -0.10 0.19]`.
- **Brightness-dependent type-A edge.** Let b = min(BL2, GL2) above floor (log10).
  - Where b ≥ `syn_typeA_min_above_floor` = 2, type A is BL2 − GL2 below the first edge (−0.40).
  - For dimmer Syn the edge moves left by `syn_typeA_dim_slope` = 0.4 per log10 of shortfall: below −0.50 at b = 1.75, below −0.60 at b = 1.5.
  - Events between that edge and the first edge stay B.
  - Warm-water Syn (T > 21) have a broad low tail below −0.40, 2–7% of all Syn, with no separate mode. It is noise in the ratio of dim Syn. The share below −0.40 falls with brightness, min(BL2, GL2) above floor: 39–64% of events at 1–1.25, 1–4% at 2–2.25, about 0 at 2.5 and above.
  - Summer Syn are mostly 1.5–2.25 above floor. Real type A keeps a flat share at any brightness above about 2 (AE2426 4–5%).
  - Share of events below −0.40 called A, and false A as a share of all Syn in summer files without a mode:

  | rule | AE2426 | AR77 | AR95 | AR92 inshore station | summer false A |
  |---|---|---|---|---|---|
  | none | 1 | 1 | 1 | 1 | 2.4–7% |
  | PE:chl > 0.05 (tried first) | 0.99 | 0.93 | 0.92 | 0.41 | 0.5–1.0% |
  | **min(BL2, GL2) ≥ floor + 2** | **1.00** | **0.99** | **1.00** | **0.74** | **0.2–1.1%** |
  | min(BL2, GL2) ≥ floor + 2.25 | 1.00 | 0.95 | 0.99 | 0.59 | 0–0.7% |

  - **Sloped edge vs a fixed limit.** The table above groups files by a separate type-A mode found with all events. That grouping misses type A where it is a minority or dim, so it overstated false A.
    - Regrouped by the share of bright Syn (b ≥ 2.25, where the ratio is precise) below −0.40: type A where that share is over 2%, none where it is under 0.5%, the rest left out.
    - Full-event histograms of ratio against brightness, 10 cruises (about 2,000 files):

    | rule | type A kept | false A (share of all Syn, files without type A) |
    |---|---|---|
    | fixed b ≥ 2 | 0.77 | 0.28% |
    | fixed b ≥ 1.75 | 0.87 | 0.84% |
    | fixed b ≥ 1.5 | 0.96 | 1.75% |
    | **sloped edge, 0.4** | **0.90** | **0.64%** |
    | sloped edge, 0.2 | 0.96 | 1.43% |
    | none | 1 | 2.9% |

    - Summer type A is dim and more common than the mode test showed: EN706, EN720 and the AR92 station keep only 25%, 75% and 70% of it with the fixed limit. AE2426, AR77 and AR95 keep over 99% under every rule.
    - The sloped edge (user's choice, 26 Nov) beats a fixed limit of 1.75 on both counts.
  - Why not PE:chl: type A's PE:chl is 0.3 above B in AE2426 and AR77, but only 0.25 above it in AR95 and 0–0.07 above it in summer type A. Preserved samples show no excess at all.
  - The rule applies to classes 2, 5 and 8 alike, since it concerns the Syn PE signal itself.
- **Summer type A.** It occurs at an inshore low-salinity station (41.32 N, 70.58 W): EN687 C26, EN720 C24 and AR92 C26 hold 7–19% of Syn in preserved runs and 9–21% in the underway files at those casts. AR92 shows a separated mode at −0.49. The 150-file summer scans missed these patches.
- **Preserved samples** (1,100 preserved CHL/PE runs from 13 cruises, classified as in process_cruise preserved mode): preservation shifts BL2 − GL2 up by about 0.06.

  | | type A peak | A/B trough | B peak | B/C trough | C peak | C/D trough |
  |---|---|---|---|---|---|---|
  | live | −0.48 | −0.40 | −0.23 | −0.10 | 0.00 | 0.19 |
  | preserved | ~−0.43 | −0.31 to −0.37 | −0.13 to −0.21 | −0.01 to −0.07 | +0.05 to +0.09 | ~0.25 |

  - `discrete_config` sets the edges 0.06 higher: [−0.34 −0.04 0.25].
  - Preserved Syn are about 0.3 brighter in PE, so the brightness rule keeps nearly all of them. Runs without a type-A mode put 0.1–0.8% of Syn below −0.34.
  - `assign_class_v2` returns `info.syn_type_edges`, `info.syn_typeA_min_above_floor` and `info.syn_typeA_dim_slope`, so the review draws the rules a run was classified with.
- **Old class files.** Files written with 3 types (`class_meta.cfg.sig.syn_type_edges` has 2 edges) are read through `syn_type_current`, which shifts their codes 1/2/3 to 2/3/4. They have no type A: its events are counted in their B.
  - `generate_attune_table` and `preserved_table` use this, and add `typeD` columns.
- **Other updates.** `score_agreement` (`new_synD`), `plot_class_cytograms`, the review frame, QC and viewer, and `merge_discrete_tables` now handle the fourth type. In the review, B/C/D keep the former A/B/C colours and A is green.
  - The review has an 11th panel, BL2 − GL2 vs PE:chl, with the edges as vertical lines. The GL2 vs BL2 panel draws the type-A boundary, including its sloped part for dim Syn.
- **Caveats.**
  - About 10% of summer type A (dim and close to the first edge) is still counted as B.
  - In preserved runs the edges assume one shift, 0.06, for all cruises (observed +0.04 to +0.10 per cruise).

### Syn types from a per-file mixture (26 Nov 2026)

The fixed windows (with the sloped type-A edge) are now only the fallback for files with fewer than 500 typed class-2 Syn. In other files `syn_type_mixture` decides which types the file has and types each event from a fit to that file. Prototype and figures: `output/syn_type_d` (`mix_proto4.m`, `mix_examples4.png`, `mix_summary.m`).

- **Why.**
  - The windows count the noisy tails of dim Syn as neighbouring types. In summer files 5–7% of Syn fall below −0.40 and were "type A", and B's spread fills the C window.
  - A shoulder or tail can't be told from a real type by a fixed edge. Mode detection on bright events missed minority types, and summer files have few bright Syn.
- **Not an offset.** Subtracting channel offsets (−60 to +120 raw counts) before the log never narrowed the dim spread in 10 files. Beads (2026 settings evaluation, BL2 330 vs 380 V) put a BL2 offset at about +35 counts, which is +0.05 at most in dim Syn, and in the wrong direction. The spread is measurement noise (`output/syn_type_d/offset_test.m`, `bead_offset.m`, `spill_test.m`).
- **The model** (`cfg.sig.syn_mix`), fitted to class-2 events binned by BL2 − GL2 (0.02) and brightness band (0.25; brightness = mean of BL2 and GL2 above their floors, which is independent of the ratio's noise):
  - types A–D are Gaussians along the ratio, at −0.48 / −0.23 / 0 / +0.45, each free to move ±0.10–0.12 per file; one position per file, so bright events place it;
  - total width = intrinsic sd (0.03–0.08 per type and file) plus band noise. The band noise is shared by all files and was fitted jointly over 1,858 files (total sd 0.21 at brightness 1.25–1.75, 0.13 at 2–2.25, 0.085 at 2.5–2.75, 0.065 above 3);
  - weights per band, pulled toward the file's overall weights (200 pseudo-events), because dim Syn can be another population than bright Syn;
  - a broad component, capped at 5% per band, for outliers;
  - backward elimination: drop types while that costs less than 30 · log n in log-likelihood. Types under 1%, and types under 10% pinned at their limit next to a kept type (a shoulder), go first.
- **Codes.** Each event takes the most probable kept type if its posterior is at least 0.8; otherwise it gets the new code **5, type ambiguous**. Classes 5 and 8 take posteriors from the class-2 fit.
  - `info.syn_mix` / `class_meta.syn_mix` hold the method, kept types, positions, widths, elimination path, and soft counts (summed posteriors, classes 2/5/8 × A–D, all events).
  - `generate_attune_table` adds `*_type_ambiguous_count` and `Syn_typeA..D_soft` (class 2, NaN where the windows were used); `preserved_table` and `merge_discrete_tables` add `syn_type_ambiguous`.
- **What changed (pen 30; scan subset, 10 cruises).**
  - Type A kept: AE2426 194 of 301 files (share 2.7%, as with the windows), AR77 85/301, AR92 117/245 (the inshore summer type A), EN706 49, EN720 22; summer EN655 6, EN668 3, EN685 1, EN687 5. The summer type-A share fell from 5–7% to 0.1–0.6%.
  - Ambiguous: 3–6% of Syn in cruises with bright Syn, 22–34% in dim summer cruises, 1–2% in preserved runs.
  - The AR92 example 15Aug2025(2) S1 keeps B, C and D: the hump at 0 (8%, sd 0.08) costs 45 · log n to drop.
  - Classes are unchanged, and the fit adds no measurable time.
- **Preserved runs:** `discrete_config` shifts the positions by +0.06, like the edges; the band noise is the live one.
- **Review.** In mixture files the type panels draw the kept types' positions in their colours instead of fixed edges; ambiguous events are drawn with the untyped ones. The side panel shows the method and kept types, and there is an "ambiguous-type Syn share" metric.
- **Caveats.**
  - Gaussian shapes: a broad dim bulk between types (EN687) is shared out among B, C and D fairly arbitrarily, though the soft counts stay consistent.
  - Whether a type is kept depends on the penalty for files near the threshold.
  - Class files written before this change have no code 5 and no `syn_mix`, and keep their window types until classified again with `process_cruise`.

### Trigger pile vs cells; Pro in runs without Pro settings; class 10 (26 Nov 2026)

Found while following an AR77 low-chl Euk branch (shared notes: `ar77_dim_euk.md`) and the AE2426 alternating Euk/Pro runs.

- **Telling the trigger-noise pile from cells** (`edge_pop_is_pile`). In about 500 underway files from every cruise:
  - the pre-2019 pile has chl pulses 0.24–0.33 (log10 A/H) shorter than the pulse reference, a chl spread (robust SD) of 0.02–0.03, 90–95% of events within 0.1 of the chl edge, and no size–chl relation;
  - every cell population near the edge (Syn, Pro, small Euks, the AR77 branch, flag 4 and 5 groups after 2019) has pulses like the reference, a chl spread of 0.07–0.4, and chl rising with size.
  - A population is the pile when its median pulse is more than `pile_AH_below_ref` (0.15) below the reference. Without a pulse-width column it is the pile when its chl spread is below 0.06, more than 70% of it lies within 0.1 of the edge, and |rho(size, chl)| is below 0.2.
- **Pro in runs without Pro settings** (`pro_in_euk_runs`).
  - The chl trigger cuts the Pro cluster, and the part above it was left unassigned. In AE2426 08Nov Euk runs it sat below the Euk chl gap (flag 2); about 30% of the neighbouring Pro runs' Pro lie above the Euk run's edge.
  - Unassigned no-PE events at or above the chl edge become class 7 when they are cells. Their size runs from Syn − 1.1 to Syn + 0.25 and their chl:size is below 1.4. If they are the pile they stay unassigned, with flag 3. `info.pro_lower_bound` is set and `class_meta.pro_lower_bound` saves it.
  - Euks are not taken. In summer cruises the small Euks are one cluster reaching Syn size, and fixed limits only sliced it: a first version took 10–45% of the Euks in TN368, EN655, EN668, EN685 and EN687.
  - **AE2426 8–10 Nov, 460 Euk runs:** Euk-run Pro divided by the neighbouring Pro runs' Pro above the Euk-run edge has a median of 0.99. Euk + Pro divided by the Pro runs' Euk + Pro is 0.96–1.03 (10–90%), and Euk is unchanged.
  - **Underway survey (508 files):** Pro appears in 207 of 444 files without Pro settings, about 193,000 events, all previously unassigned. AR39B (94,600) and EN657 (23,600) gain most. EN657's old classes called many of these Euk.
  - Pre-2019 files where the pile dominates (`edge_pile_capped`) get no Pro.
- **Class 10, low-chl Euk-like cells** (`dim_euk_class`, identity uncertain).
  - **The AR77 branch (11–13 Oct):** a compact no-PE group 0.35–0.6 above Syn in size, at chl:size about 0.6 against about 1.4 for bright Euks, with narrow chl pulses like Syn's. It is 40–80% of the Euk-like events while present, and almost absent from the matching preserved sample.
  - **The flipping it fixes:** the Euk chl:size cut's pile gate called the branch Euk in one file and unassigned (flag 5) in the next.
  - **The rule:** fixed limits relative to Syn, so the count follows the cells rather than what else is in the file. No PE, size Syn + 0.25 to Syn + 1.0, chl:size between the detritus cut and `dim_euk_max_chlSSC` (0.8), at least 200 of them, and cells rather than the pile.
  - **Alternatives tried:**
    - a limit of 1.0 also took 10–30% of summer Euks, whose small-Euk cluster peaks at chl:size 1.0–1.4;
    - a test for a separate low mode made the AR77 count flip with the branch's shape (33 jumps over 3× in 65 files).
  - **AR77, 65 files:** Euk jumps over 3× fall from 8 to 4; class 10 jumps 11 times, against 6 for the same events counted regardless of class. Summer cruises lose 2–8% of Euks to class 10.
- **Downstream:** class 10 is in `generate_attune_table` (Euk_like_low_chl count, biovolume and carbon; not in the Euk columns), `process_cruise` (codes and the class log), the review config, QC and viewer, and `plot_class_cytograms`. Pro_lower_bound_count in files without Pro settings is class 7 + flag 4.
- **Preserved runs:** `discrete_config` turns both switches off. The discrete wrapper assigns Pro in CHL runs itself, and class 10 is not yet tested against the hand gates.
- **Class files written before this change** (e.g. AR77, classified 25 Nov) are out of step with v2 until classified again with `process_cruise`.

### Euk with Syn PE floor tied to the main Syn cluster (28 Nov 2026)

Found in AR43 11Mar2020a Sample(112), where 99 of 244 Euk with Syn sat below the Syn median PE.

- **The problem:** the class 5 floor (`euksyn_min_syn_PE_prc`, the 1st percentile of the file's Syn GL2) is meant as "at least one Syn's PE". In that file 18% of the Syn form a dim group, about 2 log below the main Syn in both PE and chl. The 1st percentile (GL2 −1.77) fell inside it, against a main Syn mode of 0.39, so Euk with Syn down to 2 log below a typical Syn passed. Their spectrum is Syn-like (GL2:GL3 0.39 against Syn 0.42 and Euk −0.55).
- **The rule** (`euksyn_min_core_sd` = 3): the floor is also at least the main Syn mode minus 3 robust SDs. The mode is the peak of a kernel density of Syn GL2 (bandwidth 0.1); the SD is from the Syn within 0.5 of it. The higher of the two floors applies, so files with one Syn cluster are unchanged. `info.euksyn_syn_mode` and `info.euksyn_syn_sd` record them. AR43 Sample(112): floor −0.17, 31 events move to Euk.
- **Held-out check** (`euksyn_floor_test.m`, applied to final classes):

  | | class 5 moved | hand Euk | hand PE euk | hand none |
  |---|---|---|---|---|
  | preserved, all class 5 (8,476) | | 41% | 10% | 47% |
  | preserved, 3 SDs | 445 (5.3%) in 22 runs | 54% | 2% | 24% |
  | preserved, 2 SDs | 899 (10.6%) in 84 runs | 47% | 2.8% | 29% |

  MVCO 2019 MV418N02 and MV418N01 hold over half of the 3-SD total. Live (240 files, 36,674 class 5): 3 SDs moves 112 in 15 files (mostly AR43 11Mar2020a and HRS2609); the old classes called 41% of them Euk, against 2% of all class 5. 2 SDs moves 800 in 133 files, two thirds of them crypto in the old classes, so 3 SDs was chosen.
- Preserved runs use the same setting (`discrete_config` starts from `classifier_config`). `classifier_version` tags for v2 and discrete are now 2026-11-28a.
