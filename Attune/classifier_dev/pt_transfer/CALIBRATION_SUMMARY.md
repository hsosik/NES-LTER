# Attune side scatter → cell volume: calibration summary (Nov 2026)

How volumes are now computed for every NES-LTER Attune cruise, why, and how sure we are.
Details and intermediate checks are in [README.md](README.md). Figures are in
`classifier_dev/output/pt_transfer/`.

## 1. Configurations at a glance

| Period | Cruises | Scatter channels | Bead reference | Volume calibration |
|---|---|---|---|---|
| Jan 2018 – Jan 2019, Jul 2019 | EN608, AR28B, AR29, EN617, AR31A, TN368 | SSC 488 nm with OD2 (one channel) | TN368: FCB beads. Others: PT performance tests transferred to FCB units | 488 nm, linear: `log10 V = 1.4066·log10(SSC/bead) + 1.2109` |
| Feb – May 2019 | EN627, AR34B, RB1904 | SSC 488 nm without OD2 (160–200 V) | AR34B, RB1904: FCB. EN627: PT conversion | same 488 nm line, **after decompressing** large-particle SSC. EN627 SSC flagged (degraded optics) |
| Aug 2019 → | EN644 … AR104 | SSC 488 nm (small cells) + GL1 532 nm with OD2 (large cells), merged in bead units | FCB beads | `log10 V = 1.1962·log10(merged/bead) + 1.0826` |

Everything is in **1 µm FCB bead units**: each channel is divided by that day's 1 µm bead peak on
the same channel, at the sample's voltage (signal ∝ hv^8.5).

## 2. Instrument effects that shaped the calibration

**Low-signal excess with OD2.** Below ~100 raw counts, scatter through the OD2 filter reads
too high, and the excess grows as the signal falls:

| Raw signal (counts) | Excess (log10) |
|---|---|
| 10^2.0 | +0.04 to +0.09 |
| 10^1.8 | +0.15 to +0.22 |
| 10^1.6 | up to +0.95 |

Seen on GL1 in 15 field cruises, and directly on SSC in the 2019 Syn culture. Syn sits right at
this knee. In the current layout Syn is therefore sized on SSC, and the merge to GL1 starts only
at 0.5× bead.

![GL1 vs SSC](../output/pt_transfer/gl1_lowsignal.png)

**Compression of SSC without OD2.** Large particles read too low on SSC-H without OD2, and more
so at lower voltage. From the same cultures run with and without OD2 (2019), in bead units:

- Linear below ~3× bead, rolling off to almost flat above ~12× bead.
- A smooth fit has rmse 0.057 log.
- The knee moves with SSC voltage (from the bead voltage sweeps): ≈ 1.00 log at 160 V,
  1.09 at 200 V, 1.085 at 220 V.
- A 6 µm bead at 220 V (2020) falls on the curve within 0.04.

![compression](../output/pt_transfer/ssc_compression.png)

**Consequence for the merge window.** The SSC→GL1 switch is now set in bead units
(**0.5–2× the 1 µm bead on GL1**), not raw counts. Both channels are linear there at the old
settings (SSC/GL1 220/220 V) and at the new ones (AR104: 240/180 V). The old 500–1000 count
window was 0.45–0.9× bead at 220 V but 2.6–5.2× bead on AR104, where SSC is already compressed.

## 3. Current-layout calibration (Aug 2019 → )

- **Cultures:** 12 cultures from July/Aug 2026 (one row per culture, as in
  `size_calibration_fit_jul2026.m`; T. weissflogii kept, because field samples contain diatoms).
- **Attune modes:** taken in log space (KDE peak of log10 signal). The earlier linear-bin
  histogram modes read 0.03–0.19 log low for broad populations.
- **Coulter modes:** per log diameter (`compile_cc_histograms.m` now uses count × diameter).
- **Result:** `log10 V = 1.1962·log10(merged) + 1.0826`; 95% CI slope 1.01–1.38, intercept
  0.88–1.29; rmse 0.26 log (r² 0.95). Now in `Attune_SC2vol('SSC-H')`.
- **Change from the 11 Sep 2026 fit** (1.1696x + 1.1592): about 0.08 log (17%) lower volume for
  small cells, nearly unchanged for large ones.

![2026 calibration](../output/pt_transfer/vol_cal_SSC-H_merge_14Nov2026.png)

Modes, linear-bin vs log space:

![attune modes](../output/pt_transfer/modes_attune_2026.png)
![coulter modes](../output/pt_transfer/modes_coulter_2026.png)

## 4. Pre-Aug 2019 calibration

- **One 488 nm calibration.** The 12 cultures of March/April 2019 run with OD2 on SSC (linear)
  give `log10 V = 1.4066·log10(SSC/bead) + 1.2109`; 95% CI slope 1.23–1.59; rmse 0.23 log. All
  modes in log space.
  - Where SSC without OD2 is still linear, bead-normalized SSC with and without OD2 agree within
    0.02 log. So this one line applies to both groups of cruises.
- **No-OD2 cruises** (AR34B, RB1904 at 160 V; EN627 at 200 V): SSC is decompressed first
  (`ssc_decompress`), then the same line is applied.
  - Beyond the invertible range (~7× bead at 160 V, ≈ 460 µm³; ~9× at 200 V, ≈ 620 µm³), the
    volume is set to that limit (a lower bound) and `vol_flag` is true.
  - `generate_attune_table` counts flagged cells as > 10 µm (out of the ≤ 10 µm bins, still in
    ≤ 20 µm); their biovolume stays at the lower bound, so large-cell biovolume is underestimated.
- **Bead references for EN608, AR28B, AR29, EN617, AR31A** come from the Attune performance tests.
  - Validated against paired FCB runs: residual RSD 0.005–0.012 log.
  - Installed as `beadstat.mat`; the old PT conversions for EN617 and AR31A are kept as
    `beadstat_PTconversion.mat`.
  - With these, Syn SSC on all OD2-on-SSC cruises agrees with TN368 (real FCB beads).

![2019 calibrations](../output/pt_transfer/final_cal_2019.png)
![bead transfer validation](../output/pt_transfer/validation_residuals.png)

## 5. Why the 488 nm and 532 nm calibrations differ

For large cells the 2019 (488 nm) line is steeper than the 2026 merged line (large cells end at
532 nm on GL1). OD2 does not cause this: it is neutral at both wavelengths. A Mie model
(`classifier_dev/mie_scatter/`) shows the cause is **pigment absorption**:

- **Bead ratios:** with a ±30–50° collection cone the model reproduces the measured 0.5/1 µm ratio.
- **Wavelength alone:** for non-absorbing cells, 488 and 532 nm agree within 0.03 log in bead units.
- **With realistic absorption** (k(440) ≈ 0.01), green and carotenoid-rich cells scatter
  relatively more at 532 nm, with the difference growing with size:
  - +0.1–0.15 at 50–60 µm³ (observed on the same 2026 cells: +0.09–0.15)
  - +0.3–0.45 at 1000–3000 µm³ (the gap between the calibrations: +0.3–0.4)
- **Cryptophytes** (phycoerythrin absorbs at 532 nm) go the other way.

Cells in the merge window (≈ 5–30 µm³) are too small for this to matter (model +0.03–0.07; field
+0.03), so the SSC/GL1 merge is consistent. Composition matters most for large cells, and differs
between the two wavelengths.

![Mie](../output/pt_transfer/mie_488_532.png)

The same pigment physics shows up more cleanly in the **chlorophyll excitation ratio
GL4/BL3** (ex 532/488). In the 2026 cultures it orders the groups as follows, largely
independent of size:

| Group | log10(GL4/BL3) |
|---|---|
| Green | −0.9 to −1.1 |
| Haptophytes | −0.8 |
| Diatom and cryptophytes | −0.6 to −0.75 |
| Dinoflagellate (Heterocapsa) and coccolithophore (SE62) | −0.4 |
| Cyanobacteria | −0.1 to −0.3 |

![pigment ratios](../output/pt_transfer/culture_pigment_ratios.png)

**Field eukaryotes** (15 current-layout cruises, 8 files each, ~490,000 eukaryote events;
`field_pigment_map.m`, rerun 15 Nov 2026 with the current classifier). This uses the excitation ratio in bead units (GL4 − BL3, which removes
drift in the two lasers' relative power), so the values are offset from the raw ratios above.
On this scale the culture groups sit at: green −0.42 to −0.58, haptophytes −0.27 to −0.34,
diatom and cryptophytes −0.07 to −0.12, Heterocapsa and SE62 +0.1, Syn +0.16 to +0.43.

**2019 cultures on the same scale** (squares in the figure). These are the March/April 2019 runs
(old layout, OD2 on SSC, green chl 680/40 on file channel GL3-H). Their volume comes from the 488 nm
calibration. They land where the 2026 cultures do, after a filter-layout change and seven years:

| Group | 2019 | 2026 |
|---|---|---|
| Green (Micromonas, Nano, Dunaliella) | −0.45 to −0.55 | −0.42 to −0.58 |
| Haptophytes (Pavlova, I. galbana, Prymnesium) | −0.17 to −0.27 | −0.27 to −0.34 |
| Cryptophyte | −0.08 | −0.07 to −0.12 |
| Syn | +0.36 | +0.16 to +0.43 |

- **Dinoflagellates differ:** the 2019 dinoflagellates span −0.38 (*Alexandrium minutum*,
  AMBOPO, ~3800 µm³) to −0.18 (*H. triquetra*), and the raphidophyte (*H. akashiwo*, −0.17) sits
  with the haptophytes, while the 2026 *Heterocapsa* sits at +0.12. Dinoflagellates are not one
  pigment group on this ratio (range about 0.5 log).
- **SE62 is a coccolithophore**, so a haptophyte, yet it sits at +0.09, well above the other
  haptophytes. On this ratio alone it would pass for a dinoflagellate.
- **Pptp02** (identity not confirmed here; ~4200 µm³) sits at −0.39.

- Small eukaryotes (≤ 3 µm³) sit at about −0.4 to −0.5, the green (prasinophyte) values, in
  both seasons.
- The ratio rises with size towards the haptophyte and diatom values (about −0.3 to −0.1 above
  100 µm³). That is consistent with the shift in community composition expected with size.
- Class 1 (eukaryotes without PE) has an upper edge near −0.12, and cryptophytes (classes 3/4)
  fill the space above it. Run on the cultures, the classifier puts Heterocapsa and SE62 99–100%
  in class 1, T. weissflogii 87–90% (6–7% in class 4), and the real cryptophytes 98–100% in
  class 4. So the edge is not diatoms or dinoflagellates being pushed into the cryptophyte
  classes, at least for these cultures.
- **The summer "stripes" of low-ratio cryptophytes in the first version of this map were an
  artifact.** They came from EN644 alone (17,246 class 4 events at a green-like −0.83). Its
  `FCSfileinfo.mat` lists 19 `Dockwater_filter_config_change` files from 6–15 Aug 2019, before the
  green-laser layout change of 16 Aug (4 with QC_flag 1).
  - The classifier code of 14 Nov read them with the current channel map. It took file GL2-H (then
    620/15, which picks up chlorophyll from green cells) as PE, and GL1-H (then PE) as scatter. That
    turned green cells into class 4 and gave the fixed-volume stripes.
  - The map now skips those files. EN644 then has 167 class 4 events at +0.20, like the other
    cruises, and nothing else changed.
  - The EN644 products never included these files: the old wrapper excludes `Dockwater`. The
    redesigned wrapper has been asked to keep that exclusion.

![field pigment map](../output/pt_transfer/field_pigment_map.png)

## 6. Uncertainty

- **Composite curve:** constrained to roughly ±25–60% across most of the size range (95% CI of the
  fitted line). A single culture, i.e. a given taxon, scatters ~×4 either way (rmse 0.23–0.26 log).
  Diatoms lie above the line (more volume per unit scatter), and carotenoid-rich cells lie below
  it, most of all at 488 nm.
- **Coulter modes:** the noisiest small cultures (Micromonas, Syn 8113) and flat-topped
  distributions (Pptp02, 2019) are uncertain by ~0.05–0.2 log in volume.
- **Bead references:** ≤ 0.02 log within a configuration. EN608 and AR28B used BL3 and green-chl
  695/40 filters; their fluorescence (not scatter) references rest on an assumption.
- **Not correctable:** EN627 scatter (optics degraded in Feb 2019).

## 7. What changed in the code

| File | Change |
|---|---|
| `Attune_SC2vol.m` | SSC-H: 1.1962 / 1.0826 (11 Sep values in comments). New case `SSC-H_488_2019`: 1.4066 / 1.2109. SSC-A not refit |
| `batch_mergedSSC2volume.m` | Merge window in bead units (`cfg.gl1_merge_bead`). For `p.OD2setting` ≠ `'GL1'`: no merge, `bead_reference` normalization, decompression when 'None', 488 nm calibration, `vol_flag`, `ssc_linear_bdnorm`, notes |
| `generate_attune_table.m` | Flagged cells (`vol_flag`) are counted as > 10 µm in the size bins (the lower bound is ≈ 9.6 µm ESD at 160 V); their biovolume and carbon keep the lower-bound volume |
| `ssc_decompress.m` (new) | Inverts SSC compression without OD2; voltage-dependent knee; flag and lower bound |
| `classifier_dev/classifier_config.m`, `assign_class_v2.m` | `cfg.gl1_merge_bead = [0.5 2]` replaces raw-count `cfg.gl1_merge` |
| `Calibration Curve/compile_cc_histograms.m` | Coulter modes per log diameter (`CCSummary_saved_14Nov2026.mat`) |
| Cruise data | `beadstat.mat` installed for EN608, AR28B, AR29, EN617, AR31A; `FCSfileinfo.mat` and `beadstat.mat` made for AR104 |

**Not done yet:**
- Class-file volumes have not been regenerated for any cruise (rerun wrapper steps 5–6).
- `compile_attune_histograms.m` still makes linear-bin modes.
- SSC-A coefficients are not refit.
- Possible 4th bead size for the Mie collection-cone fit (not used yet): a 9 µm bead run 11 min
  after an FCB mix at the same settings (21 Jan 2020, current layout, OD2 on GL1, SSC 220; also
  SSC 180 and two runs after a filter swap). `\sosiknas1\lab_data\Attune\ATX_FSC_misc_experiments\
  bead_hv_experiment2\filter_swap_experiment\` ("filter experiment neutral.fcs" then "neutral 9um
  beads.fcs"). Bead material/refractive index unknown; check before using.
- Don't rerun wrapper step 4 on EN608, AR28B, AR29, EN617 or AR31A: it would overwrite the
  installed bead tables. Rebuild them with `classifier_dev/pt_transfer/make_perf_beadstat.m`.
