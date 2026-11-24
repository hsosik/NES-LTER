# PT performance tests as a bead transfer standard (Nov 2026)

Goal: better FCB-equivalent bead references for the early cruises without FCB bead runs
(EN608, AR28B, AR29, EN617, AR31A, EN627), so that bead-normalized SSC (for the culture
SSC-to-volume calibration) and the fluorescence channels are consistent with later cruises.

## Data
- `\\sosiknas1\Lab_data\Attune\performance_test_fcs`: 230 Attune performance tests,
  25 Jan 2018 – 28 Aug 2019 (215 pass QC). In each test the Attune sets every detector voltage
  so that the PT bead reaches a fixed target (SSC-H ≈ 10^5.48), so the **voltage records the
  channel's sensitivity**. The date is in `#PTDATE` (these files have no `$DATE`).
- FCB bead runs (`beadstat.mat`) on AR34B and RB1904 (no OD2), TN368 (OD2 on SSC) and EN644
  (new layout). Nearly all are within 24 h of a performance test.

## Method
`log10 FCB peak = PT perf peak + k·log10(hv_run / hv_perf) + offset(channel, bead, optical state)`

Optical state = filter layout + OD2 placement. The offset is the median over FCB runs paired with
the nearest performance test in the same state. For the early cruises, each performance test
during the cruise becomes one row of a beadstat-format table, at the test's own voltages.
`bead_reference` scales those rows to each sample's voltages, as for any cruise.

| Script | Purpose |
|---|---|
| `pt_perf_table.m` | extract the performance tests → `output/pt_perf.mat` |
| `validate_perf_transfer.m` | pair FCB runs with performance tests, fit offsets, residuals → `output/perf_transfer_validation.mat` |
| `make_perf_beadstat.m` | build + install one cruise's PT-transfer `beadstat.mat` (wrapper entry point; backup, no overwrite of non-transfer tables) |
| `perf_beadstats.m` | beadstat tables for the early cruises → `output/pt_transfer/<cruise>_beadstat.mat` (**not** written to `bead_calibrated`) |
| `compare_perf_refs.m` | compare new and current references, plus a Syn SSC check → `output/pt_transfer/compare_perf_refs.mat` |
| `plot_perf_transfer.m` | figures in `output/pt_transfer/` |

## Results

### 1. The performance tests transfer bead calibration very well
![validation](../output/pt_transfer/validation_residuals.png)

The residuals (observed FCB peak − prediction) are as tight as repeat bead runs:

| Bead | Residual RSD (log10) |
|---|---|
| SSC (1 µm) | 0.006 |
| FSC (1 µm) | 0.012 |
| BL1 (1 µm) | 0.006 |
| BL2 (0.5 µm) | 0.007 |
| BL3 (1 µm) | 0.008 |
| GL2 (0.5 µm) | 0.005 |
| GL3 (6 µm) | 0.010 |
| GL4 (0.5 µm) | 0.008 |

The two no-OD2 cruises (AR34B in April, RB1904 in May) give offsets within ≤ 0.02 of each other.

The offsets differ between optical states. For SSC, FCB − PT perf is:
- 1 µm: −0.81 without OD2, −0.94 with OD2 on SSC
- 6 µm: +0.20 without OD2, +0.29 with OD2 on SSC

So **OD2 changes the FCB/PT ratio on SSC** (the filter-swap experiment found the same), and the
current PT conversion's assumption that "OD2 does not change bead ratios" is wrong by about
0.1–0.14 log for EN617 and AR31A.

### 2. Early cruises: bead-normalized Syn SSC
![syn](../output/pt_transfer/syn_ssc_early_cruises.png)

Median normalized Syn SSC (log10, Syn − 1 µm bead):

| Cruise | OD2 | Performance tests used | Perf transfer | Current PT conversion |
|---|---|---|---|---|
| EN608 (Feb 2018) | SSC | 12 (filters 695/40) | −1.06 | none |
| AR28B (Apr 2018) | SSC | 8 (695/40) | −1.05 | none |
| AR29 (Apr 2018) | SSC | 16 | −1.08 | none |
| EN617 (Jul 2018) | SSC | 7 | −1.13 | −1.23 |
| AR31A (Oct 2018) | SSC | 3 | −1.07 | −1.21 |
| TN368 (Jul 2019, **real FCB beads**) | SSC | – | −1.04 | – |
| EN627 (Feb 2019) | none | 9 | −0.95 | −0.85 |
| FCB cruises without OD2 | none | – | Feb ≈ −1.18, Jul ≈ −1.31, Oct ≈ −1.27 | – |

- **All OD2-on-SSC cruises now agree with each other and with TN368**, the one OD2 cruise that
  has real FCB beads (−1.04 to −1.13). The current PT conversion only *looked* typical of the
  season for EN617 and AR31A: its OD2 error (~0.1–0.14) happened to cancel the real offset
  between the OD2 and no-OD2 states.
- **SSC with OD2 reads Syn about 0.15–0.25 log higher relative to the 1 µm bead than SSC
  without OD2.** TN368 shows this with real beads, so it is a measurement property of the OD2
  state (Syn sits near the SSC floor under OD2, and the response is size dependent), not a bead
  calibration error. Bead normalization alone can't remove it. It needs an explicit
  OD2-to-neutral SSC correction (next step below).
- **EN627 is still ~0.23 log high** (it was ~0.36). The performance tests confirm that SSC and
  FSC sensitivity were abnormally low and noisy through Feb 2019 (SSC needed 259–280 V vs ~227 V
  from 28 Feb). Fluorescence was barely affected, which points to the flow cell or scatter path.
  A single-size transfer can't correct the part that depends on particle size. EN627 SSC should
  be flagged.

### 3. Filter configuration found in the performance tests
BL3 and the green chlorophyll channel need ~20 V more with 695/40 than with 680/40. The tests
show:
- **695/40:** 25 Jan – 1 Mar 2018 (EN608), 1–12 Apr 2018 (AR28B), early Jul 2018
- **680/40:** 6–27 Mar 2018, from 16 Apr 2018 (**AR29**, although the settings history lists
  695/40 with a "?"), and from 16 Jul 2018 on

For EN608 and AR28B, the BL3 and GL4 offsets (measured with 680/40) assume the FCB/PT ratio is
the same through 695/40. This is untested and is flagged in `meta.notes` and the `filterset`
column.

## Caveats
- k = 8.5 is used to cover the gap between the voltages of the performance tests (SSC ~380–400 V
  with OD2, ~227 V without) and of the bead runs. Within a state that error is absorbed by the
  offset, and the 2018 test voltages are close to TN368's (381 vs 401 V on SSC).
- Within the OD2 state the offsets come from one cruise (TN368), so how well they hold over time
  is tested only in the no-OD2 state (AR34B → RB1904, ≤ 0.02).
- FSC references move by −0.11 to −0.22 against the current PT conversion (FSC test voltages
  jump around more). Not investigated further.

## Recommendations
1. **Use these tables for the OD2-on-SSC cruises.** Replace the PT conversion for EN617 and
   AR31A, and add EN608, AR28B and AR29, which have no bead reference now. This means copying
   `output/pt_transfer/<cruise>_beadstat.mat` to `<cruise>/bead_calibrated/beadstat.mat`
   (keeping the current file as a backup). That needs your OK.
2. **Treat OD2-on vs neutral SSC as a separate, explicit correction** before SSC-to-volume. Build
   it from the March 2019 cultures run with OD2 (SSC 260) and without (SSC 200), together with
   their Coulter sizes. That calibration also gives the size dependence that beads alone can't.
3. **Flag EN627 SSC** (degraded scatter optics). Fluorescence channels can use the transfer.
4. Optional: export the post-Aug 2019 performance tests too, for QC and dating events (for
   example, the misaligned blue laser on AR70B in Sep–Nov 2022).

## Installed (14 Nov 2026)
New tables installed as `<cruise>/bead_calibrated/beadstat.mat` for EN608 (12 tests), AR28B (8),
AR29 (16), EN617 (7) and AR31A (3); `meta.source = "PT performance test transfer"`. The previous
EN617 and AR31A tables (PT conversion) are kept as `beadstat_PTconversion.mat`. AR28B and AR29 had
no `bead_calibrated` folder; one was created. EN627 is unchanged (PT conversion), with SSC flagged.
**Do not rerun step 4 of the processing wrapper for these cruises.** It would rebuild
`beadstat.mat` from PT bead runs (`p.beadtype = 'PT'`) or find no FCB beads. Rebuild with
`make_perf_beadstat` instead.

**Rebuilding and installing: `make_perf_beadstat(cruisedir, opts)`** (one call per cruise; the
wrapper's bead step for cruises where `attune_config(cruise).fcb_beads` is false). It loads or
rebuilds `pt_perf.mat` and the offsets, runs `perf_beadstats` for that cruise, and installs
`bead_calibrated/beadstat.mat`:
- replaces an earlier PT-transfer table after backing it up (`beadstat_backup_<time>.mat`)
- refuses to replace any other table (FCB pipeline, old PT conversion) unless `opts.force`
- refuses cruises that `attune_config` says have FCB beads unless `opts.force`
- `opts.install = false` builds only; `opts.rebuild = true` re-reads all performance-test FCS

Tested 15 Nov 2026 with installs to a scratch folder: EN608 rebuilt identical to the installed table;
a rerun backed up and replaced; EN617 with its PT-conversion table refused, then replaced with
force; TN368 refused; AR28B 8, AR29 16, AR31A 3, EN627 9 tests built. EN627 keeps its PT conversion
(not forced) because its SSC is unreliable anyway.

## 2019 volume calibrations for the pre-Aug 2019 cruises (sizecal_2019.m)
![2019](../output/pt_transfer/sizecal_2019.png)

The same 12 cultures as the 2019 fit (Syn = strain 1820i, 9 Apr; the others 20 Mar), all modes in
log space:
- **Attune:** KDE peak of log10 SSC-H, events gated on fluorescence only (the same gate in both
  runs).
- **1 µm bead:** `find_beads` on the same day's bead runs. OD2: 10^3.002 (20 Mar) and 10^2.994
  (9 Apr). No OD2: 10^4.090 (9 Apr). The 20 Mar no-OD2 run has too few beads, so 9 Apr is used,
  as in the 2019 analysis.
- **Coulter:** mode per log diameter above `CCsize_threshold`
  (`modes_coulter_2019.png`). It moves Nano (+0.08), Micromonas (+0.06) and Pavlova (+0.06),
  where low-size noise ramps up just above the threshold. Pptp02 has a flat top: the mode is
  18.7 µm (2019) or 21.4 µm (per log d), +0.17 log. Treat these volumes as ±0.1–0.2 uncertain.

| Configuration (cruises) | Fit | 95% CI p1 / p2 | rmse | n |
|---|---|---|---|---|
| **OD2 on SSC, 260 V** (EN608, AR28B, AR29, EN617, AR31A, TN368) | **1.4066x + 1.2109** | 1.23–1.59 / 1.03–1.40 | 0.233 | 12 |
| No OD2, 200 V, raw SSC < 10^4.8 (AR34B, RB1904) | 1.2637x + 1.1465 | 0.63–1.90 / 0.77–1.53 | 0.250 | 5 |
| No OD2, 200 V, all cultures (compressed above ~5× bead) | 1.6847x + 1.3041 | 1.31–2.06 / 1.00–1.61 | 0.395 | 12 |

- **OD2 on SSC** is steeper than the 2026 SSC/GL1 merge (1.196). Small cells agree; large cells
  have less scatter per volume at 488 nm through OD2 than on GL1 at 532 nm. Syn (raw 10^2.03 with
  OD2) includes its low-signal excess (≈ +0.08), which partly matches the same excess in field Syn
  on those cruises.
- **No OD2:** SSC at 200 V is linear only up to raw ≈ 10^4.8, or ~5× the bead (≈ 50–100 µm³,
  Pavlova). Above that, cultures of very different volume (190–3500 µm³) fall within 0.25 log of
  SSC. The linear-range fit has only 5 points and a wide CI. The all-culture fit maps the
  compression empirically, but it depends on cell type (rmse 0.40). At sea (220 V) compression
  starts at about 3–5× bead.
- **`batch_mergedSSC2volume` must not run on these cruises as it stands.** It merges `GL1-H`,
  which in the pre-Aug 2019 layout is PE fluorescence (575/36), not side scatter.

## Why 488 nm and 532 nm calibrations differ for large cells: Mie check (../mie_scatter/)
![mie](../output/pt_transfer/mie_488_532.png)

`bhmie.m` (Bohren & Huffman, checked against reference values), `ssc_mie.m` (unpolarized side
scatter over a cone around 90°), `mie_488_532_check.m`.
- **Beads:** a cone of ±30–50° gives the measured 0.5/1 µm ratio (−0.77 to −0.85; measured −0.83
  to −0.86); ±64° does not (−0.6). Mie 6/1 µm is 1.37–1.50. The linear measurements (SSC with
  OD2 1.26, GL1 1.28) are close to that; SSC without OD2 (0.89–0.98) is far below it, which
  confirms that SSC without OD2 compresses the 6 µm bead.
- **Wavelength alone does not explain it.** Non-absorbing cells (m = 1.04) give |532 − 488|
  < 0.03 log in bead units at every size.
- **Pigment absorption does.** With k = a·λ/4π and k(440 nm) = 0.010, cells with chlorophyll
  a/b or many carotenoids scatter relatively more at 532 nm:
  - +0.05–0.1 below ~10 µm³
  - +0.1–0.15 at 50–60 µm³ (observed on the same 2026 cells: I. galbana +0.09, Isochrysis
    +0.15, Pavlova +0.14)
  - +0.3–0.45 at 1000–3000 µm³ (the gap between the 2026 GL1 and 2019 OD2-on-SSC calibrations:
    +0.3–0.4)

  The size of the effect scales with k(440). 0.004 gives about half of it; 0.008–0.012 matches.
  Cryptophytes (phycoerythrin absorbs at 532 nm) go the other way (−0.1 to −0.2 for large
  cells), consistent with Crypto being the 2026/2019 outlier at only +0.06.
- **Consequence:** the 488 nm and 532 nm calibrations differ for a physical reason (absorption
  that grows with cell size and pigment content). Separate calibrations per configuration are
  needed, and the composition uncertainty is larger at 488 nm and for large cells.

## One 488 nm calibration for all pre-Aug 2019 cruises
For cultures where SSC without OD2 is still linear, the bead-normalized SSC with and without
OD2 agree to within 0.02 log: Syn −0.966/−0.965, Micromonas −0.562/−0.553, Nano −0.340/−0.355,
I. galbana 0.162/0.163, Pavlova 0.638/0.656. Both are the same 488 nm measurement, so the
OD2-on-SSC fit (1.4066x + 1.2109, 12 cultures) is the linear 488 nm calibration for both groups.
The straight-line fit through all cultures without OD2 (1.68x + 1.30) is not a meaningful form
across the transition from linear to compressed. For the no-OD2 cruises the open question is
how to handle the compressed range: correct for it with a curve from the paired culture runs, or
flag.

## Compression correction for SSC without OD2 (ssc_compression_fit.m, ssc_decompress.m)
![compression](../output/pt_transfer/ssc_compression.png)

- **Data:** the same 13 cultures run without OD2 (200 V) and with OD2 on SSC (260 V), gated on
  fluorescence only. SSC-H quantiles 10–90% are paired across the two runs, each normalized to
  its own 1 µm bead (t = true from OD2, s = measured without OD2). Pairs with raw OD2 < 10^2.4
  are left out.
- **Model:** `s = t − c/b·log10(1 + 10^(b(t − t0)))`. Slope 1 below the knee, 1 − c above, with
  b setting how gradual the transition is.
  - Fit: t0 = 1.09 (12× bead), b = 1.70, slope above the knee 0.065, i.e. nearly flat. rmse
    0.057 log, n = 193.
  - The culture peaks fall within ±0.05. The off-curve tails (Htriq, H. akashiwo, Prymn pav
    10–40%) are debris.
  - Check: the 6 µm bead of the 2020 filter swap at **220 V** (true 1.257) measures 0.964,
    against 1.001 predicted. The curve in bead units carries over to 220 V.
- **Voltage:** AR34B and RB1904 sampled at **SSC 160 V** (EN627 at 200 V), where compression is
  worse. `ssc_decompress(s, hv)` moves the knee by the change in the 6/1 µm ratio in the SSC
  voltage sweeps: knee 1.00 at 160 V, 1.09 at 200 V, 1.085 at 220 V.
- **Inversion:** `ssc_decompress` returns linear t (then use the OD2-on-SSC calibration,
  1.4066x + 1.2109). It returns NaN and flags where the curve is too flat (slope < 0.5, so
  errors in s would be amplified more than 2×):
  - at 160 V: up to measured 7× bead (true 11×, ≈ 460 µm³)
  - at 200 V: up to 8.7× bead (true 13×, ≈ 620 µm³)

  Larger cells are flagged.

## Pigment signals in the 2026 cultures (culture_pigment_ratios.m)
![pigments](../output/pt_transfer/culture_pigment_ratios.png)

- **Chl fluorescence excitation ratio, log10(GL4/BL3)** (ex 532 / ex 488, both 680/40,
  voltage-corrected). The groups separate largely independent of size:

  | Group | log10(GL4/BL3) |
  |---|---|
  | Green | −0.92 to −1.09 |
  | Haptophytes | −0.77 to −0.84 |
  | Diatom (T. weissflogii) and cryptophytes | −0.57 to −0.75 |
  | Dinoflagellate (Heterocapsa) and coccolithophore (SE62) | −0.39 to −0.41 |
  | Cyanobacteria | −0.08 to −0.33 |

- **Scatter ratio, 532/488** (GL1 vs SSC decompressed at 220 V, bead units). It mostly tracks
  size, as Mie predicts: ~0 for cells < 10 µm³, +0.11–0.16 for haptophytes at 50–60 µm³.
  - At 190–260 µm³ the ordering matches Mie: Dunaliella (green) +0.40 > Rhodomonas +0.25 >
    Crypto +0.17.
  - At 1300–1700 µm³, SSC is at or beyond the limit of the compression correction (Heterocapsa
    52% of events flagged; SE62 near the limit). The values (+0.35 to +0.69) mix pigment with
    leftover compression.
- **Conclusion:** the excitation ratio is the cleaner pigment signal (no compression, little size
  dependence). The scatter ratio adds a pigment-related spread of ~0.2 at 100–600 µm³ but is
  confounded by size and, above that, by compression.
- **Caveats:** one diatom and one or two dinoflagellates. Cultures are not acclimated like field
  cells. In large cells the excitation ratio can also be affected by self-shading of pigment
  (packaging).

## Low-signal nonlinearity under OD2 (gl1_lowsignal_check.m, apply_lowsignal_curve.m)
![gl1](../output/pt_transfer/gl1_lowsignal.png)

On 15 cruises with the current layout (Aug 2019 – 2025, 12 files each), every particle is
measured on SSC (488 nm, no OD2) and on GL1 (532 nm through OD2). Each channel is normalized to
its own 1 µm bead.

- **Plateau:** for raw GL1 between ~10^2.1 and 10^3.8, GL1 and SSC agree to within +0.03
  (the 488/532 nm difference), across all cruises.
- **Low end:** below raw GL1 ≈ 10^2 (100 counts), GL1 reads increasingly **high**: about
  +0.04 at 10^2.0, +0.15 at 10^1.8 and +0.95 at 10^1.6. That is low-signal nonlinearity under OD2.
- **Syn:** it sits right at the knee (raw GL1 10^1.8–2.0) and reads +0.10 to +0.27 log high on
  GL1 relative to SSC (median +0.16), on every cruise.
- **High end:** above raw GL1 ≈ 10^3.8 (SSC_norm > +0.5) the difference rises again. That is
  SSC without OD2 compressing large particles, consistent with the voltage-sweep finding.
- **The merge window** for the current layout (`cfg.gl1_merge` = 500–1000 raw GL1 counts,
  10^2.7–3.0) lies on the plateau, which supports the operational SSC/GL1 scheme.

**Does this explain the OD2-on-SSC cruises?** Partly. Their Syn sits at raw SSC 10^1.9–2.07.
Applying the GL1 curve at those levels removes only 0.02–0.07 log of their 0.15–0.25 offset
(TN368: −1.04 → −1.11). Two things limit the transfer:
- The curve is very steep, and where the knee falls depends on each channel's noise and
  baseline. SSC with OD2 is a different detector, at different voltages from GL1.
- Syn's own GL1 − SSC excess (+0.16) is larger than the all-particle curve at the same raw level
  (+0.04–0.10). Some of it may be spectral: phycoerythrin absorbs at 532 nm.

The direct test for SSC with OD2 is the March 2019 culture set: the same cultures, including Syn
and small eukaryotes, run without OD2 (SSC 200 V) and with OD2 on SSC (260 V). The OD2 loss per
culture, plotted against the raw level under OD2, gives the low-signal curve for that channel
itself.

## Direct test on SSC: March 2019 cultures with and without OD2 (culture_od2_check.m, culture_od2_fluorgate.m, syn_od2_correction.m)
![cultures](../output/pt_transfer/culture_od2.png)

Each culture was run minutes apart without OD2 (SSC 200 V) and with OD2 on SSC (260 V); all
other voltages were identical. The figure shows OD2 loss = Δlog SSC-H − 8.5·log10(260/200), for
quantiles 5–90% of each culture.
- **Gating:** each culture is gated on fluorescence only (chl > 10^2.5, plus PE > 10^3.5 for
  Syn), identically in both runs, so event counts match (Syn 1447 vs 1451). The existing gated
  OD2 Syn file has only 484 events against 1127 without OD2; its gate cuts the low-SSC tail.
- **Plateau:** −2.04 (Micromonas, Nano, I. galbana, SA2, Pavlova and the 0.5/1 µm beads all
  within ±0.05).
- **Syn reads high with OD2, more so at lower signal:** +0.22 at raw 10^1.8, +0.15 at 10^1.9,
  +0.09 at 10^2.0, +0.05 at 10^2.1, 0 by 10^2.5.
  **Micromonas at the same raw level (10^2.0) shows none**, so the excess is not a function of
  signal height alone. Pulse shape or width for the smallest cells is a candidate cause.
- **The 1 µm bead reads 0.03 low with OD2** (loss −2.07).
- **High end:** large cells are compressed in the run *without* OD2 (200 V). The excess is
  +0.05–0.1 from raw 10^4.6–4.8 (no OD2), +0.1–0.2 for H. akashiwo, +0.13 for the 6 µm bead, and
  **+0.45 to +0.7 for PPTP02 and AMBOPO (raw 10^5.2)**. At the same raw level H. akashiwo is
  compressed more than Dunaliella, so this also depends on more than signal height.
  **The SSC-to-volume calibration for large cells must come from the OD2 (or GL1) data, not
  from neutral SSC at 200 V.**

**Does it explain the OD2-cruise Syn offsets?** Most of it. Predicted offset = Syn excess at the
cruise's raw Syn level + 0.03 (bead term). Observed offset = normalized Syn − median of the FCB
cruises without OD2 in the same month.

| Cruise | Raw Syn (OD2) | Observed | Predicted | Residual |
|---|---|---|---|---|
| EN608 | 2.07 | 0.11 | 0.09 | 0.02 |
| EN617 | 1.98 | 0.18 | 0.13 | 0.05 |
| AR29 | 2.03 | 0.17 | 0.10 | 0.07 |
| AR31A | 2.02 | 0.20 | 0.11 | 0.09 |
| AR28B | 2.05 | 0.20 | 0.09 | 0.10 |
| TN368 | 1.91 | 0.28 | 0.18 | 0.10 |

April has no FCB cruise; its seasonal reference is the mean of Feb and Jul. The Syn
culture's low-signal curve accounts for half to most of each offset, leaving 0.02–0.10. Possible
contributors to the rest:
- natural Syn differing from the culture (Micromonas shows the effect is cell-type dependent)
- the seasonal reference
- uncertainty in the transfer offsets

## July/Aug 2026 culture calibration vs March 2019 (compare_sizecal_2019_2026.m)
![sizecal](../output/pt_transfer/sizecal_2019_2026.png)

2026 uses the standard configuration: SSC without OD2 at 220 V, GL1 (532 nm side scatter through
OD2) at 220 V. Signals are normalized to the same day's 1 µm FCB bead on each channel. Volumes
come from Coulter modes in both years (2026: `CCSummary_saved_10Sep2026.mat`, CCMode; modes are
0.01–0.04 log smaller than the means, which changes none of the conclusions).

- **SSC without OD2 compresses all cells above ~200 µm³ in 2026.** Crypto, Dunaliella,
  Rhodomonas, T. weissflogii, Heterocapsa and SE62 span 200–1800 µm³ (0.95 log) but only
  6.1–8.7× bead in SSC-H (0.15 log), at raw SSC-H 10^5.5–5.7. GL1 spreads them over 9–55× bead
  (0.77 log). SSC-A is compressed as well.
  - The 1 µm bead sits higher in raw SSC in 2026 (10^4.74 at 220 V vs 10^4.09 at 200 V in 2019),
    so more cells fall in the compressed range than in 2019. That makes the SSC calibration look
    different.
- **Small and mid-sized cells:** GL1 and SSC agree for Micromonas, Nano and Syn 7335 (within
  0.04). I. galbana and Isochrysis are +0.09 to +0.15; Pavlova (raw SSC 10^5.36) is +0.14.
- **Syn 8113** (raw GL1 10^2.09) reads +0.17 high on GL1: the low-signal excess under OD2 again.
- **2026 GL1 vs 2019 SSC-with-OD2:** at the same volume, 2026 GL1 sits ~0.1–0.6 log higher for
  mid-sized and large cells (Dunaliella +0.6, I. galbana +0.2, Pavlova +0.1, Crypto +0.06).
  GL1 is at 532 nm, the 2019 data at 488 nm. Normalizing to a bead at each wavelength doesn't
  make cells and beads scale the same way (Mie size parameter), and pigments absorb differently
  at the two wavelengths (carotenoids at 488 nm, phycoerythrin at 532 nm). **The 2019 SSC-with-OD2
  calibration can't be applied to GL1 directly; GL1 needs its own calibration, which is what
  the 2026 cultures provide.**

## Scatter voltages by cruise, and AR104's new settings (scatter_hv_survey.m)
- Every cruise with the current layout (EN644 – HRS2609) ran **SSC 220 / GL1 220 V**. The
  exceptions are a few files on EN644 and DY131 (SSC 260). **AR104 (Oct 2026) runs SSC 240 /
  GL1 180 V** for both samples and bead runs (the "new" settings of the 2026 evaluation).
- For AR104, `bead_calibrated/FCSfileinfo.mat` (FCS_DateTimeList: 2856 files, all QC good) and
  `beadstat.mat` (process_beads_PT_adjust_3 unchanged: 7 runs, all QC good except one 6 µm) were
  made this session.
  - 1 µm bead: SSC-H 10^5.07, GL1-H 10^2.29, just above the GL1 low-signal knee.
  - 6 µm bead: SSC-H 10^5.95 (near saturation; 6/1 µm on SSC 0.88, on GL1 1.28).
- **SSC vs GL1 on AR104 samples** (15 files, `output/pt_transfer/gl1_lowsignal_AR104.mat`):
  - The two channels agree (|Δ| < 0.05) only for SSC_norm −0.4 to +0.6 log bead units
    (raw GL1 10^1.96–2.9, raw SSC 10^4.7–5.6).
  - Below that, GL1 reads high: Syn at raw GL1 10^1.69 is **+0.58**.
  - Above it, SSC compresses: +0.31 at raw SSC 10^5.9. At 240 V, SSC stays linear to about
    10^5.6, higher than at 220 V (~10^5.2–5.3), as the voltage sweeps predicted.
- **The merge window should be in bead units, not raw GL1 counts.** At GL1 180 V the current
  `cfg.gl1_merge` (500–1000 raw counts) sits at SSC_norm +0.4 to +0.75, where SSC is already
  compressing. The range where both channels are linear on both settings is about
  **0.5–2 × the 1 µm bead** (log −0.3 to +0.3). At 220/220 that is raw GL1 ≈ 560–2200, close to
  the current window. Changing `cfg.gl1_merge` is a classifier change, so it needs your OK.

## Merged-scatter volume calibration, July/Aug 2026 cultures (sizecal_merged_2026.m)
![merged](../output/pt_transfer/sizecal_merged_2026.png)

- **Merged scale, per event:** SSC/bead below 0.5× the 1 µm bead (on GL1), GL1/bead above 2×,
  and the mean of the two in between. Each value uses the same day's FCB 1 µm bead on its own
  channel. A culture's value is the KDE mode of log merged scatter; volume is the Coulter mode.
- **Fit without T. weissflogii:** `log10 vol = 1.179·log10(merged) + 1.027`, rms 0.17 log.
  - With T. weissflogii: 1.220x + 1.051, rms 0.23. The diatom scatters little for its volume.
  - Largest residuals: Pavlova −0.35, Dunaliella −0.3, I. galbana/Isochrysis +0.2. These are
    culture-specific (Dunaliella is high on GL1 at 532 nm).
- **At 220/220 the window choice hardly matters.** The raw-count window (500–1000) gives
  1.172x + 1.042. Only I. galbana and Isochrysis, which straddle the window, move (by ≤ 0.07).
- **The current `Attune_SC2vol` SSC-H coefficients (1.1696, 1.1592; 11 Sep 2026) have the same
  slope but an intercept 0.13 log higher.** Against these culture modes, on FCB-bead-normalized
  merged scatter, they overestimate volume by ~35% (rms 0.21 vs 0.17). The cause is still to
  check: a different bead normalization (e.g. the F8819 1 µm bead, which reads 0.08 log lower
  than the FCB 1 µm on SSC-H), mean vs mode, or a different merge.
- **Why the 11 Sep fit differs (reproduced exactly from `size_calibration_fit_jul2026.m`;
  same FCB 1 µm bead normalization):**
  1. T. weissflogii is in the fit. `~contains(Row,'T_weissflogii_CCMP')` removes the single
     runs but keeps `28-Jul-2026_T_weissflogii_pooled`. Dropping it gives 1.130x + 1.114.
  2. The culture modes in `AttuneMode` (`compile_attune_histograms.m`) are modes of
     linear-bin histograms. For a broad, roughly log-normal population that mode is lower than
     the mode in log space by about ln(10)·σ² (σ = log10 SD). Measured: Syn 8113 4075 vs KDE log
     mode 5622 (−0.14 log, σ 0.27); Micromonas −0.14 (σ 0.23); Nano −0.07; Crypto ≈ 0 (σ 0.08).
     The 1 µm beads are narrow, so their two modes agree. The bias depends on how broad each
     culture is, not on its size, and it pushes the small cells left, which raises the intercept.
     The field pipeline works in log space (beadstat KDE peaks; per-event log scatter), so a
     log-space mode (or median) matches it.
  Each effect on its own, using the 11 Sep script's rows and culture-level switch (SSC when
  SSC ≤ 1× bead, else GL1; I. galbana/Isochrysis averaged; Nano pooled):

  | Fit | Modes | T. weissflogii | Result |
  |---|---|---|---|
  | A | linear-bin histogram | in | 1.1696x + 1.1592 (= the 11 Sep fit) |
  | B | log-space (KDE) | in | 1.2069x + 1.0583 |
  | C | linear-bin histogram | out | 1.1299x + 1.1139 |
  | D | log-space (KDE) | out | 1.1663x + 1.0307 |

  `sizecal_merged_2026.m` gives 1.179x + 1.027. It also uses the per-event 0.5–2× blended merge,
  keeps I. galbana, Isochrysis, both Nano runs and C2 Micromonas as separate points, and drops
  T. weissflogii. It is within 0.01 of D.
- **Coulter modes** (`compile_cc_histograms.m`): the Coulter channels are evenly spaced in
  *diameter* (0.0357 µm steps), so the CCMode values are also linear-bin modes. The cultures are
  narrow, so the bias is small: modes in log space are larger by −0.01 to +0.04 log10 volume
  (Micromonas +0.036, T. weissflogii +0.035, Syn 8113 +0.033, Crypto +0.029; most others ≤0.014).
  Fit B (11 Sep rows, log-space Attune modes, diatom kept) becomes **1.2022x + 1.0780** (rmse
  0.267), against 1.2069x + 1.0583 with the current CCMode. Simplest fix in that script: take
  the mode of `count .* diameter_micron` (density per log diameter) instead of `count`.
- **Coulter modes redone (14 Nov).** `compile_cc_histograms.m` now takes the mode of
  `count .* diameter` (per log diameter). It was run with `printNow = 0` (so the NAS plots
  were not overwritten), and the result is saved as `CCSummary_saved_14Nov2026.mat` (the 10 Sep
  file is unchanged).
  - Volume shifts: Micromonas +0.058, C2 Micromonas +0.036, Syn 8113 +0.031, Syn 7335 +0.023,
    I. galbana +0.022, Crypto +0.020; the rest 0 to +0.014 log.
  - Plots: `modes_attune_2026.png` and `modes_coulter_2026.png` (plot_mode_comparison_2026.m).
  - **Fit, 11 Sep rows, log-space Attune modes, diatom kept, 14 Nov Coulter modes:
    log10 vol = 1.1962·log10(merged) + 1.0826** (95% CI p1 1.009–1.383, p2 0.879–1.286;
    rmse 0.261, n = 12). With the 10 Sep Coulter modes it is 1.2065x + 1.0591.
    **Now in `Attune_SC2vol.m` (SSC-H)**; plot `vol_cal_SSC-H_merge_14Nov2026.png`.
- **Classifier updated:** `cfg.gl1_merge_bead = [0.5 2]` (bead units) replaces
  `cfg.gl1_merge = [500 1000]` (raw counts) in `classifier_config.m` and `assign_class_v2.m`.
  - The old raw window was 0.45–0.9× bead at 220/220 and 2.6–5.2× bead on AR104.
  - Test on 3 AR92 and 3 AR104 files: 0.06–0.3% of events change class (1.6% on one AR104 file).
  - `batch_mergedSSC2volume.m` (the production volume step) now uses the same bead-unit
    window (`cfg.gl1_merge_bead`, applied to `gl1_bdnorm`). The window is saved in `mergeTable`
    (`GL1_merge_bead`) and written to each class file's `vol_notes`.

## Related findings (same session)
- OD2 filter swap (21 Jan 2020): OD2 on SSC removes 2.08 log from 1 µm beads but only 1.81 from
  6 µm and 1.32 from 9 µm beads. Bead ratios with OD2 reproduce TN368's anomaly (6/1 µm 1.26 vs
  1.23).
- SSC voltage sweeps (Dec 2019, Jul 2020): with the same sample, the 6/1 µm ratio rises from 0.89
  at 150 V to 0.98 at 230 V, so the detector response to large particles depends on voltage
  (detector nonlinearity). For the current layout (SSC for small cells, GL1 with OD2 for large
  cells) no general response model seems needed.
- The Mie geometry fit originally planned was set aside: the bead-ratio anomalies are explained
  by OD2 and detector effects, not by changes in scatter geometry.

## Field pigment map update (15 Nov 2026)
- The summer low-ratio cryptophyte stripes were EN644 `Dockwater_filter_config_change` files
  (6–15 Aug 2019, old layout) read with the current channel map. They are now excluded in
  `field_pigment_map`. See CALIBRATION_SUMMARY section 5.
- The March/April 2019 cultures were added on the same scale (local `cultures_2019`; bead
  ref from find_beads on the same day's FCB run; volume `SSC-H_488_2019`). Green, haptophyte,
  cryptophyte and Syn values match 2026 within ~0.1.
- SE62 is a coccolithophore (haptophyte); relabelled in `culture_pigment_ratios` and the map.
- AMBOPO is *Alexandrium minutum*; labelled dinoflagellate in the map (−0.38). Pptp02 stays
  unconfirmed.
- `field_pigment_map('replot')` redraws from the saved .mat.

## Bead ratio history (bead_ratio_history.m, Nov 2026)
`T = bead_ratio_history()` collects per-run 1/0.5 and 6/1 µm FCB bead peak ratios (FSC-H, SSC-H,
GL1-H for the current layout) from every cruise's `beadstat.mat`, and plots them against the Thermo
service dates. Output: `output/pt_transfer/bead_ratio_history.png` and `.mat`.
FSC steps by ±0.1–0.18 at the service events. SSC 6/1 drifts −0.1 from 2020 to 2026. GL1 is stable
to ±0.03. See CALIBRATION_SUMMARY section 5.
