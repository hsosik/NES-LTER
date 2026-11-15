# Bead stats: review and simplification (handoff)

Starting point for reviewing and simplifying how bead statistics are computed and saved. It was written from the classifier redesign work in `classifier_dev/` (Nov 2026).

## Current pipeline

- **[process_beads_PT_adjust_3.m](../process_beads_PT_adjust_3.m)** is called from the processing wrapper with `(basepath, FCSfileinfo, beadfiles2include, beadtype, OD2setting)`.
  - It picks bead files by name match, clusters the 0.5, 1 and 6 µm FCB beads (classes 1, 2, 3), and sets a QC flag.
  - It writes `<cruise>/bead_calibrated/beadstat_2026.mat`. That file holds the variables `beadstat_2026`, `notes`, `beadtype`, `bead_ch_names` and `beadfiles2include`.
  - It writes figures to `bead_calibrated/bead_plots_2026/`.
  - There is a separate PT-bead branch with 2.4 and 3.2 µm beads.
  - The SSC-A/GL1-A 1 µm ratio (`mean1micron_SSCA2GL1A`) is also computed, for the OD2 size merge.
- **Uncommitted local edit, keep it:** the `mode1micron_SSCA2GL1A` line (around line 490) is commented out.
- An older generation, `beadstat_2021`, may still be read elsewhere (`process_wrapper_*`, the `generate_attune_table` path). Check before renaming anything.

## What the classifier prototype reads

The new classifier reads the bead table only through `classifier_dev/bead_reference.m`. A simplified table must still provide, or `bead_reference` must be adapted to:

| Field | Use |
|---|---|
| `time` | bead runs within ±`cfg.bead_window_days` (1.5 d) of the sample; nearest run if none |
| `QC_flag` | logical; only flagged-good runs are used |
| `median1`, `median2`, `median3` | median peak **height** of the 0.5, 1 and 6 µm populations, indexed by channel name `'<ch>-H'`, for FSC, SSC, GL1, BL1, BL2, BL3, GL2, GL3, GL4 |
| `hv` | detector voltage per channel for each run |
| `parname` | channel names that line up with the `hv` columns (`parname(row, 2:end)`) |

- **Which bead each channel uses** (`cfg.bead_pop`):
  - 0.5 µm: BL2, GL2, GL4
  - 1 µm: FSC, SSC, GL1, BL1, BL3
  - 6 µm: GL3
- `attune_features.m` brings the bead position to each sample's voltage, with signal proportional to hv^k (`cfg.k`: 8.5, or 7.1 for SSC).
- Areas, means, counts, the 2.4/3.2 µm PT beads and the OD2 ratio columns are not used by the classifier.

## Findings from the classifier work

These are relevant to what is worth keeping.

- **Voltage scaling.**
  - Detector response follows a power law, signal ∝ hv^8.5, fitted across 332 QC'd bead runs from 39 cruises (`collect_bead_runs.m`, `fit_gain_curves.m`).
  - Syn in samples alternating between Pro and standard settings gives the same exponent (8.4 ± 0.6).
  - Correcting to a reference voltage removes the settings offset (`analyze_normalization.m`).
- **Scatter becomes cruise-independent once bead-relative.**
  - Between-cruise SD of the Syn position: SSC 0.39 → 0.07, GL1 0.74 → 0.06.
  - Fluorescence brightness does not become cruise-independent: the remaining spread is biological, tied to temperature.
- **0.5 vs 1 µm on SSC.**
  - Both peaks are equally repeatable, with run-to-run SD 0.005 log10 within a cruise.
  - They sit a constant 0.84 ± 0.01 log10 apart across 20 cruises.
  - The classifier uses the 1 µm bead.
- **GL1 vs SSC.** SSC and GL1 (OD2 scatter), each on its own 1 µm bead scale, are treated as equivalent for the size merge, as in `batch_mergedSSC2volume`. In the linear range they differ by −0.02 to −0.11.
- **OD2 placement.** The GL1 size merge assumes OD2 is on GL1. Cruises with OD2 on SSC, or with no OD2, need that recorded per cruise; `OD2setting` already exists as an input.
- **Classifier README.** Details are in `classifier_dev/README.md`, sections "Test 1" and "SSC bead reference".

## Possible simplification directions (to discuss)

- One tidy table per cruise: one row per bead run per population per channel, holding time, file, population, channel, peak height median, robust SD, n, hv and QC. The wide `median1/2/3` cell layout would go.
- Keep only what downstream code reads, so the clustering diagnostics live in the plots rather than the table.
- Store the voltage exponent and the reference voltage with the table, or keep them only in config.
- A single bead-reference function shared by the classifier and the existing processing chain.
- Decide whether the PT-bead branch is still needed.

## Constraints

- **Classifier contract:** the classifier prototype must keep working. Rerun `classifier_dev` scoring (`score_agreement` on 15 files each from AR92, AR95, EN727 and HRS2609) after any change to the table.
- **Branch and uncommitted edits:** the work is on branch `classifier-redesign`. Other uncommitted edits belong to the user and must be left alone: `Discrete_Samples/*`, `assign_class_summer.m`, `attune_lter_moviemaker.m`, `process_wrapper_2024b.m`, and the IFCB files.
- **Commits:** don't commit unless the user asks.
- **Data paths:** data is under `\\sosiknas1\lab_data\Attune\cruise_data\<cruise>\bead_calibrated\`. In bash heredocs, write the UNC path as `['\' '\sosiknas1...']`.

## Implemented (Nov 2026)

`process_beads_PT_adjust_3.m` was rewritten. It now writes `bead_calibrated/beadstat.mat` (variables `beadstat`, `meta`) and plots in `bead_calibrated/beadstat_plots/`. The old `beadstat_2026.mat` and `beadstat_2021.mat` files are left in place but nothing reads them any more. Read the new file with `load_beadstat(bp)`.

- **Table:** one row per bead run: `filename`, `time`, `hv` (table by physical channel), and per population P = 05, 1, 6: `nP`, `QCP` (true = good), `whyP`, `peakP` (KDE mode, linear), `medP` (median, linear), `rsdP` (robust SD of log10). Channels are physical channels named in the current filter layout (`FSC_H` … `GL4_H`, `SSC_A`, `GL1_A`).
- **Code (classifier_dev):** `attune_config` (OD2, filter layout and FCB-bead availability per cruise, from the configuration history sheet), `bead_signature` (reference bead positions, gating channels, QC limits; k from `classifier_config`), `find_beads` (gating and QC), `bead_peak` (KDE mode), `measure_beads` (two passes per cruise and the table), `plot_bead_gates`, `generate_beadstats` (batch over all or listed cruises).
- **Consumers updated:** `bead_reference.m` (QC and peak per population), `batch_mergedSSC2volume.m` (1 µm peak, QC'd runs only), `collect_bead_runs.m`, `extract_cluster_centers.m`, and 27 dev scripts (`load_beadstat`).
- **Regenerated:** 42 cruises, 762 runs. QC good: 0.5 µm 736, 1 µm 752, 6 µm 677. EN608 to EN627 have no FCB bead runs.
- **Classifier check** (`score_agreement`, 15 files each; `output/beads_baseline` vs `output/beads_new`): mean agreement AR92 0.9833 → 0.9833, AR95 0.9513 → 0.9490, EN727 0.9226 → 0.9216, HRS2609 0.9778 → 0.9774. The bead references move by ≤ 0.011 log, the largest being GL1, where the median sat low. One AR95 AltPro file drops 0.961 → 0.905 because 800 events are added instead of 218.
- **Design and test scripts:** in `bead_design/` (`run_bead_prototype`, `compare_bead_prototype`, `plot_bead_disagreements`; outputs in `output/bead_proto*`). The pipeline does not need them.

## Cruises with only PT beads (Nov 2026)

`generate_pt_beadstats.m` (with `find_pt_beads.m`) writes `beadstat.mat` for EN617, AR31A and EN627 in the same format as the FCB cruises, with `meta.source = "PT conversion"`. EN608, AR28B and AR29 have no PT bead data files.

- **PT beads:** one size, three intensity levels (~1.4 log apart) per fluorescence channel, repeatable to about 0.02 log. The PT runs are in `cruise_data/beads/PT_beads_seawater_settings` and `PT_beads_bead_settings`.
- **Conversion:** FCB peak = PT level + (FCB − PT) offset. The ratio of two beads on one detector does not depend on gain or laser power, so the offset carries over between dates with the same filters.
- **Calibration:** FCB runs of 20 Mar – 12 Apr 2019 against PT runs of 28 Feb – 5 Mar 2019, at identical voltages.
- **Cross-check:** PT runs at bead settings on the same days as the FCB runs (scaled with k) give offsets within 0.01–0.02 log (fluorescence) and 0.015–0.03 (FSC, SSC), with run-to-run SD 0.005–0.027.
- **AR31A:** only bead-settings PT runs exist. They are brought to seawater settings with the July 2018 same-day bead→seawater offsets (pair SD ≤ 0.03). Its 25–26 Oct PT runs show no PT levels on the green channels (laser or delay off?), so GL2–GL4 are NaN. A 16 Nov 2018 run looks normal but falls outside the cruise window.
- **EN627:** the PT scatter (FSC, SSC) moves between runs by up to ~0.5 log, which looks like a real change in instrument state (cleaning). The synthetic FSC and SSC therefore vary too (SD 0.19, 0.09); fluorescence varies 0.02–0.06.
- **GL1:** no green scatter channel before Aug 2019, so there is no GL1 reference. `bead_reference` now warns and returns NaN for a channel without QC'd runs.
- **Outputs:** calibration and per-cruise tables in `output/pt_conversion.mat`.
