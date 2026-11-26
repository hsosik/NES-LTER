# NES-LTER Attune NxT: instrument and processing history (2018–2026)

The underway Attune NxT (serial 2AFC216190617; blue 488 nm, 50 mW and green 532 nm, 100 mW) has
run on NES-LTER and partner cruises since January 2018. This page brings together what is recorded
about how the instrument was set up and serviced, and how its data have been processed. Part 1 is a
narrative history; Part 2 is a per-cruise reference table; Part 3 lists every *Dunaliella* culture run.

**Sources.**
- Settings history sheet: `\\sosiknas1\lab_data\Attune\Attune_Settings_configuration_history_examplecopy_20221130.xlsx`. It ends in 2022.
- Lab logs:
  - `Sosik_lab_Attune_documents\old\Attune_LOG_EndOct9_2020.rtf` (Jan 2018 – Sep 2020)
  - `Sosik_lab_Attune_documents\Attune_Log_4_14_2021_backup_Nov2025.docx` (Apr 2021 – Nov 2025)
- Thermo Fisher service reports in `\\sosiknas1\lab_data\Attune\Service reports\`.
- From each cruise's `bead_calibrated\` folder: `FCSfileinfo.mat` (triggers and sample voltages) and `beadstat.mat` (bead runs).
- In the repository: the processing wrappers, `classifier_dev/README.md`, `classifier_dev/pt_transfer/CALIBRATION_SUMMARY.md`, and the git history.

**Gaps.**
- No lab log covers Sep 2020 – Apr 2021 (EN657, EN661) or anything after Nov 2025.
- Where the sources disagree, the conflict is listed in Part 1, section 7.

Compiled Nov 2026.

---

## Part 1. Narrative

### 1. Optical configuration: five eras

**Era A: SSC through OD2, no green scatter (Jan 2018 – Oct 2018)**
- Cruises: EN608, AR28B, AR29, EN617, AR31A.
- Blue laser channels:
  - SSC 488/10 with an OD2 neutral-density filter
  - BL1 530/30
  - BL2 590/40, then 575/40 from EN617
  - BL3 695/40, then 680/40 from EN617
- Green laser channels:
  - GL1 575/36 (PE)
  - GL2 620/15
  - GL3 695/40, later relabelled 680/40
  - GL4 780/60
- There is no green side-scatter channel.
- Trigger: OR on BL3 (chl) and GL1 (PE), at 100 and 200.
- Sample voltages, from the FCS headers (set by Alexi on 30 Jan 2018): FSC 160, SSC 260, BL1 280, BL2 320, BL3 300, GL1 (PE) 300, GL2 400, GL3 320, GL4 400.
  - EN617: GL1 (PE) → 340 (Syn were clipped at the threshold) and BL3 → 340. AR31A kept these.
  - The performance-test (PT) and bead-settings runs used other voltages (e.g. SSC about 380).
- Pulse-width threshold: 1 until AR31A, then 0.1. This improved size resolution.
- Intake mesh: 80 µm on EN608, then 150 µm.

**Era B: OD2 removed from SSC (Feb – May 2019)**
- Cruises: EN627, AR34B, RB1904.
- On 1 Feb 2019, at the start of EN627, the SSC filter was changed to a plain 488/10 with no OD.
- SSC voltage:
  - EN627: tried at 260 on 1 Feb; all sample files were run at 200. BL3 320 (450 files at 340); GL1 (PE) 340.
  - AR34B: 160. Without OD2, the EN627 settings had put large particles off scale.
  - RB1904: switched to 220 on 17 May and back to 160 on 20 May. At 220, SSC-A only covered sizes up to about 20 µm.
- Without OD2, SSC-H compresses for large particles. The 2026 calibration undoes this with `ssc_decompress`; above about 7× (160 V) or 9× (200 V) the 1 µm bead, volumes are only a lower bound.
- EN627 scatter is considered uncorrectable because the optics were degraded in Feb 2019.

**Era C: OD2 back on SSC (Jul 2019)**
- Cruise: TN368.
- The original OD2 SSC filter was reinstalled on 4 Jul 2019.
- Sample voltages: FSC 280, SSC 260, BL2 360, BL3 340, GL1 (PE) 420.
- Higher PE threshold to reduce coincidence: BL3 300 OR GL1 (PE) 1100.
- "Exclude coincident events" was left off, because Attune counts matched IFCB without it.
- The green laser failed after the cruise.

**Era D: current filter layout (Aug 2019 →)**
- Cruises: EN644 onward.
- Before EN644 (6 Aug 2019), the green laser was replaced, the blue laser realigned and the flow cell cleaned.
- On 16 Aug 2019 a new green-laser filter set went in:
  - GL1 530/30 with the OD2 glued on top: green side scatter
  - GL2 575/36 (PE)
  - GL3 620/15
  - GL4 680/40
- SSC is a plain 488/10 again. The blue channels keep their roles.
- Standard sample voltages from EN644 to HRS2609: FSC 280, SSC 220, BL1 280, BL2 360, BL3 340, GL1 220, GL2 360, GL3 320, GL4 260.
  - BL2 and GL2 went to 380 during EN655.
  - GL2 was 400 on the alternating templates (EN688, EN706, part of AR78).
  - Pro-measuring files raise only BL3, to 370–420.
  - The scatter voltages (FSC, SSC, GL1) never changed in this era.
- Size now comes from SSC for small cells and GL1 (through OD2) for large ones. The two are merged in 1 µm bead units, switching between 0.5× and 2× the bead on GL1.

**Era E: new settings (Oct 2026)**
- Cruise: AR104.
- New voltages from the `2026_settings_evaluation` bead checks (Sep 2026): FSC 250, SSC 240, GL1 180, BL2 330, BL3 410, GL4 230.
- Standard samples now trigger on SSC AND BL3 (900 AND 1000 or 1200). Pro runs use SSC 500 AND BL3 300.
- The merge window in bead units still holds at these voltages.

### 2. Triggers and alternating settings

| Period | Cruises | Trigger (threshold in counts) |
|---|---|---|
| Jan 2018 – May 2019 | EN608 → RB1904 | BL3 100 OR GL1 (PE) 200 |
| Jul 2019 | TN368 | BL3 300 OR GL1 (PE) 1100, to reduce coincidence |
| Aug 2019 – Nov 2021 | EN644 → AR63 | BL3 300 OR GL2 (PE) 300 |
| Feb 2022 → | AT46 onward | Standard "winter" setting: BL3 500 OR GL2 2000 |
| Jul–Aug 2022 | EN685, EN687 | Summer: BL3 200 / GL2 300 (EN685), BL3 300 / GL2 300 (EN687) |
| Aug 2022 | EN688 (OTZ) | First alternating underway run: SSC AND GL2 500 (Syn) with SSC AND BL3 300 (Euk). GL2 hv 400 on that template |
| Aug 2023 | EN706 | Alternating SSC AND BL3 / SSC AND GL2. On 08 Aug all files ran BL3 500 OR GL2 2000 (template not copied) |
| Oct 2023 | AR77 | BL3 200 OR GL2 1000 for most of the cruise |
| Nov 2023 → | AR78, AR79, AR82A, AE2426, AR87A/B, AR92, AR95, AR96, AR98A, HRS2609 | Pro-measuring runs alternate with standard ones: SSC AND BL3 (SSC 100–500, BL3 200–400), BL3 hv 370–420, sometimes 150 µl |
| Summer/fall | EN720, AR87A/B, AR92, AR95, AR96, HRS2609 | GL2 threshold lowered to 500–1000 so dim Syn are not cut off |
| Oct 2026 | AR104 | Standard: SSC 900 AND BL3 1000 (or 1200); Pro: SSC 500 AND BL3 300 |

- **Whether a file measured Pro** is now flagged per file (`pro_measured` in `FCSfileinfo`). The rule scales each trigger threshold to the file's voltage (signal ∝ hv^8.5).
- **Detector response** with voltage follows about hv^8.5 for SSC, GL1 and the fluorescence channels, and hv^2.1 for FSC. This comes from bead voltage sweeps and the 2026 settings evaluation (`hv_sweep_findings`).
- **Other voltage changes:**
  - EN655, 26 Jul 2020: BL2 and GL2 360 → 380.
  - EN688 alternating template: GL2 400.

### 3. Beads

- **2018: performance-test (PT) beads only.**
  - EN608, AR28B and AR29 had no beads run as samples.
  - From Jul 2018 (EN617) PT beads were also run as samples, after the engineer advised it.
  - EN608–AR31A and EN627 have no FCB runs. Their bead references are transferred from the PT runs to FCB units (`pt_transfer/make_perf_beadstat`; `meta.source` "PT performance test transfer"). Validated against paired runs: residual 0.005–0.012 log.
- **FCB bead mix (0.5, 1 and 6 µm), made Feb–Mar 2019 and run routinely from AR34B (Apr 2019).**
  - Recipe in `\\sosiknas1\Lab_data\Attune\chemicals and reagents\FCB bead mix recipe.docx`.
  - New working mixes were made in Jun 2019, on 3 Aug 2022 (extra 1 µm) and on 11 Dec 2024 (less 0.5 µm).
  - The Dec 2024 mix was a new dilution of the same stock bottles; the logs don't say whether the earlier mixes used new lots.
- **Bead tables:**
  - Three generations: `beadstat_2021`, then `beadstat_2026`, then `beadstat.mat` (Nov 2026).
  - `beadstat.mat` is now the only one read. It has one row per run, per-population QC, KDE peak positions, and review exclusions (`bead_review_decisions.csv`).
  - Everything is normalized to the 1 µm bead on the same channel at the sample's voltage.
- **Bead peak doubling:**
  - Peak doubling or broadening at low signal began between AE2426 and EN727 (Nov 2024 – Jan 2025).
  - The GL1 6/1 µm ratio dropped from 1.32 to 1.25 between 17 Dec 2024 and 14 Jan 2025. That was the instrument, not the bead mix.
  - The Feb 2025 flow-cell replacement tightened the peaks but did not restore the ratio.

### 4. Service and optics events

The full dated list, with evidence from the bead ratios, is in the shared session notes
(`attune_optics_changes.md`). The events most likely to shift scatter or fluorescence scales:

| Date | Event | Between cruises |
|---|---|---|
| Jan 2018 | Install; software 2.6; blocker bar aligned | before EN608 |
| 3 Jul 2018 | Gear pump, syringe and filters replaced; blocker bar adjusted; software reinstalled | AR29 / EN617 |
| 30 Oct 2018 | Gear pump and filters replaced, flow cell flushed by hand; software 2.7.0 | AR31A / EN627 |
| 24 Jan 2019 | Service call for low or unstable SSC: Contrad soak of the flow cell, blocker bar adjusted, alignment verified; bead scatter then "lowest ever" | before EN627 |
| Feb 2019 | 25% Hellmanex and Contrad soaks (Contrad later said to etch) | EN627 / AR34B |
| 6 Aug 2019 | Green laser replaced; blue realigned; flow cell exterior cleaned; blocker bar adjusted "to maximize FSC" | TN368 / EN644 |
| 23 Oct 2019 | Preventive maintenance: alignment, pinholes, blocker bar verified | AR39B / EN649 |
| 7 Jul 2020 | Gear pump replaced | AR43 / EN655 |
| Aug 2020 | Rotary valve, syringe pump; bent bracket gave slight SIP misalignment; Windows 10 and software 4.2 | EN655 / EN657 |
| 22 Feb 2021 | Contamination cleared from the flow cell; blocker bar adjusted. **FSC 6/1 +0.10, SSC 6/1 −0.06** | EN661 / EN668 |
| Apr 2021 | New computer (data desync); beta FPGA 20.3.4.12 (pulse-width limit) | EN661 / EN668 |
| 13 May 2022 | Preventive maintenance | AR66B / EN685 |
| 22 Sep 2022 | Blue laser replaced; flow cell flushed with Hellmanex; gear pump replaced | EN688 / AR70B |
| 29 Nov 2022 | Blue laser realigned (it had been run misaligned on AR70B). **FSC 6/1 −0.18 across Sep–Nov** | AR70B / EN695 |
| 5 Oct 2023 | Software 5.3; new FPGA firmware with a 375 pulse-width limit | EN706 / AR77 |
| 17 Jan 2024 | Gear pump replaced | AR79 / EN712 |
| 21 Jun 2024 | Blue laser alignment adjusted; rotary valve head replaced | EN715 / EN720 |
| 2–8 Jul 2024 | Interrogation assembly (flow cell) and tube lifter replaced; laser path optics adjusted; blocker bar and sled adjusted. **FSC 6/1 +0.16** | EN715 / EN720 |
| Dec 2024 – Jan 2025 | Bead doubling at low signal. **GL1 6/1 −0.07, no service visit** | AE2426 / EN727 |
| 11–14 Feb 2025 | Interrogation assembly (flow cell) and rotary valve replaced; PZT power adjusted; software 6.2.2 (FPGA 20.10.12.8); preventive maintenance with alignment checks on 14 Feb | EN727 / AR87A |
| 12 Feb 2026 | Preventive maintenance | AR99 / HRS2601 |
| 13 Jul 2026 | Sheath diaphragm pump replaced | HRS2601 / HRS2609 |

Software versions: 2.6 → 2.7.0 (Oct 2018) → 3.1.2 (Jul 2019) → 4.2 (Aug 2020) → beta FPGA
20.3.4.12 (Apr 2021) → 5.3 (Oct 2023) → 6.2.2 (Feb 2025).

### 5. Processing generations

| Period | Main code | What changed |
|---|---|---|
| 2018 | Mira Santos' QC scripts; summary-statistics exports matched to underway data | Counts from Attune software gates |
| Nov 2018 – 2020 | `assign_class_<cruise>.m` per cruise (first for AR29/SPIROPA); `FCSfileinfo` (2019) | Polygon gates in MATLAB on FCS events; bead clustering first done Jun 2020 |
| 2021–2023 (Bethany Fowler) | `process_wrapper_2021.m` | Modular steps: file table, classes, beads (`beadstat_2021`), SSC-A volume calibration with GL1 merge, Attune table, movies, underway match, volume tables. First EDI tables 2022–2023. QC flag for garbled files (Dec 2023, after EN688) |
| Feb 2024 – 2026 | `process_wrapper_2024.m` / `2024b` | Defaults coded so only the cruise name is needed; per-cruise classifiers (Emily Peacock) for 2022–2026 cruises |
| 2025–Sep 2026 | same wrappers | `pro_measured` flag (Nov 2025); revised SSC/GL1 merge and culture size calibration (Mar–Sep 2026, switched to pulse height); `beadstat_2026`; Mahalanobis-refined gates (Sep 2026) |
| Nov 2026 | `classifier_dev/` | One classifier for all cruises (`assign_class_v2`); bead table redesign (`beadstat.mat`); detector exponents from beads; PT transfer for 2018 cruises; new volume calibrations for both eras (`CALIBRATION_SUMMARY.md`); new wrapper `process_cruise.m` (in progress) |

**How v2 differs from the per-cruise classifiers:**
- Every channel is put on a bead-relative scale at a reference voltage before classifying.
- Population identity comes from spectral ratios (PE:chl, GL2:GL3, BL2:GL3), which hold across cruises and seasons.
- Size comes from bead-normalized scatter.
- Clusters are located per file with a Gaussian mixture.
- New classes and labels:
  - Euk with Syn (5) and Syn on particle (8)
  - dim and bright crypto (3, 4) separated by spectrum
  - Syn spectral type A/B/C
  - flags for excluded events

See `classifier_dev/README.md`.

**Size calibration** (`CALIBRATION_SUMMARY.md`). Volume calibrations by era:
- **Current layout (Era D):**
  - Uses the merged SSC/GL1 signal in bead units.
  - 12 cultures, Jul/Aug 2026: `log10 V = 1.1962·x + 1.0826`.
- **Pre-Aug 2019 (Eras A–C):**
  - Uses SSC at 488 nm in bead units.
  - 12 cultures, Mar/Apr 2019: `log10 V = 1.4066·x + 1.2109`.
  - Large-particle SSC is decompressed first on the no-OD2 cruises.
- The class files' volumes have not yet been regenerated with these calibrations.

### 6. Known problems by cruise

- **Time stamps:**
  - EN608 started on Eastern time.
  - EN627, AR34B and RB1904 file times were local, not UTC.
  - AR70B: the Attune clock was 1 h 27 min and 5 days behind. Original files are kept in `copy_FCS_original_badtimestamp`.
  - AR92: dead clock battery.
- **Sample volume:** EN644 templates were set to 200 µl, not 400 µl, until 21 Aug 2019 02:35.
- **Wrong settings in alternating runs:**
  - EN688 10Aug2022(3)
  - EN706, all of 08Aug2023 (no alternation)
  - AR87A 3Apr2025(2) (111 extra samples on Pro settings; "ARA" file names, 3–6 Apr)
  - AR95 AltPro(28–31)
  - AR98A 31Oct(3) and 6Nov(13)
- **Bad or missing data:**
  - AR29 late cruise: spread fluorescence and air from 23–26 Apr 2018; about 70 blank files deleted.
  - AR31A: software froze on 23 Oct 2018 and part of a transect was missed.
  - EN655 ended early: samples with no cells, then sample-pump failure (28 Jul 2020).
  - EN688: overwrite errors and garbled files. This led to the `QC_dataintegrity` check.
  - EN695: samples 2–17 of 16 Jan 2023 not recorded.
  - EN715 4May2024(2): empty, because the tubing clamp was still on.
  - EN720 06Sep2024_winter(12): mostly bad. The instrument overheated and gave GreenLaser error 0108.
- **Optics:**
  - EN627 scatter degraded.
  - AR70B blue laser misaligned.
- **Particle-associated signals:** EN661 has a laser-delay split of Syn into blue- and green-only events (LASER2DELAY 305).
- **Fluidics and underway supply:**
  - air in the underway line: AR29, AR34B, AE2426
  - sheath ran out: TN368, EN668, AR61B, HRS2303, EN720
  - syringe leak: AR78
  - underway plumbing changed to bypass the debubbler on AR95 (19 Oct 2025)

### 7. Conflicts between sources

- **TN368 OD2.**
  - The settings sheet, the processing input table, `attune_config`, the lab log and the calibration summary all say OD2 was on SSC.
  - `process_wrapper_2024.m` and `2024b.m` default TN368 to GL1. The new `process_cruise.m` takes OD2 placement from `attune_config`.
- **EN685 beads.** The settings sheet says FCB bead files were not saved, but the lab log records bead runs on 13 and 17 Jul 2022.
- **Hellmanex.** The settings sheet says cleaning with Hellmanex began at AR34. The log shows first use at EN627 (Feb 2019).
- **EN655 voltages.** The settings sheet lists only BL2 360 → 380; the log says GL2 changed too.
- **HRS2303 GL2 threshold.** The log says "200 to 100" and then "1.0 x 1000". The FCS headers show 1000 and 2000.
- **Jul 2024 visit.** The lab log says "the whole PZT assembly including the flow cell" was replaced on 2 and 9 Jul. The service report (2–8 Jul) lists the interrogation assembly and tube lifter, plus laser-path, blocker-bar and sled adjustments.
- **24 Jan 2019 visit.** The lab log records only "a thorough cleaning". The service report shows a Contrad soak of the flow cell and a blocker-bar adjustment, one week before EN627.

---

## Part 2. Per-cruise reference

**Column guide.**
- **Dates** are the first and last underway sample files in `FCSfileinfo`.
- **Era** follows Part 1, section 1: A = OD2 on SSC, v2019 filters; B = no OD2; C = OD2 on SSC again; D = current layout (OD2 on GL1); E = AR104 settings.
- **hv** is the standard sample voltages, in the order FSC / SSC / GL1 scatter · BL3 · PE channel. The PE channel is GL1 575/36 in eras A–C and GL2 575/36 from D.
- **Trigger** is the main setting, threshold in counts. "PE" means the PE channel of that era.
- **Pro files** counts files with Pro-measuring settings (`pro_measured`).
- **Beads** gives the reference type and number of bead runs. PT = transferred from performance-test runs; FCB = 0.5/1/6 µm mix.
- Voltages and triggers come from the FCS headers (Nov 2026). The notes come from the logs and service reports.

| Cruise | Dates | Era | hv FSC/SSC/GL1 · BL3 · PE | Trigger | Pro files | Beads | Notes |
|---|---|---|---|---|---|---|---|
| EN608 | 2018-01-30 – 02-05 | A | 160/260/– · 300 · 300 | BL3 100 OR PE 200 | 0 | PT, 12 | Clock on Eastern time at start; first mouse-macro automation; 80 µm mesh; pulse-width threshold 1 |
| AR28B | 2018-03-29 – 04-13 | A | 160/260/– · 300 · 300 | BL3 100 OR PE 200 | 0 | PT, 8 | Software closed itself after ~200 samples (macro capped at 200) |
| AR29 | 2018-04-16 – 04-29 | A | 160/260/– · 300 · 300 | BL3 100 OR PE 200 | 0 | PT, 16 | Spread fluorescence and air in underway line from 23–26 Apr; ~70 blank files deleted |
| EN617 | 2018-07-20 – 07-25 | A | 160/260/– · 300→340 · 300→340 | BL3 100 OR PE 200 | 0 | PT, 7 | Syn clipped at the threshold → PE and BL3 raised; water leak onto computer 22 Jul |
| AR31A | 2018-10-20 – 10-23 | A | 160/260/– · 340 · 340 | BL3 100 OR PE 200 | 0 | PT, 3 | Pulse-width threshold 0.1 from here; software froze 23 Oct (part of transect missed); green laser trouble after |
| EN627 | 2019-02-01 – 02-06 | B | 160/200/– · 320 · 340 | BL3 100 OR PE 200 | 0 | PT, 9 | OD2 removed from SSC 1 Feb; **scatter not correctable (degraded optics)**, one week after a Contrad flow-cell soak and blocker-bar adjustment (24 Jan); file times local; first Hellmanex |
| AR34B | 2019-04-14 – 04-24 | B | 160/160/– · 320 · 340 | BL3 100 OR PE 200 | 0 | FCB, 19 | First FCB bead runs; times local; laser-power error and missed transit 19 Apr |
| RB1904 | 2019-05-12 – 05-25 | B | 160/160 (220 on 17–20 May)/– · 320 · 340 | BL3 100 OR PE 200 | 0 | FCB, 90 | SSC 220 for 1,489 files, then back to 160; times local; underway flow stopped 15–16 May |
| TN368 | 2019-07-05 – 07-18 | C | 280/260/– · 340 · 420 | BL3 300 OR PE 1100 | 0 | FCB, 84 | OD2 back on SSC; sheath ran out 6 Jul; `process_wrapper_2024` wrongly defaults OD2 to GL1 |
| EN644 | 2019-08-19 – 08-25 | D | 280/220/220 · 340 · 360 | BL3 300 OR GL2 300 (some 800–1100) | 0 | FCB, 22 | New filter layout; green laser replaced 6 Aug; **200 µl samples until 21 Aug 02:35**; ship power loss 19 Aug. `FCSfileinfo` also holds test files from 6 Aug |
| AR39B | 2019-10-05 – 10-09 | D | 280/220/220 · 340 · 360 | BL3 300 OR GL2 300 | 0 | FCB, 10 | Software crashes |
| EN649 | 2020-02-01 – 02-06 | D | 280/220/220 · 340 · 360 | BL3 300 OR GL2 300 | 0 | FCB, 40 | Oily underway water at start; gap 4 Feb |
| AR43 | 2020-03-11 – 03-16 | D | 280/220/220 · 340 · 360 | BL3 300 OR GL2 300 | 0 | FCB, 9 | Macro ran only once on 14 Mar; sheath nearly empty 15 Mar |
| EN655 | 2020-07-25 – 07-28 | D | 280/220/220 · 340 · 360→380 | BL3 300 OR GL2 300 | 0 | FCB, 19 | Beta pulse-width software; BL2 and GL2 360→380 on 26 Jul; ended early 28 Jul (sample pump, COVID) |
| EN657 | 2020-10-15 – 10-18 | D | 280/220/220 · 340 · 380 | BL3 300 OR GL2 300 | 0 | FCB, 19 | Only 43 sample files in `FCSfileinfo` (rest in `all_FCS`); no lab log |
| EN661 | 2021-02-03 – 02-08 | D | 280/220/220 · 340 · 380 | BL3 300 OR GL2 300 | 0 | FCB, 13 | No lab log; laser-delay split of Syn (LASER2DELAY 305); flow-cell contamination found 22 Feb |
| EN668 | 2021-07-16 – 07-21 | D | 280/220/220 · 340 · 380 | BL3 300 OR GL2 300 | 0 | FCB, 16 | First after new PC and beta FPGA; rinse-required stops; sheath ran out 20 Jul |
| AR61B | 2021-11-08 – 11-15 | D | 280/220/220 · 340 · 380 | BL3 300 OR GL2 300 | 0 | FCB, 28 | Software crashes; sheath ran out 15 Nov |
| AR62 | 2021-11-20 – 11-23 | D | 280/220/220 · 340 · 380 | BL3 300 OR GL2 300 | 0 | FCB, 5 | |
| AR63 | 2021-11-26 – 11-28 | D | 280/220/220 · 340 · 380 | BL3 300 OR GL2 300 | 0 | FCB, 2 | Student cruise |
| AT46 | 2022-02-16 – 02-21 | D | 280/220/220 · 340 · 380 | **BL3 500 OR GL2 2000** | 0 | FCB, 10 | Raised thresholds start; 92 test files at BL3 1000 OR GL2 300 (BL1 480, FSC 300); air pocket 17 Feb |
| AR66B | 2022-04-20 – 04-28 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 37 | FCB, 22 (18 good) | Small trigger tests; Euks/Pro below BL3 threshold offshore; dirty instrument |
| EN685 | 2022-07-11 – 07-18 | D | 280/220/220 · 340 · 380 | BL3 200 OR GL2 300 | 0 | FCB, 3 | Run by URI; settings sheet says no bead files, but 3 runs exist |
| EN687 | 2022-07-29 – 08-03 | D | 280/220/220 · 340 · 380 | BL3 300 OR GL2 300 | 0 | FCB, 9 | |
| EN688 | 2022-08-05 – 08-19 | D | 280/220/220 · 340 · 400 | Alternating SSC 300 AND GL2 500 / SSC 300 AND BL3 300 | 26 | FCB, 39 | First alternating cruise; overwrite errors and garbled files (→ `QC_dataintegrity`); 10Aug(3) wrong settings; `FCSfileinfo` also holds later lab runs |
| AR70B | 2022-11-17 – 11-22 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 0 | FCB, 5 | **Blue laser misaligned** (realigned 29 Nov); **clock 1 h 27 min and 5 days behind** |
| EN695 | 2023-01-10 – 01-16 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 0 | FCB, 11 | Samples 2–17 of 16 Jan not recorded; underway off in gale 13 Jan |
| HRS2303 | 2023-04-27 – 05-05 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 1000 / 2000 | 0 | FCB, 11 | GL2 threshold lowered for offshore Syn; sheath ran out 5 May |
| EN706 | 2023-08-07 – 08-12 | D | 280/220/220 · 340 · 400 | Alternating SSC AND GL2 500 / SSC AND BL3 300 | 0 | FCB, 9 | All of 08Aug ran BL3 500 OR GL2 2000 (520 files, no alternation) |
| AR77 | 2023-10-11 – 10-15 | D | 280/220/220 · 340 · 380 | BL3 200 OR GL2 1000 | 0 | FCB, 14 | Only 59 sample files in `FCSfileinfo` (rest in `all_FCS`); software 5.3 installed just before |
| AR78 | 2023-10-23 – 11-04 | D | 280/220/220 · 340 (Pro 370/390) · 380/400 | BL3 500 OR GL2 2000; alternating until 24 Oct; Pro SSC 100 AND BL3 300 on 1–2 Nov | 482 | FCB, 37 | Alternation dropped because of crashes; syringe leaked much of the cruise |
| AR79 | 2023-11-08 – 11-12 | D | 280/220/220 · 340 (Pro 380) · 380 | BL3 500 OR GL2 2000; Pro SSC 100 AND BL3 300 | 1,260 | FCB, 12 | Mold in sheath intake (found after) |
| EN712 | 2024-02-09 – 02-13 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 0 | FCB, 9 | Uneventful |
| AR82A | 2024-04-01 – 04-10 | D | 280/220/220 · 340 (Pro 380/390) · 380 | BL3 500 OR GL2 2000; Pro SSC 300/200 AND BL3 200/300 | 374 | FCB, 21 | Pro volume briefly 150 µl |
| EN715 | 2024-05-03 – 05-08 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 0 | FCB, 7 | 4May(2) empty (clamp on); ship power loss 4 May |
| EN720 | 2024-09-06 – 09-11 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 600–2000 | 0 | FCB, 13 | First after flow-cell (PZT) replacement; **overheating 7 Sep, 06Sep2024_winter(12) bad**; sheath ran out 9 Sep |
| AE2426 | 2024-11-06 – 11-11 | D | 280/220/220 · 340 (Pro 390) · 380 | BL3 500 OR GL2 2000; Pro SSC 300 AND BL3 200 | 465 | FCB, 10 | Heavy air in underway supply; PT failures; 09Nov(4) "messed up" |
| EN727 | 2025-01-24 – 01-28 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 0 | FCB, 8 | **GL1 bead 6/1 ratio steps down** (low-signal doubling); intake leak 24 Jan exp(1) |
| AR87A | 2025-03-28 – 04-08 | D | 280/220/220 · 340 (Pro 390) · 380 | BL3 500 OR GL2 2000/1000; Pro SSC 300 AND BL3 200 | 3,502 | FCB, 19 | First after flow-cell replacement (Feb 2025); file names "ARA" 3–6 Apr; 3Apr(2) 111 extra Pro-setting samples |
| AR87B | 2025-04-10 – 04-20 | D | 280/220/220 · 340 (Pro 390) · 380 | BL3 500 OR GL2 1000/2000; Pro SSC 300 AND BL3 200 | 3,239 | FCB, 41 | |
| AR88 | 2025-04-24 – 04-29 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 0 | FCB, 11 | Software closed itself 28 Apr |
| AR92 | 2025-08-14 – 08-19 | D | 280/220/220 · 340 (Pro 400) · 380 | BL3 500 OR GL2 500; Pro SSC 300 AND BL3 200 | 1,094 | FCB, 9 | Dead clock battery; Pro runs 150 µl |
| AR95 | 2025-10-18 – 10-23 | D | 280/220/220 · 340 (Pro 400) · 380 | BL3 500 OR GL2 500; Pro SSC 300 AND BL3 200 | 1,647 | FCB, 8 | Underway plumbing changed 19 Oct (debubbler bypassed); AltPro(28–31) all on Pro settings |
| AR96 | 2025-10-24 – 10-26 | D | 280/220/220 · 340 (Pro 400) · 380 | BL3 500 OR GL2 500; Pro SSC 300 AND BL3 200 | 736 | FCB, 2 | Student cruise |
| AR98A | 2025-10-30 – 11-09 | D | 280/220/220 · 340 (Pro 400) · 380 | BL3 500 OR GL2 2000; Pro SSC 500/400 AND BL3 400/300 | 3,306 | FCB, 13 | Three sensitivity levels; 31Oct(3) no alternation, 6Nov(13) all Pro; memory corruption 4 Nov |
| AR99 | 2026-01-12 – 01-18 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 0 | FCB, 17 | No lab log after Nov 2025 |
| HRS2601 | 2026-04-22 – 04-27 | D | 280/220/220 · 340 · 380 | BL3 500 OR GL2 2000 | 0 | FCB, 10 | After Feb 2026 preventive maintenance |
| HRS2609 | 2026-08-14 – 08-19 | D | 280/220/220 · 340 (Pro 370–420) · 380 | BL3 500 OR GL2 500; Pro SSC 300 AND BL3 200 | 1,988 | FCB, 6 | After sheath pump replacement (Jul 2026); FCS names lack leading NESLTER |
| AR104 | 2026-10-19 – 10-24 | E | 250/240/180 · 410 · 380 | SSC 900 AND BL3 1000/1200; Pro SSC 500 AND BL3 300 | 1,280 | FCB, 7 | New voltage settings; trigger on SSC AND chl (all but 5 test files) |

**Cruise folders without an underway series.** These hold preserved or discrete samples only, or runs that were not
logged as underway, and are not in the table:
- AR32, 20190303_AL, HB1907, AR38 and SG2105 are preserved samples only.
- The `FCSfileinfo` dates for SR1812 and DY131 are 2022 lab runs of preserved or thawed samples.

---

## Part 3. Dunaliella culture runs

*Dunaliella* ("Dun") is the culture most often used to set up and check the instrument: settings choices,
filter changes, coincidence tests, dilution series and size calibrations.

**How the runs were found** (Nov 2026):
- File names containing "dun" anywhere under `\\sosiknas1\lab_data\Attune`.
- Sample and experiment names inside every `.atx` experiment archive. An `.atx` is a zip holding each sample's FCS file plus `experiment.xml` with the sample names.
- The lab logs.

**What counts as a run:**
- Dates and times come from each FCS header (`$DATE`, `$BTIM`).
- Copies of the same run in backups and re-saved workspaces are counted once.
- The FCS files in the May 2026 size checks are gated exports (one file per gated population).
- Settings are the FCS header voltages, FSC / SSC / GL1 / BL3. GL1 is the PE channel before Aug 2019 and green side scatter after.

| Date | Experiment | Dun runs | hv FSC/SSC/GL1/BL3 | Purpose | Where |
|---|---|---|---|---|---|
| 2018-01-25 | Experiment | 2 (Dunaliella; Dun and Syn mix 1) | 220 / 220, 350 / 400 / 220 | First tests at installation | `cruise_data/20180131_EN608/FCS/setup`, `ATX_FSC_misc_experiments/ATX/Experiment.atx` |
| 2018-01-30 | Exp_LTER_test1 | 8 (Dun alone; Dun spiked into Iselin dock water at 100 and 200 µl/min) | 140–220 / 260–320 / 300–340 / 300 | Choosing the EN608 settings (Alexi): final FSC 160, SSC 260 | `cruise_data/20180131_EN608/FCS/setup`, `ATX_FSC_misc_experiments/FCS` |
| 2018-07-16 | Size calibration cruise template; Size calibration bead template | 2 (cruise method; bead method) | 160 / 260 / 300 / 300; 380 / 382 / 282 / 300 | July 2018 size calibration (with Coulter counter) | `Size_calibration_July2018` |
| 2018-12-07 | Size calibration cruise template | 1 (Heterocapsa + Dun mix) | 160 / 260 / 300 / 250 | Culture mixture; log: "only see one cluster", BL3 lowered to 250 | `Size_calibration_July2018/Het_dun_mix_Dec2018` |
| 2019-02-01 | New SSC filter check | 6 (Syn + Dun mix, original OD2 vs new filter; SSC 200, 210, 230, 250) | 160 / 200–250 / 340 / 340 | OD2 removed from SSC at the start of EN627 | `ATX_FSC_misc_experiments/ATX/New SSC filter check.atx`; FCS in `Backups_before_win10_install/Settings decisions` |
| 2019-02-26 to 02-28 | Dunaliella_Isochrisis continuous | 115 (114 consecutive on 26 Feb, one "two days later") | 160 / 200 / 340 / 320 | All-day continuous run of Dun + *I. galbana* on salt-water sheath (stability; samples 38–40 are DI water with carry-over) | `ATX_FSC_misc_experiments/Dunaliella_Isochrisis` |
| 2019-03-20 | size_cal_OD2_Mar2019; Size_cal_nuetral_Mar2019 | 2 (with OD2; without OD) | 160 / 260 / 340 / 320; 160 / 200 / 340 / 320 | March 2019 size calibration, used for the pre-Aug 2019 calibration | `Size_calibration_March2019/Attune` |
| 2019-04-12 | New SSC filter check (redo) | 6 (Syn + Dun, SSC 160, 160b, 170, 180, 200, 210) | 160 / 160–210 / 340 / 340 | Setting SSC voltage at the start of AR34B; SSC 160 chosen | `Backups_before_win10_install/Settings decisions` (no `.atx` found) |
| 2019-07-02 to 07-03 | Dock water June 2019(1) | 6 (Dun and Syn separately; mixes, with and without "avoid coincidence") | 160 / 240 / 350 / 350 | Coincidence troubleshooting | `ATX_FSC_misc_experiments/Coincidence troubleshooting` |
| 2019-07-04 | Thompson Mob 3 | 10 (Sample(19)–(28) – Dun) | 260 / 260 / 340 / 340 | TN368 mobilization: Dun on Attune and IFCB; counts matched without "exclude coincident events". Six were sent to Mike Ward on 23 Sep 2019 (`Dunaliella to Mike Ward.atx`, same files) | `cruise_data/20190705_TN368/ATX` |
| 2019-08-15 | Dockwater_filter_config_change | 1 | 280 / 260 / 360 / 340 | Checking the new green-laser filter layout before EN644 | `cruise_data/20190820_EN644/FCS` |
| 2019-10-21 | Syn and Dun coincidence | 7 (Dunaliella; mixture with Syn ×3; mixture with "exclude coincident events" ×3), plus 1 Syn alone | 280 / 220 / 220 / 340 | Coincidence examples for Mike Ward (Thermo) | `ATX_FSC_misc_experiments/Syn and Dun coincidence 102119.atx` |
| 2020-07-13 | dun_mira_dilution_test, dun_mira_dilution_test(1) | 24 (two dilution series, 1X to 1:32, duplicates), plus 2 filtered-seawater blanks | 280 / 220 / 220 / 340 | Mira's dilution tests; log: "concentrations look low" | `dilution_experiments_2020`, `MVCO_Attune_Runs_Mira` |
| 2020-07-16 | dun_mira_dilution_test(2) | 12 (dilution series 1X to 1:32, duplicates) | 280 / 220 / 220 / 340 | Repeat dilution test | `MVCO_Attune_Runs_Mira` |
| 2025-03-27 | (database backup) | 1 | 280 / 220 / 220 / 340 | AR87A mobilization: Dun run to check flow after the setup | `data_base_backup_Nov2025/Userdata/admin/…/Dun.fcs` |
| 2026-05-14 | size check | 2 (Dun; *I. galbana* + Dun + FCB beads) | 280 / 220 / 220 / 340 | Culture and bead size checks | `size_check_May2026_beads_culture`, `size_calibration_July2026/Attune/Gated_FCS/14May_TIso_Dun_beads` |
| 2026-05-15 | size check | 7 (Dun 1x, 0.5x and 0.1x, each twice; Dun + FCB and F8819 beads) | 280 / 220 / 220 / 340 | Concentration series and bead comparison | `size_check_May2026_beads_culture` |
| 2026-05-19 | cultures_may19 | 2 (underway settings; discrete settings) | 280 / 220 / 220 / 340 | Same culture at underway vs discrete-sample settings | `size_check_May2026_beads_culture/cultures_may19` |
| 2026-07-23 | size_calibration_July2026 | 3 (Dun, Dun(1); Dun(2) at BL3 250) | 280 / 220 / 220 / 340 (250) | July 2026 culture size calibration (with Coulter counter and Mastersizer) | `size_calibration_July2026/Attune` |
| 2026-09-21 | culture_check_21Sept2026 | 2 (Dun(1), Dun(3)) with an FCB bead check | 280 / 220 / 220 / 340 | Culture and bead check after the calibration | `size_calibration_July2026/Attune/ATX/culture_check_21Sept2026.atx` |

**Total: about 219 Dun runs on 20 days**, 115 of them in the Feb 2019 continuous run.

**Log mentions without runs:**
- 25 Sep 2019: planned Syn and Dun coincidence runs, not done because the beads were dirty. They were done on 21 Oct.
- 28 Apr 2021: Dun was considered for a counting test, but *Isochrysis* was used instead (Dun is off scale in scatter).
