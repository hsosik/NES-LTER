# Cruise class review: user guide

The class review page lets you step through every sample file of a cruise and look at
how `assign_class_v2` classified it. Each file is shown as a set of cytograms with
its class counts and QC numbers, and files that look unlike their neighbours are
flagged as **unusual**. You can then record a decision for each file.

There is a separate page, the **bead review** (`bead_calibrated\bead_review\`), for the
bead runs or performance tests. Those are covered at the end.

## Building the page

The page is built automatically by `process_cruise` during the classify step, unless
you pass `Review=false`.

```matlab
process_cruise("EN727")                        % core steps, review included
process_cruise("EN727", Steps="classify")      % reclassify and rebuild the review
review_cruise('EN727')                         % review only: classifies, writes no class files
review_cruise('EN727', struct('files', 200))   % quick look: 200 files spread over the cruise
```

The output goes in `<cruise>\bead_calibrated\review\`:

| File | Contents |
|---|---|
| `index.html` | The page |
| `review_data.js` | Per-file QC numbers and earlier decisions |
| `frames\` | One cytogram image per file |
| `review_qc.mat` | The same QC numbers, for MATLAB |

Preserved (discrete) samples have their own page in
`<cruise>\preserved\bead_calibrated\review\`. SYBR runs have one in `...\review\hbac\`.

The page's HTML is `review_viewer_template.html` with a version stamp added. A template
change therefore reaches an existing page by copying the template over its `index.html`,
with no reclassification.

## Opening the page and saving decisions

Use **Chrome or Edge**. Other browsers can show the page but cannot save to a file.

Saving needs Chrome's folder-access feature, which Chrome blocks for pages opened from
a `\\server\share` path. Open the page in one of these ways:

- **From a mapped drive.** Map `\\sosiknas1\Lab_data` as a drive letter, say `L:`, and
  open `L:\Attune\cruise_data\<cruise>\bead_calibrated\review\index.html`.
- **Through a local web server.** For example, run this in the review folder:

  ```
  python -m http.server 8000
  ```

  Then open `http://127.0.0.1:8000/index.html`.

If you open it from a network path, the save line says so.

Saving the first time:

1. Mark a file. The page then asks for a folder, or you can click **Choose folder** on
   the save line.
2. Pick the cruise's `bead_calibrated` folder, the one above `review`.
3. Allow editing when Chrome asks.

Decisions save to `bead_calibrated\review_decisions.csv` as you work. Next time,
**Reconnect** restores saving with one click.

If saving fails:
- The red text on the save line gives the reason.
- **CSV** downloads the decisions, and you copy the file into `bead_calibrated`
  yourself.
- Until they are saved, decisions are kept in the browser.

Decisions from several reviewers or browsers are merged by file, and the newest change
wins.

## The layout

- **Header.** It has the following controls:
  - ◀ ▶ step through files, and **◀ unusual / unusual ▶** jump to flagged files.
  - The slider moves anywhere in the cruise.
  - The badge shows "ok" or the number of reasons a file is unusual, plus its review
    status.
  - **Skip reviewed** makes the unusual jumps pass over files that already have a
    status.
  - The run selector shows all runs, Pro-setting runs only (SSC trigger, low
    threshold) or non-Pro runs only. Your choice is remembered per cruise.
- **Class counts strip.** It shows the counts of each class through the cruise, on a
  log scale.
  - Red ticks mark unusual files. Ticks above them show the review status.
  - Click the strip to jump to a file. The tick boxes show or hide classes, and the
    x axis can be sample order or time.
- **QC metric strip.** It shows one metric through the cruise, with the band its
  neighbours allow (see below). Pick the metric from the menu.
- **Frame.** These are the cytograms for the current file, with fixed axes in
  bead-relative log10 units, so frames are comparable across files and cruises.
  - Size vs chl, and chl vs PE.
  - Chl:size vs GL2:GL3, which shows the Pro colour.
  - Chl pulse shape, which shows the Pro wide-pulse cut and trigger-floor debris.
  - The Syn type view, PE (GL2) vs PE (BL2), with Syn coloured by type A/B/C/D (in order of BL2 − GL2; A = lowest PUB:PEB). It shows brightness, which explains untyped events; the type edges are the dashed diagonals.
  - The Syn type A/B view, BL2 − GL2 (the excitation ratio the types are cut on) vs PE:chl, with the same colours. The type edges are the dashed vertical lines. The ratio of dim Syn is noisy, so for Syn whose BL2 or GL2 is less than `syn_typeA_min_above_floor` (2) above its floor, the type-A edge moves left by `syn_typeA_dim_slope` (0.4) per log10 of shortfall; dim events between that edge and the first edge stay type B. The type-A boundary is drawn in the type-A colour in the GL2 vs BL2 view. Pages built before 26 Nov 2026 have 10 panels and no type A/B view.
  - Since the per-file mixture (26 Nov 2026), most files are typed by a fit to the file, not by fixed edges. Both type panels then draw the positions of the types the file kept, as solid lines in the type colours, and no edges. Events whose type is not clear (code 5, ambiguous) are drawn dark grey with the untyped ones. Files with fewer than 500 typed Syn still use the edges and type-A boundary above. The side panel's "Syn typing" line says which, with the kept types and the ambiguous count.
  - Unassigned events are grey. Flagged unassigned events are coloured by flag.
  - Lines mark the trigger edges and the rules this file used.
- **Side panel.** It shows:
  - the review buttons, note, reviewer name and save line;
  - class counts (old and v2 when old classes exist);
  - why the file is unusual;
  - the QC table.
- **Unusual list.** Below the frame is every unusual file in the current view, with
  its reasons and status. Click a row to open that file.

**Keys:**

| Key | Action |
|---|---|
| ← → | Step one file |
| Shift + ← → | Step 10 files |
| n / p | Next or previous unusual file |
| Home / End | First or last file |
| o | Mark reviewed – OK |
| x | Mark problem |
| e | Mark exclude |
| ? | Mark needs another look |
| m | Type a note (Esc to leave the note box) |

Pressing a status key a second time clears that status. A URL ending in `#123` opens
file 123.

## How "unusual" is decided

Every metric is compared with the file's neighbours:

- **Neighbours** are up to 10 files on each side that have the same Pro setting and
  are within 6 hours of the file.
  - With fewer than 4 neighbours, a file gets no score for that metric. This happens
    at the start and end of the cruise and after long gaps.
- **Score:** z = (value − neighbours' median) / scale.
  - The scale is the neighbours' spread (MAD × 1.4826).
  - It is never below a minimum set for each metric. This stops tiny natural
    differences being flagged in very steady water.
- **Flagged** when |z| exceeds the threshold, 4 by default.
- **Two metrics also have fixed limits:** a catch-all Euk share or an untyped-Syn
  share above 0.25 is always flagged.
- **Gap changes:** a Euk chl gap or Syn PE gap that is found in a file but in fewer
  than 30% of its neighbours, or missing but found in more than 70%, is also a reason.

The QC metric strip draws this band, median ± z·scale. Points outside it are flagged.

**Controls.** Detection runs in the page, so changes apply immediately and don't change
the page on the NAS.

- **flag at |z| >**: raise it, to 6 say, to keep only the large jumps.
- **count gap changes (and below-gap / chl-edge pile flags) as unusual**: off by
  default.
  The gaps often switch on and off between neighbouring files with little effect on
  Syn or Euk counts, so these reasons are shown only as grey notes unless this is
  ticked. Metrics marked \* in the list below behave this way.
- **metrics that count as unusual**: one checkbox per metric, with **all** and
  **none** buttons.
  - A file flagged only by unticked metrics is no longer unusual. The reason is still
    shown, greyed, as a note.
  - The choice is remembered per cruise in this browser.

To find the main cause on a cruise, open a few unusual files and read their side-panel
reasons. Then decide whether that metric is telling you something or just tracking
natural variation, and untick it if it isn't.

## The QC metrics

Count metrics use log10(count + 1). Cluster and reference positions are in
bead-relative log10 units: size from SSC relative to the 1 µm bead, chl from BL3, and
PE from GL2. A cluster position is the median of the class's events, and is blank when
the class has too few events (50 for Syn and Euk, 20 for each crypto class).

| Metric | Min. scale | What it is | When it fires |
|---|---|---|---|
| Euk, Syn, Pro, unassigned counts | 0.12 (Pro 0.15) | Events in each class | A real patch or front; or a class taking or losing a population (check the frame) |
| crypto dim + bright | 0.25 | Dim and bright crypto together | As above; crypto counts are small and noisy |
| Euk with Syn, Syn on particle | 0.35 | Coincident Euk + Syn events, and Syn riding on a larger particle | Usually small numbers; worth a look when Syn counts also jump |
| flag: Syn-pigment particle | 0.3 | Unassigned events with Syn pigment and low chl (degraded cells, grazers, detritus) | Particle-rich water, or Syn boundary trouble |
| flag: below Euk chl gap \* | 0.3 | Unassigned events below the Euk chl gap | Follows the Euk chl gap coming and going |
| flag: chl-edge noise pile \* | 0.35 | Unassigned trigger-noise pile at the chl edge (flag 3; capped files, pre-2019). Pages built before 25 Nov 2026 call it "Pro tail" | Trigger or threshold settings changed; on EN727, set only in files without a Euk chl gap |
| flag: Pro lower bound | 0.3 | Pro-like events (no PE, chl at or above the chl trigger, no larger than Syn) in files run without Pro settings (flag 4); a lower bound for Pro | A switch between Pro and non-Pro settings, or Pro abundance changing |
| flag: Euk below chl:size cut | 0.3 | Euk-sized events left out by the file's Euk chl:size cut (flag 5) | The per-file cut moved; check the Euk population in the size vs chl panel |
| catch-all Euk share of Euk | 0.04; always > 0.25 | Share of Euks added by the catch-all, which gives cell-like leftover events the class of their nearest assigned neighbours | The main rules missed part of the Euk population |
| catch-all Syn share of Syn | 0.04 | The same for Syn | Syn spread beyond the reference rules |
| untyped Syn share | 0.04; always > 0.25 | Syn-bearing events whose PE spectrum can't be typed | Dim or noisy PE, often near the trigger |
| ambiguous-type Syn share | 0.05 | Syn-bearing events no kept type claims with posterior 0.8 (per-file mixture) | Dim Syn whose ratio noise spans two types; high in summer |
| Syn cluster size, chl, PE, size spread (IQR) | 0.03, 0.08, 0.08, 0.015 | Position and size spread of the final Syn class | Syn class taking a different population (e.g. small Euks or PE particles) |
| Euk cluster size, chl | 0.06, 0.07 | Position of the final Euk class | Euk class shifted, e.g. a bloom of one size, or debris taken in |
| dim / bright crypto cluster size, chl, PE | 0.1, 0.12, 0.12 | Position of each crypto class | Crypto class shifted or split differently |
| Syn ref PE:chl, GL2:GL3, chl, size | 0.06, 0.05, 0.06, 0.06 | The Syn reference cluster the Syn rules are built from | Reference picked from another population; check the Syn type view |
| Syn PE floor | 0.08 | Lowest PE accepted as Syn, when set statistically (blank when the PE gap sets it) | Noise level or PE trigger changed |
| Euk chl gap \* | 0.08 | Chl of the gap below the Euk population (blank = no gap found) | Gap moved, or found / missing compared with neighbours |
| Syn PE gap \* | 0.08 | PE of the valley between Syn and no-PE events (blank = none) | As above |
| agreement with old classes | 0.03 | Share of events with the same class in the old (pre-v2) classification; only when old class files exist | v2 and the old gates diverge |

Flags 4 and 5 are blank on pages built before 25 Nov 2026. Flags 6–9 (split-Syn parts)
are not shown.

SYBR (hbac) pages use their own list: bacteria, HNA, LNA, Syn, Pro, Euk and unassigned
counts, HNA fraction, the HNA/LNA cut and valley dip (dip is \*), chl and PE spillover
offsets, and core events.

The full definitions are in `review_viewer_template.html` (the `METRICS_PH` and
`METRICS_HB` lists) and `review_qc.m`.

## Review statuses and what happens to them

| Status | Meaning |
|---|---|
| reviewed – OK | Classification looks right |
| problem | Something is wrong; describe it in the note |
| exclude | The file should not be used (a bad run or instrument trouble) |
| needs another look | Undecided; come back to it |

Each row of `review_decisions.csv` holds: file, status, note, reviewer, saved (UTC) and
build (the classification run you were looking at). To read the decisions in MATLAB:

```matlab
T = review_decisions('EN727');
bad = T.file(strcmp(T.status, 'exclude'));
```

**Exclude is the final post-review flag.** `AttuneTable` keeps every row and has a
`Review_exclude` column, set to 1 for files marked exclude. The EDI tables
(`generate_attune_table_EDI_*`) leave those rows out. The other statuses and the notes
are records only.

- `generate_attune_table` (the table step) sets the column from `review_decisions.csv`.
- Reviews usually come after the tables are built. After reviewing a cruise, update the
  column in the existing tables (`AttuneTable`, `AttuneTable_uw_match` and
  `AttuneVolTable`), then regenerate the EDI table:

  ```matlab
  apply_review_excludes("AR92")
  generate_attune_table_EDI_one_cruise('\\sosiknas1\Lab_data\Attune\cruise_data\20250814_AR92\bead_calibrated\')
  ```

Rebuilding the page keeps earlier decisions; each one shows the build it was made on.
`archive_cruise_outputs` carries the decisions into a new `bead_calibrated` only with
`CarryReview=true`.

## The bead review

`bead_calibrated\bead_review\index.html` shows the cruise's FCB bead runs, or, for
cruises without FCB beads (EN608 to EN627), the performance tests transferred to
FCB-equivalent peaks. It shows:

- bead peaks through the cruise, as a voltage-corrected offset or as detector voltage;
- each run's histograms;
- QC flags.

It opens and saves the same way as the class review. Decisions go to
`bead_calibrated\bead_review_decisions.csv`.

Unlike class review decisions, **bead exclusions are applied**. `load_beadstat` skips
bead runs or tests marked exclude, so the classify and volume steps no longer use them.
Rerun from classify after changing an exclusion:

```matlab
process_cruise("AR29", Steps=["classify" "volume" "table" "metadata"])
```

`archive_cruise_outputs` always carries the bead decisions forward.
