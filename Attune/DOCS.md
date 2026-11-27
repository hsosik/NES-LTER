# Attune documentation index

The written documentation for the Attune processing and the v2 classifier redesign. When you add a document, add a line here.

## Instrument and calibration

- [ATTUNE_HISTORY.md](ATTUNE_HISTORY.md): instrument and processing history, 2018–2026. Narrative of optical configurations, voltages and service; a per-cruise reference table; every *Dunaliella* culture run.
- [classifier_dev/pt_transfer/CALIBRATION_SUMMARY.md](classifier_dev/pt_transfer/CALIBRATION_SUMMARY.md): how side scatter is turned into cell volume for every cruise, why, and how certain it is.
- [classifier_dev/pt_transfer/README.md](classifier_dev/pt_transfer/README.md): performance tests used as a bead transfer standard for the early cruises without FCB bead runs; details behind the calibration summary.

## Classifier

- [classifier_dev/README.md](classifier_dev/README.md): the v2 classifier (`assign_class_v2`). Scripts, class rules and a record of every rule change and its effects.
- [classifier_dev/reports/class6/index.html](classifier_dev/reports/class6/index.html): report on class 6, a preservation shift in large Euks. Figures from the code as of Nov 2026; the rule itself is kept in the classifier README.

## Reviewing classifications

- [classifier_dev/review/CRUISE_REVIEW_GUIDE.md](classifier_dev/review/CRUISE_REVIEW_GUIDE.md): how to use the cruise class review and bead review pages. The pages themselves are built on the NAS by `process_cruise`.

## Bead statistics

- [classifier_dev/BEADSTATS_HANDOFF.md](classifier_dev/BEADSTATS_HANDOFF.md): how bead statistics are computed and saved, written as a starting point for simplifying them.
- [classifier_dev/bead_design/README.md](classifier_dev/bead_design/README.md): the scripts used to design and check the bead gating. Not needed by the pipeline.
