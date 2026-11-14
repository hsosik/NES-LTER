# Bead stats design and test scripts

These scripts were used to design and check the bead gating in `process_beads_PT_adjust_3` (Nov 2026). The pipeline does not need them.

- `run_bead_prototype.m`: runs `find_beads` on every bead-check file of every cruise and tabulates the results alongside the old `beadstat_2026` / `beadstat_2021` values. It writes `classifier_dev/output/bead_proto.mat`.
- `compare_bead_prototype.m`: QC pass rates and failure reasons per cruise; peak vs median; new vs old medians; the 1 µm peak vs the old `centend2` modes.
- `plot_bead_disagreements.m`: gating plots for runs where the new and old results disagree. It writes `classifier_dev/output/bead_proto_plots/`.

They use the functions in `classifier_dev` (`find_beads`, `bead_signature`, `bead_peak`, `attune_config`, `plot_bead_gates`), so run them with `classifier_dev` on the path. The findings are written up in `../BEADSTATS_HANDOFF.md`.
