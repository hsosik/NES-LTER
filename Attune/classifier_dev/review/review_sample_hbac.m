function q = review_sample_hbac(cls, info, cfg, outdir, fcsname, t, fcshdr, R)
% function q = review_sample_hbac(cls, info, cfg, outdir, fcsname, t, fcshdr, R)
% Review output for one SYBR-stained bacteria run: writes the frame to
% outdir/frames/<name>.png (review_frame_hbac) and returns the QC struct
% (review_qc_hbac) for review_viewer. Called after assign_class_hbac inside the
% process_cruise classify loop (preserved samples).
%
% Nov 2026, preserved-sample review (SYBR runs)

if nargin < 8 || isempty(R), R = review_config_hbac(); end
frame = ['frames/' regexprep(fcsname, '\.fcs$', ['.' R.format])];
nrow = review_frame_hbac(cls, info, fullfile(outdir, frame), R);
q = review_qc_hbac(cls, info, cfg);
q.file = fcsname; q.time = t; q.frame = frame; q.nrow = nrow;
q.pro_measured = 1;
q.trigger = sprintf('%s %g', fcshdr.tr1_par, fcshdr.tr1_level);
end
