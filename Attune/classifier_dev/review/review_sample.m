function q = review_sample(fcsdat, fcshdr, cls, info, old, outdir, fcsname, t, pro_measured, R)
% function q = review_sample(fcsdat, fcshdr, cls, info, old, outdir, fcsname, t, pro_measured, R)
% Review output for one classified sample: writes the frame image to
% outdir/frames/<name>.png and returns the QC struct for review_viewer. Meant to be
% called right after assign_class_v2, inside a parfor (process_cruise step 3 or
% review_cruise). old: previous classes for the same events, or [] (production).
%
% Nov 2026, classification review

if nargin < 10 || isempty(R), R = review_config(); end
frame = ['frames/' regexprep(fcsname, '\.fcs$', ['.' R.format])];
nrow = review_frame(fcsdat, fcshdr, cls, info, old, fullfile(outdir, frame), R);
q = review_qc(cls, info, old);
q.file = fcsname; q.time = t; q.frame = frame; q.nrow = nrow;
q.pro_measured = pro_measured;
q.trigger = sprintf('%s %g', fcshdr.tr1_par, fcshdr.tr1_level);
end
