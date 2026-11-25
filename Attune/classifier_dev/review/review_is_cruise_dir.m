function tf = review_is_cruise_dir(outdir)
% function tf = review_is_cruise_dir(outdir)
% True if outdir is a cruise's own classification review folder (...\bead_calibrated\review,
% underway or preserved). That review must match the class files, so only
% process_cruise (classify step, which writes both) may write it; review_cruise and
% review_prototype, which classify without writing class files, refuse it.
%
% Nov 2026, classification review
tf = ~isempty(regexpi(char(outdir), '[\\/]bead_calibrated[\\/]review[\\/]*$', 'once'));
end
