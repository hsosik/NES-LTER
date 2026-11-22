function old = load_old_class(P, fn)
% function old = load_old_class(P, fn)
% Old (pre-v2) class codes, first column, for FCS file fn, from P.old_class
% (P from cruise_paths). Stops if the file is a v2 class file (has class_meta,
% written by process_cruise): scoring v2 against v2 would look like perfect
% agreement. That happens if a cruise was reprocessed into bead_calibrated
% without being archived first.
%
% Nov 2026, classifier redesign (scoring after archiving)

f = fullfile(P.old_class, regexprep(fn, '\.fcs$', '.mat'));
if ismember('class_meta', who('-file', f))
    error('load_old_class:v2', ['%s is a v2 class file (class_meta), not the old classifier''s; ' ...
        'archive the cruise first (archive_cruise_outputs) or pass the archive Tag to cruise_paths'], f)
end
c = load(f, 'class');
old = c.class(:, 1);
end
