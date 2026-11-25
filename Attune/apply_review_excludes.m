function n = apply_review_excludes(cruise, opt)
% function n = apply_review_excludes(cruise, Name=Value)
% Set the Review_exclude column from the cruise's class review decisions in the tables
% already built, without rebuilding them: AttuneTable.mat, AttuneTable_uw_match.mat and
% AttuneVolTable.mat in <cruise>\bead_calibrated\ (whichever exist). Run it after
% reviewing a cruise (review_decisions.csv saved by the review page), then regenerate
% the EDI table (generate_attune_table_EDI_one_cruise), which leaves out excluded rows.
% generate_attune_table sets the column too, so a cruise reprocessed after its review
% needs no extra step.
%
% cruise   cruise name ("AR92") or a bead_calibrated folder
% Options
%   Decisions  review decisions file (default <bead_calibrated>\review_decisions.csv)
%   Root       folder of cruise folders (default from attune_cruise_settings)
% Returns n, a struct: number excluded in each table, and decisions marked exclude.
% e.g. apply_review_excludes("AR92")
%
% Nov 2026

arguments
    cruise {mustBeTextScalar}
    opt.Decisions (1,1) string = ""
    opt.Root (1,1) string = ""
end
cruise = char(cruise);
if contains(cruise, {'\' '/'})
    bp = cruise;
else
    s = attune_cruise_settings(cruise, char(opt.Root));
    bp = fullfile(s.basepath, 'bead_calibrated');
end
decfile = char(opt.Decisions); if isempty(decfile), decfile = fullfile(bp, 'review_decisions.csv'); end
if ~isfile(decfile), warning('apply_review_excludes:nodecisions', 'no review decisions file: %s (all rows set to 0)', decfile); end

tabs = {'AttuneTable.mat' 'AttuneTable'; 'AttuneTable_uw_match.mat' 'Attune_uw_match'; 'AttuneVolTable.mat' 'AttuneVolTable'};
n = struct();
for k = 1:size(tabs, 1)
    f = fullfile(bp, tabs{k, 1});
    if ~isfile(f), continue, end
    S = load(f);
    if ~isfield(S, tabs{k, 2}), warning('apply_review_excludes:var', '%s has no %s', f, tabs{k, 2}); continue, end
    T = S.(tabs{k, 2});
    [ex, ndec] = review_exclude_column(T.Filename, decfile);
    T.Review_exclude = ex;
    S.(tabs{k, 2}) = T;
    save(f, '-struct', 'S')
    n.(tabs{k, 2}) = sum(ex);
    n.decisions_exclude = ndec;
    fprintf('%s: Review_exclude = 1 for %d of %d rows\n', f, sum(ex), height(T));
end
if isempty(fieldnames(n)), warning('apply_review_excludes:notables', 'no Attune tables in %s', bp); end
end
