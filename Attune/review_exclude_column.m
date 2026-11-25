function [ex, ndec] = review_exclude_column(filenames, decfile)
% function [ex, ndec] = review_exclude_column(filenames, decfile)
% Review_exclude for an Attune table: 1 for files marked "exclude" in the class review
% (review_decisions.csv, written by the cruise review page), 0 otherwise. The final
% post-review flag: AttuneTable keeps every row, and the EDI tables leave out rows
% with Review_exclude = 1 (generate_attune_table_EDI_*).
%   filenames  FCS file names (the table's Filename column)
%   decfile    path to review_decisions.csv (normally <cruise>\bead_calibrated\)
% ex is a double column (0/1) the length of filenames; ndec counts the exclude
% decisions in the file (some may name files not in the table).
% See classifier_dev\review\CRUISE_REVIEW_GUIDE.md.
%
% Nov 2026

filenames = cellstr(filenames);
ex = zeros(numel(filenames), 1); ndec = 0;
if ~isfile(decfile), return, end
if ~exist('review_decisions', 'file')
    addpath(fullfile(fileparts(mfilename('fullpath')), 'classifier_dev', 'review'))
end
T = review_decisions(decfile);
bad = T.file(strcmp(T.status, 'exclude'));
ndec = numel(bad);
ex(ismember(filenames, bad)) = 1;
end
