function T = review_decisions(src, keep_cleared)
% function T = review_decisions(src, keep_cleared)
% Read the classification review decisions saved by the review viewer.
%   src  path to review_decisions.csv, or a cruise name (reads
%        <cruise>\bead_calibrated\review_decisions.csv via attune_cruise_settings)
%   keep_cleared  true: also return rows whose status and note were cleared in the
%        viewer (kept in the file so a cleared decision isn't revived when saves from
%        two browsers are merged); default false
% Returns a table with one row per file: file, status, note, reviewer, saved, build.
%   status: ok (reviewed, OK), problem, exclude, recheck (needs another look), or ''
%   saved: when the decision was last changed (UTC, ISO 8601)
%   build: which classification run (review build) the reviewer was looking at
% e.g. T = review_decisions('EN727'); bad = T.file(strcmp(T.status, 'exclude'));
%
% Nov 2026, classification review

if nargin < 2, keep_cleared = false; end
cols = {'file' 'status' 'note' 'reviewer' 'saved' 'build'};
src = char(src);
if ~contains(src, {'\' '/' '.csv'})
    s = attune_cruise_settings(src, '');
    src = fullfile(s.basepath, 'bead_calibrated', 'review_decisions.csv');
end
T = cell2table(cell(0, numel(cols)), 'VariableNames', cols);
if ~exist(src, 'file'), return, end
o = detectImportOptions(src, 'Delimiter', ',', 'Encoding', 'UTF-8', 'TextType', 'char', ...
    'VariableNamingRule', 'preserve', 'ReadVariableNames', true);
o = setvartype(o, 'char');
R = readtable(src, o);
if height(R) == 0, return, end
for c = cols
    if ismember(c{1}, R.Properties.VariableNames)
        v = R.(c{1}); v(cellfun(@(x) ~ischar(x), v)) = {''};
    else
        v = repmat({''}, height(R), 1);
    end
    R.(c{1}) = v;
end
T = R(:, cols);
T = T(~strcmp(T.file, ''), :);
if ~keep_cleared, T = T(~(strcmp(T.status, '') & strcmp(T.note, '')), :); end
end
