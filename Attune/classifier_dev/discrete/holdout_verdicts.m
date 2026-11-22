function T = holdout_verdicts(src)
% function T = holdout_verdicts(src)
% Read the verdicts saved by the held-out review viewer (holdout_review).
%   src  path to holdout_verdicts.csv (default output/discrete/holdout/review/holdout_verdicts.csv)
% Returns a table, one row per run: file, set, runtype, class_affected,
% verdict_blind (A_wrong/B_wrong/ambiguous/ok, blind mode only), classifier_row (A/B),
% verdict (classifier_wrong/hand_wrong/ambiguous/ok), verdict_revised (after the
% reveal, if changed), note, reviewer, saved, build.
% Use verdict_revised where set, else verdict, as the final call:
%   T = holdout_verdicts(); v = T.verdict; k = T.verdict_revised ~= ""; v(k) = T.verdict_revised(k);
%
% Nov 2026, discrete-sample classifier development

if nargin < 1 || isempty(src)
    here = fileparts(mfilename('fullpath'));
    src = fullfile(fileparts(here), 'output', 'discrete', 'holdout', 'review', 'holdout_verdicts.csv');
end
cols = {'file' 'set' 'runtype' 'class_affected' 'verdict_blind' 'classifier_row' 'verdict' 'verdict_revised' 'note' 'reviewer' 'saved' 'build'};
T = cell2table(cell(0, numel(cols)), 'VariableNames', cols);
if ~exist(src, 'file'), return, end
o = detectImportOptions(src, 'Delimiter', ',', 'Encoding', 'UTF-8', 'TextType', 'string', 'VariableNamingRule', 'preserve');
o = setvartype(o, 'string');
R = readtable(src, o);
R.Properties.VariableNames = regexprep(lower(R.Properties.VariableNames), '^\W+', '');
for c = cols
    if ~ismember(c{1}, R.Properties.VariableNames), R.(c{1}) = repmat("", height(R), 1); end
    R.(c{1}) = fillmissing(R.(c{1}), 'constant', "");
end
T = R(:, cols);
end
