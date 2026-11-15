function S = summarize_early_layout(sub)
% function S = summarize_early_layout(sub)
% sub - folder under output/ (default 'early_layout'; 'early_layout_robust' = with the
%       trigger-noise pile and rare-Syn handling)
% One row per cruise from output/early_layout/agreement_<CRUISE>_layout.mat
% (validate_early_layout): median agreement for v2 with the attune_config
% channel map and for the pre-change v2 (raw channel names), mean counts.
% Nov 2026
if nargin < 1, sub = 'early_layout'; end
d = fullfile(fileparts(mfilename('fullpath')), 'output', sub);
f = dir(fullfile(d, 'agreement_*_layout.mat'));
S = table;
for i = 1:numel(f)
    V = load(fullfile(d, f(i).name)).V;
    if isempty(V), continue, end
    r = table(V.cruise(1), height(V), 'VariableNames', {'cruise' 'nfiles'});
    for v = {'v2_agreement' 'v2_agreement_pe_merged' 'orig_agreement' 'orig_agreement_pe_merged'}
        if ismember(v{1}, V.Properties.VariableNames), r.(['med_' v{1}]) = median(V.(v{1}), 'omitnan'); r.(['min_' v{1}]) = min(V.(v{1})); else, r.(['med_' v{1}]) = NaN; r.(['min_' v{1}]) = NaN; end
    end
    if ismember('v2_error', V.Properties.VariableNames), r.n_v2_error = sum(strlength(V.v2_error) > 0); r.v2_error = V.v2_error(find(strlength(V.v2_error) > 0, 1)); if isempty(r.v2_error), r.v2_error = ""; end
    else, r.n_v2_error = 0; r.v2_error = ""; end
    c = {'v2_n1' 'v2_n2' 'v2_n3' 'v2_n4' 'v2_n5' 'v2_n7' 'v2_n8'};
    c = c(ismember(c, V.Properties.VariableNames));
    m = mean(V{:, c}, 1, 'omitnan'); r.v2_mean_n_1_2_3_4_5_7_8 = round(m);
    if ismember('old_n', V.Properties.VariableNames), r.old_mean_n_1_2_3_4_7 = round(mean(V.old_n, 1, 'omitnan')); else, r.old_mean_n_1_2_3_4_7 = nan(1, 5); end
    S = [S; r]; %#ok<AGROW>
end
save(fullfile(d, 'summary.mat'), 'S')
end
