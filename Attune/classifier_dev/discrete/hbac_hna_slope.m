function T = hbac_hna_slope(survey, slopes)
% function T = hbac_hna_slope(survey, slopes)
% Does a tilted HNA/LNA cut separate better than a flat SYBR (BL1) cut?
% For each slope a, the bacteria (class 9) of each survey run are projected
% on y = BL1 - a * (SSC + 2) (SSC centred near the LNA cluster, so y stays
% on the BL1 scale) and the deepest valley in a window around -1.27 is found
% as in hbac_hna_split. Prints the median dip and the fraction of runs with
% a dip below 0.9 / 0.7 per slope.
%
% Nov 2026, discrete-sample classifier development (HNA/LNA split)

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 2, slopes = [-0.2 0 0.1 0.2 0.3 0.4 0.5 0.6 0.8]; end
load(fullfile(od, [survey '.mat']), 'R', 'E')
p = struct('win', [-1.45 -1.10], 'flank', 0.35, 'bw', 0.04, 'min_n', 2000, 'max_dip', 0.9, 'default_cut', -1.27);
D = nan(height(R), numel(slopes)); C = D;
for r = 1:height(R)
    e = E(E.run == r & E.class == 9, :);
    if height(e) < 2000, continue, end
    for k = 1:numel(slopes)
        p.slope = slopes(k);
        [~, s] = hbac_hna_split(double(e.BL1), double(e.SSC), true(height(e), 1), p);
        D(r, k) = s.dip; C(r, k) = s.valley;
    end
end
T = table(slopes(:), median(D, 'omitnan')', mean(D < 0.9)', mean(D < 0.7)', median(C, 'omitnan')', ...
    'VariableNames', {'slope' 'median_dip' 'frac_dip_lt_0p9' 'frac_dip_lt_0p7' 'median_valley'});
disp(T)
[~, kb] = min(D, [], 2, 'omitnan');
fprintf('best slope per run: '), disp(tabulate(slopes(kb(all(isfinite(D), 2))))')
end
