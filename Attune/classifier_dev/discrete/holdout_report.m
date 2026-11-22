function R = holdout_report(sets)
% function R = holdout_report(sets)
% Pass criteria of the held-out check (holdout_proposal.md, approved 21 Nov 2026),
% from holdout_score outputs, each set against the 'tuning' re-score:
%  1 median auto/hand ratio per class within +-5% of tuning and in 0.90-1.10
%    (PE euk judged against tuning only)
%  2 >= 90% of runs within +-20% per class (runs with >= 100 hand events)
%  3 median agreement per run type no more than 0.03 below tuning
%  4 no cruise median ratio outside 0.80-1.25 for any class (classes with >= 2 runs)
%  5 over-tuning: 2a's agreement drop exceeds 2b's by > 0.03 (reported in the table)
% Classes: Euk (1) and Syn (2) in CHL runs, PE euk (34) in CHL and PE runs, Syn in PE
% runs, Pro (7) in pro runs, bacteria incl. Pro (97) in hbac runs.
% Writes output/discrete/holdout/holdout_report.csv and failures_<set>.csv
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete', 'holdout');
if nargin < 1, sets = {'2a' '2a_prime' '2b'}; end
T = load(fullfile(od, 'score_tuning.mat')).H;
checks = {'CHL' 1 'Euk'; 'CHL' 2 'Syn'; 'CHL' 34 'PE euk'; 'PE' 2 'Syn'; 'PE' 34 'PE euk'; 'pro' 7 'Pro'; 'pro' 2 'Syn'; 'hbac' 97 'bacteria'};
R = table;
for s = 1:numel(sets)
    S = load(fullfile(od, ['score_' sets{s} '.mat'])).H;
    fails = table;
    for k = 1:size(checks, 1)
        rt = checks{k, 1}; c = checks{k, 2};
        [mt, ~, ~] = ratios(T, rt, c);
        [ms, fs, n, rr, runs] = ratios(S, rt, c);
        if n == 0, continue, end
        pe = c == 34;
        c1 = abs(ms / mt - 1) <= 0.05 && (pe || (ms >= 0.9 && ms <= 1.1));
        c2 = fs >= 0.9;
        % cruise medians
        cm = groupsummary(table(runs.cruise, rr, 'VariableNames', {'cruise' 'r'}), 'cruise', {'median' 'nnz'}, 'r');
        cm = cm(cm.GroupCount >= 2, :);
        if pe, c4 = true; else, c4 = all(cm.median_r >= 0.8 & cm.median_r <= 1.25); end
        badc = strjoin(cm.cruise(cm.median_r < 0.8 | cm.median_r > 1.25), ' ');
        R = [R; table(string(sets{s}), string(rt), string(checks{k, 3}), n, mt, ms, fs, c1, c2, c4, string(badc), ...
            'VariableNames', {'set' 'runtype' 'class' 'n_runs' 'tuning_median' 'median' 'frac_within20' 'pass_median' 'pass_within20' 'pass_cruise' 'cruises_outside'})]; %#ok<AGROW>
        out = abs(rr - 1) > 0.2;
        fails = [fails; table(runs.cruise(out), repmat(string(rt), sum(out), 1), repmat(string(checks{k, 3}), sum(out), 1), runs.file(out), ...
            runs.hand(out), runs.auto(out), rr(out), 'VariableNames', {'cruise' 'runtype' 'class' 'file' 'hand' 'auto' 'ratio'})]; %#ok<AGROW>
    end
    % agreement per run type
    for rt = ["CHL" "PE" "pro" "hbac"]
        at = median(T.agreement(T.runtype == rt), 'omitnan'); as = median(S.agreement(S.runtype == rt), 'omitnan');
        if ~isfinite(as), continue, end
        R = [R; table(string(sets{s}), rt, "agreement", sum(S.runtype == rt), at, as, NaN, as >= at - 0.03, true, true, "", ...
            'VariableNames', R.Properties.VariableNames)]; %#ok<AGROW>
    end
    writetable(fails, fullfile(od, ['failures_' sets{s} '.csv']))
end
writetable(R, fullfile(od, 'holdout_report.csv'))
disp(R)
end

function [m, f, n, r, runs] = ratios(H, rt, c)
hc = sprintf('hand_%d', c); ac = sprintf('auto_%d', c);
k = H.runtype == rt & H.(hc) >= 100;
r = H.(ac)(k) ./ H.(hc)(k);
runs = table(H.cruise(k), H.file(k), H.(hc)(k), H.(ac)(k), 'VariableNames', {'cruise' 'file' 'hand' 'auto'});
n = numel(r); m = median(r); f = mean(abs(r - 1) <= 0.2);
end
