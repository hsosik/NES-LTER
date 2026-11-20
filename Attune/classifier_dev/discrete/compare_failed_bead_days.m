function T = compare_failed_bead_days()
% function T = compare_failed_bead_days()
% For lab days on which preserved samples were run and a bead population
% failed QC (missing_bead_days.csv, days with a lab bead file): is the
% nearest lab run where that population passed a fair stand-in?
% Test: the populations that passed on the failed day should sit where they
% sit in that nearest run. Positions are voltage-corrected log10 peaks
% (cfg.k, referred to 300 V) on every channel.
%
% T - one row per failed day x failed population: same-day run used (the one
%     with most good populations), nearest run with the failed population
%     good, days apart, and for the populations good in both runs the median
%     and maximum |difference| over channels, plus the difference on the
%     failed population's own reference channel(s) where the failed run still
%     measured a peak
% Saves output/discrete/failed_bead_days.csv and failed_bead_days.png
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
cfg = classifier_config();
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
o = detectImportOptions(fullfile(od, 'missing_bead_days.csv'));
o = setvartype(o, {'cruise' 'missing_pops_um' 'qc_failures'}, 'string');
M = readtable(fullfile(od, 'missing_bead_days.csv'), o);
M = M(M.lab_runs_that_day > 0, :);
lday = floor(datenum(lab.time));
P = {'05' '1' '6'}; PN = {'0.5' '1' '6'};
chans = cfg.chans;
refch = @(p) chans(cellfun(@(c) strcmp(P{cfg.bead_pop.(c)}, p), chans)); %channels this population is the reference for
T = table;
[days, ~] = unique(floor(datenum(M.run_day)));
for u = days'
    rows = M(floor(datenum(M.run_day)) == u, :);
    failed = unique(split(strjoin(rows.missing_pops_um, ' ')))';
    same = find(lday == u);
    ngood = sum([lab.QC05(same) lab.QC1(same) lab.QC6(same)], 2);
    [~, b] = max(ngood); s = same(b);
    for f = failed
        pf = P{strcmp(PN, f{1})};
        cand = find(lab.(['QC' pf]) & lday ~= u);
        [dd, j] = min(abs(datenum(lab.time(cand)) - datenum(lab.time(s))));
        n = cand(j);
        good = P(cellfun(@(p) lab.(['QC' p])(s) && lab.(['QC' p])(n), P));
        d = []; dr = [];
        for p = good
            for c = chans
                d(end + 1) = pos(lab, s, p{1}, c{1}, cfg) - pos(lab, n, p{1}, c{1}, cfg); %#ok<AGROW>
            end
            for c = refch(p{1}) %channels where this population is the normalizing reference
                dr(end + 1) = pos(lab, s, p{1}, c{1}, cfg) - pos(lab, n, p{1}, c{1}, cfg); %#ok<AGROW>
            end
        end
        % the failed population itself, on its reference channels, if a peak was measured
        rc = refch(pf); df = nan(1, numel(rc));
        for k = 1:numel(rc), df(k) = pos(lab, s, pf, rc{k}, cfg) - pos(lab, n, pf, rc{k}, cfg); end
        row = table(string(datestr(u, 1)), strjoin(unique(rows.cruise), ' '), string(f{1}), string(lab.(['why' pf])(s)), ...
            string(lab.filename(s)), string(lab.filename(n)), dd, ...
            string(strjoin(strrep(good, '05', '0.5'), ' ')), median(abs(d), 'omitnan'), max(abs(d)), max(abs(dr)), ...
            string(strjoin(rc, ' ')), mean(df, 'omitnan'), sum(isfinite(df)), ...
            'VariableNames', {'day' 'cruises' 'failed_pop_um' 'why' 'run_failed_day' 'nearest_good_run' 'days_apart' ...
            'pops_compared_um' 'median_abs_diff' 'max_abs_diff' 'max_abs_diff_ref_chans' 'failed_pop_ref_chans' 'failed_pop_diff' 'n_failed_chans_measured'});
        row.diff_by_chan = {d};
        T = [T; row]; %#ok<AGROW>
    end
end
writetable(removevars(T, 'diff_by_chan'), fullfile(od, 'failed_bead_days.csv'))
disp(T(:, {'day' 'failed_pop_um' 'days_apart' 'median_abs_diff' 'max_abs_diff' 'max_abs_diff_ref_chans' 'failed_pop_diff'}))
% reference: day-to-day spread between consecutive good lab runs (all populations good)
allg = find(lab.QC05 & lab.QC1 & lab.QC6);
dn = [];
for i = 2:numel(allg)
    a = allg(i - 1); b = allg(i);
    for p = P, for c = chans, dn(end + 1) = pos(lab, b, p{1}, c{1}, cfg) - pos(lab, a, p{1}, c{1}, cfg); end, end %#ok<AGROW>
end
fprintf('consecutive all-good lab runs: median |diff| %.3f, 95th percentile %.3f (n = %d run pairs)\n', ...
    median(abs(dn), 'omitnan'), prctile(abs(dn), 95), numel(allg) - 1)
fig = figure('visible', 'off', 'position', [50 50 1000 400]); hold on
for i = 1:height(T), plot(i + 0.1 * randn(size(T.diff_by_chan{i})), T.diff_by_chan{i}, 'k.'), end
plot(1:height(T), T.failed_pop_diff, 'r*')
yline([-1 1] * prctile(abs(dn), 95), 'b:')
set(gca, 'xtick', 1:height(T), 'xticklabel', strcat(T.day, " (", T.failed_pop_um, ")"), 'xticklabelrotation', 60, 'fontsize', 8)
ylabel('failed-day run - nearest good run (log10, voltage-corrected)'), grid on
title('black: good populations, all channels; red: failed population on its reference channel; blue: 95% of consecutive-run differences', 'fontsize', 8)
exportgraphics(fig, fullfile(od, 'failed_bead_days.png'), 'resolution', 90), close(fig)
end

function v = pos(b, i, p, c, cfg)
x = b.(['peak' p]).([c '_H'])(i); hv = b.hv.(c)(i);
if ~(x > 0 && hv > 0), v = NaN; return, end
v = log10(x) - cfg.k.(c) * log10(hv / 300);
end
