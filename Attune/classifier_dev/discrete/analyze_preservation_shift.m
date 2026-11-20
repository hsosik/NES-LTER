function [P, T] = analyze_preservation_shift(scheme)
% function [P, T] = analyze_preservation_shift(scheme)
% Does preservation shrink cells enough to need a volume correction?
% Pairs each surface bottle's preserved population with the same population in
% underway (live) files from the same water (preservation_shift2 output,
% presshift2_<cruise>.mat) and summarizes preserved - live:
%   size  - median bead-relative merged SSC/GL1 (log10, the volume step's input)
%   vol   - log10 volume = p1 * size + p2 (Attune_SC2vol 'SSC-H', p1 1.1962),
%           so a median shift in size is p1 times that in log10 volume
% Populations (live: v2 on underway files; preserved: scheme on the run used
% for that group in the discrete table):
%   Syn 2 (PE run, else CHL), Euk 1 (CHL, else PE), Pro 7 (pro run; live from
%   Pro-setting underway files only), dim crypto 3 and bright crypto 4 (CHL)
% The preserved hand gates (scheme 'hand') use old codes: Euk 1, Syn 2, Pro 4.
% Also the live null: the two nearest underway files compared with each other.
%
% scheme - 'disc' (assign_class_discrete, default) or 'hand'
% P      - one row per bottle x population pair (also csv)
% T      - per population summary (also csv)
% Figures in classifier_dev/output/discrete/preservation_correction_*.png
%
% Nov 2026, preservation-correction analysis

if nargin < 1, scheme = "disc"; end
scheme = string(scheme);
here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
p1 = 1.1962; %Attune_SC2vol SSC-H slope (14 Nov 2026)
min_n = 50; min_n_crypto = 20;
% pop, live class, preserved class (disc), preserved class (hand), preserved run order
pops = {'Syn' 2 2 2 {'PE' 'CHL'}; 'Euk' 1 1 1 {'CHL' 'PE'}; 'Pro' 7 7 4 {'pro'}; ...
    'crypto_dim' 3 3 NaN {'CHL' 'PE'}; 'crypto_bright' 4 4 NaN {'CHL' 'PE'}};
qv = {'size' 'size_p10' 'size_p25' 'size_p75' 'size_p90' 'size_sd' 'BL3' 'GL2' 'PEchl'};

f = dir(fullfile(od, 'presshift2_*.mat'));
P = table; N = table;
for i = 1:numel(f)
    S = load(fullfile(od, f(i).name)).S;
    if isempty(S), continue, end
    [bt, ~, ib] = unique(S(:, {'cast' 'niskin'}));
    for b = 1:height(bt)
        Sb = S(ib == b, :);
        for p = 1:size(pops, 1)
            pc = pops{p, 3}; if scheme == "hand", pc = pops{p, 4}; end
            if isnan(pc), continue, end
            mn = min_n; if startsWith(pops{p, 1}, 'crypto'), mn = min_n_crypto; end
            live = Sb(Sb.kind == "live" & Sb.scheme == "v2" & Sb.class == pops{p, 2} & Sb.n >= mn, :);
            if pops{p, 1} == "Pro", live = live(live.runtype == "underway_pro", :);
            else, live = live(live.runtype == "underway", :); end
            pres = table;
            for rt = pops{p, 5}
                pres = Sb(Sb.kind == "preserved" & Sb.scheme == scheme & Sb.class == pc & Sb.runtype == rt{1} & Sb.n >= mn, :);
                if ~isempty(pres), break, end
            end
            if isempty(live) || isempty(pres), continue, end
            [~, j] = max(pres.n); pres = pres(j, :);
            row = table(Sb.cruise(1), bt.cast(b), bt.niskin(b), Sb.depth_m(1), datetime(Sb.time_sampled(1), 'ConvertFrom', 'datenum'), ...
                string(pops{p, 1}), pres.runtype, pres.n, sum(live.n), height(live), ...
                pres.time_run - Sb.time_sampled(1), min(abs(live.time_run - Sb.time_sampled(1))) * 1440, pres.days_to_beads, ...
                'VariableNames', {'cruise' 'cast' 'niskin' 'depth_m' 'time_sampled' 'pop' 'pres_run' 'n_pres' 'n_live' 'n_live_files' ...
                'storage_days' 'live_minutes' 'days_to_beads'});
            for k = 1:numel(qv)
                lv = sum(live.(qv{k}) .* live.n) / sum(live.n); %event-weighted mean of file medians
                row.(['live_' qv{k}]) = lv;
                row.(['pres_' qv{k}]) = pres.(qv{k});
                row.(['d_' qv{k}]) = pres.(qv{k}) - lv;
            end
            P = [P; row]; %#ok<AGROW>
            % live null: two underway files of this bottle against each other
            if height(live) >= 2
                N = [N; table(Sb.cruise(1), string(pops{p, 1}), live.size(1) - live.size(2), ...
                    abs(live.time_run(1) - live.time_run(2)) * 1440, 'VariableNames', {'cruise' 'pop' 'd_size' 'minutes'})]; %#ok<AGROW>
            end
        end
    end
end
P.d_vol = p1 * P.d_size;
P.vol_factor = 10.^P.d_vol;
mo = month(P.time_sampled);
season = ["winter" "winter" "spring" "spring" "spring" "summer" "summer" "summer" "fall" "fall" "fall" "winter"];
P.season = season(mo)';
P.year = year(P.time_sampled);
writetable(P, fullfile(od, "preservation_correction_pairs_" + scheme + ".csv"))
writetable(N, fullfile(od, "preservation_correction_livenull.csv"))

% per population summary; cruise is the unit for the uncertainty
T = table;
for p = unique(P.pop, 'stable')'
    x = P(P.pop == p, :);
    C = groupsummary(x, 'cruise', 'median', 'd_size');
    nl = N.d_size(N.pop == p);
    lm = fitlm(table(log10(max(x.storage_days, 1)), x.live_size, x.depth_m, x.d_size, ...
        'VariableNames', {'log_storage' 'live_size' 'depth' 'd_size'}), 'd_size ~ log_storage + live_size + depth');
    % season and cruise as factors: share of variance between cruises
    a = anova1(x.d_size, x.cruise, 'off');
    pe = anova1(x.d_size, x.season, 'off');
    r = table(p, height(x), numel(unique(x.cruise)), median(x.d_size), mean(x.d_size), std(x.d_size), ...
        median(C.median_d_size), std(C.median_d_size), std(C.median_d_size) / sqrt(height(C)), ...
        prctile(C.median_d_size, 10), prctile(C.median_d_size, 90), mad(nl, 1) * 1.4826, ...
        p1 * median(C.median_d_size), 10^(p1 * median(C.median_d_size)), ...
        median(x.d_size_p10), median(x.d_size_p90), median(x.d_size_sd), ...
        lm.Coefficients.Estimate(2), lm.Coefficients.pValue(2), lm.Coefficients.Estimate(3), lm.Coefficients.pValue(3), ...
        lm.Coefficients.Estimate(4), lm.Coefficients.pValue(4), a, pe, ...
        'VariableNames', {'pop' 'n_bottles' 'n_cruises' 'median_d_size' 'mean_d_size' 'sd_d_size' ...
        'cruise_median' 'cruise_sd' 'cruise_se' 'cruise_p10' 'cruise_p90' 'live_null_sd' ...
        'd_logvol' 'vol_factor' 'median_d_p10' 'median_d_p90' 'median_d_sd' ...
        'slope_logstorage' 'p_logstorage' 'slope_livesize' 'p_livesize' 'slope_depth' 'p_depth' 'p_cruise' 'p_season'});
    T = [T; r]; %#ok<AGROW>
end
writetable(T, fullfile(od, "preservation_correction_summary_" + scheme + ".csv"))
disp(T)

% --- figures ---
col = lines(7);
popl = unique(P.pop, 'stable');
cr = unique(P(:, {'cruise'}), 'stable');
[~, o] = sort(arrayfun(@(c) min(P.time_sampled(P.cruise == c)), cr.cruise)); cr = cr.cruise(o);
fig = figure('visible', 'off', 'position', [50 50 1500 300 * numel(popl)]);
tiledlayout(numel(popl), 3, 'TileSpacing', 'compact')
for p = 1:numel(popl)
    x = P(P.pop == popl(p), :);
    [~, ic] = ismember(x.cruise, cr);
    % by cruise
    nexttile, hold on
    plot(ic + 0.15 * randn(size(ic)), x.d_size, '.', 'color', col(p, :), 'markersize', 9)
    C = groupsummary(x, 'cruise', 'median', 'd_size'); [~, jc] = ismember(C.cruise, cr);
    plot(jc, C.median_d_size, 'k_', 'markersize', 14, 'linewidth', 2)
    yline(0, 'k:'); yline(median(C.median_d_size), '-', sprintf('%.3f', median(C.median_d_size)), 'color', col(p, :));
    set(gca, 'xtick', 1:numel(cr), 'xticklabel', cr, 'XTickLabelRotation', 60, 'fontsize', 8), xlim([0.5 numel(cr) + 0.5])
    ylabel(sprintf('%s: \\Delta log_{10} size', strrep(popl(p), '_', ' '))), grid on
    title(sprintf('%s, preserved - live: cruise median %.3f (volume x%.2f)', strrep(popl(p), '_', ' '), ...
        median(C.median_d_size), 10^(p1 * median(C.median_d_size))), 'fontsize', 9)
    ylim([-0.5 0.5])
    % storage time
    nexttile, hold on
    for s = ["winter" "spring" "summer" "fall"]
        k = x.season == s;
        plot(max(x.storage_days(k), 1), x.d_size(k), '.', 'markersize', 10, 'displayname', s)
    end
    set(gca, 'xscale', 'log', 'fontsize', 8, 'XMinorGrid', 'off'), yline(0, 'k:'), grid on, xlabel('storage (days, sampling to run)'), ylim([-0.5 0.5])
    if p == 1, legend("" + ["winter" "spring" "summer" "fall"], 'location', 'best', 'fontsize', 7), end
    % cell size
    nexttile, hold on
    plot(x.live_size, x.d_size, '.', 'color', col(p, :), 'markersize', 10)
    yline(0, 'k:'), grid on, xlabel('live median size (log_{10} rel. 1 \mum bead)'), set(gca, 'fontsize', 8), ylim([-0.5 0.5])
end
exportgraphics(fig, fullfile(od, "preservation_correction_" + scheme + ".png"), 'resolution', 90)
close(fig)

% quantile shift: does the whole distribution move by the same amount?
fig = figure('visible', 'off', 'position', [50 50 1100 380]);
tiledlayout(1, 2, 'TileSpacing', 'compact')
nexttile, hold on
qs = {'d_size_p10' 'd_size_p25' 'd_size' 'd_size_p75' 'd_size_p90'};
for p = 1:numel(popl)
    x = P(P.pop == popl(p), :);
    m = cellfun(@(q) median(x.(q)), qs);
    plot([10 25 50 75 90] + p - 3, m, '-o', 'color', col(p, :), 'linewidth', 1.5, 'displayname', popl(p))
end
yline(0, 'k:'), grid on, xlabel('percentile of the size distribution'), ylabel('median \Delta log_{10} size (preserved - live)')
legend('location', 'best', 'interpreter', 'none')
nexttile, hold on
for p = 1:numel(popl)
    y = P.d_size(P.pop == popl(p)); nl = N.d_size(N.pop == popl(p));
    plot(p - 0.15 + 0.08 * randn(size(y)), y, '.', 'color', col(p, :))
    plot(p + 0.15 + 0.08 * randn(size(nl)), nl, '.', 'color', [0.6 0.6 0.6])
end
yline(0, 'k:'), set(gca, 'xtick', 1:numel(popl), 'xticklabel', popl, 'TickLabelInterpreter', 'none'), grid on
ylabel('\Delta log_{10} size'), title('preserved - live (colour) vs live - live, two underway files (grey)', 'fontsize', 9)
exportgraphics(fig, fullfile(od, "preservation_correction_quantiles_" + scheme + ".png"), 'resolution', 90)
close(fig)
end
