% analyze_normalization
% Do voltage and bead normalization make cluster positions comparable across
% cruises? Uses output of extract_cluster_centers and collect_bead_runs.
%  raw  : log10 median of the existing class, height channel
%  volt : raw corrected to a reference voltage with signal ~ hv^k
%  bead : volt-corrected, then expressed relative to the nearest bead run
%         (bead brought to the sample's voltage first)
% Spread is split into between-cruise (SD of cruise medians) and
% within-cruise parts; Pro-setting vs standard offset is a direct test of
% the voltage correction since those samples alternate on the same water.
%
% Nov 2026, classifier redesign exploration

outdir = fullfile(fileparts(mfilename('fullpath')), 'output');
load(fullfile(outdir, 'bead_runs.mat'), 'B', 'chans')
load(fullfile(outdir, 'cluster_centers.mat'), 'T')
T = T(T.nevents > 0, :);
T.month = month(datetime(T.time, 'ConvertFrom', 'datenum'));
fprintf('%d files from %d cruises\n', height(T), numel(unique(T.cruise)))

%% 1. voltage exponent from beads, and independently from Syn in alternating samples
G = fit_gain_curves(B, chans);
G.k_syn = nan(height(G), 1); G.k_syn_se = nan(height(G), 1);
for i = 1:height(G)
    ch = char(G.chan(i));
    d = table(T.(['syn_' ch]), log10(T.(['hv_' ch])), categorical(T.cruise), 'VariableNames', {'y' 'lhv' 'cruise'});
    d = d(isfinite(d.y) & isfinite(d.lhv), :);
    [g, ~] = findgroups(d.cruise);
    nv = splitapply(@(x) numel(unique(x)), d.lhv, g);
    d = d(nv(g) > 1, :);
    if height(d) < 6, continue, end
    d.cruise = removecats(d.cruise);
    m = fitlm(d, 'y ~ lhv + cruise');
    G.k_syn(i) = m.Coefficients{'lhv', 'Estimate'};
    G.k_syn_se(i) = m.Coefficients{'lhv', 'SE'};
end
kdefault = 8.5;
G.k_use = G.k;
G.k_use(isnan(G.k_use)) = G.k_syn(isnan(G.k_use));
G.k_use(isnan(G.k_use)) = kdefault;
disp(G(:, {'chan' 'k' 'k_se' 'k_syn' 'k_syn_se' 'k_use'}))

%% 2. normalized positions
methods = {'raw' 'volt' 'bead'};
cls = {'syn' 'euk' 'pro'};
for c = cls
    for i = 1:height(G)
        ch = char(G.chan(i)); k = G.k_use(i);
        hv = T.(['hv_' ch]);
        hvref = mode(hv);
        y = T.([c{1} '_' ch]);
        T.([c{1} '_' ch '_raw']) = y;
        T.([c{1} '_' ch '_volt']) = y - k * log10(hv ./ hvref);
        % pick the bead population that is best on scale in this channel
        bb = [median(T.(['bead1_' ch]), 'omitnan') median(T.(['bead2_' ch]), 'omitnan') median(T.(['bead3_' ch]), 'omitnan')];
        [~, bsel] = min(abs(bb - 4)); G.bead_used(i) = bsel;
        bead_at_sample_hv = T.(sprintf('bead%d_%s', bsel, ch)) + k * log10(hv ./ T.(['beadhv_' ch]));
        T.([c{1} '_' ch '_bead']) = y - bead_at_sample_hv;
    end
end

%% 3. spread metrics
S = table;
for c = cls
    for i = 1:height(G)
        ch = char(G.chan(i));
        for m = methods
            y = T.([c{1} '_' ch '_' m{1}]);
            ok = isfinite(y);
            if sum(ok) < 20, continue, end
            [g, gn] = findgroups(T.cruise(ok));
            cm = splitapply(@median, y(ok), g);
            within = sqrt(mean(splitapply(@(x) var(x, 1), y(ok), g)));
            % Pro-setting minus standard, paired within cruise
            pm = T.pro_measured(ok);
            off = splitapply(@(x, p) median(x(p == 1)) - median(x(p == 0)), y(ok), pm, g);
            S = [S; table(string(c{1}), string(ch), string(m{1}), sum(ok), numel(gn), std(y(ok), 'omitnan'), std(cm), within, mean(abs(off), 'omitnan'), ...
                'VariableNames', {'class' 'chan' 'method' 'nfiles' 'ncruise' 'sd_total' 'sd_between_cruise' 'sd_within_cruise' 'pro_std_offset'})]; %#ok<AGROW>
        end
    end
end
save(fullfile(outdir, 'normalization_results.mat'), 'T', 'G', 'S')
writetable(S, fullfile(outdir, 'normalization_spread.csv'))
disp(S(S.class == "syn", :))
disp(S(S.class == "euk" & ismember(S.chan, ["SSC" "BL3" "GL2"]), :))

%% 4. what is left after normalization: season / temperature?
for ch = {'GL2' 'BL3' 'SSC'}
    y = T.(['syn_' ch{1} '_bead']);
    ok = isfinite(y) & isfinite(T.temp);
    if sum(ok) > 20
        r = corr(y(ok), T.temp(ok));
        fprintf('syn %s (bead-normalized) vs temperature: r = %.2f (n = %d)\n', ch{1}, r, sum(ok))
    end
end

%% 5. figures
plotch = {'SSC' 'BL3' 'GL2' 'GL3'};
f = figure('Visible', 'off', 'Position', [50 50 1500 1000]);
tl = tiledlayout(numel(plotch), 3, 'TileSpacing', 'compact');
for i = 1:numel(plotch)
    for j = 1:3
        nexttile
        y = T.(['syn_' plotch{i} '_' methods{j}]);
        scatter(datetime(T.time, 'ConvertFrom', 'datenum'), y, 14, T.temp, 'filled', 'MarkerFaceAlpha', .7)
        hold on
        p = T.pro_measured == 1;
        scatter(datetime(T.time(p), 'ConvertFrom', 'datenum'), y(p), 22, 'k')
        if i == 1, title(methods{j}), end
        if j == 1, ylabel(['Syn ' plotch{i} '-H (log10)']), end
        grid on
    end
end
colormap(turbo); cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'temperature (°C)';
title(tl, 'Syn cluster median by file: raw, voltage-corrected, bead-relative (black rings = Pro settings)')
exportgraphics(f, fullfile(outdir, 'syn_centers_vs_time.png'), 'Resolution', 110)

f = figure('Visible', 'off', 'Position', [50 50 1400 500]);
tiledlayout(1, 3, 'TileSpacing', 'compact')
for j = 1:3
    nexttile
    scatter(T.(['syn_BL3_' methods{j}]), T.(['syn_GL2_' methods{j}]), 12, T.temp, 'filled', 'MarkerFaceAlpha', .6); hold on
    scatter(T.(['euk_BL3_' methods{j}]), T.(['euk_GL2_' methods{j}]), 12, [.5 .5 .5], '^')
    title([methods{j} ': Syn (filled) and Euk (triangles)']); xlabel('BL3-H chl'); ylabel('GL2-H PE'); grid on
end
colormap(turbo)
exportgraphics(f, fullfile(outdir, 'syn_euk_chl_pe.png'), 'Resolution', 110)
disp('figures saved')
