function P = pro_in_chl_runs(maxper)
% function P = pro_in_chl_runs(maxper)
% Could CHL runs give usable Pro counts? For samples (cast/niskin) that have
% both a CHL run and a Pro run, all cruises with Gated_Tables:
%   - hand Pro concentration in the CHL run (where the gater drew a Pro gate)
%     vs in the Pro run
%   - the Pro run's hand-gated Pro on the bead scale vs the CHL run's chl (BL3)
%     and SSC trigger edges: what fraction of the Pro population the CHL
%     run's settings could record at all
% maxper - samples per cruise (default 8)
% Saves output/discrete/pro_in_chl_runs.csv and .png
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 1, maxper = 8; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = discrete_config();
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
d = dir(fullfile(root, '*_*', 'preserved', 'outputs', 'Gated_Table.mat'));
P = table;
for i = 1:numel(d)
    try, G = load(fullfile(d(i).folder, d(i).name)).gated_table; catch, continue, end
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    if ~all(ismember({'cast' 'niskin' 'class'}, G.Properties.VariableNames)), continue, end
    cd_ = char(extractBetween(d(i).folder, root, filesep)); cruise = extractAfter(cd_, '_');
    if strcmp(attune_config(cruise).layout, 'v2019'), continue, end
    cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
    if ismember('vol_analyzed_ml', G.Properties.VariableNames), vol = G.vol_analyzed_ml; else, vol = nan(height(G), 1); end
    ispro = contains(G.fcslist, 'pro', 'IgnoreCase', true) & contains(G.fcslist, 'phyto_CHL');
    ischl = contains(G.fcslist, 'phyto_CHL') & ~ispro;
    ok = ~cellfun(@isempty, G.class);
    [u, ~, ib] = unique([G.cast G.niskin], 'rows');
    both = find(accumarray(ib, ispro & ok) > 0 & accumarray(ib, ischl & ok) > 0);
    both = both(unique(round(linspace(1, numel(both), min(maxper, numel(both))))));
    for b = both'
        rp = find(ib == b & ispro & ok, 1); rc = find(ib == b & ischl & ok, 1);
        try
            [dp, hp] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', G.fcslist{rp}));
            [dc, hc] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', G.fcslist{rc}));
        catch
            continue
        end
        refp = discrete_bead_reference(lab, datenum([hp.date ' ' hp.starttime]), cfg);
        refc = discrete_bead_reference(lab, datenum([hc.date ' ' hc.starttime]), cfg);
        if isempty(refp) || isempty(refc), continue, end
        [Fp, ip] = attune_features(dp, hp, refp);
        [~, ic] = attune_features(dc, hc, refc);
        ec = trig_edges(ic, cfg); %CHL run trigger edges on the bead scale
        hpro = G.class{rp}(:) == 4; hchl = G.class{rc}(:) == 4;
        x = Fp(hpro, :);
        f_bl3 = mean(x.BL3 > ec.BL3); f_both = mean(x.BL3 > ec.BL3 & x.SSC > ec.SSC);
        P = [P; table(string(cruise), u(b, 1), u(b, 2), sum(hpro), sum(hchl), sum(hpro) / vol(rp), sum(hchl) / vol(rc), ...
            median(x.BL3), ec.BL3, f_bl3, f_both, 'VariableNames', {'cruise' 'cast' 'niskin' 'n_pro_run' 'n_chl_run' ...
            'conc_pro_run' 'conc_chl_run' 'pro_BL3_median' 'chl_run_BL3_edge' 'frac_above_chl_BL3_edge' 'frac_recordable_chl_run'})]; %#ok<AGROW>
    end
    fprintf('%s: %d samples\n', cruise, numel(both))
end
P.conc_ratio_chl_over_pro = P.conc_chl_run ./ P.conc_pro_run;
writetable(P, fullfile(od, 'pro_in_chl_runs.csv'))
x = P(P.n_pro_run >= 200, :);
fprintf('\nsamples with >= 200 Pro in the Pro run: %d\n', height(x))
fprintf('fraction of Pro-run Pro recordable with CHL-run triggers: median %.2f (10-90%%: %.2f-%.2f)\n', median(x.frac_recordable_chl_run), prctile(x.frac_recordable_chl_run, [10 90]))
y = x(x.n_chl_run > 0, :);
fprintf('CHL runs with a hand Pro gate: %d; CHL/Pro concentration ratio median %.2f (10-90%%: %.2f-%.2f)\n', height(y), median(y.conc_ratio_chl_over_pro, 'omitnan'), prctile(y.conc_ratio_chl_over_pro, [10 90]))
fprintf('  same ratio vs recordable fraction: r = %.2f\n', corr(y.conc_ratio_chl_over_pro, y.frac_recordable_chl_run, 'rows', 'complete'))
fig = figure('visible', 'off', 'position', [50 50 1000 420]);
subplot(1, 2, 1), histogram(x.frac_recordable_chl_run, 0:0.05:1), xlabel('fraction of Pro-run Pro above CHL-run trigger edges'), ylabel('samples'), grid on
subplot(1, 2, 2), plot(y.frac_recordable_chl_run, y.conc_ratio_chl_over_pro, 'o'), hold on, plot([0 1], [0 1], 'k:')
xlabel('recordable fraction'), ylabel('hand Pro conc, CHL run / Pro run'), grid on, set(gca, 'yscale', 'log')
exportgraphics(fig, fullfile(od, 'pro_in_chl_runs.png'), 'resolution', 90), close(fig)
end

function e = trig_edges(info, cfg)
e = struct('BL3', -Inf, 'SSC', -Inf);
trigs = [info.trigger info.trigger2]; levels = [info.trigger_level info.trigger2_level];
for t = 1:2
    ch = char(regexp(trigs(t), '(?<=_)\w+$', 'match', 'once'));
    if isfield(e, ch), e.(ch) = info.floor.(ch) + log10(levels(t)) - log10(cfg.floor_raw); end
end
end
