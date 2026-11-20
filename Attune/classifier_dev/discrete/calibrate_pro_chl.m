function C = calibrate_pro_chl(maxper)
% function C = calibrate_pro_chl(maxper)
% Calibrates the lower-bound flag for Pro counted in CHL runs. For samples with
% both a CHL run and a Pro run (all cruises with Gated_Tables, current filter
% layout), classifies both with assign_class_discrete and compares the CHL
% run's Pro concentration with the Pro run's, against the CHL run's Pro
% completeness metrics (info.pro: peak_above_edge, frac_near_edge).
% Concentrations are counts / $VOL of the whole run (both runs alike).
% maxper - samples per cruise (default 10)
% Saves output/discrete/calibrate_pro_chl.csv and .png
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 1, maxper = 10; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = discrete_config();
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
d = dir(fullfile(root, '*_*', 'preserved', 'outputs', 'Gated_Table.mat'));
C = table;
for i = 1:numel(d)
    try, G = load(fullfile(d(i).folder, d(i).name)).gated_table; catch, continue, end
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    if ~all(ismember({'cast' 'niskin'}, G.Properties.VariableNames)), continue, end
    cd_ = char(extractBetween(d(i).folder, root, filesep)); cruise = extractAfter(cd_, '_');
    if strcmp(attune_config(cruise).layout, 'v2019'), continue, end
    cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
    ispro = contains(G.fcslist, 'pro', 'IgnoreCase', true) & contains(G.fcslist, 'phyto_CHL') & ~contains(G.fcslist, 'skip');
    ischl = contains(G.fcslist, 'phyto_CHL') & ~ispro & ~contains(G.fcslist, 'skip');
    [u, ~, ib] = unique([G.cast G.niskin], 'rows');
    both = find(accumarray(ib, ispro) > 0 & accumarray(ib, ischl) > 0);
    both = both(unique(round(linspace(1, numel(both), min(maxper, numel(both))))));
    for b = both'
        rp = find(ib == b & ispro, 1); rc = find(ib == b & ischl, 1);
        try
            [dp, hp] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', G.fcslist{rp}));
            [dc, hc] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', G.fcslist{rc}));
        catch
            continue
        end
        refp = discrete_bead_reference(lab, datenum([hp.date ' ' hp.starttime]), cfg);
        refc = discrete_bead_reference(lab, datenum([hc.date ' ' hc.starttime]), cfg);
        if isempty(refp) || isempty(refc), continue, end
        [cp, ip] = assign_class_discrete(dp, hp, refp, cfg, 'pro');
        [cc, ic] = assign_class_discrete(dc, hc, refc, cfg, 'CHL');
        depth = NaN; if ismember('depth_m', G.Properties.VariableNames), depth = G.depth_m(rc); end
        C = [C; table(string(cruise), u(b, 1), u(b, 2), depth, sum(cp == 7), sum(cc == 7), ...
            sum(cp == 7) / hp.VOL, sum(cc == 7) / hc.VOL, ip.pro.peak_above_edge, ic.pro.peak_above_edge, ...
            ic.pro.frac_near_edge, ic.pro.lower_bound, ...
            'VariableNames', {'cruise' 'cast' 'niskin' 'depth_m' 'n_pro_run' 'n_chl_run' 'conc_pro_run' 'conc_chl_run' ...
            'pro_run_peak_above_edge' 'chl_run_peak_above_edge' 'chl_run_frac_near_edge' 'chl_run_lower_bound'})]; %#ok<AGROW>
    end
    fprintf('%s: %d samples\n', cruise, numel(both))
end
C.ratio = C.conc_chl_run ./ C.conc_pro_run;
writetable(C, fullfile(od, 'calibrate_pro_chl.csv'))
x = C(C.n_pro_run >= 200 & C.n_chl_run >= cfg.sig.pro_edge_min_n, :);
for lb = [false true]
    y = x(x.chl_run_lower_bound == lb, :);
    fprintf('lower_bound %d: %d samples, CHL/Pro ratio median %.2f (10-90%%: %.2f-%.2f)\n', lb, height(y), median(y.ratio), prctile(y.ratio, [10 90]))
end
fig = figure('visible', 'off', 'position', [50 50 1100 420]);
subplot(1, 2, 1), scatter(x.chl_run_peak_above_edge, x.ratio, 20, x.chl_run_frac_near_edge, 'filled'), colorbar
set(gca, 'yscale', 'log'), yline(1, 'k:'), xline(cfg.sig.pro_lb_min_peak_above_edge, 'r:'), grid on
xlabel('CHL run: Pro chl mode above trigger edge'), ylabel('Pro conc, CHL run / Pro run'), title('colour: fraction at edge', 'fontsize', 8)
subplot(1, 2, 2), plot(x.chl_run_frac_near_edge, x.ratio, 'o'), set(gca, 'yscale', 'log'), yline(1, 'k:'), xline(cfg.sig.pro_lb_max_frac_near_edge, 'r:'), grid on
xlabel('CHL run: fraction of Pro within 0.1 of the edge'), ylabel('Pro conc, CHL run / Pro run')
exportgraphics(fig, fullfile(od, 'calibrate_pro_chl.png'), 'resolution', 90), close(fig)
end
