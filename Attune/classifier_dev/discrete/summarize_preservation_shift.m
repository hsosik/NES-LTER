function D = summarize_preservation_shift(cruises)
% function D = summarize_preservation_shift(cruises)
% Preserved minus live, per surface bottle, from preservation_shift output:
%   Syn: hand Syn (2) in the PE run (else the CHL run) vs live v2 Syn (2)
%   Euk: hand Euk (1) in the CHL run vs live v2 Euk (1)
%   Pro: hand Pro (4) in the pro run vs live v2 Pro (7)
% Also how v2, run unchanged on preserved runs, compares with the hand gates
% (count ratio v2/hand for the same population and run).
% Writes output/discrete/preservation_shift_summary.csv and a figure.
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 1, cruises = {'EN720' 'EN695' 'EN706' 'EN715'}; end
feat = {'size' 'SSC' 'FSC' 'BL1' 'BL2' 'BL3' 'GL2' 'GL3' 'GL4' 'PEchl' 'GL2GL3' 'BL2GL2' 'chlSSC'};
pops = {'Syn' 2 2 {'PE' 'CHL'}; 'Euk' 1 1 {'CHL'}; 'Pro' 4 7 {'pro'}};
D = table;
for c = cruises
    f = fullfile(od, ['presshift_' c{1} '.mat']);
    if ~exist(f, 'file'), continue, end
    S = load(f).S;
    [bt, ~, ib] = unique(S(:, {'cast' 'niskin'}));
    for b = 1:height(bt)
        Sb = S(ib == b, :);
        for p = 1:size(pops, 1)
            live = Sb(Sb.kind == "live" & Sb.scheme == "v2" & Sb.class == pops{p, 3}, :);
            pres = table;
            for rt = pops{p, 4}
                pres = Sb(Sb.kind == "preserved" & Sb.scheme == "hand" & Sb.class == pops{p, 2} & Sb.runtype == rt{1}, :);
                if ~isempty(pres), break, end
            end
            if isempty(live) || isempty(pres), continue, end
            [~, i] = max(pres.n); pres = pres(i, :);
            % the same population by v2 on that preserved run
            v2p = Sb(Sb.kind == "preserved" & Sb.scheme == "v2" & Sb.file == pres.file & Sb.class == pops{p, 3}, :);
            row = table(string(c{1}), bt.cast(b), bt.niskin(b), Sb.depth_m(1), string(pops{p, 1}), pres.runtype, pres.n, sum(live.n), ...
                sum(v2p.n) / pres.n, pres.days_to_beads, 'VariableNames', ...
                {'cruise' 'cast' 'niskin' 'depth_m' 'pop' 'pres_run' 'n_pres' 'n_live' 'v2_over_hand_pres' 'days_to_beads'});
            for k = 1:numel(feat)
                lv = sum(live.(feat{k}) .* live.n) / sum(live.n); %event-weighted mean of file medians
                row.(['d_' feat{k}]) = pres.(feat{k}) - lv;
                row.(['live_' feat{k}]) = lv;
            end
            D = [D; row]; %#ok<AGROW>
        end
    end
end
writetable(D, fullfile(od, 'preservation_shift_summary.csv'))
% overview: median preserved - live per population and cruise
G = groupsummary(D, {'pop' 'cruise'}, 'median', [{'v2_over_hand_pres'} strcat('d_', feat)]);
disp(G(:, [{'pop' 'cruise' 'GroupCount' 'median_v2_over_hand_pres'} strcat('median_d_', feat)]))
fig = figure('visible', 'off', 'position', [50 50 1400 700]);
P = unique(D.pop, 'stable');
for k = 1:numel(feat)
    subplot(3, 5, k), hold on
    for p = 1:numel(P)
        x = D.(['d_' feat{k}])(D.pop == P(p));
        plot(p + 0.15 * randn(size(x)), x, '.', 'markersize', 10)
    end
    yline(0, 'k:'), xlim([0.5 numel(P) + 0.5]), set(gca, 'xtick', 1:numel(P), 'xticklabel', P)
    title(['preserved - live: ' feat{k}], 'fontsize', 9), grid on
end
subplot(3, 5, 15), hold on
for p = 1:numel(P), x = D.v2_over_hand_pres(D.pop == P(p)); plot(p + 0.15 * randn(size(x)), x, '.', 'markersize', 10), end
yline(1, 'k:'), set(gca, 'xtick', 1:numel(P), 'xticklabel', P, 'yscale', 'log'), title('v2 / hand count, preserved run', 'fontsize', 9), grid on
exportgraphics(fig, fullfile(od, 'preservation_shift.png'), 'resolution', 90)
close(fig)
end
