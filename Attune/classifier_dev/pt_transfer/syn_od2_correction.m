function T = syn_od2_correction()
% function T = syn_od2_correction()
% Does the low-signal excess measured directly on SSC with OD2 (March 2019
% Syn culture, culture_od2_fluorgate) explain the high bead-normalized Syn
% SSC of the OD2-on-SSC cruises?
%   excess(raw) = Syn OD2 loss at quantile - plateau loss, interpolated in raw
%                 OD2 SSC (quantile pairs 5-90%), 0 above the highest point
%   bead term   = 1 um bead loss - plateau (culture_od2_check): the bead
%                 itself reads slightly low under OD2
%   predicted offset vs no-OD2 cruises = excess(Syn raw peak) - bead term
% compared with (observed normalized Syn) - (median of FCB no-OD2 cruises in
% the same month; April: mean of Feb and Jul, no FCB cruise in April).
% Also draws output/pt_transfer/culture_od2.png.
%
% Nov 2026, PT transfer / OD2 low-signal check

od = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'output');
s = load(fullfile(od, 'pt_transfer', 'culture_od2.mat'), 'T'); C = s.T;
s = load(fullfile(od, 'pt_transfer', 'culture_od2_fluorgate.mat'), 'Q'); Q = s.Q;
plateau = C.plateau;
bead = C.pk.loss(C.pk.pop == "beads_1um") - plateau;
p = [5 10 25 50 75 90];
q = Q(Q.pop == "Syn", :);
x = arrayfun(@(k) q.(sprintf('raw_od2_p%d', k)), p);
y = arrayfun(@(k) q.(sprintf('loss_p%d', k)), p) - plateau;
[x, o] = sort(x); y = max(y(o), 0);
f = @(r) interp1(x, y, min(max(r, x(1)), x(end)), 'linear');
fprintf('plateau %.3f, 1 um bead term %.3f\n', plateau, bead)
fprintf('Syn excess at raw 1.8/1.9/2.0/2.1: %.2f %.2f %.2f %.2f\n', f(1.8), f(1.9), f(2.0), f(2.1))

s = load(fullfile(od, 'syn_ssc_check.mat'), 'S'); S = s.S;
ref = groupsummary(S(S.source == "FCB" & S.cruise ~= "TN368", :), 'month', 'median', 'norm_ssc');
seas = @(m) interp1([ref.month; 4], [ref.median_norm_ssc; mean(ref.median_norm_ssc(ismember(ref.month, [2 7])))], m);
s = load(fullfile(od, 'pt_transfer', 'compare_perf_refs.mat'), 'C'); E = s.C;
E = E(E.cruise ~= "EN627" & E.norm_new < -0.5, :);
A = table(E.cruise, E.month, E.syn_ssc, E.norm_new, 'VariableNames', {'cruise', 'month', 'raw', 'norm'});
St = S(S.cruise == "TN368" & S.norm_ssc > -1.6, :);
A = [A; table(St.cruise, St.month, St.syn_ssc, St.norm_ssc, 'VariableNames', A.Properties.VariableNames)];
T = groupsummary(A, 'cruise', 'median', {'month', 'raw', 'norm'});
T.Properties.VariableNames = erase(T.Properties.VariableNames, 'median_');
T.seasonal = seas(T.month);
T.observed_offset = T.norm - T.seasonal;
T.predicted_offset = f(T.raw) - bead;
T.residual = T.observed_offset - T.predicted_offset;
disp(T(:, {'cruise', 'month', 'raw', 'norm', 'seasonal', 'observed_offset', 'predicted_offset', 'residual'}))

% figure: OD2 loss - plateau vs raw OD2 SSC, all cultures (quantiles) and beads
fg = figure('visible', 'off', 'position', [0 0 1000 600]); hold on
col = [lines(7); 0.5 0.5 0.5; 0 0 0; 0.6 0.4 0.2; 0.9 0.5 0.8]; col = col(1:height(Q), :);
for i = 1:height(Q)
    xr = arrayfun(@(k) Q.(sprintf('raw_od2_p%d', k))(i), p);
    yr = arrayfun(@(k) Q.(sprintf('loss_p%d', k))(i), p) - plateau;
    plot(xr, yr, 'o-', 'color', col(i, :), 'markerfacecolor', col(i, :), 'DisplayName', strrep(Q.pop(i), '_', ' '))
end
b = startsWith(C.pk.pop, "beads");
plot(C.pk.raw_od2(b), C.pk.loss(b) - plateau, 'ks', 'markerfacecolor', 'y', 'markersize', 9, 'DisplayName', 'FCB beads (0.5, 1, 6 um)')
yline(0, 'k:', 'HandleVisibility', 'off')
xlabel('raw SSC-H with OD2 (log10), culture quantiles 5-90%')
ylabel('OD2 loss - plateau (log10): > 0 reads high with OD2 or compressed without')
title(sprintf('Same cultures without OD2 (SSC 200 V) and with OD2 on SSC (260 V); plateau loss %.2f (k = %.1f)', plateau, C.k))
legend('location', 'eastoutside'); grid on; ylim([-0.4 0.8])
exportgraphics(fg, fullfile(od, 'pt_transfer', 'culture_od2.png'), 'Resolution', 110)
end
