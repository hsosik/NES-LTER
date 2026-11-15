function T = apply_lowsignal_curve()
% function T = apply_lowsignal_curve()
% Test: does the low-signal excess seen on GL1 (OD2) in the current layout
% (gl1_lowsignal_check) explain the high bead-normalized Syn SSC of the
% OD2-on-SSC cruises? Excess(raw) = median GL1norm - SSCnorm in bins of raw
% GL1 minus its plateau (raw 2.3-3.5); applied at each file's raw Syn SSC
% peak on the OD2-on-SSC cruises (perf-transfer references, and TN368 with
% FCB beads).
%
% Nov 2026, PT transfer / OD2 low-signal check
od = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'output');
s = load(fullfile(od, 'pt_transfer', 'gl1_lowsignal.mat'), 'G'); B = s.G.bins;
e = 1.4:0.05:4; [~, ~, bi] = histcounts(B.raw_gl1, e);
m = accumarray(bi(bi > 0), B.dGL1(bi > 0), [numel(e) - 1 1], @median, NaN);
x = e(1:end-1)' + 0.025;
plateau = median(m(x >= 2.3 & x <= 3.5), 'omitnan');
exc = max(m - plateau, 0); exc(x > 2.3) = 0;
ok = isfinite(exc);
f = @(r) interp1(x(ok), exc(ok), min(max(r, x(find(ok, 1))), 4), 'linear', 0);
fprintf('plateau %.3f; excess at raw 1.6/1.8/2.0/2.2: %.2f %.2f %.2f %.2f\n', plateau, f(1.6), f(1.8), f(2.0), f(2.2))
s = load(fullfile(od, 'pt_transfer', 'compare_perf_refs.mat'), 'C'); C = s.C;
C = C(C.cruise ~= "EN627" & C.norm_new < -0.5, :);
T = table(C.cruise, C.month, C.syn_ssc, C.norm_new, 'VariableNames', {'cruise', 'month', 'raw', 'norm'});
s = load(fullfile(od, 'syn_ssc_check.mat'), 'S'); S = s.S(s.S.cruise == "TN368" & s.S.norm_ssc > -1.6, :);
T = [T; table(S.cruise, S.month, S.syn_ssc, S.norm_ssc, 'VariableNames', T.Properties.VariableNames)];
T.excess = f(T.raw); T.norm_corr = T.norm - T.excess;
disp(groupsummary(T, 'cruise', 'median', {'month', 'raw', 'norm', 'excess', 'norm_corr'}))
end
