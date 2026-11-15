function [t, flag, t0, tmax] = ssc_decompress(s, hv, par, min_slope)
% function [t, flag, t0, tmax] = ssc_decompress(s, hv, par, min_slope)
% Undo the compression of SSC-H without the OD2 filter for large particles
% (cruises before Aug 2019 with no OD2: AR34B, RB1904, EN627).
% s = measured log10 SSC-H in 1 um bead units; t = log10 linear SSC in bead
% units, equivalent to SSC with OD2, so the 488 nm volume calibration
% (Attune_SC2vol 'SSC-H_488_2019') applies.
% Model: s = t - c/b log10(1 + 10^(b (t - t0))), par = [t0 b c], fitted on the
% March/April 2019 cultures run with and without OD2 at SSC 200 V
% (classifier_dev/pt_transfer/ssc_compression_fit.m; default par below).
% hv (optional, SSC voltage): compression is worse at low voltage, so the knee
% t0 is moved by the change in the 6 um bead's apparent ratio to the 1 um bead
% in the SSC voltage sweeps (bead_hv_experiment2, Dec 2019 + Jul 2020; true
% 6/1 um = 1.257 log from the 2020 filter swap with OD2), relative to 200 V.
% Where the curve is too flat to invert reliably (local slope ds/dt <
% min_slope, default 0.5: errors in s amplified more than 2x), t is NaN and
% flag is true; tmax is the largest t that can be recovered (a lower bound
% for flagged particles).
%
% Nov 2026, classifier_dev/pt_transfer/README.md

if nargin < 3 || isempty(par), par = [1.089743 1.697809 0.934884]; end
if nargin < 4 || isempty(min_slope), min_slope = 0.5; end
t0 = par(1); b = par(2); c = par(3);
if nargin > 1 && ~isempty(hv)
    % 6/1 um SSC-H ratio vs SSC voltage, mean of the two sweeps
    hvs = [150 160 170 180 190 200 210 220 230];
    r61 = [0.891 0.9125 0.927 0.947 0.960 0.977 0.970 0.974 0.983];
    t6 = 1.257;
    knee = @(s6) t6 - log10(10.^((t6 - s6) * b / c) - 1) / b;
    t0 = t0 + knee(interp1(hvs, r61, min(max(hv, hvs(1)), hvs(end)))) - knee(r61(hvs == 200));
end
model = @(t) t - c / b * log10(1 + 10.^(b * (t - t0)));
sigmax = (1 - min_slope) / c; % largest t with slope 1 - c*sig >= min_slope
if sigmax >= 1, tmax = Inf; else, tmax = t0 + log10(sigmax / (1 - sigmax)) / b; end
smax = model(tmax);
tg = linspace(-4, min(tmax, 4), 4000);
sg = model(tg);
t = interp1(sg, tg, s, 'linear', NaN);
lo = s < sg(1);
t(lo) = s(lo); % far below the knee: linear
flag = s > smax;
t(flag) = NaN;
end
