function [pk, med, rsd] = bead_peak(y)
% function [pk, med, rsd] = bead_peak(y)
% Peak (mode), median and robust SD of one bead population on one channel.
% y - log10 signal of the gated bead events
%
% The mode is the maximum of a Gaussian kernel density estimate in log10
% space, so it targets the bead peak and ignores skewed tails or a minor
% second mode (the purpose of the old centend2). Bandwidth scales with the
% population's robust SD and n^(-1/7) (the rate suited to mode estimation),
% so smoothing is the same in relative terms for every bead and channel.
%
% Nov 2026, bead stats redesign prototype

y = y(isfinite(y));
pk = NaN; med = NaN; rsd = NaN;
if numel(y) < 10, return, end
med = median(y);
rsd = 1.4826 * mad(y, 1);
h = max(0.003, 1.06 * rsd * numel(y) ^ (-1/7));
bw = h / 5;
lo = med - max(6 * rsd, 10 * h); hi = med + max(6 * rsd, 10 * h);
e = lo : bw : hi + bw;
n = histcounts(y, e);
k = exp(-0.5 * ((-25:25) / 5) .^ 2);
d = conv(n, k / sum(k), 'same');
m = max(d);
j = find(d == m);
pk = mean(e(j)) + bw / 2;
end
