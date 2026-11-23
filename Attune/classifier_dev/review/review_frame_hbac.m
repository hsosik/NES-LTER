function nrow = review_frame_hbac(cls, info, outfile, R)
% function nrow = review_frame_hbac(cls, info, outfile, R)
% Review image for one SYBR-stained bacteria run (assign_class_hbac output), one
% row of density panels drawn into pixels on the fixed axes of
% review_config_hbac. Bacteria are drawn as LNA / HNA (info.bact_type), then
% Euk, Syn and Pro on top. The histogram panel shows the bacteria on the tilted
% SYBR axis of the HNA/LNA split. Axes, labels and the class boundaries are drawn
% by the viewer over the image (rule lines from review_qc_hbac).
%
% Nov 2026, preserved-sample review (SYBR runs)

if nargin < 4 || isempty(R), R = review_config_hbac(); end
F = info.features;
cls = double(cls(:));
bt = zeros(size(cls)); if isfield(info, 'bact_type'), bt = double(info.bact_type(:)); end
c = cls; c(cls == 9 & bt == 1) = 91; c(cls == 9 & bt == 2) = 92; c(cls == 9 & ~ismember(bt, [1 2])) = 91;
slope = 0.35; if isfield(info, 'hna') && isfield(info.hna, 'slope'), slope = info.hna.slope; end
F.hnay = F.BL1 - slope * (F.SSC + 2);
nv = size(R.views, 1); W = R.px; g = R.gap;
img = ones(W, nv * W + (nv - 1) * g, 3);
order = [91 92 1 2 7];
for v = 1:nv
    P = ones(W, W, 3);
    lx = R.lim.(R.views{v, 1}); ly = R.lim.(R.views{v, 2});
    if strcmp(R.views{v, 2}, 'hist')
        P = hist_panel(P, F.hnay, c, lx, ly, R);
    else
        [ix, ok] = binxy(F.(R.views{v, 1}), F.(R.views{v, 2}), lx, ly, W);
        P = dens(P, ix, ok & c == 0, R.col0, .5);
        for k = order
            P = dens(P, ix, ok & c == k, R.cols(R.codes == k, :), .9);
        end
    end
    x0 = (v - 1) * (W + g);
    img(1:W, x0 + (1:W), :) = P;
end
nrow = 1;
img = uint8(round(img * 255));
if strcmpi(R.format, 'jpg')
    imwrite(img, outfile, 'Quality', 90)
else
    [X, map] = rgb2ind(img, R.ncolors, 'nodither');
    imwrite(X, map, outfile)
end
end

function P = hist_panel(P, y, c, lx, ly, R)
% bacteria on the tilted SYBR axis: one column per pixel, height = count / max,
% coloured LNA or HNA by the majority in the column
W = size(P, 1);
isb = ismember(c, [91 92]) & isfinite(y);
col = 1 + floor((y(isb) - lx(1)) / diff(lx) * W); in = col >= 1 & col <= W;
cb = c(isb); col = col(in); cb = cb(in);
nl = accumarray(col, cb == 91, [W 1]); nh = accumarray(col, cb == 92, [W 1]);
n = nl + nh;
if ~any(n), return, end
n = movmean(n, 3); nl = movmean(nl, 3); nh = movmean(nh, 3);
h = n / max(n); % relative density
for k = 1:W
    if n(k) <= 0, continue, end
    top = max(1, W - round((h(k) - ly(1)) / diff(ly) * W) + 1);
    if nh(k) > nl(k), cc = R.cols(R.codes == 92, :); else, cc = R.cols(R.codes == 91, :); end
    for ch = 1:3, P(top:W, k, ch) = cc(ch); end
end
end

function [ix, ok] = binxy(x, y, lx, ly, W)
ok = isfinite(x) & isfinite(y);
c = 1 + floor((x - lx(1)) / diff(lx) * W); c = min(max(c, 1), W);
r = W - floor((y - ly(1)) / diff(ly) * W); r = min(max(r, 1), W);
ix = r + (c - 1) * W;
ix(~ok) = 1;
end

function P = dens(P, ix, m, col, amax, amin)
if ~any(m), return, end
if nargin < 6, amin = 0.3; end
W = size(P, 1);
n = accumarray(ix(m), 1, [W * W 1]);
a = zeros(W * W, 1); k = n > 0;
a(k) = amax * (amin + (1 - amin) * log1p(n(k)) / log1p(max(n)));
a = reshape(a, W, W);
for ch = 1:3
    P(:, :, ch) = P(:, :, ch) .* (1 - a) + col(ch) * a;
end
end
