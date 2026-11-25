function nrow = review_frame(fcsdat, fcshdr, cls, info, old, outfile, R)
% function nrow = review_frame(fcsdat, fcshdr, cls, info, old, outfile, R)
% Write one review image for a sample: a grid of cytogram panels drawn straight into
% pixels (no MATLAB figure), on the fixed axes of review_config. Abundant classes are
% 2D-histogram density images, as are flagged events (class 0, by flag); rare classes
% (dim/bright crypto, Euk with Syn, Syn on particle) are drawn as outlined points on top.
%   default: one row, assign_class_v2 classes
%   R.compare_rows true and old given: three rows, old classes / v2 classes / disagreements
% Axes, labels and rule lines are not in the image: the viewer draws them over it.
% Returns the number of rows drawn.
%
% Nov 2026, classification review prototype

if nargin < 7 || isempty(R), R = review_config(); end
F = info.features;
pn = {fcshdr.par.name};
j = find(strcmp(pn, 'FSC-W'), 1);
if isempty(j), F.FSCW = nan(height(F), 1); else, F.FSCW = log10(max(fcsdat(:, j), 1)); end
cls = double(cls(:));
eb = review_euk_line(F, cls);
F.dEuk = F.BL3 - (eb(1) + eb(2) * F.size); % chl relative to the Euk chl-size line (NaN without a fit)
% Pro-colour panels (chl:size, GL2:GL3, chl pulse shape): fill any feature an older attune_features lacks
if ~ismember('chlSSC', F.Properties.VariableNames), F.chlSSC = F.BL3 - F.size; end
if ~ismember('GL2GL3', F.Properties.VariableNames), F.GL2GL3 = F.GL2 - F.GL3; end
if ~ismember('BL3_AH', F.Properties.VariableNames), F.BL3_AH = nan(height(F), 1); end
st = zeros(size(cls)); if isfield(info, 'syn_type'), st = double(info.syn_type(:)); end
fl = zeros(size(cls)); if isfield(info, 'flag'), fl = double(info.flag(:)); end

nv = size(R.views, 1); W = R.px; g = R.gap;
if isempty(old) || ~R.compare_rows, rows = {'v2'}; else, rows = {'old' 'v2' 'diff'}; end
nrow = numel(rows);
img = ones(nrow * W + (nrow - 1) * g, nv * W + (nv - 1) * g, 3);
if ~isempty(old)
    old = double(old(:)); old(~ismember(old, [0 R.codes])) = 0;
end
for r = 1:nrow
    for v = 1:nv
        lx = R.lim.(R.views{v, 1}); ly = R.lim.(R.views{v, 2});
        if isfield(R, 'view_lim') && ~isempty(R.view_lim{v}), lx = R.view_lim{v}(1, :); ly = R.view_lim{v}(2, :); end
        [ix, ok] = binxy(F.(R.views{v, 1}), F.(R.views{v, 2}), lx, ly, W);
        P = ones(W, W, 3);
        switch rows{r}
            case {'old' 'v2'}
                if strcmp(rows{r}, 'old'), c = old; else, c = cls; end
                if strcmp(rows{r}, 'v2') && v == nv
                    % Syn type view: all Syn-bearing events (Syn, Euk with Syn, Syn on
                    % particle) coloured by type, everything else light grey beneath
                    synb = ismember(c, [2 5 8]);
                    P = dens(P, ix, ok & ~synb, [.6 .6 .6], .45);
                    P = dens(P, ix, ok & synb & ~ismember(st, 1:3), R.syn_untyped_col, .9, .5);
                    [~, o] = sort(arrayfun(@(t) sum(synb & st == t), 1:3), 'descend');
                    for t = o, P = dens(P, ix, ok & synb & st == t, R.syn_type_cols(t, :), 1, .6); end
                    y0 = (r - 1) * (W + g); x0 = (v - 1) * (W + g);
                    img(y0 + (1:W), x0 + (1:W), :) = P;
                    continue
                end
                if strcmp(rows{r}, 'v2') && isfield(R, 'show_flags') && R.show_flags
                    % flagged events (class 0) as their own density layers
                    P = dens(P, ix, ok & c == 0 & fl == 0, R.col0, .55);
                    for t = 1:size(R.flag_cols, 1), P = dens(P, ix, ok & c == 0 & fl == t, R.flag_cols(t, :), .8); end
                else
                    P = dens(P, ix, ok & c == 0, R.col0, .55);
                end
                ab = setdiff(R.codes, R.rare, 'stable');
                [~, o] = sort(arrayfun(@(k) sum(c == k), ab), 'descend');
                for k = ab(o)
                    P = dens(P, ix, ok & c == k, R.cols(R.codes == k, :), .9);
                end
                for k = R.rare
                    P = pts(P, ix, ok & c == k, R.cols(R.codes == k, :), R.rare_rad, R.rare_outline);
                end
            case 'diff'
                P = dens(P, ix, ok & ~(old == cls & old ~= 0), [.6 .6 .6], .35);
                P = dens(P, ix, ok & old == cls & old ~= 0, [.35 .35 .35], .6);
                d = {old == 0 & cls ~= 0, old ~= 0 & cls == 0, old ~= 0 & cls ~= 0 & old ~= cls};
                for t = 1:3, P = pts(P, ix, ok & d{t}, R.diff_cols(t, :), 1, false); end
        end
        y0 = (r - 1) * (W + g); x0 = (v - 1) * (W + g);
        img(y0 + (1:W), x0 + (1:W), :) = P;
    end
end
img = uint8(round(img * 255));
if strcmpi(R.format, 'jpg')
    imwrite(img, outfile, 'Quality', 90)
else
    % indexed colour (no dither) keeps PNGs small; densities need few distinct colours
    [X, map] = rgb2ind(img, R.ncolors, 'nodither');
    imwrite(X, map, outfile)
end
end

function [ix, ok] = binxy(x, y, lx, ly, W)
% pixel index (column-major, row 1 at top) of each event; values outside the axis
% range are clamped to the edge pixel
ok = isfinite(x) & isfinite(y);
c = 1 + floor((x - lx(1)) / diff(lx) * W); c = min(max(c, 1), W);
r = W - floor((y - ly(1)) / diff(ly) * W); r = min(max(r, 1), W);
ix = r + (c - 1) * W;
ix(~ok) = 1;
end

function P = dens(P, ix, m, col, amax, amin)
% composite one class as a density image: opacity rises with log count per pixel,
% from amin * amax for a single event to amax at the densest pixel
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

function P = pts(P, ix, m, col, rad, outline)
% stamp each event as a (2*rad+1)^2 square, with a dark one-pixel outline if asked
if ~any(m), return, end
W = size(P, 1);
[r, c] = ind2sub([W W], unique(ix(m)));
if outline, P = stamp(P, r, c, rad + 1, [.1 .1 .1]); end
P = stamp(P, r, c, rad, col);
end

function P = stamp(P, r, c, rad, col)
W = size(P, 1);
[dr, dc] = meshgrid(-rad:rad);
rr = min(max(r + dr(:)', 1), W); cc = min(max(c + dc(:)', 1), W);
k = rr(:) + (cc(:) - 1) * W;
for ch = 1:3, P(k + (ch - 1) * W * W) = col(ch); end
end
