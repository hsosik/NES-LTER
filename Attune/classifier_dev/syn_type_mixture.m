function [st, M] = syn_type_mixture(ex, bm, fit, mix)
% function [st, M] = syn_type_mixture(ex, bm, fit, mix)
% Syn spectral type of each event from a per-file mixture on the excitation ratio
% ex = BL2-GL2, with the noise width set by brightness bm (mean of BL2 and GL2 above
% their floors). Fitted to the events in fit (logical, class 2 Syn); every event in ex
% then gets posteriors from that fit.
%   Types A-D are Gaussian along ex: position within mix.win of mix.mu0 (one per file,
%   shared by all brightness bands, so bright events place it), intrinsic sd in
%   mix.sgrid (per type and file), widened by band noise sn(bm)^2 = band_sd^2 - s0ref^2
%   (band_sd fitted jointly over 1,858 files, 10 cruises; output/syn_type_d/mix_proto).
%   Weights per brightness band (dim Syn can be another population than bright Syn),
%   pulled toward the file's overall weights by mix.alpha pseudo-events. A broad
%   component (mean mix.bg_mu, sd mix.bg_sd plus noise, at most mix.bg_max per band)
%   takes outliers.
%   Types are dropped one at a time (backward elimination) while dropping costs less
%   than mix.pen * log(n) in log-likelihood; dropped first, regardless of cost: types
%   with share below mix.wmin, and types under mix.shoulder_max whose position is
%   pinned at the limit facing a kept neighbour (a shoulder of that neighbour).
% st: 1-4 = type A-D where the most probable kept type has posterior >= mix.pmin,
%   5 = ambiguous (no type that probable).
% M: method, n_fit, keep (1x4), mu, s, w (overall weights), post (n x 4 posteriors of
%   every event), path (dropped type and cost / log n at each elimination step).
%
% Nov 2026, classifier_dev/output/syn_type_d (mix_proto4.m is the prototype)

ex = ex(:); bm = bm(:); fit = fit(:);
xe = mix.xbins(:); xc = (xe(1:end-1) + xe(2:end)) / 2;
be = [-Inf mix.band_edges(2:end) Inf]; nq = numel(be) - 1;
sn2 = max(mix.band_sd(:)'.^2 - mix.s0ref^2, 1e-6);
% binned histogram of the fit events: ratio (clamped into range) x band
xf = min(max(ex(fit), xe(1) + 1e-9), xe(end) - 1e-9);
H = histcounts2(xf, bm(fit), xe, be); % nx x nq
n = sum(H(:));
M = struct('method', 'mixture', 'n_fit', n, 'keep', true(1, 4), 'mu', nan(1, 4), 's', nan(1, 4), 'w', zeros(1, 4), 'post', [], 'path', {{}});
kk = true(1, 4); F = fitmix(H, xc, sn2, kk, mix);
while sum(kk) > 1
    best = Inf; bi = 0;
    for t = find(kk)
        k2 = kk; k2(t) = false; F2 = fitmix(H, xc, sn2, k2, mix);
        d = F.L - F2.L;
        forced = F.share(t) < mix.wmin || (F.share(t) < mix.shoulder_max && pinned(F.mu, t, kk, mix));
        if forced, d = -Inf; end
        if d < best, best = d; bi = t; B2 = F2; end
    end
    if best >= mix.pen * log(n), break, end
    tn = 'ABCD'; M.path(end + 1, :) = {tn(bi), best / log(n)};
    kk(bi) = false; F = B2;
end
M.keep = kk; M.mu = F.mu; M.mu(~kk) = NaN; M.s = F.s; M.s(~kk) = NaN; M.w = F.share;
% posteriors of every event (band weights of its own band)
q = discretize(bm, be); q(isnan(q)) = 1;
lp = -Inf(numel(ex), 4);
for t = find(kk)
    v = F.s(t)^2 + sn2(q)';
    lp(:, t) = log(max(F.W(t, q)', 1e-300)) - 0.5 * log(v) - (ex - F.mu(t)).^2 ./ (2 * v);
end
lp = lp - max(lp, [], 2); p = exp(lp); p = p ./ sum(p, 2);
p(~isfinite(p)) = 0;
M.post = p;
[pm, a] = max(p, [], 2);
st = uint8(a); st(pm < mix.pmin) = 5;
end

function F = fitmix(H, xc, sn2, kk, mix)
% EM: band weights (MAP toward the file weights), shared positions, per-type widths
nq = size(H, 2); N = sum(H(:)); Nq = sum(H, 1);
mu = mix.mu0(:)'; s = mix.s0ref * ones(1, 4); W = repmat([kk / sum(kk) * 0.98, 0.02]', 1, nq);
Lold = -Inf;
for it = 1:150
    [L, R] = estep(H, xc, sn2, W, mu, s, kk, mix);
    Nqt = squeeze(sum(H .* R, 1)); if nq == 1, Nqt = Nqt(:); else, Nqt = Nqt'; end % 5 x nq
    wbar = sum(Nqt, 2) / N; wbar(1:4) = wbar(1:4) .* kk(:);
    W = (Nqt + mix.alpha * wbar) ./ (Nq + mix.alpha); W(~[kk true], :) = 0;
    over = W(5, :) > mix.bg_max;
    if any(over), W(1:4, over) = W(1:4, over) .* (1 - mix.bg_max) ./ max(sum(W(1:4, over), 1), eps); W(5, over) = mix.bg_max; end
    for t = find(kk)
        HR = H .* R(:, :, t); if sum(HR(:)) <= 0, continue, end
        pw = HR ./ (s(t)^2 + sn2);
        mu(t) = min(max(sum(pw .* xc, 'all') / sum(pw(:)), mix.mu0(t) - mix.win(t)), mix.mu0(t) + mix.win(t));
        best = -Inf;
        for sg = mix.sgrid
            vv = sg^2 + sn2; e = sum(HR .* (-0.5 * log(vv) - (xc - mu(t)).^2 ./ (2 * vv)), 'all');
            if e > best, best = e; s(t) = sg; end
        end
    end
    if L - Lold < 1e-6 * N, break, end
    Lold = L;
end
F.L = L; F.W = W; F.mu = mu; F.s = s; F.share = sum(Nqt(1:4, :), 2)' / N;
end

function [L, R] = estep(H, xc, sn2, W, mu, s, kk, mix)
nq = size(H, 2); lp = -Inf(numel(xc), nq, 5);
for t = find(kk)
    v = s(t)^2 + sn2;
    lp(:, :, t) = log(max(W(t, :), 1e-300)) - 0.5 * log(2 * pi * v) - (xc - mu(t)).^2 ./ (2 * v);
end
vb = mix.bg_sd^2 + sn2;
lp(:, :, 5) = log(max(W(5, :), 1e-300)) - 0.5 * log(2 * pi * vb) - (xc - mix.bg_mu).^2 ./ (2 * vb);
mx = max(lp, [], 3); ss = sum(exp(lp - mx), 3); R = exp(lp - mx) ./ ss;
L = sum(H .* (mx + log(ss)), 'all');
end

function p = pinned(mu, t, kk, mix)
% position at its limit on the side of a kept neighbour
tol = 0.006; p = false;
if abs(mu(t) - (mix.mu0(t) - mix.win(t))) < tol && any(kk(1:t - 1)), p = true; end
if abs(mu(t) - (mix.mu0(t) + mix.win(t))) < tol && any(kk(t + 1:end)), p = true; end
end
