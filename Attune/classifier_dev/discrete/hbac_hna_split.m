function [bt, s] = hbac_hna_split(BL1, SSC, isbact, p)
% function [bt, s] = hbac_hna_split(BL1, SSC, isbact, p)
% HNA / LNA split of heterotrophic bacteria (class 9) in one SYBR run.
% Called by assign_class_hbac (info.bact_type, like info.syn_type);
% settings p = cfg.hbac.hna (discrete_config).
%
% Cut on tilted SYBR, y = BL1 - p.slope * (SSC + 2) (bead-relative log10;
% the HNA cloud rises with SSC, and SSC -2 is near the LNA cluster, so y
% stays on the BL1 scale): the deepest density valley of the bacteria's y
% (ksdensity, bandwidth p.bw) inside p.win, where the dip
% ratio (valley density / lower of the two flanking maxima within p.flank)
% is below p.max_dip and there are at least p.min_n bacteria. Otherwise the
% fixed cut p.default_cut. (The window keeps the cut off the trigger-edge
% noise bump near BL1 -1.85 and off the bright tail.)
%
% bt - uint8 per event: 0 not bacteria, 1 LNA, 2 HNA
% s  - cut (on y), slope, method ('valley' or 'fixed'), dip, n_lna, n_hna, hna_frac
%
% Nov 2026, discrete-sample classifier development (HNA/LNA split)

if nargin < 4 || isempty(p), c = discrete_config(); p = c.hbac.hna; end
bt = zeros(numel(BL1), 1, 'uint8');
y = double(BL1) - p.slope * (double(SSC) + 2);
x = y(isbact & isfinite(y));
s = struct('cut', p.default_cut, 'method', 'fixed', 'dip', NaN, 'valley', NaN);
if numel(x) >= p.min_n
    xg = (p.win(1) - p.flank:0.01:p.win(2) + p.flank)';
    f = ksdensity(x, xg, 'bandwidth', p.bw);
    best = Inf; v = NaN;
    for k = find(xg >= p.win(1) & xg <= p.win(2))'
        if k == 1 || k == numel(xg) || f(k) > f(k - 1) || f(k) > f(k + 1), continue, end
        L = max(f(xg >= xg(k) - p.flank & xg < xg(k)));
        R = max(f(xg > xg(k) & xg <= xg(k) + p.flank));
        d = f(k) / min(L, R);
        if d < best, best = d; v = xg(k); end
    end
    s.dip = best; s.valley = v; if isinf(best), s.dip = NaN; end %no local minimum in the window
    if best < p.max_dip, s.cut = v; s.method = 'valley'; end
end
bt(isbact & y <= s.cut) = 1;
bt(isbact & y > s.cut) = 2;
s.n_lna = sum(bt == 1); s.n_hna = sum(bt == 2);
s.hna_frac = s.n_hna / max(1, s.n_lna + s.n_hna);
s.names = {'LNA' 'HNA'}; s.slope = p.slope;
end

