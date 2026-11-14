function [class, pops] = find_beads(fcsdat, fcshdr, sig, config)
% function [class, pops] = find_beads(fcsdat, fcshdr, sig, config)
% Find the 0.5, 1 and 6 um FCB bead populations in one bead-check file.
% For each population:
%   1. predict its log10 position on every channel from the bead signature
%      (sig, from bead_signature), the file's voltages and, when set, the
%      cruise's own offsets (pop.shift, learned by a first pass)
%   2. keep events within sig.window of the prediction on the gating
%      channels, and below the limit channel (drops bead aggregates)
%   3. take the densest point (smoothed 2-D histogram) as the seed
%   4. robust ellipse: MCD covariance of events near the seed, keep events
%      within the 99.5% Mahalanobis contour, refit once
%   5. QC: enough events, tight (robust SD), close to the prediction on the
%      gating channels, and (peak) on its signature channels (pop.sigch)
%
% fcsdat     - numeric events x parameters (as from fca_readfcs)
% fcshdr     - header from fca_readfcs (par.name, par.hv)
% config     - from attune_config: OD2 placement and which file parameter
%              holds each physical channel (filter layout changed Aug 2019);
%              a char 'GL1'/'SSC'/'None' is taken as the OD2 placement with
%              the current layout
%
% class - events x 1: 0 not bead, 1 = 0.5 um, 2 = 1 um, 3 = 6 um
% pops  - struct array, one per population, with gating details and QC;
%         chan_off: median - expected (before shift) on every channel
%
% Nov 2026, bead stats redesign prototype

if nargin < 4 || isempty(config), config = 'GL1'; end
if ischar(config) || isstring(config)
    od2 = char(config); config = attune_config('');
    config.OD2 = od2;
end
names = {fcshdr.par.name};
hv = nan(1, numel(names));
for j = 1:numel(names)
    if ~isempty(fcshdr.par(j).hv), hv(j) = fcshdr.par(j).hv; end
end
class = zeros(size(fcsdat, 1), 1);
pops = struct([]);
for ip = 1:numel(sig.pop)
    p = sig.pop(ip);
    r = struct('name', p.name, 'gate', {p.gate}, 'expected', nan(1, 2), ...
        'center', nan(1, 2), 'rsd', nan(1, 2), 'offset', nan(1, 2), ...
        'nwindow', 0, 'n', 0, 'found', false, 'QC', false, 'why', "", ...
        'chan_off', struct(), 'gatech', {{}});
    % expected position on each gating channel at this file's voltages
    gate = p.gate;
    [~, e0] = expect(gate, p, names, hv, sig, config, false);
    if (any(e0 > sig.saturation - 0.2) || any(e0 < sig.floor) || any(isnan(e0))) && ~isempty(p.alt)
        gate = p.alt; r.gate = gate;
    end
    [cols, r.expected] = expect(gate, p, names, hv, sig, config, true);
    if any(cols == 0)
        r.why = "gating channel not installed";
        pops = [pops r]; %#ok<AGROW>
        continue
    end
    r.gatech = names(cols);
    X = log10(max(fcsdat(:, cols), 1));
    win = all(abs(X - r.expected) < sig.window, 2) & class == 0;
    if ~isempty(p.limit)
        [lc, le] = expect(p.limit(1), p, names, hv, sig, config, true);
        if lc > 0, win = win & log10(max(fcsdat(:, lc), 1)) < le + p.limit{2}; end
    end
    r.nwindow = sum(win);
    if r.nwindow < sig.qc.min_n / 2
        r.why = "too few events near expected position";
        pops = [pops r]; %#ok<AGROW>
        continue
    end
    % seed: densest point in the window
    bw = 0.02;
    e1 = r.expected(1) - sig.window : bw : r.expected(1) + sig.window;
    e2 = r.expected(2) - sig.window : bw : r.expected(2) + sig.window;
    H = histcounts2(X(win, 1), X(win, 2), e1, e2);
    g = exp(-((-2:2)' .^ 2 + (-2:2) .^ 2) / 2);
    H = conv2(H, g / sum(g(:)), 'same');
    [~, im] = max(H(:));
    [i1, i2] = ind2sub(size(H), im);
    seed = [e1(i1) e2(i2)] + bw / 2;
    % robust ellipse around the seed
    near = win & all(abs(X - seed) < 0.12, 2);
    if sum(near) < 20
        r.why = "no dense cluster in window";
        pops = [pops r]; %#ok<AGROW>
        continue
    end
    cut = chi2inv(0.995, 2);
    [S, mu] = robustcov(X(near, :));
    in = win & mahal_d2(X, mu, S) < cut;
    if sum(in) >= 20
        [S, mu] = robustcov(X(in, :));
        in = win & mahal_d2(X, mu, S) < cut;
    end
    r.n = sum(in);
    r.center = mu;
    r.rsd = sqrt(diag(S))';
    r.offset = mu - r.expected;
    r.found = r.n > 0;
    class(in) = ip;
    % position on every channel, relative to the unshifted signature
    bad = strings(0);
    for c = sig.chans
        [col, e] = expect(c, p, names, hv, sig, config, false);
        if col == 0, r.chan_off.(c{1}) = NaN; continue, end
        m = bead_peak(log10(max(fcsdat(in, col), 1)));
        r.chan_off.(c{1}) = m - e;
        es = e + shift_of(p, c{1});
        if ismember(c{1}, p.sigch) && es >= sig.sig_range(1) && es <= sig.sig_range(2) ...
                && abs(m - es) > sig.sig_tol
            bad(end+1) = string(c{1}); %#ok<AGROW>
        end
    end
    % QC
    why = strings(0);
    if r.n < sig.qc.min_n, why(end+1) = "n<" + sig.qc.min_n; end %#ok<AGROW>
    if any(r.rsd > p.max_rsd), why(end+1) = "broad"; end %#ok<AGROW>
    if any(abs(r.offset) > sig.qc.max_offset), why(end+1) = "off position"; end %#ok<AGROW>
    if ~isempty(bad), why(end+1) = "signature " + strjoin(bad, "/"); end %#ok<AGROW>
    r.QC = isempty(why);
    r.why = strjoin(why, ", ");
    pops = [pops r]; %#ok<AGROW>
end
end

function [cols, e] = expect(gate, p, names, hv, sig, config, shifted)
% file column and expected log10 position of each physical channel in gate;
% column 0 and NaN where that channel was not installed
cols = zeros(1, numel(gate)); e = nan(1, numel(gate));
for d = 1:numel(gate)
    ch = gate{d};
    j = find(strcmp(names, config.chan.(ch)), 1);
    if isempty(j), continue, end
    cols(d) = j;
    e(d) = p.pos.(ch) + sig.k.(ch) * log10(hv(j) / sig.hvref.(ch)) + od2_shift(ch, config.OD2);
    if shifted, e(d) = e(d) + shift_of(p, ch); end
end
end

function s = shift_of(p, ch)
s = 0;
if isstruct(p.shift) && isfield(p.shift, ch) && isfinite(p.shift.(ch)), s = p.shift.(ch); end
end

function d2 = mahal_d2(X, mu, S)
Z = X - mu;
d2 = sum((Z / S) .* Z, 2);
end

function s = od2_shift(ch, OD2setting)
% bead signature positions assume OD2 on GL1; shift if it was elsewhere
s = 0;
switch upper(OD2setting)
    case 'SSC'
        if strcmp(ch, 'GL1'), s = 2; elseif strcmp(ch, 'SSC'), s = -2; end
    case 'NONE'
        if strcmp(ch, 'GL1'), s = 2; end
end
end
