function P = pt_test_peaks(d, h, conf, chans)
% function P = pt_test_peaks(d, h, conf, chans)
% PT-bead peaks of one Attune performance test (as used by pt_perf_table and
% the bead review): singlets = central cluster on FSC-H x SSC-H (within 0.15
% log10 of the medians); per physical channel the log10 peak = median of the
% singlets (FSC, SSC: the single scatter peak) or of the brightest PT
% intensity level (fluorescence: singlets above the 70th percentile of log10
% signal), and the robust CV (IQR/1.349/median) of those events.
%
% d, h   - fca_readfcs output (d double)
% conf   - attune_config for the test's date (physical channel -> file parameter)
% chans  - physical channels (default FSC SSC BL1 BL2 BL3 GL1 GL2 GL3 GL4)
%
% P - struct: sing (logical, singlets), and per channel P.ch.(ch) with
%     col (column in d, [] if absent), hv, pk (log10), rcv, use (logical,
%     the events the peak is taken from)
%
% Nov 2026, factored out of pt_perf_table for the bead review

if nargin < 4 || isempty(chans), chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'}; end
names = {h.par.name};
fs = log10(max(d(:, strcmp(names, 'FSC-H')), 1));
ss = log10(max(d(:, strcmp(names, 'SSC-H')), 1));
g = abs(fs - median(fs)) < 0.15 & abs(ss - median(ss)) < 0.15;
P.sing = g;
for j = 1:numel(chans)
    ch = chans{j}; p = conf.chan.(ch);
    c = []; if ~isempty(p), c = find(strcmp(names, p), 1); end
    q = struct('col', c, 'hv', NaN, 'pk', NaN, 'rcv', NaN, 'use', false(size(g)));
    if ~isempty(c)
        q.hv = h.par(c).hv;
        x = d(:, c);
        if ismember(ch, {'FSC' 'SSC'})
            u = g;
        else % brightest PT intensity level
            lx = log10(max(x, 1));
            u = g & lx > prctile(lx(g), 70);
        end
        v = x(u);
        q.use = u;
        q.pk = log10(median(v));
        q.rcv = iqr(v) / 1.349 / median(v);
    end
    P.ch.(ch) = q;
end
end
