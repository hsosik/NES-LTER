function pt = find_pt_beads(fcsdat, fcshdr, conf, ref)
% function pt = find_pt_beads(fcsdat, fcshdr, conf, ref)
% Attune performance-tracking (PT) beads in one file: a single bead size with
% three fluorescence intensity levels (~1.4 log10 apart) on each channel.
%   1. gate the bead cluster: densest point on FSC-H x SSC-H among events
%      with both > 10^3.5, keep events within 0.2 log10 of it
%   2. FSC and SSC: KDE peak (bead_peak) of the gated events
%   3. fluorescence: peaks of the smoothed log10 histogram (prominence >= 8%
%      of the highest, >= 0.25 apart) between 10^1.8 and saturation; with ref,
%      each peak is assigned to the nearest reference level within 0.35
%      (generate_pt_beadstats labels levels by aligning the whole pattern)
%
% conf - attune_config (physical channel -> file parameter)
% ref  - optional struct, ref.(ch) = [L1 L2 L3] expected log10 levels at this
%        file's settings; without it only the raw peaks (pt.peaks) are returned
%
% pt.n          events in the gate
% pt.FSC, pt.SSC  log10 peak; pt.FSCrsd, pt.SSCrsd robust SD of log10
% pt.levels.(ch)  [L1 L2 L3] log10 (NaN where not found, saturated or in noise)
% pt.peaks.(ch)   all peaks found (log10, ascending)
% pt.hv.(ch)      detector voltage
%
% Nov 2026, bead stats redesign (PT-only cruises)

chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
names = {fcshdr.par.name};
L = @(nm) log10(max(fcsdat(:, strcmp(names, nm)), 1));
for c = chans
    j = strcmp(names, conf.chan.(c{1}));
    pt.hv.(c{1}) = NaN;
    if any(j) && ~isempty(fcshdr.par(j).hv), pt.hv.(c{1}) = fcshdr.par(j).hv; end
end
fs = L(conf.chan.FSC); ss = L(conf.chan.SSC);
ok = fs > 3.5 & ss > 3.5;
pt.n = 0; pt.FSC = NaN; pt.SSC = NaN; pt.FSCrsd = NaN; pt.SSCrsd = NaN;
for c = chans(3:end), pt.levels.(c{1}) = nan(1, 3); pt.peaks.(c{1}) = []; end
if sum(ok) < 200, return, end
e = 3.5:0.02:6.1;
H = histcounts2(fs(ok), ss(ok), e, e);
H = conv2(H, ones(3) / 9, 'same');
[~, im] = max(H(:)); [i1, i2] = ind2sub(size(H), im);
g = ok & abs(fs - e(i1) - 0.01) < 0.2 & abs(ss - e(i2) - 0.01) < 0.2;
pt.n = sum(g);
[pt.FSC, ~, pt.FSCrsd] = bead_peak(fs(g)); [pt.SSC, ~, pt.SSCrsd] = bead_peak(ss(g));
sat = log10(1048575) - 0.05;
ee = 1.0:0.01:6.1; xc = ee(1:end-1) + 0.005;
k = exp(-0.5 * ((-6:6) / 2) .^ 2); k = k / sum(k);
for c = chans(3:end)
    nm = conf.chan.(c{1});
    if isempty(nm), continue, end
    y = L(nm); y = y(g);
    h = conv(histcounts(y, ee), k, 'same');
    [~, pl] = findpeaks(h, xc, 'MinPeakProminence', max(h) * 0.08, 'MinPeakDistance', 0.25);
    pl = pl(pl > 1.8 & pl < sat);
    pt.peaks.(c{1}) = pl;
    if nargin > 3 && isfield(ref, c{1})
        r = ref.(c{1});
        for lv = 1:3
            if ~isfinite(r(lv)), continue, end
            [dmin, jj] = min(abs(pl - r(lv)));
            if ~isempty(dmin) && dmin < 0.35, pt.levels.(c{1})(lv) = pl(jj); end
        end
    end
end
end
