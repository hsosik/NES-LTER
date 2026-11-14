function [beadstat, meta, classes, popsall] = measure_beads(dat, hdr, files, conf, sig)
% function [beadstat, meta, classes, popsall] = measure_beads(dat, hdr, files, conf, sig)
% Bead statistics for all FCB bead-check runs of one cruise, two passes:
%   1. find_beads with a wide window around the bead_signature positions
%      (scaled to each file's voltages); each population's median offset on
%      every channel over the good runs becomes the cruise's shift
%   2. find_beads with a narrow window around the shifted positions and the
%      full QC (count, width, position, signature)
%
% dat, hdr - cell arrays: events (numeric) and header (fca_readfcs) per run
% files    - cellstr file names, same order
% conf     - attune_config for the cruise
% sig      - bead_signature (optional)
%
% beadstat - table, one row per run:
%   filename, time (datetime), hv (table of detector voltage per physical
%   channel), and per population P = 05, 1, 6 (0.5, 1, 6 um):
%   nP      events in the gate
%   QCP     true = good
%   whyP    QC failure reasons ("" when good)
%   peakP   table, peak (KDE mode, bead_peak) per channel, linear units
%   medP    table, median per channel, linear units
%   rsdP    table, robust SD of log10 signal per channel
%   Channels are physical channels named as in the current filter layout
%   (FSC_H ... GL4_H, SSC_A, GL1_A), whatever the file parameter was called
%   (attune_config maps them; NaN where a filter was not installed).
% meta     - struct: cruise configuration, signature with the learned shifts,
%   per-population number of runs used for the shift
% classes  - cell array, bead class per event for each run (0 none, 1-3)
% popsall  - cell array, find_beads population details for each run (for plots)
%
% Nov 2026, bead stats redesign

if nargin < 5 || isempty(sig), sig = bead_signature(); end
chans = [strcat(sig.chans, '_H') {'SSC_A' 'GL1_A'}];
P = {'05' '1' '6'};
nr = numel(files);
ok = ~cellfun(@isempty, hdr);

% pass 1
sig1 = setpass(sig, sig.pass1);
off = nan(nr, numel(sig.chans), numel(sig.pop));
for f = find(ok(:))'
    [~, pops] = find_beads(double(dat{f}), hdr{f}, sig1, conf);
    for ip = 1:numel(pops)
        if pops(ip).n >= sig.qc.min_n && all(pops(ip).rsd <= sig.pop(ip).max_rsd)
            off(f, :, ip) = cellfun(@(c) pops(ip).chan_off.(c), sig.chans);
        end
    end
end
sig2 = setpass(sig, sig.pass2);
nshift = zeros(1, numel(sig.pop));
for ip = 1:numel(sig.pop)
    o = off(:, :, ip); o = o(any(isfinite(o), 2), :); % uninstalled channels stay NaN
    nshift(ip) = size(o, 1);
    st = struct();
    if ~isempty(o)
        m = median(o, 1, 'omitnan');
        for c = 1:numel(sig.chans), st.(sig.chans{c}) = m(c); end
    end
    sig2.pop(ip).shift = st;
end

% pass 2
beadstat = table;
classes = cell(nr, 1); popsall = cell(nr, 1);
for f = find(ok(:))'
    fcsdat = double(dat{f}); fcshdr = hdr{f};
    names = {fcshdr.par.name};
    [cl, pops] = find_beads(fcsdat, fcshdr, sig2, conf);
    classes{f} = uint8(cl); popsall{f} = pops;
    x = table;
    x.filename = string(files{f});
    x.time = datetime([fcshdr.date ' ' fcshdr.starttime]);
    hv = array2table(nan(1, numel(sig.chans)), 'VariableNames', sig.chans);
    for c = sig.chans
        j = strcmp(names, conf.chan.(c{1}));
        if any(j) && ~isempty(fcshdr.par(j).hv), hv.(c{1}) = fcshdr.par(j).hv; end
    end
    x.hv = hv;
    for ip = 1:numel(pops)
        pk = nan(1, numel(chans)); md = pk; sd = pk;
        for k = 1:numel(chans)
            j = strcmp(names, file_param(chans{k}, conf));
            if any(j) && any(cl == ip)
                [pk(k), md(k), sd(k)] = bead_peak(log10(max(fcsdat(cl == ip, j), 1)));
            end
        end
        x.(['n' P{ip}]) = pops(ip).n;
        x.(['QC' P{ip}]) = pops(ip).QC;
        x.(['why' P{ip}]) = string(pops(ip).why);
        x.(['peak' P{ip}]) = array2table(10 .^ pk, 'VariableNames', chans);
        x.(['med' P{ip}]) = array2table(10 .^ md, 'VariableNames', chans);
        x.(['rsd' P{ip}]) = array2table(sd, 'VariableNames', chans);
    end
    beadstat = [beadstat; x]; %#ok<AGROW>
end

meta.config = conf;
meta.signature = sig2;
meta.nshift = array2table(nshift, 'VariableNames', strcat('b', P));
meta.notes = ["peakP: KDE mode of log10 signal (bead_peak), linear units; medP: median, linear; rsdP: robust SD of log10"; ...
    "QCP true = good; channels are physical channels named in the current filter layout (attune_config)"];
end

function nm = file_param(c, conf)
% physical channel name (e.g. 'GL3_H', 'SSC_A') -> this cruise's file parameter
parts = split(c, '_');
ch = parts{1}; ah = parts{2};
nm = conf.chan.(ch);
if isempty(nm), return, end
nm = [extractBefore(nm, '-') '-' ah];
end

function s = setpass(s, p)
s.window = p.window; s.qc.max_offset = p.max_offset; s.sig_tol = p.sig_tol;
end
