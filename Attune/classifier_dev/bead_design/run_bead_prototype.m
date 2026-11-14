function R = run_bead_prototype(outfile, cruises)
% function R = run_bead_prototype(outfile, cruises)
% Run find_beads on every bead-check file of every cruise and tabulate, per
% run and population: events found, gating centre, robust SD, offset from the
% expected position, QC, and for every height channel (plus SSC-A, GL1-A) the
% log10 peak (KDE mode, bead_peak), median and robust SD.
%
% Two passes per cruise:
%   1. wide window around the bead_signature position (scaled to each file's
%      voltages); the cruise's median offset of found from expected becomes
%      that population's shift
%   2. narrow window around the shifted position, with a tight position QC
% Old beadstat values (2026: medians, centend2 modes, QC; 2021: QC, counts)
% are added for comparison, matched by file name.
%
% cruises - optional cellstr of cruise folder names (default: all)
%
% Nov 2026, bead stats redesign prototype

dev = fileparts(fileparts(mfilename('fullpath'))); % classifier_dev
addpath(dev, fullfile(dev, '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
sig = bead_signature();
chans = [strcat(sig.chans, '-H') {'SSC-A' 'GL1-A'}];
if nargin < 2 || isempty(cruises)
    D = dir([root '20*_*']); D = D([D.isdir]); cruises = {D.name};
end
rows = {};
for ic = 1:numel(cruises)
    cdir = cruises{ic};
    cruise = extractAfter(cdir, '_');
    F = dir(fullfile(root, cdir, 'FCS', '*.fcs'));
    F = F(contains({F.name}, 'bead', 'IgnoreCase', true));
    if isempty(F), continue, end
    conf = attune_config(cruise);
    if ~conf.fcb_beads, fprintf('%s: no FCB bead runs (configuration history)\n', cdir); continue, end
    old = load_old(fullfile(root, cdir, 'bead_calibrated'));
    fprintf('%s: %d bead files\n', cdir, numel(F));
    % read once, keep for both passes
    dat = cell(numel(F), 1); hdr = cell(numel(F), 1);
    for f = 1:numel(F)
        try
            [d, h] = fca_readfcs(fullfile(F(f).folder, F(f).name));
            dat{f} = single(d); hdr{f} = h;
        catch err
            fprintf('  read failed %s: %s\n', F(f).name, err.message)
        end
    end
    ok = ~cellfun(@isempty, hdr);
    % pass 1: wide search, learn the cruise's offsets
    sig1 = setpass(sig, sig.pass1);
    off = nan(numel(F), numel(sig.chans), numel(sig.pop));
    for f = find(ok)'
        [~, pops] = find_beads(double(dat{f}), hdr{f}, sig1, conf);
        for ip = 1:numel(pops)
            if pops(ip).n >= sig.qc.min_n && all(pops(ip).rsd <= sig.pop(ip).max_rsd)
                off(f, :, ip) = cellfun(@(c) pops(ip).chan_off.(c), sig.chans);
            end
        end
    end
    sig2 = setpass(sig, sig.pass2);
    nshift = zeros(numel(sig.pop), 1); shift = cell(numel(sig.pop), 1);
    for ip = 1:numel(sig.pop)
        o = off(:, :, ip); o = o(any(isfinite(o), 2), :); % uninstalled channels stay NaN
        nshift(ip) = size(o, 1);
        st = struct();
        if ~isempty(o)
            m = median(o, 1, "omitnan");
            for c = 1:numel(sig.chans), st.(sig.chans{c}) = m(c); end
        end
        sig2.pop(ip).shift = st; shift{ip} = st;
    end
    % pass 2
    for f = find(ok)'
        fcsdat = double(dat{f}); fcshdr = hdr{f};
        [class, pops] = find_beads(fcsdat, fcshdr, sig2, conf);
        names = {fcshdr.par.name};
        x = struct();
        x.cruise = string(cruise);
        x.file = string(F(f).name);
        x.time = datenum(datetime([fcshdr.date ' ' fcshdr.starttime]));
        x.OD2 = string(conf.OD2); x.layout = string(conf.layout);
        x.nevents = size(fcsdat, 1);
        for c = sig.chans
            j = strcmp(names, conf.chan.(c{1}));
            x.(['hv_' c{1}]) = NaN;
            if any(j) && ~isempty(fcshdr.par(j).hv), x.(['hv_' c{1}]) = fcshdr.par(j).hv; end
        end
        for ip = 1:numel(pops)
            p = pops(ip); b = p.name;
            for c = sig.chans
                v = NaN; if isfield(shift{ip}, c{1}), v = shift{ip}.(c{1}); end
                x.([b '_shift_' c{1}]) = v;
            end
            x.([b '_nshift']) = nshift(ip);
            x.([b '_n']) = p.n; x.([b '_nwin']) = p.nwindow;
            x.([b '_gate']) = string(strjoin(p.gate, '/'));
            x.([b '_exp']) = p.expected; x.([b '_ctr']) = p.center;
            x.([b '_rsd']) = p.rsd; x.([b '_QC']) = p.QC; x.([b '_why']) = string(p.why);
            for c = chans
                % physical channel (current-layout name) -> this file's parameter
                fn = strrep(c{1}, '-', '_');
                nm = c{1};
                if endsWith(nm, '-H'), nm = conf.chan.(extractBefore(nm, '-H')); end
                if strcmp(c{1}, 'GL1-A') && isempty(conf.chan.GL1), nm = ''; end
                j = strcmp(names, nm);
                pk = NaN; md = NaN; s = NaN;
                if any(j) && any(class == ip)
                    [pk, md, s] = bead_peak(log10(max(fcsdat(class == ip, j), 1)));
                end
                x.([b '_pk_' fn]) = pk; x.([b '_md_' fn]) = md; x.([b '_rsd_' fn]) = s;
            end
        end
        x = add_old(x, old, F(f).name, sig.chans);
        rows{end+1, 1} = x; %#ok<AGROW>
    end
    R = struct2table([rows{:}]');
    if nargin > 0 && ~isempty(outfile), save(outfile, 'R', 'sig'), end
end
end

function s = setpass(s, p)
s.window = p.window; s.qc.max_offset = p.max_offset; s.sig_tol = p.sig_tol;
end

function x = add_old(x, old, fname, chans)
x.old_QC = NaN; x.old_gen = "";
for b = 1:3, x.(sprintf('old_n%d', b)) = NaN; end
for b = 1:3, for c = chans, x.(sprintf('old_b%d_%s', b, c{1})) = NaN; end, end
for c = {'SSC_A' 'SSC_H' 'GL1_A' 'GL1_H'}, x.(['old_c1_' c{1}]) = NaN; end
if isempty(old.T), return, end
r = find(strcmp(old.T.filename, fname), 1);
if isempty(r), return, end
v = old.T.Properties.VariableNames;
x.old_QC = double(old.T.QC_flag(r)); x.old_gen = old.gen;
for b = 1:3
    nm = sprintf('number%d', b);
    if ismember(nm, v), x.(sprintf('old_n%d', b)) = old.T.(nm)(r); end
    m = [];
    if ismember(sprintf('median%d', b), v), m = old.T.(sprintf('median%d', b)); end
    if istable(m)
        for c = chans
            h = [c{1} '-H'];
            if ismember(h, m.Properties.VariableNames) && m{r, h} > 0
                x.(sprintf('old_b%d_%s', b, c{1})) = log10(m{r, h});
            end
        end
    end
end
% 1 um centend2 modes used by batch_mergedSSC2volume
src = {'NoOD2centersA' 'SSC_A'; 'NoOD2centersH' 'SSC_H'; 'OD2centersA' 'GL1_A'; 'OD2centersH' 'GL1_H'};
for i = 1:size(src, 1)
    if ismember(src{i, 1}, v) && old.T.(src{i, 1})(r, 2) > 0
        x.(['old_c1_' src{i, 2}]) = log10(old.T.(src{i, 1})(r, 2));
    end
end
end

function old = load_old(bp)
old.T = []; old.gen = "";
for g = {'2026' '2021'}
    fn = fullfile(bp, ['beadstat_' g{1} '.mat']);
    if exist(fn, 'file')
        s = load(fn, ['beadstat_' g{1}]);
        old.T = s.(['beadstat_' g{1}]); old.gen = string(g{1});
        return
    end
end
end
