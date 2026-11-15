function G = gl1_lowsignal_check(cruises, nfiles, outfile)
% function G = gl1_lowsignal_check(cruises, nfiles, outfile)
% Is side scatter through the OD2 filter non-linear at low signal? On cruises
% with the current layout (Aug 2019 on), every particle is measured on SSC
% (488 nm, no OD2, linear for small cells) and on GL1 (532 nm side scatter
% through OD2). Both are normalized to the 1 um bead on the same channel
% (bead_reference). If GL1 is linear, GL1norm - SSCnorm is flat with signal
% level (up to a small 488 vs 532 nm difference); a rise at low raw GL1
% means small signals read too high under OD2 (as suspected for Syn on
% TN368 and the 2018 OD2-on-SSC cruises).
%
% Per sample file:
%   - Syn (same gate as bead_design/syn_ssc_check): KDE peaks of SSC and GL1,
%     normalized; raw GL1 and SSC peak
%   - all events with SSC and GL1 above 10 counts: median of GL1norm - SSCnorm
%     in bins of SSCnorm (the linear axis), with the median raw GL1 per bin
%
% cruises - cruise folder names (default: a spread of seasons, Aug 2019 on)
% G       - struct: syn (table, one row per file), bins (table, file x bin)
%
% Nov 2026, PT transfer / OD2 low-signal check

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
warning('off', 'bead_reference:noQC')
root = ['\' '\sosiknas1\lab_data\Attune\cruise_data\'];
if nargin < 1 || isempty(cruises)
    cruises = {'20190820_EN644' '20191005_AR39B' '20200201_EN649' '20200725_EN655' '20201013_EN657' ...
        '20210203_EN661' '20210716_EN668' '20220216_AT46' '20220806_EN688' '20230201_EN695' ...
        '20230803_EN706' '20240205_EN712' '20240830_EN720' '20250124_EN727' '20250819_AR92' ...
        '20261019_AR104'}; % AR104: new settings, SSC 240 V / GL1 180 V (others 220/220)
end
if nargin < 2 || isempty(nfiles), nfiles = 12; end
D = dir(root); D = {D([D.isdir]).name};
cfg = classifier_config();
edges = -3:0.1:1; mid = edges(1:end-1) + 0.05;
syn = {}; bins = {};
for ic = 1:numel(cruises)
    cdir = cruises{ic};
    if ~ismember(cdir, D) % tolerate a wrong date prefix
        m = D(endsWith(D, ['_' extractAfter(cdir, '_')]));
        if isempty(m), fprintf('%s: no folder\n', cdir); continue, end
        cdir = m{1};
    end
    cruise = extractAfter(cdir, '_');
    try
        bs = load_beadstat(fullfile(root, cdir, 'bead_calibrated'));
    catch
        fprintf('%s: no beadstat\n', cruise); continue
    end
    F = dir(fullfile(root, cdir, 'FCS', '*.fcs'));
    F = F(~contains({F.name}, {'bead' 'Bead' 'test' 'Test' 'Rinse' 'rinse' 'discrete' 'Discrete' 'blank' 'Blank'}) & [F.bytes] > 5e5);
    [~, o] = sort({F.name}); F = F(o);
    pick = unique(round(linspace(1, numel(F), min(nfiles, numel(F)))));
    for k = pick
        try
            [d, h] = fca_readfcs(fullfile(F(k).folder, F(k).name));
        catch
            continue
        end
        t = datenum(datetime([h.date ' ' h.starttime]));
        conf = attune_config(cruise, t);
        if conf.layout ~= "current", continue, end
        names = {h.par.name};
        col = @(ch) find(strcmp(names, conf.chan.(ch)), 1);
        if isempty(col('GL1')) || isempty(col('GL2')) || isempty(col('BL3')), continue, end
        ref = bead_reference(bs, t, cfg);
        if ~isfinite(ref.bead_log.GL1) || ~isfinite(ref.bead_log.SSC), continue, end
        hvS = h.par(col('SSC')).hv; hvG = h.par(col('GL1')).hv;
        bS = ref.bead_log.SSC + cfg.k.SSC * log10(hvS / ref.bead_hv.SSC);
        bG = ref.bead_log.GL1 + cfg.k.GL1 * log10(hvG / ref.bead_hv.GL1);
        L = @(j) log10(max(d(:, j), 1));
        ss = L(col('SSC')); g1 = L(col('GL1')); pe = L(col('GL2')); chl = L(col('BL3'));
        % all events: GL1norm - SSCnorm vs SSCnorm
        ok = ss > 1 & g1 > 1;
        xs = ss(ok) - bS; dg = (g1(ok) - bG) - xs; rg = g1(ok);
        [~, ~, bi] = histcounts(xs, edges);
        for b = 1:numel(mid)
            q = bi == b;
            if sum(q) < 30, continue, end
            bins{end+1, 1} = struct('cruise', string(cruise), 'file', string(F(k).name), 'ssc_norm', mid(b), ...
                'dGL1', median(dg(q)), 'raw_gl1', median(rg(q)), 'n', sum(q), 'hv_gl1', hvG); %#ok<AGROW>
        end
        % Syn
        cand = pe > chl + 0.3 & pe > 3 & ss > 1;
        if sum(cand) < 100, continue, end
        e1 = 1:0.05:6.1; e2 = 3:0.05:6.1;
        H = conv2(histcounts2(ss(cand), pe(cand), e1, e2), ones(3) / 9, 'same');
        [~, im] = max(H(:)); [i1, i2] = ind2sub(size(H), im);
        gg = cand & abs(ss - e1(i1) - 0.025) < 0.25 & abs(pe - e2(i2) - 0.025) < 0.3;
        if sum(gg) < 50, continue, end
        syn{end+1, 1} = struct('cruise', string(cruise), 'file', string(F(k).name), 'time', t, ...
            'month', month(datetime(t, 'ConvertFrom', 'datenum')), 'nsyn', sum(gg), ...
            'raw_ssc', bead_peak(ss(gg)), 'raw_gl1', bead_peak(g1(gg)), ...
            'ssc_norm', bead_peak(ss(gg)) - bS, 'gl1_norm', bead_peak(g1(gg)) - bG, ...
            'frac_gl1_floor', mean(g1(gg) <= 1), 'hv_ssc', hvS, 'hv_gl1', hvG, ...
            'bead_gl1_raw', bG, 'bead_ssc_raw', bS); %#ok<AGROW>
    end
    fprintf('%s done\n', cruise)
end
G.syn = struct2table([syn{:}]');
G.syn.d = G.syn.gl1_norm - G.syn.ssc_norm;
G.bins = struct2table([bins{:}]');
if nargin > 2 && ~isempty(outfile), save(outfile, 'G'), end
end
