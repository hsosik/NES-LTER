function S = syn_ssc_check(cruises, nfiles, outfile)
% function S = syn_ssc_check(cruises, nfiles, outfile)
% Bead-normalized Syn SSC across cruises, to check bead references (e.g. the
% PT-to-FCB conversion) against cruises in the same season. For nfiles sample
% files spread through each cruise:
%   - Syn: events with PE (physical GL2, 575/36) > chl (BL3) + 0.3 log10 and
%     PE > 10^3; densest point on SSC-H x PE, box +-0.25 (SSC) / +-0.3 (PE)
%   - Syn SSC = KDE peak (bead_peak) of log10 SSC-H in the box
%   - normalized = Syn SSC - (bead_log.SSC + k log10(hv_sample / bead_hv.SSC)),
%     from bead_reference (1 um bead, load_beadstat), as in attune_features
% A simple gate, independent of the classifier and of the filter layout.
%
% Nov 2026, bead stats redesign (design/test)

dev = fileparts(fileparts(mfilename('fullpath')));
addpath(dev, fullfile(dev, '..'))
warning('off', 'bead_reference:noQC')
if nargin < 2 || isempty(nfiles), nfiles = 15; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
rows = {};
for ic = 1:numel(cruises)
    cdir = cruises{ic}; cruise = extractAfter(cdir, '_');
    conf = attune_config(cruise);
    try
        bs = load_beadstat(fullfile(root, cdir, 'bead_calibrated'));
    catch
        fprintf('%s: no beadstat.mat\n', cruise); continue
    end
    F = dir(fullfile(root, cdir, 'FCS', '*.fcs'));
    F = F(~contains({F.name}, {'bead' 'Bead' 'test' 'Test' 'Rinse' 'rinse' 'discrete' 'Discrete' 'blank' 'Blank'}) & [F.bytes] > 5e5);
    [~, o] = sort([F.datenum]); F = F(o);
    pick = unique(round(linspace(1, numel(F), min(nfiles, numel(F)))));
    for k = pick
        try
            [d, h] = fca_readfcs(fullfile(F(k).folder, F(k).name));
        catch
            continue
        end
        names = {h.par.name};
        col = @(ch) find(strcmp(names, conf.chan.(ch)), 1);
        if isempty(col('GL2')) || isempty(col('BL3')), continue, end
        L = @(j) log10(max(d(:, j), 1));
        ss = L(col('SSC')); pe = L(col('GL2')); chl = L(col('BL3'));
        cand = pe > chl + 0.3 & pe > 3 & ss > 1; % PE > 10^3 clears the PE floor streak (~10^2.5)
        if sum(cand) < 100, continue, end
        e1 = 1:0.05:6.1; e2 = 3:0.05:6.1;
        H = conv2(histcounts2(ss(cand), pe(cand), e1, e2), ones(3) / 9, 'same');
        [~, im] = max(H(:)); [i1, i2] = ind2sub(size(H), im);
        g = cand & abs(ss - e1(i1) - 0.025) < 0.25 & abs(pe - e2(i2) - 0.025) < 0.3;
        if sum(g) < 50, continue, end
        t = datenum(datetime([h.date ' ' h.starttime]));
        ref = bead_reference(bs, t, cfg);
        hvs = h.par(col('SSC')).hv;
        sscpk = bead_peak(ss(g));
        norm = sscpk - (ref.bead_log.SSC + cfg.k.SSC * log10(hvs / ref.bead_hv.SSC));
        rows{end+1, 1} = struct('cruise', string(cruise), 'file', string(F(k).name), 'time', t, ...
            'month', month(datetime(t, 'ConvertFrom', 'datenum')), 'nsyn', sum(g), 'syn_ssc', sscpk, ...
            'syn_pe', bead_peak(pe(g)), 'ssc_hv', hvs, 'bead_ssc', ref.bead_log.SSC, 'bead_hv', ref.bead_hv.SSC, ...
            'norm_ssc', norm, 'source', string(beadsource(fullfile(root, cdir, 'bead_calibrated')))); %#ok<AGROW>
    end
    fprintf('%s done\n', cruise)
end
S = struct2table([rows{:}]');
if nargin > 2 && ~isempty(outfile), save(outfile, 'S'), end
end

function s = beadsource(bp)
m = load(fullfile(bp, 'beadstat.mat'), 'meta');
s = "FCB";
if isfield(m.meta, 'source'), s = m.meta.source; end
end
