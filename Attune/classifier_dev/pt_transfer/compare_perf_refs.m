function C = compare_perf_refs(cruises, nfiles, outfile)
% function C = compare_perf_refs(cruises, nfiles, outfile)
% Compare bead references from the performance-test transfer
% (perf_beadstats, output/pt_transfer/<cruise>_beadstat.mat) with the
% cruise's current beadstat.mat (PT conversion, where it exists), at nfiles
% sample files spread through each cruise:
%   - per channel: reference (log10, at the sample file's voltage) new - old
%   - Syn SSC normalized with each reference, using the same simple Syn gate
%     as bead_design/syn_ssc_check (PE > chl + 0.3 and PE > 10^3; densest
%     point on SSC x PE; KDE peak of SSC in a box around it)
%
% Nov 2026, PT transfer for early cruises

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
warning('off', 'bead_reference:noQC')
if nargin < 1 || isempty(cruises)
    cruises = {'20180131_EN608' '20180404_AR28B' '20180414_AR29' '20180720_EN617' '20181020_AR31A' '20190201_EN627'};
end
if nargin < 2 || isempty(nfiles), nfiles = 15; end
root = ['\' '\sosiknas1\lab_data\Attune\cruise_data\'];
cfg = classifier_config();
chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL2' 'GL3' 'GL4'};
rows = {};
for ic = 1:numel(cruises)
    cdir = cruises{ic}; cruise = extractAfter(cdir, '_');
    conf = attune_config(cruise);
    s = load(fullfile(dev, 'output', 'pt_transfer', [cruise '_beadstat.mat']), 'beadstat');
    bnew = s.beadstat;
    bold = [];
    try, bold = load_beadstat(fullfile(root, cdir, 'bead_calibrated')); catch, end
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
        t = datenum(datetime([h.date ' ' h.starttime]));
        r = struct('cruise', string(cruise), 'file', string(F(k).name), 'time', t, ...
            'month', month(datetime(t, 'ConvertFrom', 'datenum')));
        rn = bead_reference(bnew, t, cfg);
        if ~isempty(bold), ro = bead_reference(bold, t, cfg); end
        for c = chans
            ch = c{1}; j = col(ch);
            if isempty(j), r.(['d_' ch]) = NaN; continue, end
            hvs = h.par(j).hv;
            at = @(ref) ref.bead_log.(ch) + cfg.k.(ch) * log10(hvs / ref.bead_hv.(ch));
            r.(['new_' ch]) = at(rn);
            if isempty(bold), r.(['d_' ch]) = NaN; else, r.(['d_' ch]) = at(rn) - at(ro); end
        end
        % Syn SSC (as syn_ssc_check)
        r.syn_ssc = NaN; r.norm_new = NaN; r.norm_old = NaN; r.nsyn = 0;
        if ~isempty(col('GL2')) && ~isempty(col('BL3'))
            L = @(j) log10(max(d(:, j), 1));
            ss = L(col('SSC')); pe = L(col('GL2')); chl = L(col('BL3'));
            cand = pe > chl + 0.3 & pe > 3 & ss > 1;
            if sum(cand) >= 100
                e1 = 1:0.05:6.1; e2 = 3:0.05:6.1;
                H = conv2(histcounts2(ss(cand), pe(cand), e1, e2), ones(3) / 9, 'same');
                [~, im] = max(H(:)); [i1, i2] = ind2sub(size(H), im);
                g = cand & abs(ss - e1(i1) - 0.025) < 0.25 & abs(pe - e2(i2) - 0.025) < 0.3;
                if sum(g) >= 50
                    r.nsyn = sum(g);
                    r.syn_ssc = bead_peak(ss(g));
                    r.norm_new = r.syn_ssc - r.new_SSC;
                    r.norm_old = r.norm_new + r.d_SSC;
                end
            end
        end
        rows{end+1, 1} = r; %#ok<AGROW>
    end
    fprintf('%s done\n', cruise)
end
C = struct2table([rows{:}]');
if nargin > 2 && ~isempty(outfile), save(outfile, 'C'), end
end
