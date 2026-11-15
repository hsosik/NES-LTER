function out = generate_pt_beadstats(cruises, writefiles)
% function out = generate_pt_beadstats(cruises, writefiles)
% Bead statistics for cruises with only Attune performance-tracking (PT) bead
% runs (no FCB bead mix): FCB-equivalent peaks are derived from the PT bead
% intensity levels and written to <cruise>\bead_calibrated\beadstat.mat in
% the same format as process_beads_PT_adjust_3, so bead_reference and the
% classifier read them unchanged (meta.source = "PT conversion").
%
% Method
%   - PT beads: one size, three intensity levels per fluorescence channel
%     (find_pt_beads). The ratio of two beads on the same detector does not
%     depend on detector gain or laser power, so an FCB - PT offset measured
%     once applies at any time with the same optical filters.
%   - Calibration: FCB bead runs of 20 Mar - 12 Apr 2019 (measure_beads)
%     against PT runs of 28 Feb - 5 Mar 2019 at identical voltages
%     (FSC 160, SSC 200, BL1 280, BL2 320, BL3 320, GL1 340, GL2 400,
%     GL3 320, GL4 400; pre-Aug 2019 filter layout, no OD2). For each
%     physical channel, FCB population and PT level: offset = median FCB
%     peak - median PT level (log10). FSC and SSC use the single PT peak.
%   - Conversion of a PT run: FCB peak = PT level + offset, using the PT
%     level nearest the FCB peak that the run has. The result refers to the
%     PT run's own voltages (stored in hv), so bead_reference scales it to
%     each sample as usual.
%   - PT runs at the high "bead settings" (e.g. AR31A) are first brought to
%     the July 2018 seawater settings with the bead -> seawater offset of the
%     same-day PT pairs in July 2018 (a pure voltage effect per channel).
%
% cruises    - cellstr of cruise folder names (default: the PT-only cruises
%              with PT files: EN617, AR31A, EN627)
% writefiles - default true
% out        - struct: calibration offsets, bead->seawater offsets, per-cruise tables
%
% Nov 2026, bead stats redesign

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
warning('off', 'stats:robustcov:SingularCovariance')
if nargin < 1 || isempty(cruises), cruises = {'20180720_EN617' '20181020_AR31A' '20190201_EN627'}; end
if nargin < 2 || isempty(writefiles), writefiles = true; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
base = fullfile(root, 'beads');
chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
fl = chans(3:end);
conf = attune_config('EN627'); conf.OD2 = 'None'; % pre-Aug 2019 layout; OD2 does not change bead ratios
cfg = classifier_config();

% ---- PT inventory (both settings folders) ----
I = pt_inventory(base, conf);
canon = pt_canonical();

% ---- calibration: FCB (20 Mar - 12 Apr 2019) vs PT (28 Feb - 5 Mar 2019), same voltages ----
F = dir(fullfile(base, 'FCB_bead_mix_experiment_settings', '**', 'FCB bead check_Group_day0_*.fcs'));
fdat = {}; fhdr = {}; fnames = {};
for i = 1:numel(F)
    [d, h] = fca_readfcs(fullfile(F(i).folder, F(i).name));
    if h.par(strcmp({h.par.name}, 'SSC-H')).hv ~= 200, continue, end
    fdat{end+1, 1} = single(d); fhdr{end+1, 1} = h; fnames{end+1, 1} = F(i).name; %#ok<AGROW>
end
B = measure_beads(fdat, fhdr, fnames, attune_config('AR34B'));
calPT = I(I.time >= datetime(2019, 2, 28) & I.time <= datetime(2019, 3, 6) & I.SSChv == 200, :);
P = {'05' '1' '6'};
off = struct();
for c = chans
    for ip = 1:3
        v = log10(B.(['peak' P{ip}]).([c{1} '_H'])); v = median(v(B.(['QC' P{ip}])), 'omitnan');
        if ismember(c{1}, {'FSC' 'SSC'})
            o = v - median(calPT.(c{1}), 'omitnan'); o = [o NaN NaN]; %#ok<AGROW>
        else
            lv = vertcat(calPT.levels.(c{1})); o = v - median(lv, 1, 'omitnan');
        end
        off.(c{1})(ip, :) = o; % rows: FCB population, cols: PT level
    end
end
out.calibration = off;
out.calibration_runs = struct('FCB', {fnames}, 'PT', {calPT.file});
out.fcb_cal_QC = [sum(B.QC05) sum(B.QC1) sum(B.QC6) height(B)];

% ---- bead settings -> July 2018 seawater settings (same-day PT pairs) ----
sw = I(I.set == "seawater" & I.time < datetime(2018, 8, 1), :);
bd = I(I.set == "bead" & I.time < datetime(2018, 8, 1), :);
dd = struct();
for c = chans, dd.(c{1}) = []; end
for i = 1:height(bd)
    [dt, j] = min(abs(sw.time - bd.time(i)));
    if dt > hours(2), continue, end
    for c = chans
        if ismember(c{1}, {'FSC' 'SSC'})
            dd.(c{1})(end+1, :) = [sw.(c{1})(j) - bd.(c{1})(i) NaN NaN];
        else
            dd.(c{1})(end+1, :) = sw.levels(j).(c{1}) - bd.levels(i).(c{1});
        end
    end
end
b2s = struct(); b2s_sd = struct();
for c = chans, b2s.(c{1}) = median(dd.(c{1}), 1, 'omitnan'); b2s_sd.(c{1}) = std(dd.(c{1}), 0, 1, 'omitnan'); end
out.bead_to_seawater = b2s; out.bead_to_seawater_sd = b2s_sd; out.n_pairs = size(dd.FSC, 1);
swhv = sw.hv(1);

% ---- convert each cruise ----
for ic = 1:numel(cruises)
    cdir = cruises{ic}; cruise = extractAfter(cdir, '_');
    Fc = dir(fullfile(root, cdir, 'FCS', '*.fcs'));
    t0 = datetime(min([Fc.datenum]), 'ConvertFrom', 'datenum') - days(10);
    t1 = datetime(max([Fc.datenum]), 'ConvertFrom', 'datenum') + days(10);
    R = I(I.time >= t0 & I.time <= t1, :);
    if any(R.set == "seawater"), R = R(R.set == "seawater", :); end
    beadstat = table;
    for r = 1:height(R)
        x = table;
        x.filename = R.file(r);
        x.time = R.time(r);
        lev = R.levels(r); fsc = R.FSC(r); ssc = R.SSC(r); hv = R.hv(r);
        if R.set(r) == "bead" % bring to the July 2018 seawater settings
            for c = fl, lev.(c{1}) = lev.(c{1}) + b2s.(c{1}); end
            fsc = fsc + b2s.FSC(1); ssc = ssc + b2s.SSC(1);
            hv = swhv;
        end
        x.hv = struct2table(hv);
        for ip = 1:3
            pk = nan(1, numel(chans));
            for k = 1:numel(chans)
                c = chans{k}; o = off.(c)(ip, :);
                if strcmp(c, 'FSC'), pk(k) = fsc + o(1);
                elseif strcmp(c, 'SSC'), pk(k) = ssc + o(1);
                else
                    have = isfinite(lev.(c)) & isfinite(o);
                    if any(have)
                        [~, j] = min(abs(o(have))); % PT level nearest the FCB peak
                        hl = find(have); j = hl(j);
                        pk(k) = lev.(c)(j) + o(j);
                    end
                end
            end
            names = [strcat(chans, '_H') {'SSC_A' 'GL1_A'}];
            x.(['n' P{ip}]) = R.n(r);
            why = strings(0);
            if R.n(r) < 200, why(end+1) = "PT gate < 200 events"; end %#ok<AGROW>
            if ~(R.SSCrsd(r) <= 0.08), why(end+1) = "PT scatter broad"; end %#ok<AGROW>
            ok = isempty(why);
            x.(['QC' P{ip}]) = ok;
            x.(['why' P{ip}]) = strjoin(why, ", ");
            x.(['peak' P{ip}]) = array2table(10 .^ [pk NaN NaN], 'VariableNames', names);
            x.(['med' P{ip}]) = array2table(nan(1, numel(names)), 'VariableNames', names);
            x.(['rsd' P{ip}]) = array2table(nan(1, numel(names)), 'VariableNames', names);
        end
        beadstat = [beadstat; x]; %#ok<AGROW>
    end
    meta = struct();
    meta.source = "PT conversion";
    meta.config = attune_config(cruise);
    meta.calibration = off;
    meta.calibration_runs = out.calibration_runs;
    meta.bead_to_seawater = b2s;
    meta.pt_runs = R(:, {'set' 'file' 'time' 'n'});
    meta.cruise = string(cruise);
    meta.processed = datetime('now');
    meta.code = "generate_pt_beadstats (Nov 2026)";
    meta.notes = ["FCB-equivalent peaks from PT bead levels: peak = PT level + (FCB - PT) offset from the Mar-Apr 2019 calibration";
        "medP and rsdP are NaN (no FCB beads measured); hv is the voltage the peaks refer to"];
    out.(matlab.lang.makeValidName(cruise)) = beadstat;
    fprintf('%s: %d PT runs (%s settings) -> %d with QC good\n', cruise, height(beadstat), strjoin(unique(R.set), '/'), sum(beadstat.QC1))
    if writefiles && height(beadstat) > 0
        bp = fullfile(root, cdir, 'bead_calibrated');
        if ~exist(bp, 'dir'), mkdir(bp), end
        save(fullfile(bp, 'beadstat.mat'), 'beadstat', 'meta')
    end
end
end

function I = pt_inventory(base, conf)
% every PT run in both settings folders, with levels assigned against the
% canonical PT signature scaled to the run's voltages
canon = pt_canonical(); cfg = classifier_config();
sets = {'PT_beads_seawater_settings' 'seawater'; 'PT_beads_bead_settings' 'bead'};
rows = {};
for s = 1:2
    F = dir(fullfile(base, sets{s, 1}, '**', '*.fcs'));
    for i = 1:numel(F)
        try
            [d, h] = fca_readfcs(fullfile(F(i).folder, F(i).name));
        catch
            continue
        end
        if size(d, 1) > 2e5, continue, end % long acquisitions, not bead checks
        p0 = find_pt_beads(d, h, conf);
        % expected levels at this run's voltages: canonical + k log10(hv / hv_canon);
        % first assignment with a wide window, then the run's own levels
        ref = struct();
        for c = fieldnames(canon.levels)'
            ref.(c{1}) = canon.levels.(c{1}) + cfg.k.(c{1}) * log10(p0.hv.(c{1}) / canon.hv.(c{1}));
        end
        pt = assign(p0, ref, 0.6);
        rows{end+1, 1} = struct('set', string(sets{s, 2}), 'file', string(F(i).name), ...
            'time', datetime([h.date ' ' h.starttime]), 'n', pt.n, 'FSC', pt.FSC, 'SSC', pt.SSC, ...
            'FSCrsd', pt.FSCrsd, 'SSCrsd', pt.SSCrsd, ...
            'SSChv', pt.hv.SSC, 'levels', pt.levels, 'hv', pt.hv); %#ok<AGROW>
    end
end
I = struct2table([rows{:}]');
I = sortrows(I, 'time');
end

function pt = assign(pt, ref, w)
% label the found peaks as levels L1-L3 by aligning the whole expected pattern:
% try every shift (found peak - expected level, |shift| < w), match each
% expected level + shift to the nearest unused peak within 0.15, keep the
% shift with the most matches (then the smallest residual)
for c = fieldnames(ref)'
    pl = pt.peaks.(c{1}); r = ref.(c{1}); best = nan(1, 3); bscore = [0 Inf];
    for i = 1:numel(pl)
        for k = find(isfinite(r))
            sh = pl(i) - r(k);
            if abs(sh) >= w, continue, end
            lv = nan(1, 3); used = false(size(pl)); res = 0;
            for kk = find(isfinite(r))
                d = abs(pl - (r(kk) + sh)); d(used) = Inf;
                [dm, j] = min(d);
                if dm < 0.15, lv(kk) = pl(j); used(j) = true; res = res + dm; end
            end
            sc = [sum(isfinite(lv)) res];
            if sc(1) > bscore(1) || (sc(1) == bscore(1) && sc(2) < bscore(2)), best = lv; bscore = sc; end
        end
    end
    % a single matched peak does not confirm the pattern (e.g. noise peaks
    % when a laser was off): need two levels where two or more are expected
    if sum(isfinite(best)) < min(2, sum(isfinite(r))), best = nan(1, 3); end
    pt.levels.(c{1}) = best;
end
end

function c = pt_canonical()
% PT levels at the July 2018 seawater settings (median of 8 runs, 10-27 Jul
% 2018), by physical channel (pre-Aug 2019 file parameters in brackets:
% GL2 = file GL1 575/36, GL3 = file GL2 620/15, GL4 = file GL3 680/40)
c.hv = struct('FSC', 160, 'SSC', 260, 'BL1', 280, 'BL2', 320, 'BL3', 300, 'GL2', 300, 'GL3', 400, 'GL4', 320);
c.levels = struct('BL1', [2.12 3.42 4.13], 'BL2', [3.05 4.45 5.27], 'BL3', [2.22 3.63 4.73], ...
    'GL2', [3.46 4.88 5.69], 'GL3', [3.95 5.40 NaN], 'GL4', [2.67 4.12 5.17]);
end

function y = ifelse(c, a, b)
if c, y = a; else, y = b; end
end
