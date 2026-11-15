function V = validate_perf_transfer(cruises, maxdt_h, R)
% function V = validate_perf_transfer(cruises, maxdt_h, R)
% Test whether performance tests (pt_perf_table) predict FCB bead peaks.
% For every FCB bead run with QC'd populations, take the nearest good
% performance test in the same optical state (layout and OD2 placement) and
% predict, per physical channel and bead population:
%   log10 FCB peak = PT perf peak + k log10(hv_run / hv_perf) + offset
% offset = one constant per channel, population and optical state (median
% over all pairs). The residuals measure how well the performance tests
% transfer bead calibration between times, i.e. how well they could stand in
% for missing bead runs on the early cruises.
%
% cruises - cruise folder names with FCB beads in the performance-test era
%           (default: AR34B, RB1904, TN368, EN644)
% maxdt_h - max hours between bead run and performance test (default 24)
% R       - pt_perf_table output (default: load ../output/pt_perf.mat)
% V       - one row per bead run x channel x population with obs, pred,
%           offset, resid, dt_h, cruise, state
%
% Nov 2026, PT transfer for early cruises

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
if nargin < 1 || isempty(cruises)
    cruises = {'20190414_AR34B' '20190512_RB1904' '20190705_TN368' '20190820_EN644'};
end
if nargin < 2 || isempty(maxdt_h), maxdt_h = 24; end
if nargin < 3 || isempty(R), s = load(fullfile(dev, 'output', 'pt_perf.mat'), 'R'); R = s.R; end
R = R(R.ok, :);
R.state = R.layout + "/" + R.od2;
cfg = classifier_config();
root = ['\' '\sosiknas1\lab_data\Attune\cruise_data\'];
chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
P = {'05' '1' '6'};
rows = {};
for ic = 1:numel(cruises)
    cdir = cruises{ic}; cruise = extractAfter(cdir, '_');
    [b, meta] = load_beadstat(fullfile(root, cdir, 'bead_calibrated')); %#ok<ASGLU>
    t = b.time; if isnumeric(t), t = datetime(t, 'ConvertFrom', 'datenum'); end
    for i = 1:height(b)
        conf = attune_config(cruise, datenum(t(i)));
        state = string(conf.layout) + "/" + string(conf.OD2);
        q = find(R.state == state);
        if isempty(q), continue, end
        [dt, j] = min(abs(hours(R.time(q) - t(i))));
        if dt > maxdt_h, continue, end
        j = q(j);
        for c = 1:numel(chans)
            ch = chans{c};
            hvr = b.hv.(ch)(i); hvp = R.(['hv_' ch])(j);
            if ~(hvr > 0 && hvp > 0), continue, end
            for p = 1:3
                if ~b.(['QC' P{p}])(i), continue, end
                v = b.(['peak' P{p}]).([ch '_H'])(i);
                if ~(v > 0), continue, end
                rows{end+1, 1} = struct('cruise', string(cruise), 'file', string(b.filename(i)), ...
                    'time', t(i), 'perf_time', R.time(j), 'dt_h', dt, 'state', state, ...
                    'chan', string(ch), 'pop', string(P{p}), 'hv_run', hvr, 'hv_perf', hvp, ...
                    'obs', log10(v), 'pred0', R.(['pk_' ch])(j) + cfg.k.(ch) * log10(hvr / hvp)); %#ok<AGROW>
            end
        end
    end
end
V = struct2table([rows{:}]);
V.d = V.obs - V.pred0;
[g, gc, gp, gs] = findgroups(V.chan, V.pop, V.state);
off = splitapply(@median, V.d, g);
V.offset = off(g);
V.resid = V.d - V.offset;
% summary: robust SD of residuals and spread of cruise medians, per group
S = table(gc, gp, gs, off, splitapply(@numel, V.d, g), ...
    splitapply(@(x) 1.4826 * mad(x, 1), V.resid, g), ...
    'VariableNames', {'chan', 'pop', 'state', 'offset', 'n', 'rsd'});
cm = groupsummary(V, {'chan', 'pop', 'state', 'cruise'}, 'median', 'resid');
for r = 1:height(S)
    m = cm.median_resid(cm.chan == S.chan(r) & cm.pop == S.pop(r) & cm.state == S.state(r));
    S.cruise_range(r) = max(m) - min(m);
    S.ncruise(r) = numel(m);
end
V = addprop(V, 'summary', 'table'); V.Properties.CustomProperties.summary = S;
end
