function R = check_discrete_refs()
% function R = check_discrete_refs()
% Checks before using lab bead runs and underway matches for preserved samples:
%  A. same-day coverage: how many preserved runs (Gated_Table, all cruises)
%     have a QC'd lab bead run on the same day, per bead population
%  B. lab vs at-sea bead positions: voltage-corrected log10 bead peaks (cfg.k)
%     for lab runs and every cruise bead table, over time; lab minus the median
%     of at-sea runs within 90 days
%  C. clock: preserved Syn concentration at surface bottles vs underway Syn
%     concentration (AttuneTable, existing classes) at time offsets of -8 to +8 h
% Saves output/discrete/check_refs.mat and figures there.
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
P = {'05' '1' '6'};

%% A. same-day coverage
d = dir(fullfile(root, '*_*', 'preserved', 'outputs', 'Gated_Table.mat'));
labday = floor(datenum(lab.time));
A = table;
for i = 1:numel(d)
    try
        G = load(fullfile(d(i).folder, d(i).name)).gated_table;
    catch
        continue
    end
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    if ~ismember('date_processed', G.Properties.VariableNames), continue, end
    dp = G.date_processed;
    if ~isdatetime(dp), dp = datetime(string(dp), 'InputFormat', 'dd-MMM-yyyy'); end
    day = floor(datenum(dp));
    cr = string(extractBetween(d(i).folder, root, filesep));
    row = table(cr, height(G), numel(unique(day(isfinite(day)))), 'VariableNames', {'cruise' 'nruns' 'ndays'});
    for p = P
        q = labday(lab.(['QC' p{1}]));
        row.(['sameday' p{1}]) = sum(ismember(day, q));
        row.(['within3d' p{1}]) = sum(arrayfun(@(x) any(abs(q - x) <= 3), day));
    end
    row.first_run = datetime(min(day), 'ConvertFrom', 'datenum');
    A = [A; row]; %#ok<AGROW>
end
disp(A)

%% B. lab vs at-sea bead positions
bd = dir(fullfile(root, '*_*', 'bead_calibrated', 'beadstat.mat'));
sea = table;
for i = 1:numel(bd)
    try
        b = load_beadstat(bd(i).folder);
    catch
        continue
    end
    b.source = repmat(string(extractBetween(bd(i).folder, root, filesep)), height(b), 1);
    sea = [sea; b(:, {'filename' 'time' 'hv' 'QC05' 'QC1' 'QC6' 'peak05' 'peak1' 'peak6' 'source'})]; %#ok<AGROW>
end
sea = sea(sea.time >= datetime(2019, 8, 16), :);
lab.source = repmat("lab", height(lab), 1);
B = table;
fig = figure('visible', 'off', 'position', [50 50 1500 800]);
for k = 1:numel(cfg.chans)
    ch = cfg.chans{k}; p = P{cfg.bead_pop.(ch)};
    [ts, vs] = corrected(sea, ch, p, cfg.k.(ch));
    [tl, vl] = corrected(lab, ch, p, cfg.k.(ch));
    dlab = nan(size(tl));
    for j = 1:numel(tl)
        w = abs(datenum(ts) - datenum(tl(j))) <= 90;
        if sum(w) >= 3, dlab(j) = vl(j) - median(vs(w)); end
    end
    B = [B; table(string(ch), string(p), numel(vl), median(dlab, 'omitnan'), prctile(dlab, 10), prctile(dlab, 90), sum(isfinite(dlab)), ...
        'VariableNames', {'chan' 'pop' 'nlab' 'lab_minus_sea_median' 'p10' 'p90' 'npaired'})]; %#ok<AGROW>
    subplot(3, 3, k)
    plot(ts, vs, '.', tl, vl, 'o', 'markersize', 4)
    title(sprintf('%s, %s um bead (log10 at 300 V, k %.1f)', ch, strrep(p, '05', '0.5'), cfg.k.(ch)), 'fontsize', 8), grid on
    if k == 1, legend('at sea', 'lab', 'location', 'best'), end
end
exportgraphics(fig, fullfile(od, 'lab_vs_sea_beads.png'), 'resolution', 80); close(fig)
disp(B)

%% C. clock offset
cr = {'20240906_EN720' '20230111_EN695' '20230807_EN706' '20240503_EN715'};
offs = -8:0.25:8;
C = table;
fig = figure('visible', 'off', 'position', [50 50 1200 350]);
for i = 1:numel(cr)
    G = load(fullfile(root, cr{i}, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    ts = todatenum(G.date_sampled);
    use = G.depth_m <= 7 & isfinite(ts) & isfinite(G.syn_conc) & G.syn_conc > 0 & contains(G.fcslist, {'phyto_PE' 'phyto_CHL'}) & ~contains(G.fcslist, 'pro', 'IgnoreCase', true);
    G = G(use, :); ts = ts(use);
    % PE run preferred per bottle
    [~, o] = sort(~contains(G.fcslist, 'phyto_PE')); G = G(o, :); ts = ts(o);
    [~, ia] = unique([G.cast G.niskin], 'rows', 'stable'); G = G(ia, :); ts = ts(ia);
    U = load(fullfile(root, cr{i}, 'bead_calibrated', 'AttuneTable.mat')).AttuneTable;
    U = U(U.QC_flag == 1 & ~contains(U.Filename, {'bead' 'Bead' 'Rinses' 'discrete'}), :);
    ut = datenum(U.StartDate); uc = U.Syn_count ./ U.VolAnalyzed_ml;
    m = nan(size(offs));
    for k = 1:numel(offs)
        r = nan(height(G), 1);
        for j = 1:height(G)
            w = abs(ut - (ts(j) + offs(k) / 24)) <= 15 / 1440;
            if any(w), r(j) = log10(G.syn_conc(j) / median(uc(w))); end
        end
        m(k) = median(abs(r - median(r, 'omitnan')), 'omitnan'); % spread of pres/live ratio
        if offs(k) == 0, r0 = r; end
    end
    [~, kb] = min(m);
    C = [C; table(string(cr{i}), height(G), m(offs == 0), offs(kb), m(kb), median(r0, 'omitnan'), ...
        'VariableNames', {'cruise' 'nbottles' 'spread_at_0' 'best_offset_h' 'spread_best' 'median_log10_pres_over_live_at_0'})]; %#ok<AGROW>
    subplot(1, numel(cr), i), plot(offs, m, '.-'), xline(0, 'k:'), xlabel('offset (h)'), ylabel('MAD log10(pres/live)')
    title(cr{i}, 'interpreter', 'none', 'fontsize', 8), grid on
end
exportgraphics(fig, fullfile(od, 'clock_offset.png'), 'resolution', 80); close(fig)
disp(C)
R.A = A; R.B = B; R.C = C;
save(fullfile(od, 'check_refs.mat'), 'R')
end

function [t, v] = corrected(b, ch, p, k)
q = b.(['QC' p]);
x = b.(['peak' p]).([ch '_H']); hv = b.hv.(ch);
ok = q & x > 0 & hv > 0;
t = b.time(ok); v = log10(x(ok)) - k * log10(hv(ok) / 300);
end

function t = todatenum(d)
if isdatetime(d), t = datenum(d); return, end
d = string(d); t = nan(size(d));
for i = 1:numel(d)
    if ismissing(d(i)) || strlength(d(i)) == 0, continue, end
    s = replace(extractBefore(d(i) + "+", "+"), "T", " ");
    try, t(i) = datenum(datetime(s, 'InputFormat', 'yyyy-MM-dd HH:mm:ss')); catch, end
end
end
