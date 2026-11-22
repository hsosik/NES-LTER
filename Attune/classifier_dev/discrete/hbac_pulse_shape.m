function S = hbac_pulse_shape(nper)
% function S = hbac_pulse_shape(nper)
% Does SYBR (BL1) pulse shape separate bacteria from noise and debris in the
% preserved SYBR runs? For nper hbac runs per cruise (the 12 tuning cruises, the
% score_hbac pick), per event BL1_AH = log10(BL1-A / BL1-H) relative to the run's
% reference (median of the core hand-gated bacteria: stained, SSC < -1.4, BL1 more
% than 0.3 above the trigger edge). Groups (hand time window only):
%   core       hand bacteria (Gated_Table 3/4), BL1 > edge + 0.3
%   edge       hand bacteria within 0.3 of the BL1 trigger edge
%   near_out   not hand-gated, SSC < -1.2 (bacteria size), within 0.3 of the edge
%   unst_big   not hand-gated, SSC >= -1.2 (larger particles)
%   auto_only  classifier bacteria (9) the hand gate leaves out
%   hand_only  hand bacteria the classifier leaves out
% Also SSC_AH, as a check. Prints quantiles and the share beyond +-0.2 per group.
% Saves output/discrete/hbac_pulse_shape.mat
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 1, nper = 6; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
cruises = {'20200201_EN649' '20200725_EN655' '20210203_EN661' '20210716_EN668' '20211108_AR61B' '20220729_EN687' ...
    '20230111_EN695' '20230807_EN706' '20240503_EN715' '20240906_EN720' '20250124_EN727' '20250814_AR92'};
grp = {'core' 'edge' 'near_out' 'unst_big' 'auto_only' 'hand_only'};
S = struct('runs', table, 'ev', table);
chk = table;
for ci = 1:numel(cruises)
    G = load(fullfile(root, cruises{ci}, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    ix = find(contains(G.fcslist, 'hbac') & ~cellfun(@isempty, G.class) & ~contains(G.fcslist, {'skip' 'test' 'badly' 'virus'}));
    ix = ix(unique(round(linspace(1, numel(ix), min(nper, numel(ix))))));
    cr = extractAfter(cruises{ci}, '_');
    for i = ix'
        fn = G.fcslist{i};
        try, [dat, hdr] = fca_readfcs(fullfile(root, cruises{ci}, 'preserved', 'FCS', fn)); catch, continue, end
        pn = {hdr.par.name};
        jA = find(strcmp(pn, 'BL1-A'), 1); jH = find(strcmp(pn, 'BL1-H'), 1);
        sA = find(strcmp(pn, 'SSC-A'), 1); sH = find(strcmp(pn, 'SSC-H'), 1);
        trig = ''; if isfield(hdr, 'trigger'), trig = hdr.trigger; end
        chk = [chk; table(string(cr), string(fn), ~isempty(jA) & ~isempty(jH), ~isempty(sA) & ~isempty(sH), 'VariableNames', {'cruise' 'file' 'has_BL1_AH' 'has_SSC_AH'})]; %#ok<AGROW>
        if isempty(jA) || isempty(jH), continue, end
        cfg = discrete_config(); cfg.cruise = cr; cfg.gl1_is_od2_scatter = true;
        ref = discrete_bead_reference(lab, datenum([hdr.date ' ' hdr.starttime]), cfg);
        if isempty(ref), continue, end
        [auto, info] = assign_class_hbac(dat, hdr, ref, cfg);
        F = info.features; e = info.edges.BL1;
        AH = log10(max(dat(:, jA), 1) ./ max(dat(:, jH), 1));
        sAH = nan(size(AH)); if ~isempty(sA) && ~isempty(sH), sAH = log10(max(dat(:, sA), 1) ./ max(dat(:, sH), 1)); end
        hd = double(G.class{i}(:)); hb = ismember(hd, [3 4]);
        late = true(size(hd)); if any(hd), late = dat(:, 1) >= min(dat(hd ~= 0, 1)); end
        ab = ismember(auto(:), [7 9]);
        core = late & hb & F.SSC < -1.4 & F.BL1 > e + 0.3 & (F.BL1 - F.SSC) > 0;
        if sum(core) < 200, continue, end
        r0 = median(AH(core)); s0 = median(sAH(core), 'omitnan');
        g = strings(size(hd));
        g(late & hb & F.BL1 > e + 0.3) = "core";
        g(late & hb & F.BL1 <= e + 0.3) = "edge";
        g(late & hd == 0 & F.SSC < -1.2 & F.BL1 <= e + 0.3) = "near_out";
        g(late & hd == 0 & F.SSC >= -1.2) = "unst_big";
        ao = late & ab & ~hb; ho = late & hb & ~ab;
        % per-run summary
        row = table(string(cr), string(fn), r0, sum(core), 'VariableNames', {'cruise' 'file' 'ref_AH' 'n_core'});
        for k = 1:numel(grp)
            switch grp{k}
                case 'auto_only', m = ao;
                case 'hand_only', m = ho;
                otherwise, m = g == grp{k};
            end
            d = AH(m) - r0;
            row.(['n_' grp{k}]) = sum(m);
            row.(['lo_' grp{k}]) = mean(d < -0.2); row.(['hi_' grp{k}]) = mean(d > 0.2);
            row.(['med_' grp{k}]) = median(d);
        end
        S.runs = [S.runs; row];
        % event subsample for pooled quantiles
        for k = 1:numel(grp)
            switch grp{k}
                case 'auto_only', m = find(ao);
                case 'hand_only', m = find(ho);
                otherwise, m = find(g == grp{k});
            end
            m = m(randperm(numel(m), min(2000, numel(m))));
            S.ev = [S.ev; table(repmat(string(grp{k}), numel(m), 1), AH(m) - r0, sAH(m) - s0, F.SSC(m), F.BL1(m) - e, ...
                repmat(string(cr), numel(m), 1), 'VariableNames', {'group' 'dAH' 'dSSC_AH' 'SSC' 'BL1_above_edge' 'cruise'})]; %#ok<AGROW>
        end
        fprintf('%s %s ref %.3f core %d | lo/hi core %.3f/%.3f edge %.3f/%.3f near_out %.3f/%.3f big %.3f/%.3f\n', cr, fn, r0, sum(core), ...
            row.lo_core, row.hi_core, row.lo_edge, row.hi_edge, row.lo_near_out, row.hi_near_out, row.lo_unst_big, row.hi_unst_big);
    end
end
S.channels = chk;
save(fullfile(od, 'hbac_pulse_shape.mat'), 'S')
fprintf('\nBL1-A/H present in %d of %d SYBR runs read; SSC-A/H in %d (cruises lacking BL1-A: %s)\n', sum(chk.has_BL1_AH), height(chk), sum(chk.has_SSC_AH), ...
    strjoin(unique(chk.cruise(~chk.has_BL1_AH))', ' '));
q = @(v) sprintf('%6.2f', prctile(v, [5 25 50 75 95]));
fprintf('\npooled dAH (BL1 A/H minus run reference) 5/25/50/75/95%%, share < -0.2 and > 0.2:\n');
for k = 1:numel(grp)
    d = S.ev.dAH(S.ev.group == grp{k});
    if isempty(d), continue, end
    fprintf('%-10s n %6d | %s | lo %.3f hi %.3f\n', grp{k}, numel(d), q(d), mean(d < -0.2), mean(d > 0.2));
end
fprintf('\nper-run medians of the shares (lo = narrower than ref - 0.2, hi = wider than ref + 0.2):\n');
for k = 1:numel(grp)
    fprintf('%-10s lo %.3f hi %.3f (runs %d)\n', grp{k}, median(S.runs.(['lo_' grp{k}]), 'omitnan'), median(S.runs.(['hi_' grp{k}]), 'omitnan'), height(S.runs));
end
end
