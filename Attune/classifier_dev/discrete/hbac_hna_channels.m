function [S, P] = hbac_hna_channels(cruisedirs, nper, seed)
% function [S, P] = hbac_hna_channels(cruisedirs, nper, seed)
% Do channels other than SYBR (BL1) and SSC help split HNA from LNA bacteria?
% For nper random (seeded) SYBR runs per cruise, the bacteria (class 9) are
% split with hbac_hna_split (BL1 tilted by SSC). Each other feature X (FSC,
% GL1, BL2, BL3, GL2, GL3, GL4, BL1 and SSC pulse shape log10(A/H), FSC - SSC)
% is reduced to its residual after a linear fit on BL1 and SSC over the
% bacteria, so only what BL1 and SSC do not already explain is left
% (SYBR spills into BL2, BL3 and GL2 along fixed lines). Per run and feature:
%   floor  share of bacteria at the channel floor (feature unusable if high)
%   d      (median residual HNA - median LNA) / pooled robust SD
%   dip    best dip of y + w * residual (w = -1..1), against the tilted y alone
% S - per run: cut, method, dip; P - per run and feature
% Saves output/discrete/hbac_hna_channels.mat
%
% Nov 2026, discrete-sample classifier development (HNA/LNA split)

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 2 || isempty(nper), nper = 4; end
if nargin < 3, seed = 202; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = discrete_config(); p = cfg.hbac.hna;
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
feats = ["FSC" "GL1" "BL2" "BL3" "GL2" "GL3" "GL4" "BL1_AH" "SSC_AH" "FSCxSSC"];
W = [-1 -0.5 -0.25 0.25 0.5 1];
S = table; P = table;
for ci = 1:numel(cruisedirs)
    cd_ = cruisedirs{ci}; cruise = extractAfter(cd_, '_'); cfg.cruise = cruise;
    G = load(fullfile(root, cd_, 'preserved', 'outputs', 'Gated_Table.mat'));
    fn = fieldnames(G); G = G.(fn{find(cellfun(@(q) istable(G.(q)), fn), 1)});
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    ix = find(contains(G.fcslist, 'hbac') & ~contains(G.fcslist, {'skip' 'test' 'badly' 'virus'}));
    rng(seed + ci); ix = ix(sort(randperm(numel(ix), min(nper, numel(ix)))));
    for i = ix'
        try, [dat, hdr] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', G.fcslist{i})); catch, continue, end
        ref = discrete_bead_reference(lab, datenum([hdr.date ' ' hdr.starttime]), cfg);
        if isempty(ref), continue, end
        [cl, info] = assign_class_hbac(dat, hdr, ref, cfg);
        F = info.features; b = cl == 9;
        if sum(b) < 2000, continue, end
        pn = {hdr.par.name};
        F.BL1_AH = pulse(dat, pn, 'BL1'); F.SSC_AH = pulse(dat, pn, 'SSC');
        F.FSCxSSC = F.FSC - F.SSC;
        bt = info.bact_type; s = info.hna;
        y = F.BL1 - p.slope * (F.SSC + 2);
        S = [S; table(string(cruise), string(G.fcslist{i}), s.cut, string(s.method), s.dip, sum(b), ...
            'VariableNames', {'cruise' 'file' 'cut' 'method' 'dip' 'n'})]; %#ok<AGROW>
        r = height(S);
        X0 = [ones(sum(b), 1) F.BL1(b) F.SSC(b)];
        for f = feats
            x = F.(f)(b);
            fl = NaN;
            if ismember(f, ["FSC" "GL1" "BL2" "BL3" "GL2" "GL3" "GL4"]) && any(isfinite(x))
                fl = mean(x <= min(x) + 1e-9); %clamped to the channel floor (attune_features)
            end
            ok = isfinite(x);
            if mean(ok) < 0.5 || numel(unique(x(ok))) < 20
                P = [P; table(r, f, fl, NaN, NaN, NaN, 'VariableNames', {'run' 'feat' 'floor' 'd' 'dip_best' 'w_best'})]; %#ok<AGROW>
                continue
            end
            c = X0(ok, :) \ x(ok); res = nan(size(x)); res(ok) = x(ok) - X0(ok, :) * c;
            sd = 1.4826 * mad(res(ok), 1); if sd == 0, sd = std(res(ok)); end
            res = (res - median(res, 'omitnan')) / sd * 0.1; %0.1 log10 per robust SD, on the y scale
            lo = bt(b) == 1 & ok; hi = bt(b) == 2 & ok;
            pd = sqrt((mad(res(lo), 1)^2 + mad(res(hi), 1)^2) / 2) * 1.4826;
            d = (median(res(hi)) - median(res(lo))) / pd;
            best = s.dip; wb = 0;
            for w = W
                z = y(b) + w * res; z(~ok) = NaN;
                [~, s2] = hbac_hna_split(z, -2 * ones(size(z)), true(size(z)), setfield(p, 'slope', 0)); %#ok<SFLD>
                if s2.dip < best, best = s2.dip; wb = w; end
            end
            P = [P; table(r, f, fl, d, best, wb, 'VariableNames', {'run' 'feat' 'floor' 'd' 'dip_best' 'w_best'})]; %#ok<AGROW>
        end
        fprintf('%s %s %s dip %.2f\n', cruise, G.fcslist{i}, s.method, s.dip)
    end
end
P.dip0 = S.dip(P.run); P.method = S.method(P.run);
save(fullfile(od, 'hbac_hna_channels.mat'), 'S', 'P')
fprintf('\n%d runs (%d fixed); per feature: median floor share, median |d| (HNA-LNA residual, robust SDs), dip gain\n', height(S), sum(S.method == "fixed"))
for f = feats
    q = P(P.feat == f, :);
    g = q.dip0 - q.dip_best;
    fx = q.method == "fixed";
    fprintf('%-8s floor %.2f  |d| %.2f (d median %+.2f)  dip gain median %.3f, >0.05 in %2.0f%% of runs; fixed runs rescued (dip < %.1f) %d of %d\n', ...
        f, median(q.floor, 'omitnan'), median(abs(q.d), 'omitnan'), median(q.d, 'omitnan'), median(g, 'omitnan'), 100 * mean(g > 0.05), ...
        p.max_dip, sum(fx & q.dip_best < p.max_dip), sum(fx))
end
end

function ah = pulse(dat, pn, ch)
ja = find(strcmp(pn, [ch '-A']), 1); jh = find(strcmp(pn, [ch '-H']), 1);
ah = nan(size(dat, 1), 1);
if ~isempty(ja) && ~isempty(jh), ah = log10(max(dat(:, ja), 1) ./ max(dat(:, jh), 1)); end
end

