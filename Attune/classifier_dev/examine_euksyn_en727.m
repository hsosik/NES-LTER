function examine_euksyn_en727(tag, outdir)
% function examine_euksyn_en727(tag, outdir)
% Check #2: in EN727 the Syn-spectrum PE euks (class 5) have ~6x less chl than
% Euks of the same size. Are they a separate group (e.g. detritus/aggregates
% with Syn, between class 8 and Euk+Syn)? Pools the scored files of one cruise and
% compares class 5 with Euks (1), Syn (2) and Syn-on-particle (8): chl relative to
% the Euk chl-size line, PE in Syn equivalents, chl in Syn equivalents, FSC, FSC-W.
% Writes outdir/euksyn_<tag>.png.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
a = load(fullfile(outdir, ['agreement_' tag '.mat']), 'A'); A = a.A;
E = table;
for i = 1:height(A)
    cpath = fullfile(root, char(A.cruise(i))); fn = char(A.file(i));
    bs = load_beadstat(fullfile(cpath, 'bead_calibrated'));
    [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, A.time(i), cfg);
    [new, info] = assign_class_v2(d, h, ref, cfg, A.pro_measured(i));
    F = info.features; pn = {h.par.name};
    eu = new == 1; sy = new == 2;
    if sum(eu) < 50 || sum(sy) < 50, continue, end
    p = robustfit(F.size(eu), F.BL3(eu));
    FSCW = log10(max(d(:, strcmp(pn, 'FSC-W')), 1));
    k = find(ismember(new, [1 5 8]) | (sy & rand(height(F), 1) < 0.05));
    E = [E; table(repmat(i, numel(k), 1), new(k), F.size(k), F.BL3(k), F.BL3(k) - (p(1) + p(2) * F.size(k)), ...
        F.GL2(k) - median(F.GL2(sy)), F.BL3(k) - median(F.BL3(sy)), F.FSC(k), FSCW(k), F.GL2GL3(k) - median(F.GL2GL3(sy)), ...
        'VariableNames', {'filei' 'class' 'size' 'BL3' 'dchl' 'PEsyn' 'chlsyn' 'FSC' 'FSCW' 'dGL2GL3'})]; %#ok<AGROW>
end
fprintf('%s pooled: Euk %d, Syn-on-particle %d, Euk-with-Syn %d\n', tag, sum(E.class == 1), sum(E.class == 8), sum(E.class == 5))
nm = {1 'Euk'; 5 'Euk with Syn'; 8 'Syn on particle'; 2 'Syn (5% sample)'};
for j = 1:4
    g = E.class == nm{j, 1};
    fprintf('  %-16s size %s  chl vs Euk line %s  chl (Syn eq) %s  PE (Syn eq) %s  FSC %s  FSC-W %s\n', nm{j, 2}, ...
        mat2str(prctile(E.size(g), [10 50 90]), 2), mat2str(prctile(E.dchl(g), [10 50 90]), 2), mat2str(prctile(E.chlsyn(g), [10 50 90]), 2), ...
        mat2str(prctile(E.PEsyn(g), [10 50 90]), 2), mat2str(prctile(E.FSC(g), [10 50 90]), 2), mat2str(prctile(E.FSCW(g), [10 50 90]), 2))
end
% does class 5 split by chl relative to Euks?
g = E.class == 5;
if sum(g) > 30
    x = E.dchl(g);
    g1 = fitgmdist(x, 1); g2 = fitgmdist(x, 2, 'Replicates', 5);
    fprintf('  Euk-with-Syn chl offset: BIC 1 comp %.0f, 2 comp %.0f; 2-comp means %s weights %s\n', g1.BIC, g2.BIC, mat2str(sort(g2.mu)', 2), mat2str(g2.ComponentProportion, 2))
end
f = figure('Visible', 'off', 'Position', [20 20 1500 450]); tiledlayout(1, 3, 'TileSpacing', 'compact');
cc = {[.6 .6 .6] [1 .55 0] [.55 .55 0] [.6 .1 .1]}; ord = [1 2 5 8];
pv = {'size' 'BL3' 'log10 side scatter rel. 1 um bead' 'chl (BL3)'; 'chlsyn' 'PEsyn' 'chl, Syn equivalents (log10)' 'PE, Syn equivalents (log10)'; 'size' 'FSCW' 'log10 side scatter rel. 1 um bead' 'FSC-W'};
for v = 1:3
    nexttile; hold on
    for j = 1:4
        g = E.class == ord(j); x = find(g); x = x(1:max(1, ceil(numel(x) / 15000)):end);
        scatter(E.(pv{v, 1})(x), E.(pv{v, 2})(x), 4, cc{j}, 'filled', 'MarkerFaceAlpha', .35)
    end
    xlabel(pv{v, 3}); ylabel(pv{v, 4}); grid on; box on
    if v == 2, xline(0, 'k:'); yline(0, 'k:'); end
    if v == 1, legend({'Euk' 'Syn' 'Euk with Syn' 'Syn on particle'}, 'Location', 'northwest', 'FontSize', 7); end
end
sgtitle(sprintf('%s: Euk-with-Syn (class 5) vs Euks and Syn on particle', tag))
exportgraphics(f, fullfile(outdir, ['euksyn_' tag '.png']), 'Resolution', 90); close(f)
end
