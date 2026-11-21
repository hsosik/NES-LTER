function comp_detail(L, outname, sigover)
% function comp_detail(L, outname, sigover)
% Every mixture component of the listed files (L: n x 2 cell, cruise folder and
% file name), classified with pro_measured = 1 and cfg.sig fields overridden by
% sigover (e.g. struct('pro_max_chlSSC_below_syn', 9) to see all Pro candidates):
% weight, events, old-class and v2 Pro inside it, mean size/chl/PE, chl:scatter and
% size relative to the Syn components, height above the chl trigger, SDs, label.
% Writes output/<outname>.txt
%
% Nov 2026, classifier redesign prototype
if nargin < 3, sigover = struct(); end
sp = fullfile(fileparts(mfilename('fullpath')), 'output');
addpath(fileparts(mfilename('fullpath')), fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
R = cell(size(L, 1), 1);
parfor i = 1:size(L, 1)
    cfg = classifier_config();
    fn0 = fieldnames(sigover);
    for q = 1:numel(fn0), cfg.sig.(fn0{q}) = sigover.(fn0{q}); end
    cfg.cruise = extractAfter(L{i, 1}, '_');
    cfg.gl1_is_od2_scatter = strcmpi(attune_config(cfg.cruise).OD2, 'GL1');
    cpath = fullfile(root, L{i, 1}); bp = fullfile(cpath, 'bead_calibrated'); fn = L{i, 2};
    s = load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo'); FI = s.FCSfileinfo; k = find(strcmp(FI.fcslist, fn), 1);
    try
        [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
        ref = bead_reference(load_beadstat(bp), FI.matdate_start(k), cfg);
        [new, info] = assign_class_v2(d, h, ref, cfg, 1);
        F = info.features; m = info.comp; e = NaN; if isfield(info.trigger_edges, 'BL3'), e = info.trigger_edges.BL3; end
        cf = fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat'));
        old7 = false(size(new)); if exist(cf, 'file'), c = load(cf, 'class'); old7 = c.class(:, 1) == 7; end
        hc = info.trace.has_chl;
        comp = zeros(size(new)); comp(hc) = cluster(info.gmm, [F.size(hc) F.BL3(hc) F.PEchl(hc) F.GL2GL3(hc)]);
        syn = m.label == 2; sync = sum(m.chlSSC(syn) .* m.weight(syn)) / sum(m.weight(syn)); syns = sum(m.SSC(syn) .* m.weight(syn)) / sum(m.weight(syn));
        str = sprintf('%s  edge %.2f  synChlSSC %.2f synSize %.2f  v2Pro %d old %d  Pro-like smallnoPE events %d\n', fn, e, sync, syns, sum(new == 7), sum(old7), ...
            sum(F.size < syns & F.BL3 < cfg.sig.pro_max_BL3 & F.GL2 < -1.8 & hc));
        str = [str sprintf('   c   wt     n  old7 new7 | size  BL3   GL2  PEchl GL2GL3 dchl dsize  aboveEdge sdS  sdB  lab\n')];
        for j = 1:height(m)
            S = info.gmm.Sigma(:, :, j);
            str = [str sprintf('  %2d %.3f %6d %5d %5d | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  %.2f %.2f %d\n', j, m.weight(j), sum(comp == j), sum(old7 & comp == j), sum(new == 7 & comp == j), ...
                m.SSC(j), m.BL3(j), m.GL2(j), m.PEchl(j), m.GL2GL3(j), m.chlSSC(j) - sync, m.SSC(j) - syns, m.BL3(j) - e, sqrt(S(1, 1)), m.sd_BL3(j), m.label(j))]; %#ok<AGROW>
        end
        str = [str sprintf('  old Pro by new class: %s\n', mat2str(arrayfun(@(x) sum(old7 & new == x), 0:8)))];
        R{i} = str;
    catch ME
        R{i} = sprintf('%s FAILED %s\n', fn, ME.message);
    end
end
fid = fopen(fullfile(sp, [outname '.txt']), 'w'); fprintf(fid, '%s', R{:}); fclose(fid);
end
