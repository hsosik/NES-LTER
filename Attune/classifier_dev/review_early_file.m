function [new, info, old] = review_early_file(cruisedir, fcsname, outdir)
% function [new, info, old] = review_early_file(cruisedir, fcsname, outdir)
% Review figures for one file of a cruise before the Aug 2019 filter change, as
% the classifier examples are reviewed: plot_class_cytograms (existing classes,
% v2, disagreements) and plot_v2_rules (v2 with its rule boundaries). The cruise
% configuration comes from attune_config (cfg.cruise, GL1 merge only where GL1 is
% OD2 scatter). Writes outdir/<CRUISE>_<file>_cytograms.png and _rules.png
%
% Nov 2026, classifier redesign: pre-2019 channel layout

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cruise = extractAfter(cruisedir, '_');
cfg = classifier_config();
cfg.cruise = cruise;
cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
cpath = fullfile(root, cruisedir);
bp = fullfile(cpath, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
k = find(strcmp(FCSfileinfo.fcslist, fcsname));
pm = 0;
if ismember('pro_measured', FCSfileinfo.Properties.VariableNames), pm = FCSfileinfo.pro_measured(k); end
[fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fcsname));
ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
[new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, pm);
c = load(fullfile(bp, 'class', regexprep(fcsname, '\.fcs$', '.mat')), 'class');
old = c.class(:, 1); old(~ismember(old, [0 1 2 3 4 7])) = 0;
either = old ~= 0 | new ~= 0;
ag = sum(old == new & old ~= 0) / sum(either);
ttl = sprintf('%s   agreement %.3f   v2 adds %d, drops %d, relabels %d of %d phyto events   trigger %s %g, %s %g (file names; physical edges %s)   noise pile %.2f%s', ...
    fcsname, ag, sum(old == 0 & new ~= 0), sum(old ~= 0 & new == 0), sum(old ~= 0 & new ~= 0 & old ~= new), sum(either), ...
    fcshdr.tr1_par, fcshdr.tr1_level, fcshdr.tr2_par, fcshdr.tr2_level, strjoin(fieldnames(info.trigger_edges)', ' '), info.edge_pile_frac, ...
    repmat(' (fit capped)', 1, info.edge_pile_capped));
if ~exist(outdir, 'dir'), mkdir(outdir), end
stem = fullfile(outdir, sprintf('%s_%s', cruise, regexprep(fcsname, '\.fcs$', '')));
plot_class_cytograms(fcsdat, fcshdr, c.class(:, 1), new, info, ttl, [stem '_cytograms.png']);
plot_v2_rules(fcsdat, fcshdr, new, info, ttl, [stem '_rules.png']);
fprintf('%s\nold counts [1 2 3 4 7] %s\nv2 counts [1 2 3 4 5 7 8] %s\nsyn_rescued %d\n', ttl, mat2str(arrayfun(@(x) sum(old == x), [1 2 3 4 7])), ...
    mat2str(arrayfun(@(x) sum(new == x), [1 2 3 4 5 7 8])), info.syn_rescued);
end
