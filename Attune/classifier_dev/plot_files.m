function plot_files(L, prefix, outdir)
% function plot_files(L, prefix, outdir)
% Draw plot_class_cytograms for a list of files: L is an n x 2 cell array of
% {cruise directory, fcs file name}; images go to outdir/examples/<prefix><i>_<cruise>.png
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
for i = 1:size(L, 1)
    cpath = fullfile(root, L{i, 1}); bp = fullfile(cpath, 'bead_calibrated'); fn = L{i, 2};
    bs = load_beadstat(bp);
    load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo'); k = find(strcmp(FCSfileinfo.fcslist, fn));
    [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
    pm = 0; if ismember('pro_measured', FCSfileinfo.Properties.VariableNames), pm = FCSfileinfo.pro_measured(k); end
    [new, info] = assign_class_v2(d, h, ref, cfg, pm);
    c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class'); old = c.class(:, 1);
    o = old; o(~ismember(o, [0 1 2 3 4 7])) = 0;
    ag = sum(o == new & o ~= 0) / sum(o ~= 0 | new ~= 0);
    ttl = sprintf('%s   agreement %.3f   old low-PE %d / v2 dim crypto %d, old high-PE %d / v2 bright crypto %d, v2 Euk+Syn %d   trigger %s %g', fn, ag, sum(o == 3), sum(new == 3), ...
        sum(o == 4), sum(new == 4), sum(new == 5), h.tr1_par, h.tr1_level);
    plot_class_cytograms(d, h, old, new, info, ttl, fullfile(outdir, 'examples', sprintf('%s%d_%s.png', prefix, i, extractAfter(L{i, 1}, '_'))));
end
end
