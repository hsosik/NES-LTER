function plot_v2_vs_old(cruisedir, fcsname, outfile)
% function plot_v2_vs_old(cruisedir, fcsname, outfile)
% Side-by-side cytograms (bead-normalized units) coloured by the existing
% class and by assign_class_v2, with trigger edges marked.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
cpath = fullfile(root, cruisedir);
bp = fullfile(cpath, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
k = find(strcmp(FCSfileinfo.fcslist, fcsname));
[fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fcsname));
ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
pm = 0;
if ismember('pro_measured', FCSfileinfo.Properties.VariableNames), pm = FCSfileinfo.pro_measured(k); end
[new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, pm);
F = info.features;
c = load(fullfile(bp, 'class', regexprep(fcsname, '\.fcs$', '.mat')), 'class');
old = c.class(:, 1); old(old == 4) = 3;

cols = containers.Map({0 1 2 3 5 6 7}, {[.7 .7 .7] [0 .6 0] [.9 .4 0] [.85 0 .85] [.7 .7 .7] [.7 .7 .7] [0 .3 1]});
pairs = {'BL3' 'GL2'; 'SSC' 'BL3'; 'GL3' 'GL2'; 'SSC' 'GL2GL3'};
f = figure('Visible', 'off', 'Position', [50 50 1700 800]);
tl = tiledlayout(2, size(pairs, 1), 'TileSpacing', 'compact');
labs = {old, new}; ttl = {'existing' 'v2'};
for r = 1:2
    for p = 1:size(pairs, 1)
        nexttile; hold on
        L = labs{r};
        for u = [0 1 2 7 3] %PE euks last so the rare ones are visible
            idx = find(L == u);
            if isempty(idx), continue, end
            if u ~= 3, idx = idx(1:max(1, ceil(numel(idx) / 5000)):end); end
            sz = 2 + 6 * (u == 3);
            scatter(F.(pairs{p, 1})(idx), F.(pairs{p, 2})(idx), sz, cols(u), 'filled', 'MarkerFaceAlpha', .4)
        end
        for e = fieldnames(info.trigger_edges)'
            if strcmp(pairs{p, 1}, e{1}), xline(info.trigger_edges.(e{1}), '--r'); end
            if strcmp(pairs{p, 2}, e{1}), yline(info.trigger_edges.(e{1}), '--r'); end
        end
        xlabel(pairs{p, 1}); ylabel(pairs{p, 2}); grid on
        if p == 1, title(ttl{r}); end
    end
end
title(tl, sprintf('%s  grey=0 green=Euk orange=Syn magenta=PE euk blue=Pro; red dashed = trigger edge', fcsname), 'Interpreter', 'none')
exportgraphics(f, outfile, 'Resolution', 85)
close(f)
end
