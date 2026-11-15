function S = explore_files(files, outdir)
% function S = explore_files(files, outdir)
% Plot events from a few files in bead-normalized feature space, coloured by
% the existing class assignments, and tabulate per-class feature medians.
% files - cell array of {cruise_folder, fcs_name}
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
cols = containers.Map({0 1 2 3 4 7}, {[.75 .75 .75] [0 .6 0] [.9 .4 0] [.8 0 .8] [.5 0 .5] [0 .3 1]});
S = table;
for i = 1:size(files, 1)
    cpath = fullfile(root, files{i, 1});
    bp = fullfile(cpath, 'bead_calibrated');
    bs = load_beadstat(bp);
    load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
    k = find(strcmp(FCSfileinfo.fcslist, files{i, 2}));
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', files{i, 2}));
    ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
    [F, info] = attune_features(fcsdat, fcshdr, ref);
    c = load(fullfile(bp, 'class', regexprep(files{i, 2}, '\.fcs$', '.mat')), 'class');
    old = c.class(:, 1);
    F.PEchl = F.GL2 - F.BL3;
    F.GL2GL3 = F.GL2 - F.GL3;
    F.chlSSC = F.BL3 - F.SSC;

    for u = unique(old)'
        idx = old == u;
        r = table(string(files{i, 2}), u, sum(idx), 'VariableNames', {'file' 'old_class' 'n'});
        for v = {'SSC' 'GL1' 'BL1' 'BL3' 'GL2' 'GL3' 'PEchl' 'GL2GL3' 'chlSSC'}
            r.(v{1}) = median(F.(v{1})(idx));
        end
        S = [S; r]; %#ok<AGROW>
    end

    f = figure('Visible', 'off', 'Position', [50 50 1600 900]);
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact');
    pairs = {'BL3' 'GL2'; 'SSC' 'BL3'; 'SSC' 'GL2'; 'GL3' 'GL2'; 'SSC' 'PEchl'; 'BL1' 'BL3'};
    for p = 1:size(pairs, 1)
        nexttile; hold on
        for u = unique(old)'
            idx = find(old == u);
            idx = idx(1:max(1, ceil(numel(idx) / 4000)):end);
            scatter(F.(pairs{p, 1})(idx), F.(pairs{p, 2})(idx), 2, cols(u), 'filled', 'MarkerFaceAlpha', .35)
        end
        xlabel(pairs{p, 1}); ylabel(pairs{p, 2}); grid on
        if isfield(info.floor, pairs{p, 1}), xline(info.floor.(pairs{p, 1}), ':k'); end
        if isfield(info.floor, pairs{p, 2}), yline(info.floor.(pairs{p, 2}), ':k'); end
    end
    title(tl, sprintf('%s  (BL3 hv %d, GL2 hv %d, trigger %s %d)  grey=0 green=1 orange=2 magenta=3 purple=4 blue=7', ...
        files{i, 2}, info.hv.BL3, info.hv.GL2, info.trigger, info.trigger_level), 'Interpreter', 'none')
    exportgraphics(f, fullfile(outdir, sprintf('explore_%02d.png', i)), 'Resolution', 90)
    close(f)
end
end
