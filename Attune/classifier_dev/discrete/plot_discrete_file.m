function plot_discrete_file(cruisedir, fn, cfg, outfile)
% function plot_discrete_file(cruisedir, fn, cfg, outfile)
% Cytograms of one preserved run: hand gates (top) vs assign_class_v2 with cfg
% (bottom; assign_class_discrete), events inside the hand gates' time window only.
% e.g. plot_discrete_file('20240503_EN715', 'NESLTER_EN715_May2024_preserved(1)_phyto_CHL_SSC_C10N02.fcs')
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 3 || isempty(cfg), cfg = discrete_config(); end
if nargin < 4, outfile = fullfile(od, 'examples', [regexprep(fn, '\.fcs$', '') '.png']); end
if ~exist(fileparts(outfile), 'dir'), mkdir(fileparts(outfile)), end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cruise = extractAfter(cruisedir, '_');
cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
G = load(fullfile(root, cruisedir, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
G.Properties.VariableNames = lower(G.Properties.VariableNames);
i = find(strcmp(G.fcslist, fn));
[fcsdat, fcshdr] = fca_readfcs(fullfile(root, cruisedir, 'preserved', 'FCS', fn));
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
ref = discrete_bead_reference(lab, datenum([fcshdr.date ' ' fcshdr.starttime]), cfg);
rt = 'CHL'; if contains(fn, 'pro', 'IgnoreCase', true), rt = 'pro'; elseif contains(fn, 'phyto_PE'), rt = 'PE'; end
[auto, info] = assign_class_discrete(fcsdat, fcshdr, ref, cfg, rt);
h = double(G.class{i}(:)); h(h == 4) = 7; h(ismember(h, [5 6])) = 34; h(h == 3) = 0;
a = double(auto); a(ismember(a, [3 4 5])) = 34;
late = true(size(h)); if any(h), late = fcsdat(:, 1) >= min(fcsdat(h ~= 0, 1)); end
F = info.features(late, :); h = h(late); a = a(late);
cs = [0 1 2 7 34 8]; cols = [.75 .75 .75; 0 .6 0; 1 0 0; .4 .4 .4; 1 0 1; 1 .6 0];
pairs = {'size' 'BL3'; 'BL3' 'GL2'; 'size' 'PEchl'; 'PEchl' 'GL2GL3'; 'size' 'GL2'};
fig = figure('visible', 'off', 'position', [50 50 1500 650]);
for row = 1:2
    lab_ = h; if row == 2, lab_ = a; end
    for k = 1:size(pairs, 1)
        subplot(2, size(pairs, 1), (row - 1) * size(pairs, 1) + k), hold on
        for c = 1:numel(cs)
            j = lab_ == cs(c);
            plot(F.(pairs{k, 1})(j), F.(pairs{k, 2})(j), '.', 'color', cols(c, :), 'markersize', 3)
        end
        xlabel(pairs{k, 1}), ylabel(pairs{k, 2}), grid on
        if k == 1
            if row == 1, t = 'hand'; else, t = 'auto'; end
            title(sprintf('%s: Euk %d Syn %d Pro %d PEeuk %d SynPart %d', t, sum(lab_ == 1), sum(lab_ == 2), sum(lab_ == 7), sum(lab_ == 34), sum(lab_ == 8)), 'fontsize', 8)
        end
    end
end
annotation('textbox', [0 0.95 1 0.05], 'string', [fn '  (green Euk, red Syn, dark gray Pro, magenta PE euk, orange Syn on particle)'], ...
    'edgecolor', 'none', 'horizontalalignment', 'center', 'interpreter', 'none')
exportgraphics(fig, outfile, 'resolution', 80), close(fig)
end
