function H = examine_hbac(cruisedirs, nper)
% function H = examine_hbac(cruisedirs, nper)
% First look at SYBR-stained bacteria runs (hbac_SYBR_SSC) on the bead scale:
% hand gates (Gated_Table class 3 bacteria (+Pro), 1 Euk, 2 Syn, 4 Pro) in
% SSC, BL1 (SYBR), BL3 (chl) and GL2 (PE). Prints feature quantiles per hand
% class and draws pooled and per-run cytograms (output/discrete/hbac_*.png).
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 2, nper = 4; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = discrete_config();
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
E = table; H = table;
for ci = 1:numel(cruisedirs)
    cd_ = cruisedirs{ci}; cruise = extractAfter(cd_, '_');
    cfg.cruise = cruise;
    G = load(fullfile(root, cd_, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    ix = find(contains(G.fcslist, 'hbac') & ~cellfun(@isempty, G.class) & ~contains(G.fcslist, {'skip' 'test' 'badly'}));
    ix = ix(unique(round(linspace(1, numel(ix), min(nper, numel(ix))))));
    for i = ix'
        try, [dat, hdr] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', G.fcslist{i})); catch, continue, end
        ref = discrete_bead_reference(lab, datenum([hdr.date ' ' hdr.starttime]), cfg);
        if isempty(ref), continue, end
        [F, info] = attune_features(dat, hdr, ref);
        h = double(G.class{i}(:));
        late = true(size(h)); if any(h), late = dat(:, 1) >= min(dat(h ~= 0, 1)); end
        F = F(late, :); h = h(late);
        F.hand = h; F.cruise = repmat(string(cruise), height(F), 1);
        k = randperm(height(F), min(20000, height(F)));
        E = [E; F(k, :)]; %#ok<AGROW>
        H = [H; table(string(cruise), string(G.fcslist{i}), height(F), sum(h == 3), sum(h == 4), sum(h == 2), sum(h == 1), ...
            string(info.trigger), string(info.trigger2), 'VariableNames', {'cruise' 'file' 'n' 'bact' 'pro' 'syn' 'euk' 'trig1' 'trig2'})]; %#ok<AGROW>
    end
end
disp(H)
names = {'0' 'Euk' 'Syn' 'bacteria' 'Pro'};
for c = 0:4
    x = E(E.hand == c, :);
    if isempty(x), continue, end
    q = @(v) sprintf('%5.2f %5.2f %5.2f', prctile(v, [5 50 95]));
    fprintf('%-9s n %6d  SSC %s  BL1 %s  BL3 %s  GL2 %s  BL3-BL1 %s\n', names{c + 1}, height(x), q(x.SSC), q(x.BL1), q(x.BL3), q(x.GL2), q(x.BL3 - x.BL1))
end
cols = [.75 .75 .75; 0 .6 0; 1 0 0; 0 0 1; .3 .3 .3];
fig = figure('visible', 'off', 'position', [50 50 1400 420]);
pairs = {'SSC' 'BL1'; 'BL1' 'BL3'; 'SSC' 'BL3'; 'BL1' 'GL2'};
for k = 1:4
    subplot(1, 4, k), hold on
    for c = [0 3 4 2 1]
        j = E.hand == c; if c == 0, j = j & rand(size(j)) < 0.3; end
        plot(E.(pairs{k, 1})(j), E.(pairs{k, 2})(j), '.', 'color', cols(c + 1, :), 'markersize', 1)
    end
    xlabel(pairs{k, 1}), ylabel(pairs{k, 2}), grid on
end
annotation('textbox', [0 0.94 1 0.05], 'string', 'hbac SYBR runs, hand gates: blue bacteria, dark gray Pro, red Syn, green Euk, light gray 0', 'edgecolor', 'none', 'horizontalalignment', 'center')
exportgraphics(fig, fullfile(od, 'hbac_pooled.png'), 'resolution', 90), close(fig)
save(fullfile(od, 'hbac_events.mat'), 'E', 'H')
end
