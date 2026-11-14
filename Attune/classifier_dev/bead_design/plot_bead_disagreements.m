function L = plot_bead_disagreements(maxplots)
% function L = plot_bead_disagreements(maxplots)
% Gating plots (plot_bead_gates) for runs where the prototype and the old
% beadstat disagree: old QC pass but a new population fails, old QC fail
% but all three new populations pass, or new peak vs old median differing by
% > 0.03 log on a bead/channel pair the classifier uses.
% Writes classifier_dev/output/bead_proto_plots/<cruise>_<file>.png and returns the list.
%
% Nov 2026, bead stats redesign prototype

if nargin < 1, maxplots = 60; end
here = fileparts(fileparts(mfilename('fullpath'))); % classifier_dev
addpath(here, fullfile(here, '..'))
load(fullfile(here, 'output', 'bead_proto.mat'), 'R', 'sig')
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
outdir = fullfile(here, 'output', 'bead_proto_plots');
if ~exist(outdir, 'dir'), mkdir(outdir), end

allnew = R.b05_QC & R.b1_QC & R.b6_QC;
reason = strings(height(R), 1);
reason(R.old_QC == 1 & ~allnew) = "old pass, new fail";
reason(R.old_QC == 0 & allnew) = "old fail, new pass";
pairs = {1 'BL2'; 1 'GL2'; 1 'GL4'; 2 'FSC'; 2 'SSC'; 2 'GL1'; 2 'BL1'; 2 'BL3'; 3 'GL3'};
pn = {'b05' 'b1' 'b6'};
for i = 1:size(pairs, 1)
    b = pairs{i, 1}; c = pairs{i, 2};
    d = abs(R.([pn{b} '_pk_' c '_H']) - R.(sprintf('old_b%d_%s', b, c)));
    k = R.old_QC == 1 & R.([pn{b} '_QC']) & d > 0.03 & reason == "";
    reason(k) = sprintf("%s %s new-old > 0.03", pn{b}, c);
end
L = R(reason ~= "", {'cruise' 'file' 'old_QC' 'b05_QC' 'b1_QC' 'b6_QC' 'b05_why' 'b1_why' 'b6_why'});
L.reason = reason(reason ~= "");
fprintf('%d runs disagree\n', height(L))
% spread the plots over cruises
[~, ord] = sortrows([findgroups(L.cruise) (1:height(L))']);
pick = ord(1:min(maxplots, end));
D = dir([root '20*_*']); dnames = {D.name};
sig2 = sig; sig2.window = sig.pass2.window; sig2.qc.max_offset = sig.pass2.max_offset; sig2.sig_tol = sig.pass2.sig_tol;
fig = figure('visible', 'off', 'position', [0 0 1300 470]);
for k = pick'
    cdir = dnames{endsWith(dnames, "_" + L.cruise(k))};
    T = R(R.cruise == L.cruise(k), :);
    s = sig2;
    for ip = 1:numel(s.pop)
        st = struct();
        for c = s.chans, st.(c{1}) = T.([s.pop(ip).name '_shift_' c{1}])(1); end
        s.pop(ip).shift = st;
    end
    [d, h] = fca_readfcs(fullfile(root, cdir, 'FCS', char(L.file(k))));
    [cl, P] = find_beads(d, h, s, attune_config(L.cruise(k)));
    figure(fig), plot_bead_gates(d, h, cl, P, sprintf('%s  %s  (old QC %d): %s', L.cruise(k), L.file(k), L.old_QC(k), L.reason(k)), s.window)
    exportgraphics(fig, fullfile(outdir, regexprep(sprintf('%s_%s.png', L.cruise(k), L.file(k)), '\.fcs', '')), 'resolution', 70)
end
writetable(L, fullfile(outdir, 'disagreements.csv'))
end
