% compare_bead_prototype
% Summarize classifier_dev/output/bead_proto.mat (run_bead_prototype):
%   - QC pass rates per cruise, old vs new (per population), learned shifts
%   - QC failure reasons
%   - peak (KDE mode) vs median: offset and within-cruise repeatability on
%     the bead/channel pairs the classifier uses
%   - new peak and median vs old median; 1 um peak vs old centend2 modes
%     (the size merge inputs)
%
% Nov 2026, bead stats redesign prototype

load(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'output', 'bead_proto.mat'), 'R', 'sig')
pn = {'b05' 'b1' 'b6'};

fprintf('\n== per cruise: runs, old QC pass, new QC pass per population, cruise shift on gating channels\n')
fprintf('%-8s %4s %4s %5s %4s %4s %4s | %11s %11s %11s\n', 'cruise', 'runs', 'gen', 'oldQC', 'q05', 'q1', 'q6', 'shift05', 'shift1', 'shift6')
cr = unique(R.cruise, 'stable');
for i = 1:numel(cr)
    T = R(R.cruise == cr(i), :);
    s = cell(1, 3);
    for b = 1:3
        g = split(T.([pn{b} '_gate'])(1), '/');
        s{b} = sprintf('%5.2f/%5.2f', T.([pn{b} '_shift_' char(g(1))])(1), T.([pn{b} '_shift_' char(g(2))])(1));
    end
    fprintf('%-8s %4d %4s %5d %4d %4d %4d | %s %s %s\n', cr(i), height(T), T.old_gen(1), sum(T.old_QC == 1), ...
        sum(T.b05_QC), sum(T.b1_QC), sum(T.b6_QC), s{:})
end
fprintf('total runs %d; old QC pass %d; new QC pass 0.5 %d, 1 %d, 6 %d, all three %d\n', height(R), ...
    sum(R.old_QC == 1), sum(R.b05_QC), sum(R.b1_QC), sum(R.b6_QC), sum(R.b05_QC & R.b1_QC & R.b6_QC))

fprintf('\n== QC failure reasons (top 8 per population)\n')
for b = 1:3
    w = R.([pn{b} '_why']); w = w(~R.([pn{b} '_QC']));
    [u, ~, j] = unique(w); c = accumarray(j, 1); [~, o] = sort(c, 'descend');
    for k = o(1:min(8, end))', fprintf('%-4s %4d  %s\n', pn{b}, c(k), u(k)); end
end

pairs = {1 'BL2'; 1 'GL2'; 1 'GL4'; 2 'FSC'; 2 'SSC'; 2 'GL1'; 2 'BL1'; 2 'BL3'; 3 'GL3'};
fprintf('\n== classifier pairs, QC-passed runs (log10)\n')
fprintf('%-4s %-4s %4s | %8s %8s | %8s %8s | %8s %8s %8s\n', 'pop', 'ch', 'n', 'pk-md', '|pk-md|95', 'SDw pk', 'SDw md', 'pk-old', 'md-old', '|pk-old|95')
for i = 1:size(pairs, 1)
    b = pairs{i, 1}; c = pairs{i, 2};
    ok = R.([pn{b} '_QC']);
    pk = R.([pn{b} '_pk_' c '_H']); md = R.([pn{b} '_md_' c '_H']); od = R.(sprintf('old_b%d_%s', b, c));
    [sp, sm] = within_sd(R.cruise(ok), pk(ok), md(ok));
    d = pk - md; e = pk - od; f = md - od; oo = ok & isfinite(od);
    fprintf('%-4s %-4s %4d | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f %8.4f\n', pn{b}, c, sum(ok), median(d(ok)), ...
        prctile(abs(d(ok)), 95), sp, sm, median(e(oo)), median(f(oo)), prctile(abs(e(oo)), 95))
end
fprintf('SDw = pooled within-cruise SD over runs\n')

fprintf('\n== 1 um size-merge inputs: new peak vs old centend2 mode (runs with both)\n')
for c = {'SSC_A' 'SSC_H' 'GL1_A' 'GL1_H'}
    pk = R.(['b1_pk_' c{1}]); md = R.(['b1_md_' c{1}]); oc = R.(['old_c1_' c{1}]);
    ok = R.b1_QC & isfinite(oc);
    [sp, sm] = within_sd(R.cruise(ok), pk(ok), oc(ok));
    fprintf('%-6s n=%3d  pk-centend2 median %7.4f |.|95 %6.4f   md-centend2 %7.4f   SDw pk %6.4f centend2 %6.4f\n', c{1}, sum(ok), ...
        median(pk(ok) - oc(ok)), prctile(abs(pk(ok) - oc(ok)), 95), median(md(ok) - oc(ok)), sp, sm)
end

fprintf('\n== QC-passed runs with |peak - median| > 0.02 on a classifier pair\n')
for i = 1:size(pairs, 1)
    b = pairs{i, 1}; c = pairs{i, 2};
    d = R.([pn{b} '_pk_' c '_H']) - R.([pn{b} '_md_' c '_H']);
    k = find(R.([pn{b} '_QC']) & abs(d) > 0.02);
    if isempty(k), continue, end
    [u, ~, j] = unique(R.cruise(k));
    fprintf('%-4s %-4s %3d runs: %s\n', pn{b}, c, numel(k), strjoin(u + "(" + accumarray(j, 1) + ")", ' '))
end

function [sa, sb] = within_sd(g, a, b)
[~, ~, j] = unique(g);
ma = accumarray(j, a, [], @mean); mb = accumarray(j, b, [], @mean);
ra = a - ma(j); rb = b - mb(j);
dof = numel(a) - max(j);
sa = sqrt(sum(ra .^ 2) / dof); sb = sqrt(sum(rb .^ 2) / dof);
end
