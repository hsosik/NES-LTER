function plot_class_cytograms(fcsdat, fcshdr, old, new, info, ttl, outfile)
% function plot_class_cytograms(fcsdat, fcshdr, old, new, info, ttl, outfile)
% One sample, three rows x six views on the common bead-relative scale:
%   row 1 existing classes, row 2 assign_class_v2 classes,
%   row 3 only the events where the two disagree, coloured by kind of change
% Views (each chosen to separate particular groups):
%   size vs chl      Pro / Syn / Euks / detritus (chl per size)
%   chl vs PE        Syn and PE euks vs Euks without PE
%   GL3 vs GL2       PE spectrum: Syn vs cryptophyte-type PE
%   size vs PE:chl   PE euks (larger, less PE per chl) vs Syn
%   chl vs BL2       broadband fluorescence: detritus / noise vs pigment
%   GL3 vs BL2       PE spectrum without GL2 (on scale where GL2 saturates):
%                    cryptophyte-type PE (dim and bright cryptos) sits below the Syn line
%   size vs chl relative to the Euk chl-size line (fitted to v2 Euks): the Euk band
%                    and its bottom (chl gap), the Euk-with-Syn cut 0.5 below the line
%                    (below it: unassigned, flag 1), the Euk minimum size
%   GL2 vs BL2       Syn spectral types A/B/C (v2 row coloured by type)
% Dashed lines mark trigger edges. Class 0 events are drawn in light grey.
%
% Nov 2026, classifier redesign prototype

F = info.features;
pn = {fcshdr.par.name};
% chl relative to the Euk chl-size line, fitted to the v2 Euks (same line in all rows)
eb = [NaN NaN];
if sum(new == 1) >= 50, eb = robustfit(F.size(new == 1), F.BL3(new == 1))'; end
F.dEuk = F.BL3 - (eb(1) + eb(2) * F.size);
old(~ismember(old, [0 1 2 3 4 5 7 8])) = 0;
codes = [1 2 3 4 5 7 8];
names = {'Euk' 'Syn' 'dim crypto' 'bright crypto' 'Euk with Syn' 'Pro' 'Syn on particle'};
cols = [0 .6 0; 1 .55 0; 0 .75 .9; .95 0 .75; .45 .1 .75; 0 .35 1; .85 0 0]; %Euk green, Syn orange, dim crypto cyan, bright crypto magenta, Euk with Syn purple, Pro blue, Syn-on-particle red
% rare classes are drawn last, larger, opaque and outlined so they stand out
% over the abundant ones; unassigned events are drawn in mid grey beneath
rare = ismember(codes, [3 4 5 8]);
views = {'size' 'BL3' 'log10 side scatter rel. 1 um bead (SSC; GL1 when large)' 'chl (BL3)'; 'BL3' 'GL2' 'chl (BL3)' 'PE (GL2)'; 'GL3' 'GL2' 'GL3' 'PE (GL2)'; ...
    'size' 'PEchl' 'log10 side scatter rel. 1 um bead (SSC; GL1 when large)' 'PE:chl (GL2-BL3)'; 'BL3' 'BL2' 'chl (BL3)' 'BL2 (green fl., blue laser)'; 'GL3' 'BL2' 'GL3' 'BL2 (PE, blue laser)'; 'size' 'dEuk' 'log10 side scatter rel. 1 um bead (SSC; GL1 when large)' 'chl minus Euk chl-size line'; ...
    'GL2' 'BL2' 'PE (GL2, green laser)' 'PE (BL2, blue laser)'};
% last view: Syn spectral subtypes (v2 row: Syn coloured by type A/B/C/D; dashed lines = type boundaries)
stcols = [0 .45 .25; .5 .2 0; 1 .55 0; .95 .85 0]; vst = size(views, 1);
edge = info.trigger_edges;

% axis limits from all events, clipped to the populated range
lim = zeros(size(views, 1), 4);
for v = 1:size(views, 1)
    x = F.(views{v, 1}); y = F.(views{v, 2});
    if all(~isfinite(y)), y = zeros(size(y)); end
    lim(v, :) = [prctile(x, 0.2) - 0.1, max(x) + 0.1, prctile(y, 0.2) - 0.1, max(y) + 0.1];
    if strcmp(views{v, 2}, 'dEuk'), lim(v, 3:4) = [max(lim(v, 3), -3) min(lim(v, 4), 2)]; end
end
% subsample class-0 background for speed
bg = find(old == 0 & new == 0);
bg = bg(1:max(1, ceil(numel(bg) / 40000)):end);

f = figure('Visible', 'off', 'Position', [20 20 2700 1050]);
tl = tiledlayout(3, size(views, 1), 'TileSpacing', 'tight', 'Padding', 'compact');
lab = {old, new};
rowname = {'existing' 'v2'};
for r = 1:2
    c = lab{r};
    % draw in order of decreasing abundance so rare classes stay visible
    [~, ord] = sort(arrayfun(@(k) sum(c == k), codes), 'descend');
    for v = 1:size(views, 1)
        nexttile; hold on
        z = find(c == 0); z = z(1:max(1, ceil(numel(z) / 40000)):end);
        scatter(F.(views{v, 1})(z), F.(views{v, 2})(z), 3, [.55 .55 .55], 'filled', 'MarkerFaceAlpha', .35)
        for o = ord(~rare(ord))
            i = c == codes(o);
            scatter(F.(views{v, 1})(i), F.(views{v, 2})(i), 3, cols(o, :), 'filled', 'MarkerFaceAlpha', .3)
        end
        for o = ord(rare(ord))
            i = c == codes(o);
            scatter(F.(views{v, 1})(i), F.(views{v, 2})(i), 14, cols(o, :), 'filled', 'MarkerFaceAlpha', .9, ...
                'MarkerEdgeColor', 'k', 'LineWidth', .25)
        end
        if r == 2 && v == vst && isfield(info, 'syn_type')
            % Syn by spectral subtype, drawn over the class colours
            st = info.syn_type;
            for t = 1:4
                i = c == 2 & st == t;
                scatter(F.GL2(i), F.BL2(i), 3, stcols(t, :), 'filled', 'MarkerFaceAlpha', .5)
            end
            nst = arrayfun(@(t) sum(c == 2 & st == t), [1 2 3 4 9]);
            text(0.02, 0.98, sprintf('Syn type A %d\nB %d\nC %d\nD %d\nuntyped %d', nst), 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7, 'BackgroundColor', [1 1 1 .7])
        end
        decorate(v)
        if v == 1
            ylabel(sprintf('%s\n%s', rowname{r}, views{v, 4}), 'FontWeight', 'bold')
            n = arrayfun(@(k) sum(c == k), codes);
            txt = arrayfun(@(k) sprintf('%s %d', names{k}, n(k)), find(n > 0), 'UniformOutput', false);
            text(0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7, 'BackgroundColor', [1 1 1 .7])
        end
    end
end
% disagreement row
add = old == 0 & new ~= 0; drop = old ~= 0 & new == 0; swap = old ~= 0 & new ~= 0 & old ~= new;
dc = [.85 0 0; 0 .3 1; .55 0 .75];
for v = 1:size(views, 1)
    nexttile; hold on
    scatter(F.(views{v, 1})(bg), F.(views{v, 2})(bg), 2, [.85 .85 .85], 'filled', 'MarkerFaceAlpha', .2)
    agreed = find(old == new & old ~= 0); agreed = agreed(1:max(1, ceil(numel(agreed) / 20000)):end);
    scatter(F.(views{v, 1})(agreed), F.(views{v, 2})(agreed), 2, [.55 .55 .55], 'filled', 'MarkerFaceAlpha', .25)
    g = {add, drop, swap};
    for j = 1:3
        scatter(F.(views{v, 1})(g{j}), F.(views{v, 2})(g{j}), 4, dc(j, :), 'filled', 'MarkerFaceAlpha', .6)
    end
    decorate(v)
    xlabel(views{v, 3})
    if v == 1
        ylabel(sprintf('disagreements\n%s', views{v, 4}), 'FontWeight', 'bold')
        % most common transitions
        d = find(old ~= new);
        [u, ~, k] = unique([old(d) new(d)], 'rows');
        cnt = accumarray(k, 1); [cnt, s] = sort(cnt, 'descend'); u = u(s, :);
        txt = arrayfun(@(i) sprintf('%s → %s  %d', cname(u(i, 1)), cname(u(i, 2)), cnt(i)), 1:min(5, numel(cnt)), 'UniformOutput', false);
        text(0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7, 'BackgroundColor', [1 1 1 .7])
    end
    if v == size(views, 1)
        h = arrayfun(@(j) scatter(NaN, NaN, 20, dc(j, :), 'filled'), 1:3);
        ha = scatter(NaN, NaN, 20, [.55 .55 .55], 'filled');
        legend([ha h], {'agreed phyto' 'v2 adds (old 0)' 'v2 drops (to 0)' 'relabelled'}, 'Location', 'southeast', 'FontSize', 7)
    end
end
% class legend on the first row
nexttile(size(views, 1)); hold on
h = arrayfun(@(o) scatter(NaN, NaN, 20 + 20 * rare(o), cols(o, :), 'filled', 'MarkerEdgeColor', tern(rare(o), 'k', 'none')), 1:numel(codes));
h0 = scatter(NaN, NaN, 20, [.55 .55 .55], 'filled');
legend([h0 h], [{'unassigned / not phyto'} names], 'Location', 'southeast', 'FontSize', 7)
title(tl, ttl, 'Interpreter', 'none', 'FontSize', 11)
exportgraphics(f, outfile, 'Resolution', 80)
close(f)

    function decorate(v)
        axis(lim(v, :)); grid on; box on
        if v == vst && isfield(info, 'syn_type_edges')
            xx = lim(v, 1:2);
            for e2 = info.syn_type_edges(:)', plot(xx, xx + e2, 'k--', 'LineWidth', .8), end
        end
        if strcmp(views{v, 2}, 'dEuk') && all(isfinite(eb))
            % Euk line, Euk-with-Syn cut, Euk chl gap (a BL3 level, sloped here), Euk minimum size
            xx = linspace(lim(v, 1), lim(v, 2), 50);
            plot(lim(v, 1:2), [0 0], '-', 'Color', [0 .5 0], 'LineWidth', .8)
            if isfield(info, 'rules') && isfield(info.rules, 'euksyn_max_chl_below_euk')
                plot(lim(v, 1:2), -info.rules.euksyn_max_chl_below_euk * [1 1], '--', 'Color', [.45 .1 .75], 'LineWidth', .8)
            end
            if isfield(info, 'euk_chl_gap') && isfinite(info.euk_chl_gap)
                plot(xx, info.euk_chl_gap - (eb(1) + eb(2) * xx), ':', 'Color', [0 .5 0], 'LineWidth', 1)
            end
            if isfield(info, 'rules') && isfield(info.rules, 'euk_min_size')
                xline(info.rules.euk_min_size, '-.', 'Color', [0 .5 0]);
            end
        end
        set(gca, 'FontSize', 7)
        ylabel(views{v, 4}) %first column overrides this with the row name
        % trigger edges where the view's axis is a trigger channel
        for ax = 1:2
            ch = views{v, ax};
            if isfield(edge, ch)
                if ax == 1, xline(edge.(ch), 'k--'); else, yline(edge.(ch), 'k--'); end
            end
        end
    end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function s = cname(k)
switch k
    case 0, s = 'none';
    case 1, s = 'Euk';
    case 2, s = 'Syn';
    case 3, s = 'dim crypto';
    case 4, s = 'bright crypto';
    case 5, s = 'Euk+Syn';
    case 7, s = 'Pro';
    case 8, s = 'Syn-particle';
    otherwise, s = num2str(k);
end
end
