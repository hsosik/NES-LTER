function plot_v2_rules(fcsdat, fcshdr, class, info, ttl, outfile)
% function plot_v2_rules(fcsdat, fcshdr, class, info, ttl, outfile)
% One sample classified by assign_class_v2, drawn with the rule boundaries that
% separate the classes (values as used for this file, from info.rules). Eight panels:
%   1 size vs chl         detritus and Euk chl:scatter cuts, Euk chl-size line,
%                         Euk minimum size, Syn size window, Syn-on-particle cut and
%                         Syn chl trend, Pro tail box; mixture components (2 SD ellipses)
%   2 chl vs PE           PE detection, Syn PE:chl band, Syn PE floor, Pro PE ceiling,
%                         PE-euk PE:chl ceiling
%   3 GL3 vs GL2          PE spectrum: Syn GL2:GL3 band, crypto GL2:GL3 ceiling
%   4 GL3 vs BL2          PE spectrum read from BL2: Syn line, crypto boundary
%   5 size vs PE:chl      Syn box, PE-euk region, Syn-on-particle box
%   6 crypto split        PE (Syn equivalents) vs spectrum offset for PE euks:
%                         dim / bright groups (2-group mixture ellipses) or PE floor
%   7 size vs chl offset  chl relative to the Euk chl-size line: Euk band, Euk with
%                         Syn -> Syn on particle cut, catch-all Euk floor
%   8 key numbers         counts and the thresholds used
% Class 0 events are drawn in light grey.
%
% Nov 2026, classifier redesign prototype

F = info.features; R = info.rules; e = info.trigger_edges;
pn = {fcshdr.par.name};
codes = [1 2 3 4 5 7 8];
names = {'Euk' 'Syn' 'dim crypto' 'bright crypto' 'Euk with Syn' 'Pro' 'Syn on particle'};
cols = [0 .6 0; 1 .55 0; 0 .8 .8; .9 0 .75; .55 .55 0; 0 .35 1; .6 .1 .1];
ccol = @(k) cols(codes == k, :);
lw = 1.1;

% derived per-event quantities for the crypto panel and the Euk-offset panel
sat = @(ch) fcsdat(:, strcmp(pn, [ch '-H'])) >= R.sat_frac * max(fcsdat(:, strcmp(pn, [ch '-H'])));
if isfield(R, 'syn_BL2GL3')
    dB = (F.BL2 - F.GL3) - R.syn_BL2GL3;
    PE = F.GL2; k = sat('GL2'); PE(k) = F.BL2(k) - R.syn_BL2GL2;
    PEs = PE - R.syn_GL2;
else
    dB = nan(height(F), 1); PEs = nan(height(F), 1);
end
F.dB = dB; F.PEs = PEs;
if isfield(R, 'euk_line')
    F.dEuk = F.BL3 - (R.euk_line(1) + R.euk_line(2) * F.size);
else
    F.dEuk = nan(height(F), 1);
end

xs = 'log10 side scatter rel. 1 um bead';
P = {'size' 'BL3' xs 'chl (BL3)' 'Size vs chl: detritus, Euk, Pro, Syn size, Syn on particle'; ...
     'BL3' 'GL2' 'chl (BL3)' 'PE (GL2)' 'Chl vs PE: PE detection, Syn PE:chl, Pro PE'; ...
     'GL3' 'GL2' 'GL3' 'PE (GL2)' 'PE spectrum (green laser): Syn vs crypto'; ...
     'GL3' 'BL2' 'GL3' 'BL2 (PE, blue laser)' 'PE spectrum (blue laser): crypto boundary'; ...
     'size' 'PEchl' xs 'PE:chl (GL2-BL3)' 'Size vs PE:chl: Syn, PE euks, Syn on particle'; ...
     'PEs' 'dB' 'PE, log10 Syn equivalents (BL2 where GL2 saturates)' 'BL2:GL3 minus Syn' 'PE euks: dim / bright crypto, Euk with Syn'; ...
     'size' 'dEuk' xs 'chl minus Euk chl-size line' 'Chl relative to Euks: Euk band, catch-all floor'};

f = figure('Visible', 'off', 'Position', [20 20 2300 1150]);
tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for v = 1:size(P, 1)
    ax = nexttile; hold on
    x = F.(P{v, 1}); y = F.(P{v, 2});
    show = true(size(class));
    if v == 6, show = ismember(class, [3 4 5]); end
    draw_points(x, y, show)
    if v == 6 && ~any(show), text(.5, .5, 'no PE euks', 'Units', 'normalized', 'HorizontalAlignment', 'center'), end
    L = limits(x(show), y(show), v);
    axis(L); grid on; box on; set(ax, 'FontSize', 7, 'Layer', 'top')
    xlabel(P{v, 3}); ylabel(P{v, 4}); title(P{v, 5}, 'FontSize', 8)
    h = boundaries(v, L);
    h = h(isgraphics(h));
    if ~isempty(h), legend(h, 'Location', 'southeast', 'FontSize', 6, 'Box', 'on', 'AutoUpdate', 'off'), end
end
% panel 8: class legend and numbers
ax = nexttile; axis off; hold on
hc = arrayfun(@(o) scatter(NaN, NaN, 25, cols(o, :), 'filled'), 1:numel(codes));
h0 = scatter(NaN, NaN, 25, [.75 .75 .75], 'filled');
n = arrayfun(@(k) sum(class == k), codes);
legend([h0 hc], [{sprintf('unassigned %d', sum(class == 0))} arrayfun(@(i) sprintf('%s %d', names{i}, n(i)), 1:numel(codes), 'UniformOutput', false)], ...
    'Location', 'northwest', 'FontSize', 8)
text(0, .52, keytext(), 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7, 'FontName', 'FixedWidth', 'Interpreter', 'none')
title(tl, ttl, 'Interpreter', 'none', 'FontSize', 11)
exportgraphics(f, outfile, 'Resolution', 80)
close(f)

    function draw_points(x, y, show)
        z = find(class == 0 & show); z = z(1:max(1, ceil(numel(z) / 40000)):end);
        scatter(x(z), y(z), 2, [.78 .78 .78], 'filled', 'MarkerFaceAlpha', .3, 'HandleVisibility', 'off')
        [~, ord] = sort(arrayfun(@(k) sum(class == k), codes), 'descend');
        for o = ord
            i = find(class == codes(o) & show); i = i(1:max(1, ceil(numel(i) / 30000)):end);
            scatter(x(i), y(i), 3, cols(o, :), 'filled', 'MarkerFaceAlpha', .4, 'HandleVisibility', 'off')
        end
    end

    function L = limits(x, y, v)
        x = x(isfinite(x)); y = y(isfinite(y));
        if isempty(x), L = [-1 1 -1 1]; return, end
        L = [prctile(x, 0.2) - 0.1, max(x) + 0.1, prctile(y, 0.2) - 0.1, max(y) + 0.1];
        if v == 6, L = [min(x) - 0.2, max(x) + 0.2, min(y) - 0.2, max(y) + 0.2]; end
        if v == 7, L(3) = max(L(3), -3); end
    end

    function h = ln(xx, yy, c, st, name)
        h = plot(xx, yy, st, 'Color', c, 'LineWidth', lw, 'DisplayName', name);
    end

    function h = boundaries(v, L)
        h = gobjects(0); xx = L(1:2); sc = R.syn_check;
        switch v
            case 1 % size vs chl
                ellipses(1, 2);
                if isfield(e, 'BL3'), h(end+1) = ln(xx, e.BL3 * [1 1], 'k', '--', 'chl trigger edge'); end
                h(end+1) = ln(xx, xx + R.detritus_chlSSC, [.4 .4 .4], ':', 'detritus cut (chl:scatter)');
                if isfield(R, 'euk_chlSSC_cut'), h(end+1) = ln(xx, xx + R.euk_chlSSC_cut, ccol(1), '--', 'Euk chl:scatter cut'); end
                if isfield(R, 'euk_line'), h(end+1) = ln(xx, R.euk_line(1) + R.euk_line(2) * xx, ccol(1), '-', 'Euk chl-size line'); end
                if isfield(info, 'catchall_euk_floor')
                    b = info.catchall_euk_floor; h(end+1) = ln(xx, b(1) + b(2) * xx + b(3), ccol(1), ':', 'catch-all Euk floor');
                end
                h(end+1) = ln(R.euk_min_size * [1 1], L(3:4), ccol(1), '-.', 'Euk min size (Syn - 0.5)');
                h(end+1) = ln([sc.size_lo sc.size_lo NaN sc.size_hi sc.size_hi], [L(3:4) NaN L(3:4)], ccol(2), '--', 'Syn size window');
                if isfield(R, 'synpart_cut')
                    h(end+1) = ln(R.synpart_cut * [1 1], L(3:4), ccol(8), '-.', 'Syn on particle: min size');
                    sx = [R.synpart_cut L(2)];
                    h(end+1) = ln(sx, R.syn_chl_trend(1) + R.syn_chl_trend(2) * sx - R.synpart_below_trend, ccol(8), '--', 'Syn on particle: below Syn chl trend');
                end
                if isfield(R, 'pro_tail')
                    h(end+1) = ln([L(1) R.pro_tail(2) R.pro_tail(2)], R.pro_tail(1) * [1 1 0] + [0 0 L(3)], ccol(7), ':', 'Pro tail at chl trigger (left 0)');
                end
            case 2 % chl vs PE
                if isfield(e, 'GL2'), h(end+1) = ln(xx, e.GL2 * [1 1], 'k', '--', 'PE trigger edge'); end
                if isfield(e, 'BL3'), xline(e.BL3, 'k--', 'HandleVisibility', 'off'); end
                h(end+1) = ln(xx, R.pe_min * [1 1], [.4 .4 .4], ':', 'PE detected');
                h(end+1) = ln([xx NaN xx], [xx + sc.PEchl_lo NaN xx + sc.PEchl_hi], ccol(2), '--', 'Syn PE:chl band');
                if isfinite(info.syn_PE_floor), h(end+1) = ln(xx, info.syn_PE_floor * [1 1], ccol(2), ':', 'Syn PE floor'); end
                if isfield(info, 'pro_GL2_max'), h(end+1) = ln(xx, info.pro_GL2_max * [1 1], ccol(7), '-.', 'Pro PE ceiling'); end
                h(end+1) = ln(xx, xx + R.peeuk_PEchl_max, ccol(4), '--', 'PE euk: PE:chl below this');
            case 3 % GL3 vs GL2
                if isfield(e, 'GL2'), h(end+1) = ln(xx, e.GL2 * [1 1], 'k', '--', 'PE trigger edge'); end
                h(end+1) = ln([xx NaN xx], [xx + sc.GL2GL3_lo NaN xx + sc.GL2GL3_hi], ccol(2), '--', 'Syn GL2:GL3 band');
                h(end+1) = ln(xx, xx + R.peeuk_max_GL2GL3, ccol(4), '-.', 'crypto-type below (GL2:GL3)');
                h(end+1) = ln(sc.gl3ok_min * [1 1], L(3:4), [.4 .4 .4], ':', 'spectrum used right of this');
            case 4 % GL3 vs BL2
                if isfield(R, 'syn_BL2GL3')
                    h(end+1) = ln(xx, xx + R.syn_BL2GL3, ccol(2), '-', 'Syn BL2:GL3');
                    h(end+1) = ln(xx, xx + R.syn_BL2GL3 + R.crypto_dB, ccol(4), '--', 'crypto-type below');
                end
            case 5 % size vs PE:chl
                h(end+1) = ln([sc.size_lo sc.size_hi sc.size_hi sc.size_lo sc.size_lo], [sc.PEchl_lo sc.PEchl_lo sc.PEchl_hi sc.PEchl_hi sc.PEchl_lo], ccol(2), '--', 'Syn check box');
                h(end+1) = ln([R.peeuk_min_size L(2) NaN R.peeuk_min_size R.peeuk_min_size], [R.peeuk_PEchl_max R.peeuk_PEchl_max NaN L(3) R.peeuk_PEchl_max], ccol(4), '-.', 'PE euk: larger, less PE:chl');
                if isfield(R, 'synpart_cut')
                    h(end+1) = ln([R.synpart_cut L(2) L(2) R.synpart_cut R.synpart_cut], R.synpart_PEchl([1 1 2 2 1]), ccol(8), '--', 'Syn on particle box');
                end
            case 6 % crypto split
                yline(R.crypto_dB, '--', 'Color', ccol(5), 'LineWidth', lw, 'HandleVisibility', 'off');
                h(end+1) = ln(NaN, NaN, ccol(5), '--', 'crypto spectrum below / Euk with Syn above');
                cs = info.crypto_split;
                if isfield(info, 'crypto_gmm')
                    g = info.crypto_gmm; o = info.crypto_gmm_order;
                    st = '--'; if cs.accepted, st = '-'; end
                    for j = 1:2
                        [ex, ey] = ellipse_xy(g.mu(o(j), :), g.Sigma(:, :, o(j)), 2);
                        hh = plot(ex, ey, st, 'Color', ccol(2 + j), 'LineWidth', lw, 'DisplayName', sprintf('%s group (2 SD)%s', names{2 + j}, tern(cs.accepted, '', ', rejected')));
                        h(end+1) = hh; %#ok<AGROW>
                    end
                end
                if ~cs.accepted
                    h(end+1) = ln(R.crypto_min_PE_above_syn * [1 1], L(3:4), ccol(4), '-.', 'PE floor: dim left / bright right');
                end
            case 7 % chl relative to Euk line
                if isfield(R, 'euk_line')
                    h(end+1) = ln(xx, [0 0], ccol(1), '-', 'Euk chl-size line');
                    h(end+1) = ln(xx, -R.euksyn_max_chl_below_euk * [1 1], ccol(8), '--', 'Euk with Syn below this -> Syn on particle');
                    if isfield(info, 'catchall_euk_floor'), h(end+1) = ln(xx, info.catchall_euk_floor(3) * [1 1], ccol(1), ':', 'catch-all Euk floor'); end
                    sx = linspace(L(1), L(2), 50); base = R.euk_line(1) + R.euk_line(2) * sx;
                    if isfield(R, 'euk_chlSSC_cut'), h(end+1) = ln(sx, sx + R.euk_chlSSC_cut - base, ccol(1), '--', 'Euk chl:scatter cut'); end
                    h(end+1) = ln(sx, sx + R.detritus_chlSSC - base, [.4 .4 .4], ':', 'detritus cut');
                    h(end+1) = ln(R.euk_min_size * [1 1], L(3:4), ccol(1), '-.', 'Euk min size');
                end
        end
    end

    function ellipses(i1, i2)
        % 2 SD ellipses of the mixture components in two of the fit dimensions
        % [size, BL3, PE:chl, GL2:GL3], coloured by the component's name
        if ~isfield(info, 'gmm') || isempty(info.gmm), return, end
        g = info.gmm; lab = info.comp.label;
        for c = 1:g.NumComponents
            [ex, ey] = ellipse_xy(g.mu(c, [i1 i2]), g.Sigma([i1 i2], [i1 i2], c), 2);
            if lab(c) == 0, cc = [.45 .45 .45]; else, cc = ccol(lab(c)); end
            plot(ex, ey, '-', 'Color', [cc .55], 'LineWidth', .8, 'HandleVisibility', 'off')
        end
        plot(NaN, NaN, '-', 'Color', [.45 .45 .45], 'DisplayName', 'mixture components (2 SD)');
    end

    function s = keytext()
        sc = R.syn_check;
        s = {sprintf('Syn: size %.2f [%.2f %.2f], PE:chl %.2f [%.2f %.2f]', sc.size, sc.size_lo, sc.size_hi, sc.PEchl, sc.PEchl_lo, sc.PEchl_hi), ...
             sprintf('     GL2:GL3 [%.2f %.2f], PE floor %.2f', sc.GL2GL3_lo, sc.GL2GL3_hi, info.syn_PE_floor)};
        if isfield(R, 'euk_chlSSC_cut'), s{end+1} = sprintf('Euk: chl:scatter cut %.2f, min size %.2f', R.euk_chlSSC_cut, R.euk_min_size); end
        if isfield(info, 'pro_GL2_max'), s{end+1} = sprintf('Pro: PE ceiling %.2f', info.pro_GL2_max); end
        s{end+1} = sprintf('PE euk: size > %.2f, PE:chl < %.2f or crypto spectrum', R.peeuk_min_size, R.peeuk_PEchl_max);
        cs = info.crypto_split;
        if cs.n > 0
            if cs.accepted, m = sprintf('two groups accepted (sep %.1f)', cs.sep);
            elseif isfield(cs, 'single_group') && cs.single_group, m = 'one compact group, labelled by its PE';
            else, m = sprintf('diffuse: PE floor %.1f Syn eq', R.crypto_min_PE_above_syn); end
            s{end+1} = sprintf('Crypto split (%d crypto-type): %s', cs.n, m);
        end
        if isfield(R, 'synpart_cut'), s{end+1} = sprintf('Syn on particle: size > %.2f, chl %.1f below Syn trend', R.synpart_cut, R.synpart_below_trend); end
        s{end+1} = sprintf('Mixture: %d components; catch-all assigned %d', info.ncomp, info.n_catchall);
    end
end

function [x, y] = ellipse_xy(m, S, k)
t = linspace(0, 2 * pi, 80);
[V, D] = eig((S + S') / 2);
p = V * sqrt(max(D, 0)) * [cos(t); sin(t)] * k;
x = m(1) + p(1, :); y = m(2) + p(2, :);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
