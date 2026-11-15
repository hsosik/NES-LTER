function M = mie_488_532_check()
% function M = mie_488_532_check()
% Can pigment absorption explain why large cells have more side scatter (in
% 1 um bead units) at 532 nm (GL1, 2026 layout) than at 488 nm (SSC, 2019
% calibration with OD2 and 2026 SSC), while small cells agree?
%   - beads: polystyrene, n 1.598 (488) / 1.595 (532); water 1.337 / 1.335
%   - cells: homogeneous spheres, real index 1.04 relative to water,
%     imaginary index from pigment absorption, k = a lambda / (4 pi) with
%     a = absorption coefficient of the cell material at that wavelength.
%     Three pigment types, scaled to k(440 nm) = 0.010 (0.004-0.012 tried) (k/1.337 relative):
%       green (chl a/b):        488/440 0.55, 532/440 0.15
%       carotenoid-rich (golden, Dunaliella): 0.85, 0.30
%       cryptophyte (PE):       0.55, 0.85
%     plus non-absorbing cells (wavelength/Mie effect only)
%   - side scatter over a cone around 90 deg, half-angle 30, 50 and 64 deg
%     (64 = NA ~1.2 in water)
% Checks the bead ratios 0.5/1 and 6/1 um against the measured ones, then
% compares log10(cell/bead) at 532 minus 488 with (a) the 2026 GL1 - SSC on
% the same cultures where SSC is linear, and (b) the gap between the 2026
% GL1 calibration (Attune_SC2vol 14 Nov) and the 2019 OD2-on-SSC calibration
% (sizecal_2019), both as a function of volume.
% Draws ../output/pt_transfer/mie_488_532.png.
%
% Nov 2026, classifier_dev/mie_scatter

od = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'output', 'pt_transfer');
lam = [0.488 0.532]; nw = [1.337 1.335]; nps = [1.598 1.595];
cones = [30 50 64];
d = logspace(log10(0.6), log10(25), 45);
k440 = 0.010; % 0.004-0.012 tried; ~0.008-0.012 matches the observations
types = {'non-absorbing' [0 0]; 'green (chl a/b)' [0.55 0.15]; 'carotenoid-rich' [0.85 0.30]; 'cryptophyte (PE)' [0.55 0.85]};
M = struct();
for ic = 1:numel(cones)
    ca = cones(ic);
    for il = 1:2
        bd = ssc_mie([0.5 1 6], nps(il) / nw(il), lam(il), nw(il), ca);
        M.bead(ic, il, :) = log10(bd / bd(2));
        for it = 1:size(types, 1)
            k = k440 * types{it, 2}(il) / nw(il);
            c = ssc_mie(d, 1.04 + 1i * k, lam(il), nw(il), ca);
            M.cell(ic, il, it, :) = log10(c / bd(2));
        end
    end
    fprintf('cone %2d deg: bead 0.5/1 um %.2f / %.2f, 6/1 um %.2f / %.2f (488 / 532)\n', ca, ...
        M.bead(ic, 1, 1), M.bead(ic, 2, 1), M.bead(ic, 1, 3), M.bead(ic, 2, 3));
end
fprintf('measured: SSC 0.5/1 -0.83 to -0.86, 6/1 0.89-0.98 (no OD2) / 1.26 (OD2); GL1 6/1 ~1.28\n')
M.d = d; M.vol = pi / 6 * d.^3; M.cones = cones; M.types = types(:, 1);
save(fullfile(od, 'mie_488_532.mat'), 'M')

% observations
lv = linspace(-0.3, 3.6, 50);
gap = (lv - 1.0826) / 1.1962 - (lv - 1.2109) / 1.4066; % 2026 GL1 cal - 2019 OD2-SSC cal, in log bead units
obs = {'Micromonas' 2.0 0.02; 'Nano' 10.3 0.03; 'Syn 7335' 3.5 0.00; 'I galbana' 47 0.09; 'Isochrysis' 50 0.15; 'Pavlova' 60 0.14};

f = figure('visible', 'off', 'position', [0 0 1500 520]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
col = lines(4);
for ic = 1:numel(cones)
    nexttile; hold on
    for it = 1:size(types, 1)
        dlt = squeeze(M.cell(ic, 2, it, :) - M.cell(ic, 1, it, :));
        plot(log10(M.vol), movmean(dlt, 3), '-', 'color', col(it, :), 'linewidth', 1.5, 'DisplayName', types{it, 1})
    end
    plot(lv, gap, 'k--', 'linewidth', 2, 'DisplayName', '2026 GL1 cal - 2019 OD2-SSC cal')
    plot(log10([obs{:, 2}]), [obs{:, 3}], 'ks', 'markerfacecolor', 'y', 'DisplayName', '2026 GL1 - SSC (same cells, SSC linear)')
    yline(0, 'k:', 'HandleVisibility', 'off')
    xlabel('log10 volume (\mum^3)'); ylabel('log10(cell/bead) at 532 - at 488 nm')
    title(sprintf('collection cone \\pm%d deg around 90 deg', cones(ic)))
    grid on; ylim([-0.4 0.8]); xlim([-0.3 3.6])
    if ic == 1, legend('location', 'northwest', 'fontsize', 7), end
end
title(tl, sprintf('Mie: side scatter in 1 \\mum bead units, 532 vs 488 nm; cells m = 1.04 + ik, k(440) = %.3f', k440))
exportgraphics(f, fullfile(od, 'mie_488_532.png'), 'Resolution', 110)
end
