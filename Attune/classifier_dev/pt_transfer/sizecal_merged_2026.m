function R = sizecal_merged_2026(win)
% function R = sizecal_merged_2026(win)
% SSC-to-volume calibration on the merged side-scatter scale in 1 um bead
% units, from the July/Aug 2026 cultures (standard configuration).
% Per event, with s = SSC-H / 1 um bead SSC-H and g = GL1-H / 1 um bead GL1-H
% (same day's FCB bead runs):
%   merged = s           where g <= win(1)
%          = (s + g) / 2 where win(1) < g <= win(2)
%          = g           where g > win(2)
% (bead-unit window; default win = [0.5 2]). For comparison the same with the
% raw-count window used so far (GL1-H 500-1000 counts).
% Culture value = KDE mode of log10 merged (bead_peak); volume = Coulter mode
% (CCSummary_saved_10Sep2026.mat). Fits log10 volume = p1 log10 merged + p2,
% with all cultures and without T. weissflogii (diatom: low scatter per
% volume), and shows the current Attune_SC2vol SSC-H curve.
% Draws output/pt_transfer/sizecal_merged_2026.png.
%
% Nov 2026, PT transfer / volume calibration

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
if nargin < 1 || isempty(win), win = [0.5 2]; end
p26 = ['\' '\sosiknas1\Lab_data\Attune\size_calibration_July2026\'];
gdir = fullfile(p26, 'Attune', 'Gated_FCS');
s = load(fullfile(p26, 'AttuneSummary_saved_11Sep2026.mat'), 'AttuneMode'); A = s.AttuneMode;
A.label = string(A.label);
s = load(fullfile(p26, 'CCSummary_saved_10Sep2026.mat'), 'CCMode'); CC = s.CCMode;
rows = {};
for i = 1:height(A)
    lab = A.label(i); stem = A.Properties.RowNames{i};
    if ~ismember(lab, CC.Properties.RowNames) || contains(stem, 'pooled'), continue, end
    fn = fullfile(gdir, [stem '.fcs']);
    if ~isfile(fn), fprintf('missing %s\n', stem); continue, end
    b = A(A.run_date == A.run_date(i) & A.label == "1_um", :);
    bS = median(b.("SSC-H")); bG = median(b.("GL1-H"));
    [d, h] = fca_readfcs(fn);
    n = {h.par.name};
    sh = d(:, strcmp(n, 'SSC-H')) / bS; gh = d(:, strcmp(n, 'GL1-H')) / bG;
    graw = d(:, strcmp(n, 'GL1-H'));
    m = merge(sh, gh, gh, win);
    mraw = merge(sh, gh, graw, [500 1000]);
    pk = @(x) bead_peak(log10(x(x > 0)));
    rows{end+1, 1} = struct('culture', lab, 'file', string(stem), 'n', numel(sh), ...
        'vol', pi / 6 * CC{lab, 'diameter_mode'}^3, 'ssc', pk(sh), 'gl1', pk(gh), ...
        'merged', pk(m), 'merged_raw', pk(mraw), 'frac_gl1', mean(gh > win(2)), ...
        'frac_blend', mean(gh > win(1) & gh <= win(2))); %#ok<AGROW>
end
R = struct2table([rows{:}]');
R = groupsummary(R, 'culture', 'median', {'vol', 'ssc', 'gl1', 'merged', 'merged_raw', 'frac_gl1', 'frac_blend'});
R.Properties.VariableNames = erase(R.Properties.VariableNames, 'median_');
R.lvol = log10(R.vol);
disp(R(:, {'culture', 'GroupCount', 'vol', 'ssc', 'gl1', 'merged', 'merged_raw', 'frac_gl1', 'frac_blend'}))

tw = R.culture == "T_weissflogii";
f1 = polyfit(R.merged, R.lvol, 1);
f2 = polyfit(R.merged(~tw), R.lvol(~tw), 1);
f3 = polyfit(R.merged_raw(~tw), R.lvol(~tw), 1);
rms = @(f, x, y) sqrt(mean((polyval(f, x) - y).^2));
fprintf('bead-unit window [%.2g %.2g]: all      p1 %.4f p2 %.4f  rms %.3f\n', win, f1, rms(f1, R.merged, R.lvol))
fprintf('bead-unit window: without T. weissflogii p1 %.4f p2 %.4f  rms %.3f\n', f2, rms(f2, R.merged(~tw), R.lvol(~tw)))
fprintf('raw window 500-1000: without T. weiss.   p1 %.4f p2 %.4f  rms %.3f\n', f3, rms(f3, R.merged_raw(~tw), R.lvol(~tw)))
fprintf('Attune_SC2vol SSC-H (11 Sep 2026):        p1 1.1696 p2 1.1592  rms on bead-unit merged %.3f\n', ...
    rms([1.1696 1.1592], R.merged(~tw), R.lvol(~tw)))
R.Properties.UserData = struct('win', win, 'fit_all', f1, 'fit_noTW', f2, 'fit_rawwin_noTW', f3);

f = figure('visible', 'off', 'position', [0 0 900 650]); hold on
x = linspace(-1.4, 2, 50);
plot(R.merged, R.lvol, 'ko', 'markerfacecolor', 'k', 'DisplayName', 'merged, bead-unit window')
plot(R.merged_raw, R.lvol, 'r+', 'DisplayName', 'merged, raw window 500-1000')
plot(R.ssc, R.lvol, 'b.', 'markersize', 12, 'DisplayName', 'SSC only')
for i = 1:height(R), text(R.merged(i) + 0.04, R.lvol(i), strrep(R.culture(i), '_', ' '), 'fontsize', 8), end
plot(x, polyval(f2, x), 'k-', 'DisplayName', sprintf('fit (no T. weissflogii): %.3f x + %.3f', f2))
plot(x, polyval(f1, x), 'k:', 'DisplayName', sprintf('fit (all): %.3f x + %.3f', f1))
plot(x, 1.1696 * x + 1.1592, 'g--', 'DisplayName', 'Attune\_SC2vol SSC-H (11 Sep 2026)')
xline(log10(win), ':', 'HandleVisibility', 'off')
xlabel('merged side scatter, 1 \mum bead units (log10)'); ylabel('Coulter mode volume (log10 \mum^3)')
legend('location', 'northwest', 'fontsize', 8); grid on
title(sprintf('July/Aug 2026 cultures: merged SSC/GL1 (window %.2g-%.2g x bead)', win))
exportgraphics(f, fullfile(dev, 'output', 'pt_transfer', 'sizecal_merged_2026.png'), 'Resolution', 110)
end

function m = merge(s, g, gsel, w)
% s, g in bead units; gsel = the quantity the window applies to
m = s;
hi = gsel > w(2); mid = gsel > w(1) & ~hi;
m(hi) = g(hi);
m(mid) = (s(mid) + g(mid)) / 2;
end
