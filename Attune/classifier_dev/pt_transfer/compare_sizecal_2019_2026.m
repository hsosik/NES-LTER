function T = compare_sizecal_2019_2026()
% function T = compare_sizecal_2019_2026()
% Culture size calibrations: March 2019 (old layout; SSC without OD2 at
% 200 V and with OD2 at 260 V) vs July/Aug 2026 (current layout: SSC
% without OD2 at 220 V, GL1 = 532 nm side scatter through OD2 at 220 V).
% All signals normalized to the 1 um FCB bead on the same channel, same day.
%   2019: ModeTable in FitResults_saved_SSC-H_01Apr2026.mat (Coulter mode
%         volume; noOD / OD2 / merged bead-normalized modes)
%   2026: AttuneMode in AttuneSummary_saved_11Sep2026.mat (SSC-H, GL1-H,
%         SSC-A modes of gated files) and Coulter mode diameters from
%         CCSummary_saved_10Sep2026.mat (CCMode)
% Draws output/pt_transfer/sizecal_2019_2026.png.
%
% Nov 2026, PT transfer / OD2 low-signal check

od = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'output', 'pt_transfer');
p19 = ['\' '\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\FitResults_saved_SSC-H_01Apr2026.mat'];
p26 = ['\' '\sosiknas1\Lab_data\Attune\size_calibration_July2026\'];
s = load(p19, 'ModeTable'); M = s.ModeTable;
s = load(fullfile(p26, 'AttuneSummary_saved_11Sep2026.mat'), 'AttuneMode'); A = s.AttuneMode;
% Coulter mode diameters (um), July/Aug 2026 (CCSummary_saved_10Sep2026.mat, CCMode),
% as the 2019 calibration uses Coulter modes
s = load(fullfile(p26, 'CCSummary_saved_10Sep2026.mat'), 'CCMode');
cc = struct();
for i = 1:height(s.CCMode), cc.(s.CCMode.Properties.RowNames{i}) = s.CCMode.diameter_mode(i); end
A.label = string(A.label);
rows = {};
for d = unique(A.run_date)'
    b = A(A.run_date == d & A.label == "1_um", :);
    if isempty(b), continue, end
    bS = median(b.("SSC-H")); bG = median(b.("GL1-H")); bA = median(b.("SSC-A"));
    a = A(A.run_date == d, :);
    for i = 1:height(a)
        lab = a.label(i);
        if ~isfield(cc, lab), continue, end
        if contains(a.Properties.RowNames{i}, 'pooled'), continue, end
        dia = cc.(lab);
        rows{end+1, 1} = struct('year', 2026, 'culture', lab, 'vol', pi / 6 * dia^3, ...
            'ssc', a.("SSC-H")(i) / bS, 'gl1', a.("GL1-H")(i) / bG, 'sscA', a.("SSC-A")(i) / bA, ...
            'raw_ssc', a.("SSC-H")(i), 'raw_gl1', a.("GL1-H")(i)); %#ok<AGROW>
    end
end
T = struct2table([rows{:}]');
T = groupsummary(T, {'year', 'culture'}, 'median', {'vol', 'ssc', 'gl1', 'sscA', 'raw_ssc', 'raw_gl1'});
T.Properties.VariableNames = erase(T.Properties.VariableNames, 'median_');
T = T(:, {'year', 'culture', 'vol', 'ssc', 'gl1', 'sscA', 'raw_ssc', 'raw_gl1'});
% 2019
M19 = table(repmat(2019, height(M), 1), string(M.files), M.cc_mode_vol, M.noOD_mode_bd_norm, M.OD2_mode_bd_norm, ...
    nan(height(M), 1), M.noOD_mode, M.OD2_mode, 'VariableNames', T.Properties.VariableNames);
T = [T; M19];
disp(T)

f = figure('visible', 'off', 'position', [0 0 1300 560]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on
y19 = T.year == 2019; y26 = T.year == 2026;
plot(log10(T.ssc(y19)), log10(T.vol(y19)), 'bo', 'markerfacecolor', 'b', 'DisplayName', '2019 SSC, no OD2 (200 V)')
plot(log10(T.gl1(y19)), log10(T.vol(y19)), 'bs', 'DisplayName', '2019 SSC with OD2 (260 V)')
plot(log10(T.ssc(y26)), log10(T.vol(y26)), 'ro', 'markerfacecolor', 'r', 'DisplayName', '2026 SSC, no OD2 (220 V)')
plot(log10(T.gl1(y26)), log10(T.vol(y26)), 'rs', 'DisplayName', '2026 GL1 (532 nm, OD2, 220 V)')
for i = find(y26)', text(log10(T.gl1(i)) + 0.03, log10(T.vol(i)), T.culture(i), 'fontsize', 7, 'color', 'r'), end
for i = find(y19)', text(log10(T.gl1(i)) + 0.03, log10(T.vol(i)), T.culture(i), 'fontsize', 7, 'color', 'b'), end
xlabel('scatter / 1 \mum bead on the same channel (log10)'); ylabel('Coulter volume (log10 \mum^3)')
legend('location', 'northwest'); grid on
title('Culture volume vs bead-normalized side scatter')
nexttile; hold on
plot(log10(T.raw_ssc(y26)), log10(T.gl1(y26) ./ T.ssc(y26)), 'ro', 'markerfacecolor', 'r', 'DisplayName', '2026: GL1_{norm} - SSC_{norm}')
plot(log10(T.raw_ssc(y19)), log10(T.gl1(y19) ./ T.ssc(y19)), 'bo', 'markerfacecolor', 'b', 'DisplayName', '2019: OD2_{norm} - noOD_{norm}')
for i = find(y26)', text(log10(T.raw_ssc(i)) + 0.02, log10(T.gl1(i) / T.ssc(i)), T.culture(i), 'fontsize', 7, 'color', 'r'), end
for i = find(y19)', text(log10(T.raw_ssc(i)) + 0.02, log10(T.gl1(i) / T.ssc(i)), T.culture(i), 'fontsize', 7, 'color', 'b'), end
yline(0, 'k:', 'HandleVisibility', 'off')
xlabel('raw SSC-H without OD2 (log10)'); ylabel('OD2 channel - no-OD2 channel, bead-normalized (log10)')
legend('location', 'northwest'); grid on
title('Where the two scatter measurements disagree')
exportgraphics(f, fullfile(od, 'sizecal_2019_2026.png'), 'Resolution', 110)
end
