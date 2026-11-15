function Q = culture_od2_fluorgate(cultures, k)
% function Q = culture_od2_fluorgate(cultures, k)
% Like culture_od2_check, but from the ungated size-calibration files with
% one fluorescence-only gate applied identically to the no-OD2 and OD2 runs
% (fluorescence voltages are the same in both), so the SSC distributions
% are not truncated by an SSC-based gate. Gate: chl (BL3) > 10^2.5, above
% the triggered noise; Syn also PE > 10^3.5.
%
% Q - one row per culture: n, SSC-H peak and 10/50/90th percentiles in both
%     runs, loss at each (OD2 - no OD2 - k log10(260/200))
%
% Nov 2026, PT transfer / OD2 low-signal check

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
if nargin < 1 || isempty(cultures), cultures = {'Syn' 'Micromonas' 'Nano' 'I_galbana' 'Pavlova'}; end
if nargin < 2 || isempty(k), k = 8.5; end
d0 = ['\' '\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\FCS\'];
dv = k * log10(260 / 200);
rows = {};
for i = 1:numel(cultures)
    fa = dir(fullfile(d0, ['Size_cal_nuetral_Mar2019_Group_day0_' cultures{i} '.fcs']));
    fb = dir(fullfile(d0, ['size_cal_OD2_Mar2019_Group_day0_' cultures{i} '.fcs']));
    if isempty(fa) || isempty(fb), continue, end
    [a, ha] = fca_readfcs(fullfile(d0, fa(1).name));
    [b, hb] = fca_readfcs(fullfile(d0, fb(1).name));
    L = @(x, h, n) log10(max(x(:, strcmp({h.par.name}, n)), 1));
    % cells: chl above the BL3 noise (triggered events sit at ~10^2.0-2.1);
    % Syn also PE (old-layout GL1-H, 575/36) > 10^3.5. Same gate on both runs.
    gate = @(x, h) L(x, h, 'BL3-H') > 2.5 & (~strcmp(cultures{i}, 'Syn') | L(x, h, 'GL1-H') > 3.5);
    ga = gate(a, ha); gb = gate(b, hb);
    sa = L(a, ha, 'SSC-H'); sb = L(b, hb, 'SSC-H');
    sa = sa(ga); sb = sb(gb);
    r = struct('pop', string(cultures{i}), 'n_neu', numel(sa), 'n_od2', numel(sb), ...
        'frac_od2_floor', mean(sb <= 1), 'raw_neu', bead_peak(sa(sa > 0)), 'raw_od2', bead_peak(sb(sb > 0)));
    r.loss_pk = r.raw_od2 - r.raw_neu - dv;
    for p = [5 10 25 50 75 90]
        qa = prctile(sa, p); qb = prctile(sb, p);
        r.(sprintf('raw_od2_p%d', p)) = qb;
        r.(sprintf('loss_p%d', p)) = qb - qa - dv;
    end
    rows{end+1, 1} = r; %#ok<AGROW>
end
Q = struct2table([rows{:}]');
end
