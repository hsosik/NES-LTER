function T = culture_od2_check(k, outfile)
% function T = culture_od2_check(k, outfile)
% Response of SSC with the OD2 filter vs without, from the March/April 2019
% size-calibration cultures: the same culture run minutes apart without OD2
% (SSC 200 V) and with OD2 on SSC (SSC 260 V), all other voltages equal.
% Gated files from Gated_FCS_March2026 (noOD_* / OD2_*).
%
% For each culture (and FCB bead population):
%   loss = log10 SSC-H(OD2) - log10 SSC-H(no OD2) - k log10(260/200)
% i.e. the apparent OD2 attenuation after removing the voltage change. A
% linear system gives the same loss for every particle. Loss smaller in magnitude at
% low raw OD2 signal = low-signal excess under OD2. Loss smaller in magnitude for
% bright particles = compression of SSC without OD2 at the high end.
% Also a quantile-quantile mapping per culture (5-95%), which follows the
% loss across each culture's own range of signal.
%
% k       - SSC voltage exponent between 200 and 260 V (default 8.5)
% T       - struct: pk (one row per population), qq (quantile pairs)
%
% Nov 2026, PT transfer / OD2 low-signal check

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
if nargin < 1 || isempty(k), k = 8.5; end
base = ['\' '\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\Gated_FCS_March2026\'];
pairs = {'noOD_20Mar2019' 'OD2_20Mar2019'; 'noOD_9Apr2019' 'OD2_9Apr2019'};
dv = k * log10(260 / 200);
pk = {}; qq = {};
for ip = 1:size(pairs, 1)
    F = dir(fullfile(base, pairs{ip, 1}, '*.fcs'));
    for i = 1:numel(F)
        nm = erase(F(i).name, '.fcs');
        stem = regexprep(nm, '_\d(_R\d)?$', ''); % 1820i_2_R2 -> 1820i (not beads_1um)
        G = dir(fullfile(base, pairs{ip, 2}, [stem '*.fcs']));
        if isempty(G), continue, end
        [a, ha] = fca_readfcs(fullfile(F(i).folder, F(i).name));
        [b, hb] = fca_readfcs(fullfile(G(1).folder, G(1).name));
        ca = strcmp({ha.par.name}, 'SSC-H'); cb = strcmp({hb.par.name}, 'SSC-H');
        xa = log10(max(a(:, ca), 1)); xb = log10(max(b(:, cb), 1));
        xa = xa(xa > 0.5); xb = xb(xb > 0.5);
        if numel(xa) < 100 || numel(xb) < 100, continue, end
        pa = bead_peak(xa); pb = bead_peak(xb);
        aa = strcmp({ha.par.name}, 'SSC-A'); ab = strcmp({hb.par.name}, 'SSC-A');
        ya = log10(max(a(:, aa), 1)); yb = log10(max(b(:, ab), 1));
        lossA = bead_peak(yb(yb > 0.5)) - bead_peak(ya(ya > 0.5)) - dv;
        pk{end+1, 1} = struct('set', string(pairs{ip, 2}), 'pop', string(stem), ...
            'n_neu', numel(xa), 'n_od2', numel(xb), 'raw_neu', pa, 'raw_od2', pb, ...
            'loss', pb - pa - dv, 'lossA', lossA, 'hv_neu', ha.par(ca).hv, 'hv_od2', hb.par(cb).hv); %#ok<AGROW>
        q = (5:5:95)';
        qa = prctile(xa, q); qb = prctile(xb, q);
        for j = 1:numel(q)
            qq{end+1, 1} = struct('pop', string(stem), 'q', q(j), 'raw_neu', qa(j), ...
                'raw_od2', qb(j), 'loss', qb(j) - qa(j) - dv); %#ok<AGROW>
        end
    end
end
T.pk = struct2table([pk{:}]');
T.qq = struct2table([qq{:}]');
T.k = k;
% plateau: populations with raw OD2 peak between 10^2.3 and 10^3.8
mid = T.pk.raw_od2 > 2.3 & T.pk.raw_od2 < 3.8;
T.plateau = median(T.pk.loss(mid));
T.pk.excess = T.pk.loss - T.plateau;
T.qq.excess = T.qq.loss - T.plateau;
if nargin > 1 && ~isempty(outfile), save(outfile, 'T'), end
end
