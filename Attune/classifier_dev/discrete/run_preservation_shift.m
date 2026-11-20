% run preservation_shift on cruises with bottle metadata in their Gated_Table
% (summer, winter, spring, summer); Nov 2026
s = load(fullfile(fileparts(mfilename('fullpath')), '..', 'output', 'discrete', 'lab_beadstat.mat'));
fprintf('lab bead runs: %d, %s to %s\n', height(s.beadstat), string(min(s.beadstat.time)), string(max(s.beadstat.time)))
if ~exist('cruises', 'var')
    cruises = {'20240906_EN720' '20230111_EN695' '20230807_EN706' '20240503_EN715' '20250814_AR92' '20250124_EN727'};
end
for c = cruises
    try
        S = preservation_shift(c{1});
        fprintf('%s: %d rows\n', c{1}, height(S))
    catch err
        fprintf('%s failed: %s\n', c{1}, err.message)
        disp(err.stack(1))
    end
end
