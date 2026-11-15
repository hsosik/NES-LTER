B = '\\sosiknas1\lab_data\Attune\cruise_data\20220806_EN688\bead_calibrated\';
s = load([B 'beadstat_2021.mat']); f = fieldnames(s); disp(f)
x = s.beadstat_2021; disp(class(x)); disp(size(x))
if istable(x), disp(x.Properties.VariableNames'); end
s = load([B 'beadstat_2026.mat']); x = s.beadstat_2026;
m = x.median2; disp(class(m)); disp(size(m));
if istable(m), disp(m.Properties.VariableNames(1:20)); disp(m(:, 1:20)); end
disp(size(x.hv)); disp(x.parname(1,1:20)); disp(x.hv(1,1:19))
% lat/lon names in uw_match across cruises
D = dir('\\sosiknas1\lab_data\Attune\cruise_data\20*_*');
for i = 1:numel(D)
    fn = fullfile(D(i).folder, D(i).name, 'bead_calibrated', 'AttuneTable_uw_match.mat');
    if exist(fn, 'file')
        t = load(fn); t = t.(char(fieldnames(t)));
        v = t.Properties.VariableNames;
        disp([D(i).name ': ' strjoin(v(contains(v, {'lat','lon','sbe','temp','tsg','sal'}, 'IgnoreCase', true)), ' ')])
    end
end
