B = '\\sosiknas1\lab_data\Attune\cruise_data\20250814_AR92\bead_calibrated\';
s = load([B 'beadstat_2026.mat']); f = fieldnames(s); disp(f)
for i = 1:numel(f)
    x = s.(f{i}); disp(f{i}); disp(class(x)); disp(size(x))
    if istable(x), disp(x.Properties.VariableNames'); disp(x(1:2,:)); end
end
s2 = load([B 'FCSfileinfo.mat']); disp(s2.FCSfileinfo.Properties.VariableNames')
disp(s2.FCSfileinfo(1:3,:))
s3 = load([B 'AttuneTable_uw_match.mat']); f3 = fieldnames(s3); disp(f3)
t = s3.(f3{1}); if istable(t), disp(t.Properties.VariableNames'); end
d = dir([B 'class\*.mat']); disp(numel(d))
c = load(fullfile(d(500).folder, d(500).name)); disp(d(500).name); disp(fieldnames(c)); disp(c.notes); disp(size(c.class))
tabulate(c.class(:,1))
