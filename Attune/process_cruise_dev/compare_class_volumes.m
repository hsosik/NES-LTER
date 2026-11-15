function [R, B] = compare_class_volumes(cruise, newclassdir, nfiles, c)
% function [R, B] = compare_class_volumes(cruise, newclassdir, nfiles, c)
% Separate the size-calibration change from the class-membership change
% between the existing class files (bead_calibrated\class, variable volume)
% and new ones from process_cruise (volume_cubic_microns), for one class.
%
% cruise      - e.g. 'EN727'
% newclassdir - class folder written by process_cruise (e.g. a scratch run
%               from compare_attune_tables)
% nfiles      - files spread through the cruise (default 40)
% c           - class (default 1, Euk)
%
% R - per-file median of: biovolume (<= 20 um ESD) of events in class c in both
%     versions, new / old (calibration only); the table biovolume new / old;
%     share of old biovolume in events only the old classifier put in c; share
%     of new biovolume in events only the new one put in c
% B - median new / old volume of the same events by old ESD bin
%
% EN727, 15 Nov 2026 (Euk): same events 0.88, table 0.91, old-only 4%, new-only
% 8%; new/old volume 1.18 below 1 um falling to 0.81 at 10-20 um: the drop in
% Euk biovolume comes from the flatter new SSC-to-volume curve, not the classifier.
%
% Nov 2026, processing wrapper redesign (development/testing script)

if nargin < 3 || isempty(nfiles), nfiles = 40; end
if nargin < 4 || isempty(c), c = 1; end
addpath(fileparts(fileparts(mfilename('fullpath'))))
s = attune_cruise_settings(cruise);
olddir = fullfile(s.basepath, 'bead_calibrated', 'class');
d = dir(fullfile(newclassdir, '*.mat'));
k = unique(round(linspace(1, numel(d), nfiles)));
esd = @(v) (v*3/4/pi).^(1/3)*2;
R = []; E = [];
for i = k
    n = load(fullfile(d(i).folder, d(i).name));
    of = fullfile(olddir, d(i).name);
    if ~exist(of, 'file'), continue, end
    o = load(of);
    if numel(o.class) ~= numel(n.class) || ~isfield(o, 'volume') || ~isfield(n, 'volume_cubic_microns'), continue, end
    vn = n.volume_cubic_microns; vo = real(o.volume);
    dn = esd(vn) <= 20; do = esd(vo) <= 20;
    both = o.class == c & n.class == c;
    R = [R; sum(vn(both & dn))/sum(vo(both & do)), sum(vn(n.class == c & dn))/sum(vo(o.class == c & do)), ...
        sum(vo(o.class == c & n.class ~= c & do))/sum(vo(o.class == c & do)), ...
        sum(vn(n.class == c & o.class ~= c & dn))/sum(vn(n.class == c & dn))]; %#ok<AGROW>
    E = [E; vo(both) vn(both)]; %#ok<AGROW>
end
R = array2table(median(R, 'omitnan'), 'VariableNames', {'same_events_biovol_ratio' 'table_biovol_ratio' 'old_only_share_of_old' 'new_only_share_of_new'});
disp(R)
de = esd(E(:, 1)); ed = [0 1 2 3 5 10 20 100];
B = table;
for j = 1:numel(ed)-1
    q = de > ed(j) & de <= ed(j+1);
    B = [B; table(ed(j), ed(j+1), sum(q), median(E(q, 2)./E(q, 1), 'omitnan'), ...
        'VariableNames', {'esd_min' 'esd_max' 'n' 'median_new_over_old'})]; %#ok<AGROW>
end
disp(B)
end
