function test_process_cruise(cruise, outdir, nfiles, steps)
% function test_process_cruise(cruise, outdir, nfiles, steps)
% Run process_cruise on a few files of a cruise into a scratch folder, so
% nothing in the cruise folder on the NAS is written. Copies FCSfileinfo.mat
% and beadstat.mat from the cruise's bead_calibrated folder first, then prints
% the class log, the new AttuneTable columns and one class file.
%
% e.g. test_process_cruise('EN727', 'C:\scratch\EN727_test', 8)
%
% Nov 2026, processing wrapper redesign (development/testing script)

if nargin < 3 || isempty(nfiles), nfiles = 8; end
if nargin < 4 || isempty(steps), steps = [3 4 5 6]; end
addpath(fileparts(fileparts(mfilename('fullpath'))))
s = attune_cruise_settings(cruise);
if ~exist(outdir, 'dir'), mkdir(outdir), end
bc = fullfile(s.basepath, 'bead_calibrated');
copyfile(fullfile(bc, 'FCSfileinfo.mat'), outdir)
copyfile(fullfile(bc, 'beadstat.mat'), outdir)

process_cruise(cruise, steps, struct('outpath', outdir, 'files', nfiles, 'parallel', nfiles > 4));

if ismember(3, steps)
    L = load(fullfile(outdir, 'class_log.mat')); disp(L.L)
end
if ismember(5, steps)
    T = load(fullfile(outdir, 'AttuneTable.mat'));
    A = T.AttuneTable;
    disp(rows2vars(A(:, {'Syn_count' 'Crypto_dim_count' 'Crypto_bright_count' 'Euk_with_Syn_count' 'Syn_on_particle_count' ...
        'Syn_typeA_count' 'Syn_typeB_count' 'Syn_typeC_count' 'Flag_lowchl_Syn_pigment_count' 'Flag_below_Euk_chl_gap_count' 'Flag_Pro_tail_count'})))
    disp(T.table_metadata)
end
d = dir(fullfile(outdir, 'class', '*.mat'));
c = load(fullfile(d(1).folder, d(1).name));
disp(c); disp(c.notes); disp(c.class_meta)
end
