function R = compare_attune_tables(cruise, outdir, run)
% function R = compare_attune_tables(cruise, outdir, run)
% Process a whole cruise with process_cruise into a scratch folder (classify to metadata,
% parallel) and compare the new AttuneTable with the existing one in the
% cruise's bead_calibrated folder: new/old ratios of counts, biovolume and
% carbon per file, median per-ml values of the new classes, the files with the
% largest differences, and a time-series figure. Nothing in the cruise folder
% is written.
%
% cruise - e.g. 'EN727'
% outdir - scratch folder for the new class files and tables
% run    - true (default): run process_cruise first; false: compare only
%
% Writes <outdir>\compare_<cruise>.png and compare_<cruise>.mat.
% EN727, 15 Nov 2026: 2230 files in 27 min on 26 workers, no failures; Syn and
% Euk counts +2%, PE euks -15%, Euk biovolume -17% (size calibration: see
% compare_class_volumes).
%
% Nov 2026, processing wrapper redesign (development/testing script)

if nargin < 3, run = true; end
addpath(fileparts(fileparts(mfilename('fullpath'))))
s = attune_cruise_settings(cruise);
bc = fullfile(s.basepath, 'bead_calibrated');
if run
    if ~exist(outdir, 'dir'), mkdir(outdir), end
    copyfile(fullfile(bc, 'FCSfileinfo.mat'), outdir)
    copyfile(fullfile(bc, 'beadstat.mat'), outdir)
    tic
    process_cruise(cruise, Steps=["classify" "volume" "table" "metadata"], OutPath=outdir, Parallel=true);
    fprintf('process_cruise: %.1f min\n', toc/60)
end

O = load(fullfile(bc, 'AttuneTable.mat')); O = O.AttuneTable;
N = load(fullfile(outdir, 'AttuneTable.mat')); N = N.AttuneTable;
[~, io, in] = intersect(O.Filename, N.Filename);
O = O(io, :); N = N(in, :);
fprintf('files: %d in common\n', numel(io))
v = {'Syn_count' 'Euk_without_PE_leq2um_count' 'Euk_without_PE_leq3um_count' 'Euk_without_PE_leq5um_count' ...
    'Euk_without_PE_leq10um_count' 'Euk_without_PE_leq20um_count' 'Euk_w_PE_leq20um_count' 'Pro_count' ...
    'Syn_biovolume' 'Euk_without_PE_leq20um_biovolume' 'Euk_w_PE_leq20um_biovolume' 'Syn_carbon' 'Euk_without_PE_leq20um_carbon'};
R = table;
for i = 1:numel(v)
    o = O.(v{i}); n = N.(v{i});
    ok = isfinite(o) & isfinite(n) & o >= 20;
    r = n(ok) ./ o(ok);
    if sum(ok) > 2, lr = corrcoef(log10(o(ok)), log10(max(n(ok), 1))); lr = lr(1, 2); else, lr = NaN; end
    R = [R; table(string(v{i}), sum(ok), sum(n(ok))/sum(o(ok)), median(r), prctile(r, 10), prctile(r, 90), lr, ...
        'VariableNames', {'variable' 'n_files' 'sum_ratio' 'median_ratio' 'p10' 'p90' 'r_log'})]; %#ok<AGROW>
end
disp(R)

vol = N.VolAnalyzed_ml;
newv = {'Crypto_dim_count' 'Crypto_bright_count' 'Euk_with_Syn_count' 'Syn_on_particle_count' 'Syn_typeA_count' ...
    'Syn_typeB_count' 'Syn_typeC_count' 'Syn_type_unusable_count' 'Flag_lowchl_Syn_pigment_count' ...
    'Flag_below_Euk_chl_gap_count' 'Flag_Pro_tail_count'};
for i = 1:numel(newv)
    fprintf('%-32s median per ml %8.1f\n', newv{i}, median(N.(newv{i})./vol, 'omitnan'))
end
for vv = {'Syn_count' 'Euk_without_PE_leq20um_count' 'Euk_w_PE_leq20um_count'}
    o = O.(vv{1}); n = N.(vv{1}); d = abs(log2((n+10)./(o+10)));
    [~, k] = sort(d, 'descend'); k = k(1:min(5, end));
    fprintf('\n%s largest differences:\n', vv{1})
    disp(table(N.Filename(k), N.StartDate(k), o(k), n(k), 'VariableNames', {'file' 'time' 'old' 'new'}))
end

f = figure('Visible', 'off', 'Position', [50 50 1400 1000]);
pv = {'Syn_count' 'Euk_without_PE_leq20um_count' 'Euk_w_PE_leq20um_count' 'Euk_without_PE_leq2um_count' ...
    'Syn_biovolume' 'Euk_without_PE_leq20um_biovolume'};
for i = 1:numel(pv)
    subplot(3, 2, i)
    semilogy(O.StartDate, O.(pv{i})./O.VolAnalyzed_ml, '.', N.StartDate, N.(pv{i})./vol, '.')
    title(strrep(pv{i}, '_', ' ')); ylabel('per ml')
    if i == 1, legend('existing table', 'process\_cruise', 'Location', 'best'); end
end
exportgraphics(f, fullfile(outdir, ['compare_' cruise '.png']))
save(fullfile(outdir, ['compare_' cruise '.mat']), 'R', 'O', 'N')
end
