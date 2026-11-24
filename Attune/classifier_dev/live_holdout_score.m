function H = live_holdout_score(sets, outdir)
% function H = live_holdout_score(sets, outdir)
% Score assign_class_v2 on the live held-out lists (live_holdout_lists) against
% the old classes (archived, via cruise_paths). One row per file: counts per class
% (old/v2), agreement as in score_agreement (agreement_pe_merged: PE euk 3/4/5 as
% one class), the v2 flag counts, and pass/fail.
%
% Pass marks, fixed before scoring (the discrete check's "within 20%"):
%   - agreement_pe_merged >= 0.85
%   - for Euk, Syn, PE euk and Pro, wherever old or v2 has >= 100 events:
%     v2/old within 20% (0.8 to 1.2)
% Set level (printed): median agreement and class ratios, share of files passing,
% compared with set T (the tuning files, same code). Over-tuning shows as L1 well
% below L2 and T.
%
% sets   - default {'T' 'L1' 'L2'} (R is the untouched reserve)
% outdir - default output/live_holdout
% Saves outdir/score_<set>.mat (H) and outdir/failures_<set>.csv (file, cruisedir,
% class, old, v2, ratio) for live_holdout_review.
%
% Nov 2026, classifier redesign (held-out check, live samples)

here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here))
if nargin < 1 || isempty(sets), sets = {'T' 'L1' 'L2'}; end
if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'output', 'live_holdout'); end
L = readtable(fullfile(outdir, 'holdout_lists.csv'), 'TextType', 'string');
[~, commit] = system(['git -C "' here '" rev-parse --short HEAD']); commit = strtrim(commit);
S = table;
for s = 1:numel(sets)
    Ls = L(L.set == sets{s}, :);
    rows = cell(height(Ls), 1);
    parfor i = 1:height(Ls)
        rows{i} = score_one(Ls(i, :));
    end
    H = vertcat(rows{:});
    H.commit = repmat(string(commit), height(H), 1);
    save(fullfile(outdir, ['score_' sets{s} '.mat']), 'H')
    % failures, one row per failed criterion
    Fl = table;
    for i = 1:height(H)
        if ~(H.agreement_pe_merged(i) >= 0.85)
            Fl = [Fl; table(H.file(i), H.cruisedir(i), "agreement", NaN, NaN, H.agreement_pe_merged(i), 'VariableNames', {'file' 'cruisedir' 'class' 'old' 'v2' 'ratio'})]; %#ok<AGROW>
        end
        for c = ["euk" "syn" "peeuk" "pro"]
            o = H.("old_" + c)(i); n = H.("new_" + c)(i);
            if max(o, n) >= 100 && ~(n / o >= 0.8 && n / o <= 1.2)
                Fl = [Fl; table(H.file(i), H.cruisedir(i), c, o, n, n / o, 'VariableNames', {'file' 'cruisedir' 'class' 'old' 'v2' 'ratio'})]; %#ok<AGROW>
            end
        end
    end
    if isempty(Fl), Fl = table(strings(0, 1), strings(0, 1), strings(0, 1), [], [], [], 'VariableNames', {'file' 'cruisedir' 'class' 'old' 'v2' 'ratio'}); end
    writetable(Fl, fullfile(outdir, ['failures_' sets{s} '.csv']));
    H.pass = ~ismember(H.file, Fl.file);
    save(fullfile(outdir, ['score_' sets{s} '.mat']), 'H')
    S = [S; summarize(H, sets{s})]; %#ok<AGROW>
end
disp(S)
writetable(S, fullfile(outdir, 'summary.csv'));
end

function r = score_one(Li)
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cd0 = char(Li.cruisedir); fn = char(Li.file);
P = cruise_paths(cd0, Root=root);
cfg = classifier_config(); cfg.cruise = extractAfter(cd0, '_');
cfg.gl1_is_od2_scatter = strcmpi(attune_config(cfg.cruise).OD2, 'GL1');
[dat, hdr] = fca_readfcs(fullfile(P.basepath, 'FCS', fn));
ref = bead_reference(load_beadstat(P.bead_dir), Li.time, cfg);
[new, info] = assign_class_v2(dat, hdr, ref, cfg, Li.pro_measured);
old = load_old_class(P, fn); old(~ismember(old, [0 1 2 3 4 7])) = 0;
o34 = old; o34(ismember(o34, [4 5])) = 3; n34 = new; n34(ismember(n34, [4 5])) = 3;
either = old ~= 0 | new ~= 0;
r = table(Li.set, Li.cruisedir, Li.file, Li.time, Li.pro_measured, numel(new), sum(either), ...
    sum(old == new & old ~= 0) / max(1, sum(either)), sum(o34 == n34 & o34 ~= 0) / max(1, sum(either)), ...
    sum(old == 1), sum(new == 1), sum(old == 2), sum(new == 2), sum(ismember(old, [3 4])), sum(ismember(new, [3 4 5])), ...
    sum(old == 7), sum(new == 7), sum(new == 8), sum(old == 0 & new ~= 0), sum(old ~= 0 & new == 0), ...
    sum(info.flag == 1), sum(info.flag == 2), sum(info.flag == 3), ...
    'VariableNames', {'set' 'cruisedir' 'file' 'time' 'pro_measured' 'n' 'n_phyto' 'agreement' 'agreement_pe_merged' ...
    'old_euk' 'new_euk' 'old_syn' 'new_syn' 'old_peeuk' 'new_peeuk' 'old_pro' 'new_pro' 'new_synpart' 'n_added' 'n_dropped' ...
    'flag_synpig' 'flag_below_gap' 'flag_edge_pile'});
end

function s = summarize(H, name)
rat = @(o, n) median(n(max(o, n) >= 100) ./ o(max(o, n) >= 100), 'omitnan');
w20 = @(o, n) mean(abs(n(max(o, n) >= 100) ./ o(max(o, n) >= 100) - 1) <= 0.2);
s = table(string(name), height(H), mean(H.pass), median(H.agreement_pe_merged), min(H.agreement_pe_merged), ...
    rat(H.old_euk, H.new_euk), rat(H.old_syn, H.new_syn), rat(H.old_peeuk, H.new_peeuk), rat(H.old_pro, H.new_pro), ...
    w20(H.old_euk, H.new_euk), w20(H.old_syn, H.new_syn), w20(H.old_peeuk, H.new_peeuk), w20(H.old_pro, H.new_pro), ...
    'VariableNames', {'set' 'files' 'frac_pass' 'agree_median' 'agree_min' 'euk_ratio' 'syn_ratio' 'peeuk_ratio' 'pro_ratio' ...
    'euk_within20' 'syn_within20' 'peeuk_within20' 'pro_within20'});
end
