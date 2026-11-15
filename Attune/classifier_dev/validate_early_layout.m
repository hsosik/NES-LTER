function V = validate_early_layout(cruisedir, nfiles, outdir, classdir, compare_orig)
% function V = validate_early_layout(cruisedir, nfiles, outdir, classdir, compare_orig)
% score_agreement for cruises before the Aug 2019 filter change, with the
% channel map from attune_config (bead_reference sets ref.chan from the date;
% cfg.cruise gives ref.OD2, cfg.gl1_is_od2_scatter from attune_config.OD2).
% Event-level agreement as in score_agreement (classes 1, 2, 3, 4, 7 kept;
% agreement_pe_merged counts 3/4/5 as one).
%
% classdir     - folder with the old class files and FCSfileinfo.mat (default
%                <cruise>/bead_calibrated; EN608: <cruise>/Trial_202008). '' and
%                no class files: v2 counts only (no agreement)
% compare_orig - function handle @(fcsdat, fcshdr, ref, cfg, pro) for a second
%                classifier to score on the same files (e.g. the pre-change v2
%                reading raw channel names), or [] for none
% Saves outdir/agreement_<CRUISE>_layout.mat (table V)
%
% Nov 2026, classifier redesign: pre-2019 channel layout

if nargin < 4, classdir = ''; end
if nargin < 5, compare_orig = []; end
addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cruise = extractAfter(cruisedir, '_');
cpath = fullfile(root, cruisedir);
bp = fullfile(cpath, 'bead_calibrated');
if isempty(classdir), classdir = bp; end
conf = attune_config(cruise);
cfg = classifier_config();
cfg.cruise = cruise;
cfg.gl1_is_od2_scatter = strcmpi(conf.OD2, 'GL1');
bs = load_beadstat(bp);
haveclass = ~isempty(dir(fullfile(classdir, 'class', '*.mat'))); %(EN608: class folder holds only plots)
if isfile(fullfile(classdir, 'FCSfileinfo.mat'))
    load(fullfile(classdir, 'FCSfileinfo.mat'), 'FCSfileinfo')
    Fi = FCSfileinfo;
    if ismember('QC_flag', Fi.Properties.VariableNames), Fi = Fi(Fi.QC_flag == 1, :); end
else
    d = dir(fullfile(cpath, 'FCS', '*.fcs'));
    Fi = table({d.name}', 'VariableNames', {'fcslist'});
    Fi.matdate_start = nan(height(Fi), 1);
end
Fi = Fi(~contains(Fi.fcslist, {'bead' 'Bead' 'test' 'Test' 'Rinse' 'discrete' 'Culture' 'Dock' 'PT'}), :);
if haveclass
    Fi = Fi(cellfun(@(c) exist(fullfile(classdir, 'class', regexprep(c, '\.fcs$', '.mat')), 'file') == 2, Fi.fcslist), :);
end
if all(isnan(Fi.matdate_start)), Fi = sortrows(Fi, 'fcslist'); else, Fi = sortrows(Fi, 'matdate_start'); end
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
pick = unique(round(linspace(1, height(Fi), min(nfiles, height(Fi)))));
rows = {};
for k = pick
    fn = Fi.fcslist{k};
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    t = Fi.matdate_start(k);
    if isnan(t), t = datenum([fcshdr.date ' ' fcshdr.starttime]); end %#ok<DATNM>
    ref = bead_reference(bs, t, cfg);
    row = table(string(cruise), string(fn), t, size(fcsdat, 1), 'VariableNames', {'cruise' 'file' 'time' 'n_events'});
    if haveclass
        c = load(fullfile(classdir, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class');
        old = c.class(:, 1); old(~ismember(old, [0 1 2 3 4 7])) = 0;
    else
        old = [];
    end
    try
        [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, Fi.pro_measured(k));
        row = [row score(old, new, 'v2')]; %#ok<AGROW>
        row.size_from_gl1 = info.size_from_gl1;
        row.ncomp = info.ncomp;
        row.trigger_edges = string(strjoin(fieldnames(info.trigger_edges)', ' '));
    catch e
        row.v2_error = string(e.message);
    end
    if ~isempty(compare_orig)
        try
            o = compare_orig(fcsdat, fcshdr, ref, cfg, Fi.pro_measured(k));
            row = [row score(old, o, 'orig')]; %#ok<AGROW>
        catch e
            row.orig_error = string(e.message);
        end
    end
    rows{end + 1} = row; %#ok<AGROW>
    disp(row)
end
V = stack_rows(rows);
save(fullfile(outdir, ['agreement_' cruise '_layout.mat']), 'V')
end

function V = stack_rows(rows)
% stack tables whose variables differ (failed runs lack the score columns)
names = {};
for i = 1:numel(rows), names = union(names, rows{i}.Properties.VariableNames, 'stable'); end
V = table;
for i = 1:numel(rows)
    r = rows{i};
    miss = setdiff(names, r.Properties.VariableNames);
    for v = reshape(miss, 1, [])
        if contains(v{1}, 'error') || strcmp(v{1}, 'trigger_edges'), r.(v{1}) = "";
        elseif strcmp(v{1}, 'old_n'), r.old_n = nan(1, 5);
        else, r.(v{1}) = NaN; end
    end
    V = [V; r(:, names)]; %#ok<AGROW>
end
end

function T = score(old, new, p)
n = arrayfun(@(x) sum(new == x), [1 2 3 4 5 7 8]);
T = array2table(n, 'VariableNames', strcat(p, '_n', {'1' '2' '3' '4' '5' '7' '8'}));
if isempty(old), return, end
either = old ~= 0 | new ~= 0;
o34 = old; o34(ismember(o34, [4 5])) = 3; n34 = new; n34(ismember(n34, [4 5])) = 3;
T.([p '_agreement']) = sum(old == new & old ~= 0) / sum(either);
T.([p '_agreement_pe_merged']) = sum(o34 == n34 & o34 ~= 0) / sum(either);
T.([p '_added']) = sum(old == 0 & new ~= 0);
T.([p '_dropped']) = sum(old ~= 0 & new == 0);
T.([p '_relabelled']) = sum(old ~= 0 & new ~= 0 & old ~= new);
if strcmp(p, 'v2'), T.old_n = arrayfun(@(x) sum(old == x), [1 2 3 4 7]); end
end
