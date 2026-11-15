function run_early_layout_validation(which, nfiles, origdir)
% function run_early_layout_validation(which, nfiles, origdir)
% Driver for the pre-Aug-2019 channel-map validation (Nov 2026).
% which   - 'regress': EN727 and AR92 (current layout), new vs pre-change v2,
%           event by event (must be identical); 'early': agreement on the early
%           cruises with old class files and v2 counts on those without
% origdir - folder holding assign_class_v2_orig / attune_features_orig (the
%           pre-change code), optional; adds the "before" scores
% Results: output/early_layout/

if nargin < 2, nfiles = 12; end
dev = fileparts(mfilename('fullpath'));
addpath(dev), addpath(fullfile(dev, '..'))
out = fullfile(dev, 'output', 'early_layout');
if ~isfolder(out), mkdir(out), end
orig = [];
if nargin > 2 && ~isempty(origdir)
    addpath(origdir)
    orig = @(d, h, r, c, p) assign_class_v2_orig(d, h, r, c, p);
end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
switch which
    case 'regress'
        R = table;
        for cd = {'20250124_EN727' '20250814_AR92'}
            if ~isfolder(fullfile(root, cd{1})), fprintf('missing %s\n', cd{1}); continue, end
            cruise = extractAfter(cd{1}, '_');
            bp = fullfile(root, cd{1}, 'bead_calibrated');
            bs = load_beadstat(bp);
            load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
            Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & ~contains(FCSfileinfo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
            Fi = sortrows(Fi, 'matdate_start');
            if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
            for k = unique(round(linspace(1, height(Fi), min(nfiles, height(Fi)))))
                [d, h] = fca_readfcs(fullfile(root, cd{1}, 'FCS', Fi.fcslist{k}));
                cfg = classifier_config(); %as score_agreement (no cruise given)
                ref = bead_reference(bs, Fi.matdate_start(k), cfg);
                [a, ia] = assign_class_v2(d, h, ref, cfg, Fi.pro_measured(k));
                [b, ib] = assign_class_v2_orig(d, h, ref, cfg, Fi.pro_measured(k));
                cfg.cruise = cruise; cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
                ref2 = bead_reference(bs, Fi.matdate_start(k), cfg);
                c2 = assign_class_v2(d, h, ref2, cfg, Fi.pro_measured(k));
                R = [R; table(string(cruise), string(Fi.fcslist{k}), numel(a), sum(a ~= b), sum(ia.syn_type ~= ib.syn_type), ...
                    sum(ia.flag ~= ib.flag), sum(c2 ~= b), 'VariableNames', {'cruise' 'file' 'n' 'class_diff' 'syntype_diff' 'flag_diff' 'class_diff_with_cruise_cfg'})]; %#ok<AGROW>
                disp(R(end, :))
            end
        end
        save(fullfile(out, 'regress_current_layout.mat'), 'R')
    case 'early'
        % with old class files
        C = {'20190705_TN368' '' ; '20190414_AR34B' ''; '20190512_RB1904' ''; '20190201_EN627' ''; ...
            '20180131_EN608' fullfile(root, '20180131_EN608', 'Trial_202008')};
        for i = 1:size(C, 1)
            fprintf('==== %s\n', C{i, 1});
            try, validate_early_layout(C{i, 1}, nfiles, out, C{i, 2}, orig); catch e, fprintf('%s failed: %s\n', C{i, 1}, e.message); end
        end
        % no old class files: v2 counts only
        for cd = {'20180404_AR28B' '20180414_AR29' '20180720_EN617' '20181020_AR31A'}
            fprintf('==== %s\n', cd{1});
            try, validate_early_layout(cd{1}, 6, out, '', orig); catch e, fprintf('%s failed: %s\n', cd{1}, e.message); end
        end
end
end
