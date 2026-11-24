function L = live_holdout_lists(nper, outdir)
% function L = live_holdout_lists(nper, outdir)
% Held-out file lists for checking assign_class_v2 on underway (live) samples
% (step 6 of the held-out plan; the discrete version is discrete/holdout_*).
% There are no hand gates for live samples, so v2 is scored against the old
% classes in each cruise's bead_calibrated_archive_* folder (cruise_paths).
%
% Sets (fixed here, before anything is scored):
%   T   the 60 tuning files (15 each from AR92, AR95, EN727, HRS2609, as picked by
%       score_agreement), rescored with the same code as a baseline for the others
%   L1  unseen cruises: two per season, one 2019-2021 and one 2022-2024, plus AR78
%       (fall 2023, the one unseen cruise with Pro-capable files), none used
%       in tuning or discussed in the README: nper random files each
%   L2  unseen files from the four tuning cruises and the Pro-survey cruises AR79,
%       AR82A and AE2426: nper random files each
% In cruises with Pro-capable files (pro_measured) half of each draw is Pro-capable.
%   R   reserve: four more unseen cruises, listed but not scored, kept for
%       confirming any fixes the L1/L2 review prompts
% Files are drawn at random with a fixed seed from QC-good sample files that have
% an old class file, leaving out the files tuning already looked at: the 15
% score_agreement files and the 8 files per cruise of the Nov 2026 underway survey
% (syn_comp_survey, evenly spaced through each cruise's old FCSfileinfo), and the\n% 12 Pro files per cruise of survey_pro.
% Left out as cruises: AR104 (Pro survey; no old classes), AR29 (no old
% classes), and EN627, EN688, AR88, AR34B, RB1904, TN368, AT46, EN655 (single
% files examined during tuning).
%
% nper   - files per cruise in L1 and L2 (default 15)
% outdir - default output/live_holdout
% Writes outdir/holdout_lists.csv (set, cruisedir, file, time, pro_measured) and
% outdir/holdout_lists_meta.mat (git commit, date).
%
% Nov 2026, classifier redesign (held-out check, live samples)

here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here))
if nargin < 1 || isempty(nper), nper = 15; end
if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'output', 'live_holdout'); end
if ~isfolder(outdir), mkdir(outdir), end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';

tune = {'20250814_AR92' '20251018_AR95' '20250124_EN727' '20260814_HRS2609'};
L1 = {'20200311_AR43' '20230429_HRS2303' ...   % spring
      '20190820_EN644' '20230807_EN706' ...    % summer
      '20191005_AR39B' '20221117_AR70B' '20231023_AR78' ... % fall (AR78: Pro-capable files)
      '20200201_EN649' '20240209_EN712'};      % winter
R = {'20210716_EN668' '20211108_AR61B' '20210203_EN661' '20260422_HRS2601'};

L = table;
add = @(L, s, cd, Fi) [L; table(repmat(string(s), height(Fi), 1), repmat(string(cd), height(Fi), 1), ...
    string(Fi.fcslist), Fi.matdate_start, Fi.pro_measured, 'VariableNames', {'set' 'cruisedir' 'file' 'time' 'pro_measured'})];
for s = {'T' 'L1' 'L2' 'R'}
    switch s{1}
        case 'T', cr = tune; case 'L1', cr = L1; case 'L2', cr = [tune {'20231108_AR79' '20240401_AR82A' '20241106_AE2426'}]; case 'R', cr = R;
    end
    for c = 1:numel(cr)
        [Fi, seen, tuned] = pool(root, cr{c});
        if s{1} == "T"
            Fi = Fi(ismember(Fi.fcslist, tuned), :);
        else
            Fi = Fi(~ismember(Fi.fcslist, seen), :);
            rng(sum(double(cr{c})) + 7 * strlength(s{1}), 'twister');
            % half the draw from Pro-capable files where a cruise has them
            ip = find(Fi.pro_measured == 1); in = find(Fi.pro_measured ~= 1);
            np = min(numel(ip), ceil(nper / 2) * ~isempty(ip)); nn = min(numel(in), nper - np);
            k = [ip(randperm(numel(ip), np)); in(randperm(numel(in), nn))];
            Fi = sortrows(Fi(k, :), 'matdate_start');
        end
        L = add(L, s{1}, cr{c}, Fi);
        fprintf('%-3s %-18s %3d files\n', s{1}, cr{c}, height(Fi));
    end
end
writetable(L, fullfile(outdir, 'holdout_lists.csv'));
[~, commit] = system(['git -C "' here '" rev-parse --short HEAD']);
meta = struct('commit', strtrim(commit), 'made', datestr(now, 'yyyy-mm-dd HH:MM'), 'nper', nper); %#ok<NASGU>
save(fullfile(outdir, 'holdout_lists_meta.mat'), 'meta')
fprintf('%d files listed: %s\n', height(L), fullfile(outdir, 'holdout_lists.csv'));
end

function [Fi, seen, tuned] = pool(root, cd)
% QC-good sample files with an old class file; seen = files tuning looked at
P = cruise_paths(cd, Root=root);
if isempty(P.old_class), error('live_holdout_lists:noold', '%s: %s', cd, P.old_note), end
% files tuning looked at, picked from the FCSfileinfo of the time (now in the archive)
q = load(fullfile(P.old, 'FCSfileinfo.mat'), 'FCSfileinfo'); Fo = q.FCSfileinfo;
Fo = Fo(~contains(Fo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
if ismember('QC_flag', Fo.Properties.VariableNames), Fo = Fo(Fo.QC_flag == 1, :); end
Fo = sortrows(Fo, 'matdate_start');
survey = Fo.fcslist(unique(round(linspace(1, height(Fo), min(8, height(Fo))))));
Fs = Fo(cellfun(@(c) exist(fullfile(P.old_class, regexprep(c, '\.fcs$', '.mat')), 'file') == 2, Fo.fcslist), :);
tuned = Fs.fcslist(unique(round(linspace(1, height(Fs), min(15, height(Fs))))));
% Pro survey (survey_pro, Nov 2026): 12 Pro-named Pro-capable files spread through the cruise
Fp = Fo(contains(Fo.fcslist, {'Pro' 'pro'}) & ~contains(Fo.fcslist, 'PT'), :);
q = load(fullfile(P.info_dir, 'FCSfileinfo.mat'), 'FCSfileinfo'); Fn = q.FCSfileinfo;
if ismember('pro_measured', Fn.Properties.VariableNames)
    Fp = Fp(ismember(Fp.fcslist, Fn.fcslist(Fn.pro_measured == 1)), :);
    prosurvey = Fp.fcslist(unique(round(linspace(1, height(Fp), min(12, height(Fp))))));
else, prosurvey = {};
end
seen = [survey; tuned; prosurvey(:)];
% the pool, from the rebuilt FCSfileinfo
q = load(fullfile(P.info_dir, 'FCSfileinfo.mat'), 'FCSfileinfo'); Fi = q.FCSfileinfo;
Fi = Fi(~contains(Fi.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
if ismember('QC_flag', Fi.Properties.VariableNames), Fi = Fi(Fi.QC_flag == 1, :); end
Fi = Fi(cellfun(@(c) exist(fullfile(P.old_class, regexprep(c, '\.fcs$', '.mat')), 'file') == 2, Fi.fcslist), :);
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi = attune_pro_measured(Fi, fullfile(P.basepath, 'FCS')); end
Fi = sortrows(Fi, 'matdate_start');
end
