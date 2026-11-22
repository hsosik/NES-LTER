function L = holdout_draw()
% function L = holdout_draw()
% Seeded draw of the held-out runs for the over-tuning check (frozen code: tag
% v2-holdout-test, 10f93ee). Sets (holdout_proposal.md in the shared notes, approved 21 Nov 2026):
%   2a       unseen cruises, 6 runs per run type (CHL, PE, pro, hbac) per cruise, rng(2026);
%            CHL and hbac split between <= 20 m and deeper where possible
%   2a_prime partly seen cruises (HB1907, AR38, AR39B, EN688), same rule, samples
%            examined during tuning excluded
%   2b       unseen runs from the 12 tuning cruises, 4 per run type, rng(2027),
%            excluding the tuning score sample, samples examined during tuning,
%            EN727 no-real-Pro runs, AR92 empty-gate Pro runs, EN687 preserved(2) reruns
%   reserve  2 per run type per cruise from the 2a and tuning cruises, rng(2028),
%            disjoint from the sets above; not scored until fixes need confirming
% Writes output/discrete/holdout/holdout_lists.csv
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete', 'holdout');
if ~exist(od, 'dir'), mkdir(od), end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
unseen = {'20180720_EN617' '20190201_EN627' '20190705_TN368' '20190820_EN644' '20200311_AR43' '20210512_SG2105' ...
    '20220216_AT46' '20240209_EN712' '20260113_AR99' 'MVCO'};
partly = {'20190725_HB1907' '20190921_AR38' '20191005_AR39B' '20220806_EN688'};
tuning = {'20200201_EN649' '20200725_EN655' '20210203_EN661' '20210716_EN668' '20211108_AR61B' '20220729_EN687' ...
    '20230111_EN695' '20230807_EN706' '20240503_EN715' '20240906_EN720' '20250124_EN727' '20250814_AR92'};
% samples (cruise_CxxNyy) examined individually during tuning: all their runs are left out
seen = ["HB1907_C15N03" "HB1907_C11N02" "HB1907_C11N04" "HB1907_C04N04" "HB1907_C15N05" "AR38_C11N13" "AR38_C10N07" "AR38_C12N18" ...
    "AR39B_C06N14" "AR39B_C14N06" "AR39B_C11N08" "AR39B_C10N09" "AR39B_C04N05" "EN688_C01N18" "EN688_C08N19" "EN688_C10N18" ...
    "EN688_C11N18" "EN688_C02N19" "EN688_C04N17" "EN688_C12N18" "EN688_C05N20" ...
    "EN655_C15N09" "EN655_C15N11" "EN687_C05N08" "EN687_C10N10" "EN687_C18N10" "EN687_C07N05" "EN687_C18N17" "EN687_C26N05" ...
    "EN706_C02N10" "EN706_C06N08" "EN706_C14N10" "EN706_C09N15" "EN706_C14N03" "EN720_C11N09" "EN720_C16N05" "EN720_C24N02" ...
    "EN720_C07N03" "EN720_C19N05" "EN720_C11N07" "EN727_C19N13" "EN727_C19N02" "EN727_C09N01" "EN727_C06N06" "EN727_C10N08" ...
    "AR92_C20N10" "AR92_C20N07" "AR92_C12N02" "AR92_C21N05" "AR92_C26N05" "EN668_C001N01" "EN668_C002N17" "EN661_C12N03" ...
    "EN661_C14N07" "EN715_C10N02" "AR61B_C09N07" "EN649_C01N02"];
ar92_empty = ["_C01N15" "_C02N01" "_C02N10" "_C06N03" "_C07N04" "_C07N08" "_C12N12" "_C14N02" "_C15N02" "_C16N02" ...
    "_C17N14" "_C21N08" "_C21N12" "_C22N05" "_C22N09" "_C22N13" "_C25N04"];
rts = ["CHL" "PE" "pro" "hbac"];
L = table;
% tuning score sample (score_discrete / score_hbac pick 6 per run type with linspace)
tun = strings(0, 1);
for ci = 1:numel(tuning)
    [f, rt] = runs(root, tuning{ci});
    for r = rts
        ix = find(rt == r); ix = ix(unique(round(linspace(1, numel(ix), min(6, numel(ix))))));
        tun = [tun; f(ix)]; %#ok<AGROW>
    end
end
taken = tun;
L = [L; drawset('2a', unseen, 6, 2026)];
L = [L; drawset('2a_prime', partly, 6, 2026)];
L = [L; drawset('2b', tuning, 4, 2027)];
L = [L; drawset('reserve', [unseen tuning], 2, 2028)];
writetable(L, fullfile(od, 'holdout_lists.csv'))
disp(groupsummary(L, {'set' 'runtype'}))

    function T = drawset(name, cruises, n, seed)
        rng(seed, 'twister');
        T = table;
        for c = 1:numel(cruises)
            cdir = cruises{c};
            [f, rt, dep, cn] = runs(root, cdir);
            cr = string(extractAfter(cdir, '_')); if cr == "", cr = string(cdir); end
            ok = ~ismember(f, taken) & ~ismember(cr + "_" + cn, seen);
            if cr == "EN687", ok = ok & ~contains(f, "preserved(2)"); end
            if cr == "AR92", ok = ok & ~(rt == "pro" & contains(f, ar92_empty)); end
            for r = rts
                ix = find(rt == r & ok);
                if isempty(ix), continue, end
                if ismember(r, ["CHL" "hbac"]) && any(isfinite(dep(ix)))
                    sh = ix(dep(ix) <= 20); dp = ix(dep(ix) > 20 | ~isfinite(dep(ix)));
                    ns = min(numel(sh), ceil(n / 2)); nd = min(numel(dp), n - ns); ns = min(numel(sh), n - nd);
                    pick = [sh(randperm(numel(sh), ns)); dp(randperm(numel(dp), nd))];
                else
                    pick = ix(randperm(numel(ix), min(n, numel(ix))));
                end
                T = [T; table(repmat(string(name), numel(pick), 1), repmat(string(cdir), numel(pick), 1), f(pick), ...
                    repmat(r, numel(pick), 1), dep(pick), 'VariableNames', {'set' 'cruisedir' 'file' 'runtype' 'depth_m'})]; %#ok<AGROW>
                taken = [taken; f(pick)]; %#ok<AGROW>
            end
        end
    end
end

function [f, rt, dep, cn] = runs(root, cdir)
G = load(fullfile(root, cdir, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
G.Properties.VariableNames = lower(G.Properties.VariableNames);
G = G(~cellfun(@isempty, G.class), :);
f = string(G.fcslist(:));
keep = ~contains(f, ["skip" "test" "badly" "virus"]);
G = G(keep, :); f = f(keep);
rt = repmat("other", numel(f), 1);
rt(contains(f, 'phyto_CHL')) = "CHL";
rt(contains(f, 'phyto_CHL') & contains(f, 'pro', 'IgnoreCase', true)) = "pro";
rt(contains(f, 'phyto_PE')) = "PE";
rt(contains(f, 'hbac')) = "hbac";
dep = NaN(numel(f), 1);
if ismember('depth_m', G.Properties.VariableNames), dep = double(G.depth_m(:));
elseif ismember('depth', G.Properties.VariableNames), dep = double(G.depth(:)); end
cn = regexp(f, 'C\d+N\d+', 'match', 'once');
% normalise leading zeros, e.g. C014N03 -> C14N03 (as in the 'seen' list)
cn = regexprep(cn, 'C0*(\d\d+)N', 'C$1N'); cn = regexprep(cn, 'C0*(\d)N', 'C0$1N');
end
