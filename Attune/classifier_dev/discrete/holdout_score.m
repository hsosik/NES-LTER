function H = holdout_score(setname, rtonly)
% function H = holdout_score(setname, rtonly)
% Score one set of the held-out check (holdout_draw) with the frozen classifier
% (tag v2-holdout-test): assign_class_discrete for phytoplankton runs (CHL, PE, pro),
% assign_class_hbac for SYBR runs, against the hand gates, inside the hand gates'
% time window. Same class mapping as score_discrete / score_hbac. All preserved
% samples were run in the lab on the current layout, so gl1_is_od2_scatter = true
% for every file (score_discrete took it from the cruise-time config, which is
% wrong for the pre-2020 cruises EN617, EN627, TN368).
% setname - '2a', '2a_prime', '2b', 'reserve' (rows of holdout_lists.csv), or
%           'tuning' (the tuning score sample, re-scored with the same metric code)
% rtonly  - optional run type ('CHL', 'PE', 'pro', 'hbac'): score only those runs,
%           saved as score_<setname>_<rtonly>.mat
% Saves output/discrete/holdout/score_<setname>.mat (H: one row per run)
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
if strcmp(setname, 'tuning')
    L = tuning_list(root);
else
    L = readtable(fullfile(od, 'holdout', 'holdout_lists.csv'), 'TextType', 'string');
    L = L(L.set == setname, :);
end
if nargin > 1 && ~isempty(rtonly), L = L(L.runtype == rtonly, :); else, rtonly = ''; end
H = table; Gc = struct();
for i = 1:height(L)
    cdir = char(L.cruisedir(i)); fn = char(L.file(i)); r = L.runtype(i);
    cr = extractAfter(cdir, '_'); if isempty(cr), cr = cdir; end
    key = matlab.lang.makeValidName(cdir);
    if ~isfield(Gc, key)
        G = load(fullfile(root, cdir, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
        G.Properties.VariableNames = lower(G.Properties.VariableNames); Gc.(key) = G;
    end
    G = Gc.(key);
    gi = find(strcmp(G.fcslist, fn), 1);
    % Pro runs gated with the CHL run's workspace (awsfilename phyto_CHL_SSC/...): the CHL
    % polygons miss the Pro run's populations (EN712 C10N06: 78 of ~1,800 Syn inside the
    % Syn gate); not a valid hand reference (workspace_mismatch.csv, 22 Nov 2026)
    if r == "pro" && contains(fn, ["EN712_Feb2024_preserved(1)_phyto_CHL_SSC_pro_C10N06" "EN627_Feb2019_preserved(1)_phyto_CHL_SSC_pro_C17N08"]), continue, end
    try
        [dat, hdr] = fca_readfcs(fullfile(root, cdir, 'preserved', 'FCS', fn));
    catch err
        fprintf('%s: read failed (%s)\n', fn, err.message); continue
    end
    cfg = discrete_config(); cfg.cruise = cr; cfg.gl1_is_od2_scatter = true;
    ref = discrete_bead_reference(lab, datenum([hdr.date ' ' hdr.starttime]), cfg);
    if isempty(ref), fprintf('%s: no same-day lab beads\n', fn); continue, end
    hand = double(G.class{gi}(:));
    % hand classes rebuilt from the gate memberships with the wider gate-name map
    % (hand_class_rebuild: SYN_all / Syn_plus Syn gates were dropped by the exact-name map)
    if all(ismember({'gate_names' 'gate_assignments'}, G.Properties.VariableNames)) && ~isempty(G.gate_assignments{gi})
        hr = hand_class_rebuild(G.gate_names{gi}, G.gate_assignments{gi}, true);
        if numel(hr) == numel(hand), hand = hr(:); end
    end
    try
        if r == "hbac"
            [auto, info] = assign_class_hbac(dat, hdr, ref, cfg);
        else
            [auto, info] = assign_class_discrete(dat, hdr, ref, cfg, char(r));
        end
    catch err
        fprintf('%s: classify failed (%s)\n', fn, err.message); continue
    end
    if numel(hand) ~= numel(auto), fprintf('%s: hand/auto length differ\n', fn); continue, end
    if r == "hbac"
        h = hand; h(ismember(h, [3 4])) = 97; h(~ismember(h, [1 2 97])) = 0;
        a = double(auto); a(ismember(a, [7 9])) = 97; a(~ismember(a, [1 2 97])) = 0;
        % EN688 and HB1907 SYBR runs (other workspace: No_SYN_Euk_Pro, MysteryPE gates): the
        % hand Euk gate is mostly the unstained-particle band (PE:chl on the bacteria's SYBR
        % spectrum), 1.3-44x the Euks per mL of the same sample's CHL run (EN688 C11N01 507 m:
        % 10,512 vs 241); Euks are left out of both labellings there (23 Nov 2026)
        if contains(fn, ["EN688" "HB1907"]), h(h == 1) = 0; a(a == 1) = 0; end
    else
        h = hand; h(h == 4) = 7; h(ismember(h, [5 6])) = 34; h(h == 3) = 0;
        if contains(fn, 'EN727') && contains(fn, ["_C19N13" "_C19N02" "_C09N01" "_C06N06"]), h(h == 7) = 0; end
        % EN687 C16N19 and C15N05 CHL (2.7-2.8 m): the hand "Pro" (size -1, BL3 0.2, GL2
        % -1.7) are the population hand-gated Euk in the same sample's pro run; real Pro
        % (BL3 ~ -1.2 to -1.4) are below the CHL trigger (user OK 23 Nov 2026)
        if contains(fn, ["EN687_Aug2022_preserved(1)_phyto_CHL_SSC_C16N19" "EN687_Aug2022_preserved(1)_phyto_CHL_SSC_C15N05"]), h(h == 7) = 0; end
        a = double(auto); a(ismember(a, [3 4 5])) = 34;
    end
    if any(h ~= 0)
        late = dat(:, 1) >= min(dat(h ~= 0, 1)); h(~late) = 0; a(~late) = 0;
    end
    either = h ~= 0 | a ~= 0;
    row = table(string(setname), string(cr), r, string(fn), L.depth_m(i), sum(either), sum(h == a & h ~= 0) / sum(either), ...
        'VariableNames', {'set' 'cruise' 'runtype' 'file' 'depth_m' 'n_either' 'agreement'});
    for c = [1 2 7 34 97]
        row.(sprintf('hand_%d', c)) = sum(h == c); row.(sprintf('auto_%d', c)) = sum(a == c);
    end
    row.auto_8 = sum(a == 8);
    row.pro_lower_bound = NaN;
    if isfield(info, 'pro') && isstruct(info.pro), row.pro_lower_bound = double(info.pro.lower_bound); end
    H = [H; row]; %#ok<AGROW>
    fprintf('%s %s %s agr %.3f Euk %d/%d Syn %d/%d Pro %d/%d PEeuk %d/%d bact %d/%d\n', setname, cr, r, row.agreement, row.hand_1, row.auto_1, ...
        row.hand_2, row.auto_2, row.hand_7, row.auto_7, row.hand_34, row.auto_34, row.hand_97, row.auto_97);
end
if isempty(rtonly), save(fullfile(od, 'holdout', ['score_' setname '.mat']), 'H')
else, save(fullfile(od, 'holdout', ['score_' setname '_' rtonly '.mat']), 'H'), end
end

function L = tuning_list(root)
% the tuning score sample: score_discrete / score_hbac pick 6 runs per run type per
% cruise with linspace (AR92 empty-gate Pro runs left out, as in score_discrete)
tuning = {'20200201_EN649' '20200725_EN655' '20210203_EN661' '20210716_EN668' '20211108_AR61B' '20220729_EN687' ...
    '20230111_EN695' '20230807_EN706' '20240503_EN715' '20240906_EN720' '20250124_EN727' '20250814_AR92'};
ar92_empty = ["_C01N15" "_C02N01" "_C02N10" "_C06N03" "_C07N04" "_C07N08" "_C12N12" "_C14N02" "_C15N02" "_C16N02" ...
    "_C17N14" "_C21N08" "_C21N12" "_C22N05" "_C22N09" "_C22N13" "_C25N04"];
L = table;
for ci = 1:numel(tuning)
    G = load(fullfile(root, tuning{ci}, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    dep = NaN(height(G), 1); if ismember('depth_m', G.Properties.VariableNames), dep = double(G.depth_m(:)); end
    f = string(G.fcslist(:)); ne = ~cellfun(@isempty, G.class);
    rt = repmat("other", numel(f), 1);
    rt(contains(f, 'phyto_CHL')) = "CHL"; rt(contains(f, 'phyto_CHL') & contains(f, 'pro', 'IgnoreCase', true)) = "pro";
    rt(contains(f, 'phyto_PE')) = "PE"; rt(contains(f, 'hbac')) = "hbac";
    for r = ["CHL" "PE" "pro" "hbac"]
        if r == "hbac", ix = find(rt == r & ne & ~contains(f, ["skip" "test" "badly" "virus"]));
        else, ix = find(rt == r & ne); end
        ix = ix(unique(round(linspace(1, numel(ix), min(6, numel(ix))))));
        if r == "pro" && contains(tuning{ci}, 'AR92'), ix = ix(~contains(f(ix), ar92_empty)); end
        L = [L; table(repmat("tuning", numel(ix), 1), repmat(string(tuning{ci}), numel(ix), 1), f(ix), repmat(r, numel(ix), 1), dep(ix), ...
            'VariableNames', {'set' 'cruisedir' 'file' 'runtype' 'depth_m'})]; %#ok<AGROW>
    end
end
end
