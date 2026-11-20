function A = score_discrete(cruisedirs, nper, cfg, tag, classify)
% function A = score_discrete(cruisedirs, nper, cfg, tag, classify)
% Event-level agreement between the hand gates (Gated_Table class) and an
% automatic classifier on preserved runs, spread through each cruise and over
% all depths. Hand classes are put on the shared scheme: 1 Euk, 2 Syn,
% 4 -> 7 Pro, 5/6 (low/high-PE euk) -> 34 "PE euk" (v2 3, 4 and 5 also -> 34).
% v2 class 8 (Syn on particle) is kept as 8 and counts as a disagreement.
% Only events inside the hand gates' time window are compared.
%
% cruisedirs - cell of cruise folders, e.g. {'20240906_EN720'}
% nper       - runs per run type (CHL, PE, pro) per cruise (default 6)
% cfg        - classifier config (default discrete_config if it exists, else classifier_config)
% tag        - name for the output file output/discrete/score_<tag>.mat
% classify   - optional handle @(fcsdat, fcshdr, ref, cfg, pro_measured) -> [class, info]
%              (default assign_class_v2)
% A          - one row per run: counts per class (hand and auto), agreement
% Also saves E: a subsample of events (features, hand class, auto class) for
% threshold design.
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 2 || isempty(nper), nper = 6; end
if nargin < 3 || isempty(cfg)
    if exist('discrete_config', 'file'), cfg = discrete_config(); else, cfg = classifier_config(); end
end
if nargin < 4 || isempty(tag), tag = 'test'; end
if nargin < 5 || isempty(classify), classify = @assign_class_v2; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
A = table; E = table;
for ci = 1:numel(cruisedirs)
    cd_ = cruisedirs{ci}; cruise = extractAfter(cd_, '_');
    cfg.cruise = cruise;
    cfg.gl1_is_od2_scatter = strcmpi(attune_config(cruise).OD2, 'GL1');
    G = load(fullfile(root, cd_, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    G = G(~cellfun(@isempty, G.class), :);
    rt = arrayfun(@(f) runtype(f), string(G.fcslist));
    for r = ["CHL" "PE" "pro"]
        ix = find(rt == r);
        ix = ix(unique(round(linspace(1, numel(ix), min(nper, numel(ix))))));
        for i = ix'
            fn = G.fcslist{i};
            try
                [fcsdat, fcshdr] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', fn));
            catch err
                warning('%s: %s', fn, err.message), continue
            end
            t = datenum([fcshdr.date ' ' fcshdr.starttime]);
            ref = discrete_bead_reference(lab, t, cfg);
            if isempty(ref), continue, end
            if contains(func2str(classify), 'discrete'), rarg = char(r); else, rarg = r == "pro"; end
            [auto, info] = classify(fcsdat, fcshdr, ref, cfg, rarg);
            hand = double(G.class{i}(:));
            if numel(hand) ~= numel(auto), continue, end
            h = hand; h(h == 4) = 7; h(ismember(h, [5 6])) = 34; h(h == 3) = 0;
            a = double(auto); a(ismember(a, [3 4 5])) = 34;
            % hand gates include the time gate (first 20% of the run dropped):
            % compare only events after the first hand-gated event
            if any(h ~= 0)
                late = fcsdat(:, 1) >= min(fcsdat(h ~= 0, 1));
                h(~late) = 0; a(~late) = 0;
            end
            either = h ~= 0 | a ~= 0;
            fb = cell2mat(struct2cell(ref.fallback));
            row = table(string(cruise), r, string(fn), t, max(fb), sum(either), sum(h == a & h ~= 0) / sum(either), ...
                'VariableNames', {'cruise' 'runtype' 'file' 'time' 'bead_fallback_days' 'n_either' 'agreement'});
            for c = [1 2 7 34]
                row.(sprintf('hand_%d', c)) = sum(h == c);
                row.(sprintf('auto_%d', c)) = sum(a == c);
            end
            row.auto_8 = sum(a == 8);
            row.hand_to_0 = sum(h ~= 0 & a == 0); row.auto_from_0 = sum(h == 0 & a ~= 0);
            A = [A; row]; %#ok<AGROW>
            % event subsample: all hand-gated events up to 2000, plus 1000 others
            F = info.features; k1 = find(h ~= 0); k0 = find(h == 0);
            k = [k1(randperm(numel(k1), min(2000, numel(k1)))); k0(randperm(numel(k0), min(1000, numel(k0))))];
            Ek = F(k, :); Ek.hand = h(k); Ek.auto = a(k);
            Ek.cruise = repmat(string(cruise), numel(k), 1); Ek.runtype = repmat(r, numel(k), 1); Ek.file = repmat(string(fn), numel(k), 1);
            E = [E; Ek]; %#ok<AGROW>
            fprintf('%s %s agreement %.3f  Euk %d/%d Syn %d/%d Pro %d/%d PEeuk %d/%d (hand/auto)\n', cruise, r, row.agreement, ...
                row.hand_1, row.auto_1, row.hand_2, row.auto_2, row.hand_7, row.auto_7, row.hand_34, row.auto_34);
        end
    end
end
save(fullfile(od, ['score_' tag '.mat']), 'A', 'E', 'cfg')
end

function rt = runtype(fn)
if contains(fn, 'hbac'), rt = "hbac";
elseif contains(fn, 'pro', 'IgnoreCase', true), rt = "pro";
elseif contains(fn, 'phyto_PE'), rt = "PE";
elseif contains(fn, 'phyto_CHL'), rt = "CHL";
else, rt = "other";
end
end
