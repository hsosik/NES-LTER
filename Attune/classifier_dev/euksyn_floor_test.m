function T = euksyn_floor_test(what, k_list, nmax)
% function T = euksyn_floor_test(what, k_list, nmax)
% Test a Euk-with-Syn (class 5) PE floor tied to the main Syn cluster instead of the
% 1st percentile of all Syn GL2 (euksyn_min_syn_PE_prc). Where a file's Syn include a
% dim group (AR43 11Mar2020a Sample(112): 18% of Syn ~2 log dimmer in PE and chl), the
% 1st percentile falls in it and Euk with Syn far below a typical Syn pass.
% Candidate floor: m - k * s, m the mode of Syn GL2 (kernel density), s the robust SD
% of Syn GL2 within 0.5 of the mode; never below the current floor. Applied to the
% final classes: class 5 below the floor -> Euk (1), as the current rule does.
%   what 'preserved': held-out preserved CHL/PE/pro runs (sets 2a, 2a_prime, 2b), with
%        the hand labels of the moved events (Euk, PE euk, Syn, none)
%        'live': live held-out files (L1, L2), with the old classes of the moved events
%   k_list  floors to try (default [2 3])
% T: one row per run: class 5 count, the floors, and per k the events moved and their
% reference labels. Saved to output/euksyn_floor_<what>.mat
%
% Nov 2026, Euk with Syn PE floor (test; v2 unchanged)

here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here), fullfile(here, 'discrete'))
if nargin < 2 || isempty(k_list), k_list = [2 3]; end
if nargin < 3 || isempty(nmax), nmax = Inf; end
warning('off', 'all')
rows = {};
switch what
    case 'preserved'
        od = fullfile(here, 'output', 'discrete');
        root = '\\sosiknas1\lab_data\Attune\cruise_data\';
        L = readtable(fullfile(od, 'holdout', 'holdout_lists.csv'), 'TextType', 'string');
        L = L(L.runtype ~= "hbac" & ismember(L.set, ["2a" "2a_prime" "2b"]), :);
        lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
        Gc = struct();
        for j = 1:min(height(L), nmax)
            cdir = char(L.cruisedir(j)); fn = char(L.file(j)); r = char(L.runtype(j));
            cr = extractAfter(cdir, '_'); if isempty(cr), cr = cdir; end
            try
                key = matlab.lang.makeValidName(cdir);
                if ~isfield(Gc, key)
                    G = load(fullfile(root, cdir, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
                    G.Properties.VariableNames = lower(G.Properties.VariableNames); Gc.(key) = G;
                end
                G = Gc.(key); gi = find(strcmp(G.fcslist, fn), 1);
                [dat, hdr] = fca_readfcs(fullfile(root, cdir, 'preserved', 'FCS', fn));
                cfg = discrete_config(); cfg.cruise = cr; cfg.gl1_is_od2_scatter = true;
                ref = discrete_bead_reference(lab, datenum([hdr.date ' ' hdr.starttime]), cfg);
                [c, info] = assign_class_discrete(dat, hdr, ref, cfg, r);
                hand = double(G.class{gi}(:));
                if all(ismember({'gate_names' 'gate_assignments'}, G.Properties.VariableNames)) && ~isempty(G.gate_assignments{gi})
                    hr = hand_class_rebuild(G.gate_names{gi}, G.gate_assignments{gi}, true);
                    if numel(hr) == numel(hand), hand = hr(:); end
                end
                % hand codes as in holdout_review: 1 Euk, 2 Syn, 4 Pro, 5/6 PE euk (34), 3 other
                h = hand; h(h == 4) = 7; h(ismember(h, [5 6])) = 34; h(h == 3) = 0;
                late = true(size(h)); if any(h ~= 0), late = dat(:, 1) >= min(dat(h ~= 0, 1)); end
                ref_lab = nan(size(c)); ref_lab(late) = h(late); % hand labels inside the hand time window only
            catch err
                fprintf('%s: %s\n', fn, err.message); continue
            end
            rows{end + 1} = one_run(fn, r, c, info, ref_lab, k_list, [1 2 34 0], {'hand_euk' 'hand_syn' 'hand_peeuk' 'hand_none'}); %#ok<AGROW>
            fprintf('%d/%d %s\n', j, height(L), fn);
        end
    case 'live'
        od = fullfile(here, 'output', 'live_holdout');
        L = readtable(fullfile(od, 'holdout_lists.csv'), 'TextType', 'string');
        L = L(ismember(L.set, ["L1" "L2"]), :);
        for j = 1:min(height(L), nmax)
            cd0 = char(L.cruisedir(j)); fn = char(L.file(j));
            try
                P = cruise_paths(cd0);
                cfg = classifier_config(); cfg.cruise = extractAfter(cd0, '_');
                cfg.gl1_is_od2_scatter = strcmpi(attune_config(cfg.cruise).OD2, 'GL1');
                [dat, hdr] = fca_readfcs(fullfile(P.basepath, 'FCS', fn));
                ref = bead_reference(load_beadstat(P.bead_dir), L.time(j), cfg);
                [c, info] = assign_class_v2(dat, hdr, ref, cfg, L.pro_measured(j));
                old = nan(size(c));
                try, o = double(load_old_class(P, fn)); if numel(o) == numel(c), old = o(:); end, catch, end %#ok<CTCH>
                old(ismember(old, [3 4])) = 34;
            catch err
                fprintf('%s: %s\n', fn, err.message); continue
            end
            rows{end + 1} = one_run(fn, 'live', c, info, old, k_list, [1 2 34 0], {'old_euk' 'old_syn' 'old_crypto' 'old_none'}); %#ok<AGROW>
            fprintf('%d/%d %s\n', j, height(L), fn);
        end
end
T = struct2table([rows{:}]);
save(fullfile(here, 'output', ['euksyn_floor_' what '.mat']), 'T', 'k_list')
fprintf('\n%s: %d runs, %d class 5 events\n', what, height(T), sum(T.n5));
for ki = 1:numel(k_list)
    kk = sprintf('k%g', k_list(ki) * 10);
    mv = T.(['moved_' kk]);
    fprintf('k = %g: %d events moved (%.1f%% of class 5) in %d runs', k_list(ki), sum(mv), 100 * sum(mv) / max(1, sum(T.n5)), sum(mv > 0));
    vn = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, [kk '_']));
    for v = vn, fprintf(', %s %d', extractAfter(v{1}, [kk '_']), sum(T.(v{1}))); end
    fprintf('\n');
end
end


function r = one_run(fn, rt, c, info, ref_lab, k_list, codes, names)
F = info.features; c = double(c(:));
s = c == 2; e = c == 5;
r = struct('file', string(fn), 'runtype', string(rt), 'n_syn', sum(s), 'n5', sum(e), ...
    'floor_now', NaN, 'syn_mode', NaN, 'syn_sd', NaN);
if isfield(info, 'euksyn_min_GL2'), r.floor_now = info.euksyn_min_GL2; end
% reference labels of the class 5 events kept by every floor (context)
for i = 1:numel(codes), r.(['all5_' names{i}]) = sum(e & ref_lab == codes(i)); end
[m, sd] = syn_core(F.GL2(s));
r.syn_mode = m; r.syn_sd = sd;
for ki = 1:numel(k_list)
    kk = sprintf('k%g', k_list(ki) * 10);
    f = m - k_list(ki) * sd;
    if isfinite(r.floor_now), f = max(f, r.floor_now); end
    if sum(s) < 50 || ~isfinite(f), f = -Inf; end
    mv = e & F.GL2 < f;
    r.(['floor_' kk]) = f; r.(['moved_' kk]) = sum(mv);
    for i = 1:numel(codes), r.([kk '_' names{i}]) = sum(mv & ref_lab == codes(i)); end
end
end


function [m, sd] = syn_core(g)
% mode of Syn GL2 and the robust SD of the Syn within 0.5 of it
m = NaN; sd = NaN; g = g(isfinite(g));
if numel(g) < 50, return, end
[d, x] = ksdensity(g, 'Bandwidth', 0.1);
[~, i] = max(d); m = x(i);
core = g(abs(g - m) < 0.5);
sd = 1.4826 * median(abs(core - median(core)));
end
