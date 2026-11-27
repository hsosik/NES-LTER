function R = pe_catchall_cut_test()
% function R = pe_catchall_cut_test()
% PE runs of the held-out check (tuning, 2a, 2a_prime, 2b): is a fixed chl:size cut for
% v2's catch-all Syn (pe_catchall_syn_min_chlSSC) stable across samples, or should it be
% relative to each run's Syn core? Each run is classified once with the cut off; the
% candidate rules are then applied to the catch-all Syn:
%   abs t   - catch-all Syn with chl:size < t (t = 0, 0.1)
%   rel k   - chl:size < core median - k * robust SD of core chl:size, core = Syn not from
%             the catch-all (k = 2, 2.5, 3)
% Per run: Syn core chl:size median and SD, and for each rule the auto/hand Syn ratio,
% hand Syn removed and non-hand removed. Saves output/discrete/holdout/pe_catchall_cut_test.mat
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete', 'holdout');
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
lab = load(fullfile(fileparts(od), 'lab_beadstat.mat')).beadstat;
L = table;
for s = {'tuning' '2a' '2a_prime' '2b'}
    H = load(fullfile(od, ['score_' s{1} '_PE.mat'])).H;
    Lh = readtable(fullfile(od, 'holdout_lists.csv'), 'TextType', 'string');
    for i = 1:height(H)
        j = find(Lh.file == H.file(i), 1);
        if isempty(j) % tuning set: cruise folder from the cruise name
            d = dir(fullfile(root, ['*_' char(H.cruise(i))])); d = d([d.isdir]); cdir = string(d(1).name);
        else, cdir = Lh.cruisedir(j); end
        L = [L; table(string(s{1}), cdir, H.file(i), 'VariableNames', {'set' 'cruisedir' 'file'})]; %#ok<AGROW>
    end
end
rules = {'abs' 0; 'abs' 0.1; 'rel' 2; 'rel' 2.5; 'rel' 3};
R = table; Gc = struct();
for i = 1:height(L)
    cdir = char(L.cruisedir(i)); fn = char(L.file(i)); cr = extractAfter(cdir, '_'); if isempty(cr), cr = cdir; end
    key = matlab.lang.makeValidName(cdir);
    if ~isfield(Gc, key)
        G = load(fullfile(root, cdir, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table; G.Properties.VariableNames = lower(G.Properties.VariableNames); Gc.(key) = G;
    end
    G = Gc.(key); gi = find(strcmp(G.fcslist, fn), 1);
    try
        [d, h] = fca_readfcs(fullfile(root, cdir, 'preserved', 'FCS', fn));
        cfg = discrete_config(); cfg.cruise = cr; cfg.gl1_is_od2_scatter = true; cfg.sig.pe_catchall_syn_min_chlSSC = -Inf;
        ref = discrete_bead_reference(lab, datenum([h.date ' ' h.starttime]), cfg);
        [c, info] = assign_class_discrete(d, h, ref, cfg, 'PE');
    catch err
        fprintf('%s: %s\n', fn, err.message); continue
    end
    hc = hand_class_rebuild(G.gate_names{gi}, G.gate_assignments{gi}, true); hc = hc(:);
    if numel(hc) ~= numel(c), continue, end
    late = true(size(c)); if any(hc), late = d(:, 1) >= min(d(hc ~= 0, 1)); end
    F = info.features; ca = info.trace.catchall(:);
    syn = c(:) == 2 & late; core = syn & ~ca; cas = syn & ca;
    cm = median(F.chlSSC(core)); csd = 1.4826 * mad(F.chlSSC(core), 1);
    hs = sum(late & hc == 2);
    row = table(L.set(i), string(cr), string(fn), hs, sum(syn), sum(cas), cm, csd, ...
        'VariableNames', {'set' 'cruise' 'file' 'hand_syn' 'auto_syn' 'catchall_syn' 'core_chlSSC' 'core_sd'});
    for r = 1:size(rules, 1)
        if strcmp(rules{r, 1}, 'abs'), t = rules{r, 2}; nm = sprintf('abs%02d', round(10 * t));
        else, t = cm - rules{r, 2} * csd; nm = sprintf('rel%02d', round(10 * rules{r, 2})); end
        rm = cas & F.chlSSC < t;
        row.(['cut_' nm]) = t;
        row.(['ratio_' nm]) = (sum(syn) - sum(rm)) / hs;
        row.(['handlost_' nm]) = sum(rm & hc == 2);
        row.(['other_' nm]) = sum(rm & hc ~= 2);
    end
    row.ratio_none = sum(syn) / hs;
    R = [R; row]; %#ok<AGROW>
    fprintf('%d/%d %s %s core %.2f sd %.2f\n', i, height(L), L.set(i), fn, cm, csd);
end
save(fullfile(od, 'pe_catchall_cut_test.mat'), 'R')
k = R.hand_syn >= 100;
names = [{'none'} cellfun(@(a, b) ternary(strcmp(a, 'abs'), sprintf('abs%02d', round(10 * b)), sprintf('rel%02d', round(10 * b))), rules(:, 1)', rules(:, 2)', 'UniformOutput', false)];
for s = ["tuning" "2a" "2a_prime" "2b"]
    m = k & R.set == s; if ~any(m), continue, end
    fprintf('\n%s (%d PE runs):\n', s, sum(m));
    for n = names
        r = R.(['ratio_' n{1}])(m);
        if strcmp(n{1}, 'none'), hl = 0; ot = 0; else, hl = sum(R.(['handlost_' n{1}])(m)); ot = sum(R.(['other_' n{1}])(m)); end
        fprintf('  %-6s median %.3f within20 %.2f min %.2f max %.2f | removed: hand Syn %6d, other %6d\n', n{1}, median(r), mean(abs(r - 1) <= 0.2), min(r), max(r), hl, ot);
    end
end
fprintf('\nSyn core chl:size by set (median of run medians, range): \n');
for s = ["tuning" "2a" "2a_prime" "2b"]
    m = R.set == s; if any(m), fprintf('  %s: core %.2f (%.2f to %.2f), sd %.2f (%.2f to %.2f)\n', s, median(R.core_chlSSC(m)), min(R.core_chlSSC(m)), max(R.core_chlSSC(m)), median(R.core_sd(m)), min(R.core_sd(m)), max(R.core_sd(m))); end
end
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end
