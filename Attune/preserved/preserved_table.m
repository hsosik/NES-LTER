function T = preserved_table(p)
% function T = preserved_table(p)
% Per-Niskin table of preserved (discrete) samples: process_cruise table step
% with Sample="preserved". One row per sample (cast, niskin, sizefrac; underway
% discretes by sample label). Each group comes from one run of that sample:
%
%   group                                          run used       if missing
%   Syn (2), Syn types, Syn on particle (8)        PE             CHL
%   Euk (1, by size), dim/bright crypto (3/4),     CHL            PE
%     Euk with Syn (5)
%   Pro (7)                                        pro            CHL, only where its Pro count is
%                                                                 not a lower bound (class_meta.pro.
%                                                                 lower_bound false); otherwise the
%                                                                 CHL count is pro_lower_bound_*
%   heterotrophic bacteria (9), Pro in SYBR runs   hbac
%     bacteria LNA / HNA (class file bact_type)
%
% Where a sample has several runs of one type (repeats), the last one is used
% and n_<type>_runs says how many there were. Counts use only events in the
% time gate (class file time_gate) and class_meta.vol_analyzed_ml (0.8 x run
% volume). Classes stay separate; bacteria_incl_pro_cells_per_ml (7 + 9 in
% the SYBR run) is an explicit sum for comparison with the older product.
% Concentrations: cells per ml; biovolume cubic microns per ml; carbon
% micrograms per litre (biovol2carbon, pg per cell). Euk size classes are
% equivalent spherical diameter from volume_cubic_microns (volume step);
% without volumes, biovolume and carbon are NaN and size classes empty.
%
% Saved as <outpath>\DiscreteTable.mat (T, table_metadata); sampling metadata
% from BottleMetadata.mat (metadata step) is joined when it exists.
%
% Nov 2026, process_cruise preserved mode

d = dir(fullfile(p.classpath, '*.mat'));
wl = warning('off', 'MATLAB:load:variableNotFound'); % bact_type: SYBR runs only
S = struct('file', {}, 'key', {}, 'runtype', {}, 'vol', {}, 'stats', {}, 'pro_lb', {}, 'hna_method', {}, 'meta', {});
for k = 1:numel(d)
    c = load(fullfile(d(k).folder, d(k).name), 'class', 'syn_type', 'bact_type', 'time_gate', 'class_meta', 'volume_cubic_microns');
    if ~isfield(c, 'class_meta') || ~isfield(c.class_meta, 'sample') || c.class_meta.sample ~= "preserved", continue, end
    m = c.class_meta;
    g = true(size(c.class)); if isfield(c, 'time_gate'), g = c.time_gate(:); end
    cl = c.class(:); cl = cl(g);
    st = zeros(size(cl), 'uint8'); if isfield(c, 'syn_type') && numel(c.syn_type) == numel(g), st = syn_type_current(c.syn_type(:), c.class_meta); st = st(g); end
    bt = []; if isfield(c, 'bact_type') && numel(c.bact_type) == numel(g), bt = c.bact_type(:); bt = bt(g); end
    hm = ""; if isfield(m, 'hbac') && isfield(m.hbac, 'hna') && isfield(m.hbac.hna, 'method'), hm = string(m.hbac.hna.method); end
    v = NaN(size(cl)); if isfield(c, 'volume_cubic_microns') && numel(c.volume_cubic_microns) == numel(g), v = c.volume_cubic_microns(:); v = v(g); end
    lb = true; if isfield(m, 'pro') && isfield(m.pro, 'lower_bound'), lb = logical(m.pro.lower_bound); end
    if m.cast > 0, key = sprintf('%d_%d_%s', m.cast, m.niskin, m.sizefrac); else, key = sprintf('0_0_%s_%s', m.sample_label, m.sizefrac); end
    meta = struct('cast', m.cast, 'niskin', m.niskin, 'sample_label', string(m.sample_label), 'sizefrac', string(m.sizefrac));
    S(end+1) = struct('file', regexprep(d(k).name, '\.mat$', '.fcs'), 'key', key, 'runtype', char(m.runtype), ...
        'vol', m.vol_analyzed_ml, 'stats', file_stats(cl, st, bt, v), 'pro_lb', lb, 'hna_method', hm, 'meta', meta); %#ok<AGROW>
end
warning(wl)
if isempty(S), error('preserved_table:none', 'no preserved class files in %s (run classify first)', p.classpath), end
F = load([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo'); F = F.FCSfileinfo;
[~, ia] = ismember({S.file}, F.fcslist);
tstart = F.matdate_start(ia);

keys = unique({S.key}, 'stable');
rows = cell(numel(keys), 1);
for i = 1:numel(keys)
    q = find(strcmp({S.key}, keys{i}));
    r = struct2table(S(q(1)).meta, 'AsArray', true);
    pick = @(rt) last_run(q, S, tstart, rt);
    [iPE, nPE] = pick('PE'); [iCHL, nCHL] = pick('CHL'); [iPRO, nPRO] = pick('pro'); [iHB, nHB] = pick('hbac');
    r.n_CHL_runs = nCHL; r.n_PE_runs = nPE; r.n_pro_runs = nPRO; r.n_hbac_runs = nHB;
    % Syn
    j = iPE; if isempty(j), j = iCHL; end
    r = add_group(r, 'syn', S, j, {'syn' 'syn_typeA' 'syn_typeB' 'syn_typeC' 'syn_typeD' 'syn_type_ambiguous' 'syn_type_unusable' 'syn_on_particle'}, {'syn'});
    % Euks
    j = iCHL; if isempty(j), j = iPE; end
    r = add_group(r, 'euk', S, j, {'euk' 'euk_leq2um' 'euk_leq3um' 'euk_leq5um' 'euk_leq10um' 'euk_leq20um' 'crypto_dim' 'crypto_bright' 'euk_with_syn' 'euk_pe_shift'}, ...
        {'euk' 'euk_leq2um' 'euk_leq3um' 'euk_leq5um' 'euk_leq10um' 'euk_leq20um' 'crypto_dim' 'crypto_bright' 'euk_with_syn' 'euk_pe_shift'});
    % Pro
    lbj = [];
    if ~isempty(iPRO), j = iPRO; src = "pro run";
    elseif ~isempty(iCHL) && ~S(iCHL).pro_lb, j = iCHL; src = "CHL run (complete)";
    else, j = []; src = ""; lbj = iCHL;
    end
    r = add_group(r, 'pro', S, j, {'pro'}, {'pro'});
    r.pro_source = src;
    if ~isempty(lbj)
        r.pro_lower_bound_cells_per_ml = S(lbj).stats.pro / S(lbj).vol;
        r.pro_lower_bound_filename = string(S(lbj).file);
    else
        r.pro_lower_bound_cells_per_ml = NaN; r.pro_lower_bound_filename = "";
    end
    % SYBR run: heterotrophic bacteria and Pro
    r = add_group(r, 'hbac', S, iHB, {'bacteria' 'bacteria_LNA' 'bacteria_HNA' 'pro_sybr' 'bacteria_incl_pro'}, {});
    r.bacteria_hna_method = ""; if ~isempty(iHB), r.bacteria_hna_method = S(iHB).hna_method; end
    rows{i} = r;
end
T = vertcat(rows{:});
T.cruise = repmat(string(p.cruise), height(T), 1);
T = movevars(T, 'cruise', 'Before', 1);
T = sortrows(T, {'cast' 'niskin'});
if exist([p.outpath 'BottleMetadata.mat'], 'file'), T = preserved_join_metadata(T, p); end

table_metadata = {"process_cruise preserved table (preserved_table)"; string(datetime('now')); p.classpath; ...
    "cells_per_ml: cells per ml; biovolume_concentration: cubic microns per ml; carbon_concentration: micrograms per litre"; ...
    "counts from events in the time gate (last 80% of acquisition) / class_meta.vol_analyzed_ml"; ...
    "syn_*: PE run (CHL if none); euk_*, crypto_*, euk_with_syn: CHL run (PE if none); euk = Euk without PE (classes 1 and 6), all sizes; euk_leqXum: classes 1 and 6 by ESD; crypto dim/bright = classes 3/4 (<= 20 um); euk_with_syn = class 5 (<= 20 um); euk_pe_shift = class 6, Euks whose PE channels are raised by preservation (all sizes; included in euk)"; ...
    "pro_*: pro run, or CHL run where its Pro cluster is complete; otherwise NaN and pro_lower_bound_cells_per_ml from the CHL run"; ...
    "bacteria = heterotrophic bacteria (class 9) in the SYBR run; pro_sybr = Pro (7) in the SYBR run; bacteria_incl_pro = 7 + 9"; ...
    "bacteria_LNA / bacteria_HNA: class 9 split by bact_type (sum = bacteria); bacteria_hna_method 'valley' = cut at the run's SYBR density valley, 'fixed' = no valley (HNA and LNA merged), split nominal"; ...
    "volumes: live size calibration applied as is (no preservation correction)"; ...
    "<group>_storage_days: run start minus collection time (days) for the run used; treat Pro from long-stored reruns as uncertain"};
save([p.outpath 'DiscreteTable.mat'], 'T', 'table_metadata')
fprintf('table: %d samples -> %sDiscreteTable.mat\n', height(T), p.outpath)
end


function s = file_stats(cl, st, bt, v)
% counts, biovolume and carbon (pg) per group for one file, time-gated events
esd = (v * 3 / 4 / pi).^(1/3) * 2;
C = real(biovol2carbon(v, 0));
grp = struct();
grp.syn = cl == 2; grp.euk = cl == 1 | cl == 6; %class 6: Euks with PE channels raised by preservation (27 Nov 2026)
grp.syn_typeA = cl == 2 & st == 1; grp.syn_typeB = cl == 2 & st == 2; grp.syn_typeC = cl == 2 & st == 3; grp.syn_typeD = cl == 2 & st == 4; grp.syn_type_ambiguous = cl == 2 & st == 5; grp.syn_type_unusable = cl == 2 & st == 9;
grp.syn_on_particle = cl == 8;
for e = [2 3 5 10 20]
    grp.(sprintf('euk_leq%dum', e)) = (cl == 1 | cl == 6) & esd <= e;
end
grp.crypto_dim = cl == 3 & ~(esd > 20); grp.crypto_bright = cl == 4 & ~(esd > 20); grp.euk_with_syn = cl == 5 & ~(esd > 20);
grp.euk_pe_shift = cl == 6; %class 6, all sizes: Euks whose PE channels are raised by preservation (matched live files show no PE); a subset of euk
grp.pro = cl == 7; grp.bacteria = cl == 9;
nobt = isempty(bt); if nobt, bt = zeros(size(cl), 'uint8'); end %class file without bact_type: LNA/HNA NaN
grp.bacteria_LNA = cl == 9 & bt == 1; grp.bacteria_HNA = cl == 9 & bt == 2; grp.pro_sybr = cl == 7; grp.bacteria_incl_pro = cl == 7 | cl == 9;
s = struct();
for f = fieldnames(grp)'
    x = grp.(f{1});
    s.(f{1}) = sum(x); s.([f{1} '_biovol']) = sum(v(x), 'omitnan'); s.([f{1} '_carbon']) = sum(C(x), 'omitnan');
    if all(isnan(v)), s.([f{1} '_biovol']) = NaN; s.([f{1} '_carbon']) = NaN; end
end
if nobt, s.bacteria_LNA = NaN; s.bacteria_HNA = NaN; end
if all(isnan(v)) % no volumes yet: no size classes
    for e = [2 3 5 10 20], s.(sprintf('euk_leq%dum', e)) = NaN; end
end
end


function [j, n] = last_run(q, S, tstart, rt)
i = q(strcmp({S(q).runtype}, rt));
n = numel(i); j = [];
if n == 0, return, end
[~, b] = max(tstart(i));
j = i(b);
end


function r = add_group(r, prefix, S, j, names, withbio)
% <name>_cells_per_ml (+ biovolume/carbon for names in withbio); <prefix>_filename, _volume_analyzed_ml
for k = 1:numel(names)
    nm = names{k};
    if isempty(j)
        r.([nm '_cells_per_ml']) = NaN;
        if ismember(names{k}, withbio), r.([nm '_biovolume_concentration']) = NaN; r.([nm '_carbon_concentration']) = NaN; end
        continue
    end
    s = S(j).stats; vol = S(j).vol;
    r.([nm '_cells_per_ml']) = s.(names{k}) / vol;
    if ismember(names{k}, withbio)
        r.([nm '_biovolume_concentration']) = s.([names{k} '_biovol']) / vol;
        r.([nm '_carbon_concentration']) = s.([names{k} '_carbon']) / vol / 1000;
    end
end
if isempty(j), r.([prefix '_filename']) = ""; r.([prefix '_volume_analyzed_ml']) = NaN;
else, r.([prefix '_filename']) = string(S(j).file); r.([prefix '_volume_analyzed_ml']) = S(j).vol;
end
end
