function q = review_qc(cls, info, old)
% function q = review_qc(cls, info, old)
% Per-file QC numbers for the review viewer, from assign_class_v2 output (class,
% info) and, if given, the old classes. All fields are scalars (NaN when not
% found or not applicable) except trig, the trigger edges struct, so a struct array
% of these converts directly to a table or JSON.
%
% Nov 2026, classification review prototype

cls = double(cls(:));
q = struct();
q.n_events = numel(cls);
for k = [0 1 2 3 4 5 7 8 10]
    q.(sprintf('n%d', k)) = sum(cls == k);
end
q.pro_lower_bound = double(fieldor(info, 'pro_lower_bound', false)); % class 7 in a file without Pro settings: a lower bound
fl = zeros(size(cls)); if isfield(info, 'flag'), fl = double(info.flag(:)); end
for t = 1:5, q.(sprintf('flag%d', t)) = sum(fl == t); end % see review_config flag_names
st = zeros(size(cls)); if isfield(info, 'syn_type'), st = double(info.syn_type(:)); end
synb = ismember(cls, [2 5 8]);
for t = 1:4, q.(sprintf('syn_type%d', t)) = sum(synb & st == t); end
q.syn_untyped = sum(synb & st == 9);
q.syn_untyped_share = q.syn_untyped / max(1, sum(synb));
q.syn_ambiguous = sum(synb & st == 5); % per-file mixture: no type with posterior >= syn_mix.pmin
q.syn_ambiguous_share = q.syn_ambiguous / max(1, sum(synb));
q.syn_mix_method = ''; q.syn_mix_keep = ''; q.syn_mix_mu = nan(1, 4);
if isfield(info, 'syn_mix')
    q.syn_mix_method = info.syn_mix.method; tn = 'ABCD'; q.syn_mix_keep = tn(info.syn_mix.keep);
    if strcmp(info.syn_mix.method, 'mixture'), q.syn_mix_mu = info.syn_mix.mu; end % NaN for dropped types
end
q.n_catchall = fieldor(info, 'n_catchall', NaN);
% catch-all additions by class (classes 1 2 3 4 5 7 8). Shares are of the class's
% final count; the total is mostly Syn in Syn-rich files, so it is not divided by
% the Euk count. Counted from the final classes of the catch-all
% events (info.trace.catchall): later steps, e.g. the per-file Euk chl:size cut, can
% unassign catch-all Euks, which info.catchall_classes (taken at the catch-all step)
% still counts (AR92 15Aug(9) Sample(32): 202 there, 3 left as Euk)
cc = nan(1, 7);
if isfield(info, 'trace') && isfield(info.trace, 'catchall') && numel(info.trace.catchall) == numel(cls)
    ca = info.trace.catchall(:); cc = arrayfun(@(x) sum(ca & cls == x), [1 2 3 4 5 7 8]);
elseif isfield(info, 'catchall_classes') && numel(info.catchall_classes) == 7
    cc = double(info.catchall_classes(:)');
end
if q.n_catchall == 0, cc = zeros(1, 7); end
q.n_catchall_euk = cc(1);
q.n_catchall_syn = cc(2);
q.catchall_euk_share = q.n_catchall_euk / max(1, q.n1);
q.catchall_syn_share = q.n_catchall_syn / max(1, q.n2);
q.euk_chl_gap = finite_or_nan(fieldor(info, 'euk_chl_gap', NaN));
q.syn_PE_gap = finite_or_nan(fieldor(info, 'syn_PE_gap', NaN));
q.syn_PE_floor = finite_or_nan(fieldor(info, 'syn_PE_floor', NaN));
q.syn_PE_noise = finite_or_nan(fieldor(info, 'syn_PE_noise', NaN));
% cluster positions of the final classes (bead-relative log units): median size, chl
% (BL3) and PE (GL2), and for Syn the size spread (IQR). The viewer compares them with
% neighbouring files, so an abrupt shift of a cluster (a class taking a different
% population) shows. NaN when the class has too few events for a stable median
q = add_cluster(q, info, cls == 2, 'syn', 50, true);
q = add_cluster(q, info, cls == 1, 'euk', 50, false);
% dim and bright crypto apart: they sit ~1.5 log apart in PE, so a median over both
% jumps with the dim:bright mix
q = add_cluster(q, info, cls == 3, 'dcry', 20, false);
q = add_cluster(q, info, cls == 4, 'bcry', 20, false);
sr = fieldor(info, 'syn_ref', struct());
for f = {'size' 'BL3' 'PEchl' 'GL2GL3'}
    q.(['syn_ref_' f{1}]) = finite_or_nan(fieldor(sr, f{1}, NaN));
end
q.edge_pile_frac = fieldor(info, 'edge_pile_frac', NaN);
q.ncomp = fieldor(info, 'ncomp', NaN);
cs = fieldor(info, 'crypto_split', struct('n', 0, 'accepted', false));
q.crypto_n = fieldor(cs, 'n', 0);
q.crypto_split_accepted = double(fieldor(cs, 'accepted', false));
q.agreement = NaN;
if ~isempty(old)
    o = double(old(:)); o(~ismember(o, [0 1 2 3 4 5 7 8 10])) = 0;
    q.agreement = sum(o == cls & o ~= 0) / max(1, sum(o ~= 0 | cls ~= 0));
    for k = [0 1 2 3 4 5 7 8 10]
        q.(sprintf('old_n%d', k)) = sum(o == k);
    end
end
% Euk chl-size line and the rules drawn on the chl-minus-Euk-line panel
q.euk_line = review_euk_line(info.features, cls);
rl = fieldor(info, 'rules', struct());
q.euksyn_cut = -fieldor(rl, 'euksyn_max_chl_below_euk', NaN);
q.euk_min_size = fieldor(rl, 'euk_min_size', NaN);
q.trig = fieldor(info, 'trigger_edges', struct());
q.pro_AH_max = fieldor(info, 'pro_AH_max', NaN); % Pro wide-pulse cut (chl pulse-shape panel)
if isfield(info, 'syn_type_edges'), q.syn_type_edges = info.syn_type_edges(:)'; else, q.syn_type_edges = [NaN NaN NaN]; end
% type A also needs min(BL2, GL2) at least this far above its floor (drawn in the
% GL2 vs BL2 view as the corner GL2 >= floor + a, BL2 >= floor + a below the first edge)
q.syn_typeA_min_above_floor = fieldor(info, 'syn_typeA_min_above_floor', NaN);
q.syn_typeA_dim_slope = fieldor(info, 'syn_typeA_dim_slope', NaN); % dimmer events: A below e1 - slope*(limit - b)
q.floor_BL2 = NaN; q.floor_GL2 = NaN;
if isfield(info, 'floor') && isfield(info.floor, 'BL2'), q.floor_BL2 = info.floor.BL2; end
if isfield(info, 'floor') && isfield(info.floor, 'GL2'), q.floor_GL2 = info.floor.GL2; end
end

function q = add_cluster(q, info, m, nm, nmin, spread)
v = nan(1, 4);
if isfield(info, 'features') && sum(m) >= nmin
    F = info.features;
    v(1:3) = [median(F.size(m), 'omitnan') median(F.BL3(m), 'omitnan') median(F.GL2(m), 'omitnan')];
    v(4) = iqr(F.size(m));
end
q.([nm '_size']) = v(1); q.([nm '_chl']) = v(2); q.([nm '_pe']) = v(3);
if spread, q.([nm '_size_iqr']) = v(4); end
end

function v = fieldor(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); if isnumeric(v) || islogical(v), v = double(v(1)); end, else, v = d; end
end

function v = finite_or_nan(v)
if ~isfinite(v), v = NaN; end
end
