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
for k = [0 1 2 3 4 5 7 8]
    q.(sprintf('n%d', k)) = sum(cls == k);
end
fl = zeros(size(cls)); if isfield(info, 'flag'), fl = double(info.flag(:)); end
for t = 1:3, q.(sprintf('flag%d', t)) = sum(fl == t); end
st = zeros(size(cls)); if isfield(info, 'syn_type'), st = double(info.syn_type(:)); end
synb = ismember(cls, [2 5 8]);
for t = 1:3, q.(sprintf('syn_type%d', t)) = sum(synb & st == t); end
q.syn_untyped = sum(synb & st == 9);
q.syn_untyped_share = q.syn_untyped / max(1, sum(synb));
q.n_catchall = fieldor(info, 'n_catchall', NaN);
q.catchall_share = q.n_catchall / max(1, q.n1);
q.euk_chl_gap = finite_or_nan(fieldor(info, 'euk_chl_gap', NaN));
q.syn_PE_gap = finite_or_nan(fieldor(info, 'syn_PE_gap', NaN));
q.syn_PE_floor = finite_or_nan(fieldor(info, 'syn_PE_floor', NaN));
q.syn_PE_noise = finite_or_nan(fieldor(info, 'syn_PE_noise', NaN));
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
    o = double(old(:)); o(~ismember(o, [0 1 2 3 4 5 7 8])) = 0;
    q.agreement = sum(o == cls & o ~= 0) / max(1, sum(o ~= 0 | cls ~= 0));
    for k = [0 1 2 3 4 5 7 8]
        q.(sprintf('old_n%d', k)) = sum(o == k);
    end
end
% Euk chl-size line and the rules drawn on the chl-minus-Euk-line panel
q.euk_line = review_euk_line(info.features, cls);
rl = fieldor(info, 'rules', struct());
q.euksyn_cut = -fieldor(rl, 'euksyn_max_chl_below_euk', NaN);
q.euk_min_size = fieldor(rl, 'euk_min_size', NaN);
q.trig = fieldor(info, 'trigger_edges', struct());
if isfield(info, 'syn_type_edges'), q.syn_type_edges = info.syn_type_edges(:)'; else, q.syn_type_edges = [NaN NaN]; end
end

function v = fieldor(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); if isnumeric(v) || islogical(v), v = double(v(1)); end, else, v = d; end
end

function v = finite_or_nan(v)
if ~isfinite(v), v = NaN; end
end
