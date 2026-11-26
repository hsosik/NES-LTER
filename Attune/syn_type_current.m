function st = syn_type_current(st, class_meta)
% function st = syn_type_current(st, class_meta)
% Per-event syn_type codes from a class file, in the current 4-type scheme
% (1/2/3/4 = types A/B/C/D in order of BL2-GL2, 9 = spectrum not usable, 0 = other).
% Class files written before 26 Nov 2026 used 3 types (1/2/3 = the present B/C/D,
% classifier setting syn_type_edges with 2 edges); their codes are shifted up by one.
% Those files have no type A: its events are inside their type 1 (now B).
%
% Nov 2026, Syn type D (classifier_dev/output/syn_type_d)

if nargin < 2 || isempty(st), return, end
e = [];
if isstruct(class_meta) && isfield(class_meta, 'cfg') && isfield(class_meta.cfg, 'sig') && isfield(class_meta.cfg.sig, 'syn_type_edges')
    e = class_meta.cfg.sig.syn_type_edges;
end
if numel(e) == 2
    m = st >= 1 & st <= 3;
    st(m) = st(m) + 1;
end
end
