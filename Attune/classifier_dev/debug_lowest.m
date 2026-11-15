% debug_lowest: why does v2 drop dim Syn in the lowest-agreement AR92 file?
addpath('..')
root = '\\sosiknas1\lab_data\Attune\cruise_data\20250814_AR92';
bp = fullfile(root, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
cfg = classifier_config();
fn = 'NESLTER_AR92_17Aug2025_ProAlt(12)_Group_day0_Euk(22).fcs';
k = find(strcmp(FCSfileinfo.fcslist, fn));
[d, h] = fca_readfcs(fullfile(root, 'FCS', fn));
ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
[new, info] = assign_class_v2(d, h, ref, cfg, FCSfileinfo.pro_measured(k));
c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class'); old = c.class(:, 1);
F = info.features; lost = old == 2 & new == 0;
disp(info.trigger_edges)
disp(info.comp)
fprintf('lost %d; median BL3 %.2f GL2 %.2f size %.2f PEchl %.2f GL2GL3 %.2f chlSSC %.2f\n', sum(lost), median(F.BL3(lost)), median(F.GL2(lost)), median(F.size(lost)), median(F.PEchl(lost)), median(F.GL2GL3(lost)), median(F.chlSSC(lost)))
r = info.syn_ref;
fprintf('syn ref size %.2f PEchl %.2f GL2GL3 %.2f BL3 %.2f\n', r.size, r.PEchl, r.GL2GL3, r.BL3)
% rerun the pipeline steps for the lost events
sig = cfg.sig;
pe_min = max(sig.pe_min_GL2, info.trigger_edges.GL2 + sig.trigger_margin);
chl_min = max(sig.chl_min, info.trigger_edges.BL3 + sig.trigger_margin);
has_pe = F.GL2 > pe_min & F.PEchl > sig.pe_min_PEchl & F.BL3 > sig.chl_min;
has_chl = (F.BL3 > chl_min | has_pe) & F.SSC > -10;
fprintf('lost: has_pe %d, has_chl(BL3 only) %d, neither %d\n', sum(lost & has_pe), sum(lost & F.BL3 > chl_min), sum(lost & ~has_pe & ~(F.BL3 > chl_min)))
fprintf('lost: GL2 below pe_min %d (pe_min %.2f), PEchl too low %d\n', sum(lost & F.GL2 <= pe_min), pe_min, sum(lost & F.PEchl <= sig.pe_min_PEchl))
