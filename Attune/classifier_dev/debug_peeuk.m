% which pe_euk condition rejects the old PE euks in one AR92 file
addpath('..')
root = '\\sosiknas1\lab_data\Attune\cruise_data\20250814_AR92\';
fn = 'NESLTER_AR92_14Aug2025(1)_Group_day0_Sample(1).fcs';
cfg = classifier_config(); sig = cfg.sig;
bs = load_beadstat([root 'bead_calibrated']);
load([root 'bead_calibrated\FCSfileinfo.mat'])
k = find(strcmp(FCSfileinfo.fcslist, fn));
[fcsdat, fcshdr] = fca_readfcs([root 'FCS\' fn]);
ref = bead_reference(bs, FCSfileinfo.matdate_start(k), cfg);
[new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, 0);
F = info.features;
c = load([root 'bead_calibrated\class\' strrep(fn, '.fcs', '.mat')]);
o = ismember(c.class(:, 1), [3 4]);
fprintf('old PE euks: %d\n', sum(o))
disp(info.trigger_edges)
syn = new == 2;
fprintf('syn PEchl %.2f size %.2f\n', median(F.PEchl(syn)), median(F.size(syn)))
disp(summary(F(o, {'size' 'SSC' 'GL1' 'BL3' 'GL2' 'PEchl' 'GL2GL3'})))
fprintf('GL2>pe_min %d, PEchl>%g %d, BL3>chl_min %d, GL2>-1 %d, size %d, PEchl below syn %d\n', ...
    sum(o & F.GL2 > info.trigger_edges.GL2 + .05), sig.pe_min_PEchl, sum(o & F.PEchl > sig.pe_min_PEchl), sum(o & F.BL3 > sig.chl_min), ...
    sum(o & F.GL2 > sig.peeuk_min_GL2), sum(o & F.size > median(F.size(syn)) + .6), sum(o & F.PEchl < median(F.PEchl(syn)) - .5))
tabulate(new(o))
pe_min = max(sig.pe_min_GL2, info.trigger_edges.GL2 + sig.trigger_margin);
has_pe = F.GL2 > pe_min & F.PEchl > sig.pe_min_PEchl & F.BL3 > sig.chl_min;
sp = median(F.PEchl(syn)); ss = median(F.size(syn));
a = has_pe & F.GL2 > sig.peeuk_min_GL2; b = F.size > ss + sig.peeuk_size_above_syn; d = F.PEchl < sp - sig.peeuk_PEchl_below_syn | F.GL2GL3 < sig.peeuk_max_GL2GL3;
fprintf('all events: a %d, a&b %d, a&b&d %d; old PE: %d %d %d; new==3 %d\n', sum(a), sum(a&b), sum(a&b&d), sum(o&a), sum(o&a&b), sum(o&a&b&d), sum(new==3))
