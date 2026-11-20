function A = score_agreement(cruisedir, nfiles, outdir)
% function A = score_agreement(cruisedir, nfiles, outdir)
% Event-level agreement between the existing class files and assign_class_v2
% for files spread through a cruise. agreement = events with the same class /
% events either method calls phytoplankton (3 dim and 4 bright cryptos and 5 Euk
% with Syn kept apart; agreement_pe_merged counts 3/4/5 as one; v2 class 8 is a disagreement wherever
% the old class differs).
% Saves output/agreement_<cruise>.mat
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
cfg.cruise = extractAfter(cruisedir, '_'); %cruise settings from attune_config, as the processing wrapper does
cfg.gl1_is_od2_scatter = strcmpi(attune_config(cfg.cruise).OD2, 'GL1');
cpath = fullfile(root, cruisedir);
bp = fullfile(cpath, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & ~contains(FCSfileinfo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
Fi = Fi(cellfun(@(c) exist(fullfile(bp, 'class', regexprep(c, '\.fcs$', '.mat')), 'file') == 2, Fi.fcslist), :);
Fi = sortrows(Fi, 'matdate_start');
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
pick = unique(round(linspace(1, height(Fi), min(nfiles, height(Fi)))));
A = table;
for k = pick
    fn = Fi.fcslist{k};
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, Fi.matdate_start(k), cfg);
    [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, Fi.pro_measured(k));
    st = info.syn_type(new == 2);
    c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class');
    old = c.class(:, 1); old(~ismember(old, [0 1 2 3 4 7])) = 0;
    o34 = old; o34(ismember(o34, [4 5])) = 3; n34 = new; n34(ismember(n34, [4 5])) = 3;
    either = old ~= 0 | new ~= 0;
    A = [A; table(string(cruisedir), string(fn), Fi.matdate_start(k), Fi.pro_measured(k), sum(either), ...
        sum(old == new & old ~= 0) / sum(either), sum(old == 0 & new ~= 0), sum(old ~= 0 & new == 0), sum(old ~= 0 & new ~= 0 & old ~= new), ...
        sum(old == 1), sum(new == 1), sum(old == 2), sum(new == 2), sum(old == 3), sum(new == 3), sum(old == 4), sum(new == 4), sum(new == 5), sum(new == 8), sum(o34 == n34 & o34 ~= 0) / sum(either), sum(st == 1), sum(st == 2), sum(st == 3), sum(st == 9), info.euk_chl_gap - getedge(info), info.n_euk_below_gap, sum(info.flag == 1), sum(info.flag == 2), sum(info.flag == 3), ...
        'VariableNames', {'cruise' 'file' 'time' 'pro_measured' 'n_phyto' 'agreement' 'n_added' 'n_dropped' 'n_relabelled' ...
        'old_euk' 'new_euk' 'old_syn' 'new_syn' 'old_peeuk' 'new_peeuk' 'old_peeuk4' 'new_peeuk4' 'new_euksyn' 'new_synpart' 'agreement_pe_merged' 'new_synA' 'new_synB' 'new_synC' 'new_syn_untyped' 'euk_gap_above_edge' 'n_euk_below_gap' 'flag_synpig' 'flag_below_gap' 'flag_edge_pile'})]; %#ok<AGROW>
    fprintf('%s agreement %.3f (added %d, dropped %d, relabelled %d)\n', fn, A.agreement(end), A.n_added(end), A.n_dropped(end), A.n_relabelled(end));
end
save(fullfile(outdir, ['agreement_' extractAfter(cruisedir, '_') '.mat']), 'A')
end

function e = getedge(info)
e = NaN; if isfield(info.trigger_edges, 'BL3'), e = info.trigger_edges.BL3; end
end
