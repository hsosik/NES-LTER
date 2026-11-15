function R = compare_v2(cruisedir, nfiles, outdir)
% function R = compare_v2(cruisedir, nfiles, outdir)
% Run assign_class_v2 on files spread through a cruise and compare with the
% existing class files: per-file counts by product group and a pooled
% event-level confusion matrix. Old classes 3 and 4 are both "Euk with PE".
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
cpath = fullfile(root, cruisedir);
bp = fullfile(cpath, 'bead_calibrated');
bs = load_beadstat(bp);
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
Fi = FCSfileinfo(FCSfileinfo.QC_flag == 1 & ~contains(FCSfileinfo.fcslist, {'bead' 'Bead' 'test' 'Rinses' 'discrete'}), :);
Fi = Fi(cellfun(@(c) exist(fullfile(bp, 'class', regexprep(c, '\.fcs$', '.mat')), 'file') == 2, Fi.fcslist), :);
Fi = sortrows(Fi, 'matdate_start');
if ~ismember('pro_measured', Fi.Properties.VariableNames), Fi.pro_measured = zeros(height(Fi), 1); end
pick = unique(round(linspace(1, height(Fi), min(nfiles, height(Fi)))));

codes = [0 1 2 3 4 5 7 8];
names = {'none' 'euk' 'syn' 'dim_crypto' 'bright_crypto' 'euk_syn' 'pro' 'syn_particle'};
conf = zeros(numel(codes));
R = table;
for k = pick
    fn = Fi.fcslist{k};
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, Fi.matdate_start(k), cfg);
    [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, Fi.pro_measured(k));
    c = load(fullfile(bp, 'class', regexprep(fn, '\.fcs$', '.mat')), 'class');
    old = c.class(:, 1);
    old(~ismember(old, codes)) = 0;
    [~, io] = ismember(old, codes); [~, in] = ismember(new, codes);
    conf = conf + accumarray([io in], 1, [numel(codes) numel(codes)]);
    r = table(string(fn), Fi.matdate_start(k), Fi.pro_measured(k), info.ncomp, 'VariableNames', {'file' 'time' 'pro_measured' 'ncomp'});
    for j = 2:numel(codes)
        r.(['old_' names{j}]) = sum(old == codes(j));
        r.(['new_' names{j}]) = sum(new == codes(j));
    end
    R = [R; r]; %#ok<AGROW>
    fprintf('%s K=%d syn %d/%d (+%d particle) euk %d/%d pe %d/%d pro %d/%d\n', fn, info.ncomp, r.old_syn, r.new_syn, r.new_syn_particle, r.old_euk, r.new_euk, r.old_pe_euk, r.new_pe_euk, r.old_pro, r.new_pro);
end
C = array2table(conf, 'VariableNames', strcat('new_', names), 'RowNames', strcat('old_', names));
disp(C)
tag = extractAfter(cruisedir, '_');
save(fullfile(outdir, ['compare_v2_' tag '.mat']), 'R', 'C')
end
