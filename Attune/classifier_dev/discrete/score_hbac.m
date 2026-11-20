function A = score_hbac(cruisedirs, nper, tag)
% function A = score_hbac(cruisedirs, nper, tag)
% assign_class_hbac vs the hand gates on SYBR runs, spread through each cruise
% (all depths), events inside the hand gates' time window only. Hand
% "bacteria" (Gated_Table class 3) included Pro, so it is compared with
% auto 9 + 7; hand Pro (4) is also put in that group.
% Saves output/discrete/score_hbac_<tag>.mat
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 2 || isempty(nper), nper = 6; end
if nargin < 3, tag = 'test'; end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = discrete_config();
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
A = table;
for ci = 1:numel(cruisedirs)
    cd_ = cruisedirs{ci}; cruise = extractAfter(cd_, '_'); cfg.cruise = cruise;
    G = load(fullfile(root, cd_, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    ix = find(contains(G.fcslist, 'hbac') & ~cellfun(@isempty, G.class) & ~contains(G.fcslist, {'skip' 'test' 'badly' 'virus'}));
    ix = ix(unique(round(linspace(1, numel(ix), min(nper, numel(ix))))));
    for i = ix'
        try, [dat, hdr] = fca_readfcs(fullfile(root, cd_, 'preserved', 'FCS', G.fcslist{i})); catch, continue, end
        ref = discrete_bead_reference(lab, datenum([hdr.date ' ' hdr.starttime]), cfg);
        if isempty(ref), continue, end
        [auto, info] = assign_class_hbac(dat, hdr, ref, cfg);
        hd = double(G.class{i}(:));
        if numel(hd) ~= numel(auto), continue, end
        hd(ismember(hd, [3 4])) = 97; a = auto; a(ismember(a, [7 9])) = 97;
        if any(hd), late = dat(:, 1) >= min(dat(hd ~= 0, 1)); else, late = true(size(hd)); end
        hd(~late) = 0; a(~late) = 0; ap = auto; ap(~late) = 0;
        either = hd ~= 0 | a ~= 0;
        A = [A; table(string(cruise), string(G.fcslist{i}), sum(either), sum(hd == a & hd ~= 0) / sum(either), ...
            sum(hd == 97), sum(a == 97), sum(ap == 7), sum(hd == 2), sum(a == 2), sum(hd == 1), sum(a == 1), ...
            info.chlx(1), info.chlx(2), info.pex(1), info.pex(2), info.n_core, ...
            'VariableNames', {'cruise' 'file' 'n_either' 'agreement' 'hand_bact' 'auto_bact_pro' 'auto_pro' 'hand_syn' 'auto_syn' ...
            'hand_euk' 'auto_euk' 'chlx_off' 'chlx_sd' 'pex_off' 'pex_sd' 'n_core'})]; %#ok<AGROW>
        fprintf('%s agreement %.3f bact %d/%d (Pro %d) Syn %d/%d Euk %d/%d\n', G.fcslist{i}, A.agreement(end), A.hand_bact(end), A.auto_bact_pro(end), ...
            A.auto_pro(end), A.hand_syn(end), A.auto_syn(end), A.hand_euk(end), A.auto_euk(end))
    end
end
save(fullfile(od, ['score_hbac_' tag '.mat']), 'A')
r = @(x, y) median(x(y > 100) ./ y(y > 100));
fprintf('\n%d runs: agreement median %.3f (10th percentile %.3f); auto/hand median: bacteria %.3f, Syn %.2f, Euk %.2f; bacteria within 10%%: %d\n', ...
    height(A), median(A.agreement), prctile(A.agreement, 10), r(A.auto_bact_pro, A.hand_bact), r(A.auto_syn, A.hand_syn), r(A.auto_euk, A.hand_euk), ...
    sum(abs(A.auto_bact_pro ./ A.hand_bact - 1) < 0.1))
end
