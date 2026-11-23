function R = review_config_hbac()
% function R = review_config_hbac()
% Review frames and viewer settings for SYBR-stained bacteria runs (preserved
% hbac_SYBR_SSC runs, assign_class_hbac). Same viewer as the phytoplankton runs
% (review_viewer, R.mode 'hbac'), with panels chosen to show how the classes
% separate and the boundaries assign_class_hbac draws:
%   1 SSC vs SYBR (BL1)    bacteria cloud; HNA/LNA tilted cut, stain line, SYBR
%                          trigger edge and edge margin, bacteria size limit
%   2 tilted SYBR histogram  y = BL1 - slope (SSC + 2) of the bacteria, LNA / HNA
%                          coloured, with the cut (valley or fixed)
%   3 SYBR vs chl (BL3)    chl spillover line of the bacteria and the chl-excess cut
%   4 SYBR vs PE (GL2)     PE spillover line and the PE-excess cut
%   5 SSC vs chl           Pro / Euk size split, chl:size minimum
%   6 chl vs PE            Syn PE:chl minimum
%   7 SSC vs PE:chl minus the bacteria's   0 = the SYBR spectrum (unstained-particle
%                          band); Syn, Euk and stained-Euk / large-Pro cuts
%   8 chl:size vs chl      large deep Pro limits (dim, low chl:size)
% Classes: bacteria (9) split by bact_type into LNA and HNA (codes 91, 92 here),
% Syn 2, Pro 7, Euk 1, unassigned 0.
%
% Nov 2026, preserved-sample review (SYBR runs)

R = review_config(); % frame size, format and colour handling as for phytoplankton runs
R.mode = 'hbac';
R.views = {'SSC' 'BL1' 'size (log10 SSC rel. 1 um bead)' 'SYBR (BL1)'; ...
    'hnay' 'hist' 'tilted SYBR: BL1 - slope (SSC + 2)' 'bacteria (rel. density)'; ...
    'BL1' 'BL3' 'SYBR (BL1)' 'chl (BL3)'; ...
    'BL1' 'GL2' 'SYBR (BL1)' 'PE (GL2)'; ...
    'SSC' 'BL3' 'size (SSC)' 'chl (BL3)'; ...
    'BL3' 'GL2' 'chl (BL3)' 'PE (GL2)'; ...
    'SSC' 'dS' 'size (SSC)' 'PE:chl minus the bacteria''s'; ...   % colour vs the SYBR spectrum: particle band ~0, Pro ~-0.6, Euks -0.6 to -0.9, Syn above
    'chlSSC' 'BL3' 'chl:size (BL3 - SSC)' 'chl (BL3)'};          % large deep Pro (dim, low chl:size) vs Euks
R.lim = struct('SSC', [-3.2 1.2], 'BL1', [-2.4 1.6], 'BL3', [-1.6 3.0], 'GL2', [-3.2 2.0], 'hnay', [-2.4 0.4], 'hist', [0 1.05], ...
    'dS', [-2.5 2.0], 'chlSSC', [-1 3]);
R.view_lim = cell(size(R.views, 1), 1);
R.codes = [91 92 1 2 7];
R.names = {'LNA bacteria' 'HNA bacteria' 'Euk' 'Syn' 'Pro'};
R.cols = [.55 .35 .85; .85 .1 .3; 0 .6 0; 1 .55 0; 0 .35 1];
R.rare = [];
R.col0 = [.55 .55 .55];
R.flag_names = {};
end
