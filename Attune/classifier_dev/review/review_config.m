function R = review_config()
% function R = review_config()
% Settings for the per-sample review frames (review_frame) and the cruise viewer
% (review_viewer). Axis limits are fixed per channel in bead-relative log10 units,
% so frames are comparable across a cruise (and across cruises).
%
% Nov 2026, classification review prototype

% views: x channel, y channel, x label, y label (same as plot_class_cytograms)
R.views = {'size' 'BL3' 'size (log10 SSC rel. 1 um bead)' 'chl (BL3)'; ...
    'BL3' 'GL2' 'chl (BL3)' 'PE (GL2)'; ...
    'GL3' 'GL2' 'GL3' 'PE (GL2)'; ...
    'size' 'PEchl' 'size' 'PE:chl (GL2-BL3)'; ...
    'BL3' 'BL2' 'chl (BL3)' 'BL2 (green fl., blue laser)'; ...
    'GL3' 'BL2' 'GL3' 'BL2 (PE, blue laser)'; ...
    'size' 'dEuk' 'size' 'chl minus Euk chl-size line'; ...
    'chlSSC' 'GL2GL3' 'chl:size (BL3 - size)' 'GL2:GL3'; ...      % Pro colour: Pro low in both, small Euks high (large Pro vs Euk)
    'BL3' 'BL3_AH' 'chl (BL3)' 'chl pulse shape log10(A/H)'; ...  % Pro wide-pulse cut; trigger-floor debris
    'GL2' 'BL2' 'PE (GL2)' 'PE (BL2)  Syn types'};                % Syn type view: keep last (review_frame, viewer)
% fixed axis limits per channel; events outside are drawn on the panel edge
R.lim = struct('size', [-2.6 2.8], 'BL3', [-1.3 3.2], 'GL2', [-3.5 1.7], 'GL3', [-2.5 2.8], ...
    'PEchl', [-4.4 1.8], 'BL2', [-2.5 2.3], 'FSCW', [1.0 2.8], 'dEuk', [-3 2], ...
    'chlSSC', [-1 3], 'GL2GL3', [-2.5 1.5], 'BL3_AH', [-0.8 0.8]);
% per-view limits that override R.lim: the Syn type view is zoomed on the Syn region
R.view_lim = cell(size(R.views, 1), 1);
R.view_lim{end} = [-1.2 1.8; -0.7 2.3]; % [x; y], equal spans so type edges (y = x + e) run at 45 degrees
R.px = 240;       % panel width and height in pixels
R.gap = 4;        % white pixels between panels
R.compare_rows = false; % true: add old-class and disagreement rows when old classes are given
R.format = 'png'; % 'png' (indexed colour) or 'jpg'
R.ncolors = 128; % colours in an indexed PNG

% classes: code, name, colour; abundant ones drawn as density, rare ones as points
R.codes = [1 2 3 4 5 7 8];
R.names = {'Euk' 'Syn' 'dim crypto' 'bright crypto' 'Euk with Syn' 'Pro' 'Syn on particle'};
R.cols = [0 .6 0; 1 .55 0; 0 .75 .9; .95 0 .75; .45 .1 .75; 0 .35 1; .85 0 0];
R.rare = [3 4 5 8];
R.col0 = [.45 .45 .45];                          % unassigned
R.rare_rad = 1;          % rare-class points: (2*rad+1)^2 pixel squares
R.rare_outline = false;  % dark outline round rare-class points
% Syn type view (last): Syn, Euk with Syn and Syn on particle coloured by type A/B/C
% (PUB:PEB low/mid/high), untyped dark grey, everything else light grey
R.syn_type_cols = [.9 .45 0; .8 0 .3; .15 .35 .85];
R.syn_untyped_col = [.2 .2 .2];
% show_flags: draw unassigned events in their flag colours below (true) or all in
% plain grey (false, default: fewer colours to read; the flag counts stay in the QC
% metrics and table)
R.show_flags = false;
% unassigned events by assign_class_v2 flag (density layers, class 0): 1 Syn-pigment
% particle, 2 below the Euk chl gap, 3 trigger-noise pile at the chl edge (pre-2019),
% 4 Pro-like event in a file without Pro settings (Pro lower bound), 5 Euk below the
% per-file chl:size cut. Flags 6-9 (split-Syn parts) are not drawn
R.flag_cols = [.6 .4 .2; .75 .68 .45; .45 .65 .85; .62 .58 .95; .6 .8 .5];
R.flag_names = {'Syn-pigment particle' 'below Euk chl gap' 'chl-edge noise pile' 'Pro lower bound' 'Euk below chl:size cut'};
R.diff_cols = [.85 0 0; 0 .3 1; .55 0 .75];        % v2 adds, v2 drops, relabelled
end
