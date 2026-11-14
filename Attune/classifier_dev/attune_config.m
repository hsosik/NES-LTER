function c = attune_config(cruise, t)
% function c = attune_config(cruise, t)
% Optical configuration of the Attune for a cruise, from the settings history
% \\sosiknas1\lab_data\Attune\Attune_Settings_configuration_history_examplecopy_20221130.xlsx
%
% cruise - cruise name, e.g. 'TN368' (used for the OD2 placement)
% t      - optional datenum of the data; decides the filter layout when given
%
% c.OD2     - where the OD2 filter sat: 'SSC', 'GL1' or 'None'
% c.fcb_beads - false for cruises without FCB bead-mix runs (EN608 to EN627:
%             performance-tracking beads only)
% c.layout  - 'v2019' (before EN644, Aug 2019) or 'current'
% c.chan    - struct: physical channel as named in the current layout
%             (bead_signature names) -> parameter name in this cruise's files,
%             '' where that filter was not installed
%
% Filter layouts (green 532 nm laser):
%   before EN644: GL1 575/36, GL2 620/15, GL3 680/40, GL4 780/60, no green scatter
%   EN644 on:     GL1 530/30+OD2 (side scatter), GL2 575/36, GL3 620/15, GL4 680/40
% Blue laser channels kept their roles (SSC 488/10, BL1 530/30, BL2 575/40 or
% 590/40, BL3 680/40 or 695/40).
%
% Nov 2026, bead stats redesign prototype

cruise = upper(char(cruise));
if ismember(cruise, {'EN608' 'AR28B' 'AR29' 'EN617' 'AR31A' 'TN368'})
    c.OD2 = 'SSC';
elseif ismember(cruise, {'EN627' 'AR34B' 'RB1904'})
    c.OD2 = 'None';
else
    c.OD2 = 'GL1';
end
c.fcb_beads = ~ismember(cruise, {'EN608' 'AR28B' 'AR29' 'EN617' 'AR31A' 'EN627'});
early = {'EN608' 'AR28B' 'AR29' 'EN617' 'AR31A' 'EN627' 'AR34B' 'RB1904' 'TN368'};
if nargin > 1 && ~isempty(t)
    isold = t < datenum(2019, 8, 19);
else
    isold = ismember(cruise, early);
end
names = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
for i = 1:numel(names), c.chan.(names{i}) = [names{i} '-H']; end
if isold
    c.layout = 'v2019';
    c.chan.GL1 = '';       % no green side scatter channel
    c.chan.GL2 = 'GL1-H';  % 575/36
    c.chan.GL3 = 'GL2-H';  % 620/15
    c.chan.GL4 = 'GL3-H';  % 680/40
else
    c.layout = 'current';
end
end
