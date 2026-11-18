function T = attune_pro_measured(T, fcs_path)
% function T = attune_pro_measured(T, fcs_path)
% Decide from the trigger settings whether each file could count
% Prochlorococcus, and set FCSfileinfo.pro_measured.
%
% T        - FCSfileinfo table (FCS_DateTimeList)
% fcs_path - optional folder of the FCS files. Rows without detector voltages
%            (hv_SSC missing or NaN, i.e. tables made before Nov 2026) get the
%            trigger settings and voltages from the FCS headers (fcs_header_info).
%
% Adds or replaces
%   trig_SSC_eq220  SSC trigger threshold scaled to SSC at 220 V (NaN if SSC is not a trigger)
%   trig_BL3_eq340  BL3 trigger threshold scaled to BL3 at 340 V (NaN if BL3 is not a trigger)
%   pro_measured    1 when the trigger lets typical Pro through, else 0
%
% Thresholds are compared in -H units scaled with signal ~ hv^8.5 (as in the
% bead voltage correction): thr_eq = thr * (hv_ref / hv)^8.5.
% A file counts Pro when SSC is one of its triggers and
%   SSC: trig_SSC_eq220 <= 500
%   BL3: trig_BL3_eq340 <= 175 (when BL3 is the other trigger)
%   any other channel as an AND trigger (PE, green scatter...): never (Pro lack it)
% Files triggered on chlorophyll or PE alone (OR_BL3 / OR_GL2) are never counted,
% including the early OR_BL3 100 runs (EN617 to RB1904, BL3 equivalent 100-167):
% they are not Pro settings and v2's Pro class has not been checked on them.
%
% Why these limits (survey of all FCS headers, Nov 2026):
% - every deliberate Pro run has BL3 equivalents of 23-146 and SSC 100-500
%   (AR98A Nov 2025 used AND_SSC 400 and 500 at 220 V)
% - the Euk halves of the EN688, EN706 and AR78 Euk/Syn alternating runs
%   (AND_SSC 300 / AND_BL3 300 at 340 V) and AR104's main runs (BL3 equivalent
%   204-244) are excluded: summer Pro have median BL3-H of about 50-160 at 340 V,
%   so these settings miss most of them (fall Pro, at 300-480, mostly pass)
% - SSC rarely limits: 80-99% of summer Pro and nearly all fall Pro have SSC-H
%   above 500 at 220 V
%
% Nov 2026, processing wrapper redesign

if nargin < 2, fcs_path = ''; end
n = height(T);
det = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
for v = {'trigger1_parameter' 'trigger2_parameter'}
    if ~ismember(v{1}, T.Properties.VariableNames), T.(v{1}) = repmat({''}, n, 1); end
    x = T.(v{1});
    if ~iscell(x), x = cellstr(string(x)); end
    x(cellfun(@(y) ~(ischar(y) || isstring(y)) || any(ismissing(string(y))), x)) = {''};
    T.(v{1}) = cellfun(@char, x, 'UniformOutput', false);
end
for v = [{'trigger1_threshhold' 'trigger2_threshhold' 'trigger1_hv' 'trigger2_hv'} strcat('hv_', det)]
    if ~ismember(v{1}, T.Properties.VariableNames), T.(v{1}) = NaN(n, 1); end
end

% backfill settings from the headers
need = find(isnan(T.hv_SSC) & isfinite(T.matdate_start));
if ~isempty(need) && ~isempty(fcs_path)
    fprintf('attune_pro_measured: reading %d FCS headers for trigger settings and voltages\n', numel(need))
    for k = need'
        try
            h = fcs_header_info(fullfile(fcs_path, T.fcslist{k}));
        catch
            continue
        end
        for i = 1:numel(det)
            if isfield(h.hv_det, det{i}), T.(['hv_' det{i}])(k) = h.hv_det.(det{i}); end
        end
        T.trigger1_parameter{k} = h.tr1_par; T.trigger1_threshhold(k) = h.tr1_level;
        T.trigger2_parameter{k} = h.tr2_par; T.trigger2_threshhold(k) = h.tr2_level;
        T.trigger1_hv(k) = hv_of(T, k, h.tr1_par);
        T.trigger2_hv(k) = hv_of(T, k, h.tr2_par);
    end
end

T.trig_SSC_eq220 = NaN(n, 1); T.trig_BL3_eq340 = NaN(n, 1);
pro = zeros(n, 1);
for k = 1:n
    tp = {T.trigger1_parameter{k} T.trigger2_parameter{k}};
    tl = [T.trigger1_threshhold(k) T.trigger2_threshhold(k)];
    use = ~cellfun(@isempty, tp);
    if ~any(use), pro(k) = NaN; continue, end
    tp = tp(use); tl = tl(use);
    ch = extractAfter(tp, '_');
    pass = false(size(tp));
    for j = 1:numel(tp)
        switch ch{j}
            case 'SSC'
                T.trig_SSC_eq220(k) = tl(j) * (220 / T.hv_SSC(k))^8.5;
                pass(j) = T.trig_SSC_eq220(k) <= 500;
            case 'BL3'
                T.trig_BL3_eq340(k) = tl(j) * (340 / T.hv_BL3(k))^8.5;
                pass(j) = T.trig_BL3_eq340(k) <= 175;
        end
    end
    issc = strcmp(ch, 'SSC');
    if ~any(issc), continue, end              % no SSC trigger: not Pro settings
    if any(startsWith(tp, 'AND')), pro(k) = all(pass); else, pro(k) = pass(issc); end
end
pro(isnan(pro)) = 0;
T.pro_measured = pro;
end

function hv = hv_of(T, k, par)
hv = NaN;
d = ['hv_' char(extractAfter(string(par), '_'))];
if ~isempty(par) && ismember(d, T.Properties.VariableNames), hv = T.(d)(k); end
end
