function T = preserved_fileinfo(T)
% function T = preserved_fileinfo(T)
% Add the preserved-sample columns to an FCSfileinfo table (FCS_DateTimeList
% on <cruise>\preserved\FCS), from the file names:
%   sample_label  last name token without the repeat suffix, e.g. 'C01N05',
%                 'C01N05A', 'Ex12345' (underway discretes), 'B01'
%   cast, niskin  from _CxxNyy (or _xxNyy); 0 when the name has none (underway
%                 or other discretes, matched to the event log by preserved_metadata)
%   repeat        n from a '(n)' right after the sample label (a rerun of the same
%                 sample and run type), 0 when none
%   runtype       'CHL' (phyto_CHL_SSC), 'PE' (phyto_PE_SSC), 'pro'
%                 (phyto_CHL_SSC_pro), 'hbac' (hbac_SYBR_SSC), or 'other'
%   sizefrac      'total', 'less5', 'less20' when in the name (DY131), else ''
%   pro_measured  1 for pro and hbac runs (SYBR-stained Pro are class 7), else 0.
%                 Replaces the voltage rule (attune_pro_measured), which gives 0
%                 for hbac runs (AND_SSC + AND_BL1) although they count Pro.
%                 The rule's value is kept as pro_measured_trigger.
%
% Nov 2026, process_cruise preserved mode

n = height(T);
name = string(T.fcslist);
stem = regexprep(name, '\.fcs$', '', 'ignorecase');

rep = regexp(stem, '\((\d+)\)$', 'tokens', 'once');
T.repeat = zeros(n, 1);
has = ~cellfun(@isempty, rep);
T.repeat(has) = cellfun(@(x) str2double(x{1}), rep(has));
stem = regexprep(stem, '\(\d+\)$', '');
lab = strings(n, 1);
for i = 1:n
    t = split(stem(i), '_'); lab(i) = t(end);
end
% the CxxNyy token can be followed by more text (EN657 ..._C12N18_medium)
cnt = regexp(stem, '(?<=_)C?\d{1,3}N\d{1,2}[A-Za-z]?(?=_|$)', 'match', 'once');
hascn = strlength(string(cnt)) > 0;
lab(hascn) = string(cnt(hascn));
T.sample_label = cellstr(lab);

cn = regexp(lab, '^C?(\d{1,3})N(\d{1,2})[A-Za-z]?$', 'tokens', 'once');
T.cast = zeros(n, 1); T.niskin = zeros(n, 1);
ok = ~cellfun(@isempty, cn);
T.cast(ok) = cellfun(@(x) str2double(x{1}), cn(ok));
T.niskin(ok) = cellfun(@(x) str2double(x{2}), cn(ok));

rt = repmat("other", n, 1);
rt(contains(name, "phyto_CHL", 'IgnoreCase', true)) = "CHL";
rt(contains(name, "phyto_PE", 'IgnoreCase', true)) = "PE";
rt(contains(name, "phyto_CHL_SSC_pro", 'IgnoreCase', true)) = "pro";
rt(contains(name, "hbac", 'IgnoreCase', true)) = "hbac";
T.runtype = cellstr(rt);

sf = strings(n, 1);
for f = ["total" "less5" "less20"]
    sf(contains(name, f, 'IgnoreCase', true)) = f;
end
T.sizefrac = cellstr(sf);

if ismember('pro_measured', T.Properties.VariableNames)
    T.pro_measured_trigger = T.pro_measured;
end
T.pro_measured = double(ismember(rt, ["pro" "hbac"]));
T = movevars(T, {'sample_label' 'cast' 'niskin' 'repeat' 'runtype' 'sizefrac'}, 'After', 'fcslist');
end
