function h = fcs_header_info(filename)
% function h = fcs_header_info(filename)
% Read only the TEXT segment of an Attune FCS file (no event data) and return
% the settings used for processing bookkeeping:
%   h.date, h.starttime, h.stoptime  - $DATE, $BTIM, $ETIM strings
%   h.total  - $TOT (events)
%   h.vol    - $VOL (as written by the Attune, microlitres x 1000)
%   h.tr1_par, h.tr1_level, h.tr2_par, h.tr2_level - #TR1/#TR2 ('' / NaN if none)
%   h.par    - cellstr of parameter names ($PnN)
%   h.hv     - detector voltage per parameter ($PnV, NaN if absent)
%   h.hv_det - struct: detector (e.g. SSC, BL3) -> voltage, from the -H parameters
%
% Much faster than fca_readfcs when only the header is needed.
%
% Nov 2026, processing wrapper redesign

fid = fopen(filename, 'r', 'b');
if fid < 0, error('fcs_header_info:open', 'Cannot open %s', filename); end
c = onCleanup(@() fclose(fid));
hdr = fread(fid, 58, '*char')';
if numel(hdr) < 26 || ~startsWith(hdr, 'FCS')
    error('fcs_header_info:format', '%s is not an FCS file', filename);
end
t0 = str2double(hdr(11:18)); t1 = str2double(hdr(19:26));
fseek(fid, t0, 'bof');
txt = fread(fid, t1 - t0 + 1, '*char')';
delim = txt(1);
kv = strsplit(txt(2:end), delim, 'CollapseDelimiters', false);
if mod(numel(kv), 2), kv(end) = []; end
keys = upper(strtrim(kv(1:2:end))); vals = kv(2:2:end);
m = containers.Map(keys, vals, 'UniformValues', false);
get = @(k) value_or(m, k);

h.date = get('$DATE'); h.starttime = get('$BTIM'); h.stoptime = get('$ETIM');
h.total = str2double(get('$TOT'));
h.vol = str2double(get('$VOL'));
[h.tr1_par, h.tr1_level] = trig(get('#TR1'));
[h.tr2_par, h.tr2_level] = trig(get('#TR2'));
n = str2double(get('$PAR'));
h.par = cell(1, n); h.hv = NaN(1, n); h.hv_det = struct;
for i = 1:n
    h.par{i} = get(sprintf('$P%dN', i));
    h.hv(i) = str2double(get(sprintf('$P%dV', i)));
    if endsWith(h.par{i}, '-H')
        h.hv_det.(matlab.lang.makeValidName(extractBefore(h.par{i}, '-H'))) = h.hv(i);
    end
end
end

function v = value_or(m, k)
if isKey(m, k), v = strtrim(m(k)); else, v = ''; end
end

function [p, l] = trig(s)
p = ''; l = NaN;
if isempty(s), return, end
s = strsplit(s, ',');
p = strtrim(s{1});
if numel(s) > 1, l = str2double(s{2}); end
end
