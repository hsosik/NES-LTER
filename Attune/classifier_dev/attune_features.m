function [F, info] = attune_features(fcsdat, fcshdr, ref)
% function [F, info] = attune_features(fcsdat, fcshdr, ref)
% Convert raw Attune height channels to a common scale:
%   log10(signal) corrected to the bead run's voltage (signal ~ hv^k) and
%   expressed relative to the bead peak, so 0 = bead peak position.
% Values <= 0 are set to the channel floor before the log.
%
% Channels are physical channels named as in the current filter layout
% (ref.chans, e.g. GL2 = 575/36 PE). Each is read from this cruise's file
% parameter through ref.chan (attune_config(...).chan, set by bead_reference):
% before Aug 2019 (layout 'v2019') PE was GL1-H, 620/15 GL2-H and 680/40 GL3-H,
% and there was no green side scatter (physical GL1 = NaN column, no floor).
% Without ref.chan the layout is taken from the FCS header date.
%
% ref  - struct from bead_reference: per channel bead log10 position, bead hv, k,
%        and chan (physical channel -> file parameter name)
% F    - table, one column per channel in ref.chans (log10 bead-relative units),
%        plus BL3_AH = log10(BL3 area / height), the chl pulse shape
% info - struct with per-channel sample hv, floor in F units, trigger info
%        (trigger names translated to physical channels, e.g. file OR_GL1 ->
%        OR_GL2 before Aug 2019; raw names in trigger_file), col (physical
%        channel -> column of fcsdat, 0 when not in the file)
%
% Nov 2026, classifier redesign prototype; physical channel map Nov 2026

if isfield(ref, 'chan')
    chan = ref.chan;
else
    t = [];
    if isfield(fcshdr, 'date') && ~isempty(fcshdr.date)
        try, t = datenum(fcshdr.date); catch, end %#ok<DATNM>
    end
    c = attune_config('', t); chan = c.chan;
end
pn = {fcshdr.par.name};
F = table;
info = struct('hv', struct(), 'floor', struct(), 'col', struct());
info.trigger_file = string(fcshdr.tr1_par);
info.trigger = phys_trigger(info.trigger_file, chan);
info.trigger_level = fcshdr.tr1_level;
info.trigger2 = "";
info.trigger2_level = NaN;
info.trigger2_file = "";
if isfield(fcshdr, 'tr2_par')
    info.trigger2_file = string(fcshdr.tr2_par);
    info.trigger2 = phys_trigger(info.trigger2_file, chan);
    info.trigger2_level = fcshdr.tr2_level;
end
for i = 1:numel(ref.chans)
    ch = ref.chans{i};
    j = [];
    if isfield(chan, ch) && ~isempty(chan.(ch)), j = find(strcmp(pn, chan.(ch)), 1); end
    if isempty(j)
        F.(ch) = nan(size(fcsdat, 1), 1);
        info.col.(ch) = 0;
        continue
    end
    info.col.(ch) = j;
    hv = fcshdr.par(j).hv;
    % bead position brought to this sample's voltage
    beadpos = ref.bead_log.(ch) + ref.k.(ch) * log10(hv / ref.bead_hv.(ch));
    x = fcsdat(:, j);
    x(x < ref.floor_raw) = ref.floor_raw;
    F.(ch) = log10(x) - beadpos;
    info.hv.(ch) = hv;
    info.floor.(ch) = log10(ref.floor_raw) - beadpos;
end
% BL3 pulse shape (Nov 2026): log10(area / height) of the chl pulse, from the
% -A column of the same physical channel. Both scale alike with voltage, so no
% bead correction. Small single cells give short pulses like Syn (about 0);
% faint signals with wide pulses (about 0.3-0.4) sit among Pro in AR95 Pro runs
% and are not Pro (assign_class_v2). NaN when the area column is missing
F.BL3_AH = nan(size(fcsdat, 1), 1);
info.col.BL3_A = 0;
if isfield(chan, 'BL3') && ~isempty(chan.BL3) && isfield(info.col, 'BL3') && info.col.BL3 > 0
    ja = find(strcmp(pn, regexprep(chan.BL3, '-H$', '-A')), 1);
    if ~isempty(ja) && ja ~= info.col.BL3
        info.col.BL3_A = ja;
        F.BL3_AH = log10(max(fcsdat(:, ja), 1) ./ max(fcsdat(:, info.col.BL3), 1));
    end
end
end

function p = phys_trigger(tr, chan)
% file trigger parameter (e.g. "OR_GL1", "BL3") -> physical channel name
p = tr;
fch = regexp(tr, '[A-Z]+\d*(?=(-[AHW])?$)', 'match', 'once');
if strlength(fch) == 0, return, end
f = fieldnames(chan);
for i = 1:numel(f)
    if strcmp(regexprep(chan.(f{i}), '-[AHW]$', ''), fch)
        p = replace(tr, fch, f{i});
        return
    end
end
end
