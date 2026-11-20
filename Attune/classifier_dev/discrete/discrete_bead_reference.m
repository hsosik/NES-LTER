function ref = discrete_bead_reference(lab, t, cfg)
% function ref = discrete_bead_reference(lab, t, cfg)
% Bead reference for a preserved (discrete) run made in the lab at time t
% (datenum): bead_reference on the lab bead runs from the same day only.
% A channel whose bead population failed QC in every same-day run takes that
% population from the nearest lab run (any day) where it passed, and is
% flagged. (Populations that passed on failed days sit within ~0.03 log10 of
% the nearest good run; compare_failed_bead_days, Nov 2026.)
%
% lab - lab bead table (lab_beadstat.mat)
% ref - as bead_reference, plus
%       ref.fallback.(ch) - 0 when from same-day runs, else days to the run used
%       ref.any_fallback  - true if any channel used a fallback
%       Empty ([]) when there is no lab bead run that day: the run is not
%       normalized (skip it).
%
% Nov 2026, discrete-sample classifier development

lday = floor(datenum(lab.time));
same = lday == floor(t);
if ~any(same), ref = []; return, end
w = warning('off', 'bead_reference:noQC');
ref = bead_reference(lab(same, :), t, cfg);
P = {'05' '1' '6'};
ref.any_fallback = false;
for c = cfg.chans
    ch = c{1};
    ref.fallback.(ch) = 0;
    if ref.ntimes.(ch) > 0, continue, end
    p = P{cfg.bead_pop.(ch)};
    good = find(lab.(['QC' p]) & lab.(['peak' p]).([ch '_H']) > 0 & lab.hv.(ch) > 0);
    if ismember('excluded', lab.Properties.VariableNames), good = good(~lab.excluded(good)); end
    if isempty(good), continue, end
    [d, j] = min(abs(datenum(lab.time(good)) - t));
    r = bead_reference(lab(good(j), :), t, cfg);
    ref.bead_log.(ch) = r.bead_log.(ch); ref.bead_hv.(ch) = r.bead_hv.(ch);
    ref.k.(ch) = r.k.(ch); ref.ntimes.(ch) = r.ntimes.(ch); ref.runs.(ch) = r.runs.(ch);
    ref.fallback.(ch) = d;
    ref.any_fallback = true;
end
warning(w)
end
