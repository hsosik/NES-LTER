function T = preserved_join_metadata(T, p)
% function T = preserved_join_metadata(T, p)
% Add sampling metadata (BottleMetadata.mat from preserved_metadata) to the
% per-Niskin preserved table, matched by sample label and size fraction (the
% metadata step may set cast/niskin of underway discretes to -9999). Adds
% <group>_storage_days for each run used (syn, euk, pro, pro_lower_bound, hbac):
% run start (FCSfileinfo) minus collection time (date_sampled), in days.
%
% Nov 2026, process_cruise preserved mode

B = load([p.outpath 'BottleMetadata.mat'], 'M'); M = B.M;
v = {'date_sampled' 'latitude' 'longitude' 'depth_m' 'potemp090c' 'salinity' 'nearest_station' 'r2r_event' 'source'};
v = v(ismember(v, M.Properties.VariableNames));
T(:, intersect([v {'metadata_source'}], T.Properties.VariableNames)) = []; % joined before: replace
kT = string(T.sample_label) + "|" + string(T.sizefrac);
kM = string(M.sample_label) + "|" + string(M.sizefrac);
[ok, j] = ismember(kT, kM);
for i = 1:numel(v)
    x = M.(v{i});
    if isdatetime(x), y = NaT(height(T), 1);
    elseif isnumeric(x), y = NaN(height(T), 1);
    else, x = string(x); y = strings(height(T), 1);
    end
    y(ok) = x(j(ok));
    T.(v{i}) = y;
end
T = renamevars(T, 'source', 'metadata_source');
T = movevars(T, intersect([v(1:end-1) {'metadata_source'}], T.Properties.VariableNames, 'stable'), 'After', 'sizefrac');
if any(~ok)
    fprintf('table: no sampling metadata for %d samples\n', sum(~ok))
end

% storage time of each run used: run start minus collection time, days (long-stored
% reruns, e.g. EN687 samples rerun in Nov 2025, can lose Pro onto the chl trigger floor)
if ismember('date_sampled', T.Properties.VariableNames)
    F = load([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo'); F = F.FCSfileinfo;
    tc = datenum(T.date_sampled);
    for g = {'syn' 'euk' 'pro' 'pro_lower_bound' 'hbac'}
        fc = [g{1} '_filename'];
        if ~ismember(fc, T.Properties.VariableNames), continue, end
        [has, jf] = ismember(cellstr(T.(fc)), F.fcslist);
        sd = NaN(height(T), 1);
        sd(has) = F.matdate_start(jf(has)) - tc(has);
        T.([g{1} '_storage_days']) = sd;
        T = movevars(T, [g{1} '_storage_days'], 'After', fc);
    end
end
end
