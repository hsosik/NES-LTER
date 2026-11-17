function ref = bead_reference(beadstat, t, cfg)
% function ref = bead_reference(beadstat, t, cfg)
% Bead peak positions to normalize a sample against: for each channel, the
% median over bead runs within cfg.bead_window_days of time t (datenum) in
% which that channel's bead population passed QC; falls back to the nearest
% such run. Each channel uses the bead population that is on scale there
% (cfg.bead_pop: 1 = 0.5, 2 = 1, 3 = 6 um) and that channel's voltage
% exponent (cfg.k); runs are brought to the first run's voltage before the
% median. Runs marked exclude in the bead review (beadstat.excluded, set by
% load_beadstat) are never used.
%
% ref.ntimes / ref.runs - per channel, how many and which bead runs (file
%            names) went into the median
% ref.chan / ref.layout - attune_config channel map for time t (read by
%            attune_features); ref.OD2 too when cfg.cruise is set
%
% beadstat - table from load_beadstat (process_beads_PT_adjust_3): per
%            population P = 05, 1, 6, QCP and peakP (KDE mode per channel,
%            linear units), and hv per channel
%
% Nov 2026, classifier redesign prototype; bead stats redesign

P = {'05' '1' '6'};
bt = datenum(beadstat.time);
ref.chans = cfg.chans;
% physical channel -> this file's parameter name. beadstat columns are already
% physical channels (measure_beads / PT transfer map them through attune_config),
% so only the sample files need the map; the layout follows from the date
% (filter change Aug 2019). cfg.cruise, when set, also gives the OD2 placement
if isfield(cfg, 'cruise'), cr = cfg.cruise; else, cr = ''; end
conf = attune_config(cr, t);
ref.chan = conf.chan;
ref.layout = conf.layout;
if ~isempty(cr), ref.OD2 = conf.OD2; end
ref.floor_raw = cfg.floor_raw;
ref.ntimes = struct();
for c = cfg.chans
    ch = c{1};
    b = P{cfg.bead_pop.(ch)};
    v = beadstat.(['peak' b]).([ch '_H']);
    hv = beadstat.hv.(ch);
    good = beadstat.(['QC' b]) & v > 0 & hv > 0;
    if ismember('excluded', beadstat.Properties.VariableNames), good = good & ~beadstat.excluded; end
    use = good & abs(bt - t) <= cfg.bead_window_days;
    if ~any(use) && any(good)
        q = find(good);
        [~, i] = min(abs(bt(q) - t)); use(q(i)) = true;
    end
    if ~any(use)
        % e.g. GL1 before Aug 2019 (no green scatter channel), or green
        % channels on AR31A (PT conversion without green signal)
        warning('bead_reference:noQC', 'no QC''d bead runs for %s', ch)
        ref.bead_log.(ch) = NaN; ref.bead_hv.(ch) = NaN; ref.k.(ch) = cfg.k.(ch); ref.ntimes.(ch) = 0; ref.runs.(ch) = strings(0, 1);
        continue
    end
    ref.runs.(ch) = string(beadstat.filename(use));
    v = v(use); hv = hv(use);
    hv0 = hv(1);
    ref.bead_log.(ch) = median(log10(v) - cfg.k.(ch) * log10(hv / hv0));
    ref.bead_hv.(ch) = hv0;
    ref.k.(ch) = cfg.k.(ch);
    ref.ntimes.(ch) = sum(use);
end
end
