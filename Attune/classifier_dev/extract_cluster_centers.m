function T = extract_cluster_centers(nPerCruise, cruiseFilter, outfile)
% function T = extract_cluster_centers(nPerCruise, cruiseFilter, outfile)
% For a sample of files from each processed cruise, record the median
% position (log10, height channels) of existing Syn / Euk / Pro classes,
% the detector voltages, the nearest-in-time FCB bead peak positions, and
% time/location. Used to test whether voltage or bead normalization makes
% cluster positions comparable across cruises.
%
% nPerCruise - number of standard files per cruise (Pro-setting files get half that)
% cruiseFilter - optional cellstr of cruise names to restrict to (e.g., {'AR92'})
% outfile - optional .mat path for results
%
% Nov 2026, classifier redesign exploration

if ~exist('nPerCruise', 'var') || isempty(nPerCruise), nPerCruise = 12; end
if ~exist('cruiseFilter', 'var'), cruiseFilter = {}; end
addpath(fullfile(fileparts(mfilename('fullpath')), '..'))  %fca_readfcs

root = '\\sosiknas1\lab_data\Attune\cruise_data\';
chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
classes = struct('name', {'syn' 'euk' 'pro'}, 'code', {2 1 7});
exclude = {'bead' 'Bead' 'Cast' 'cast' 'test' 'Dockwater' 'discrete' 'Rinses' 'Filter config' 'Cultures' 'Grazer' 'Dilution'};

D = dir([root '20*_*']);
D = D([D.isdir]);
rows = {};
for i = 1:numel(D)
    cruise = extractAfter(D(i).name, '_');
    if ~isempty(cruiseFilter) && ~ismember(cruise, cruiseFilter), continue, end
    bp = fullfile(root, D(i).name, 'bead_calibrated');
    if ~exist(fullfile(bp, 'FCSfileinfo.mat'), 'file') || ~exist(fullfile(bp, 'class'), 'dir'), continue, end
    try
        cruiseRows = do_cruise(cruise, fullfile(root, D(i).name), bp, nPerCruise, chans, classes, exclude);
        rows = [rows; cruiseRows]; %#ok<AGROW>
        fprintf('%s: %d files\n', cruise, numel(cruiseRows));
    catch ME
        fprintf('%s: FAILED %s\n', cruise, ME.message);
    end
end
T = struct2table([rows{:}]', 'AsArray', true);
if exist('outfile', 'var') && ~isempty(outfile)
    save(outfile, 'T', 'chans')
end
end

function rows = do_cruise(cruise, cpath, bp, nPerCruise, chans, classes, exclude)
load(fullfile(bp, 'FCSfileinfo.mat'), 'FCSfileinfo')
F = FCSfileinfo;
F = F(~contains(F.fcslist, exclude) & F.QC_flag == 1, :);
classfiles = regexprep(F.fcslist, '\.fcs$', '.mat');
F = F(cellfun(@(c) exist(fullfile(bp, 'class', c), 'file') == 2, classfiles), :);
F = sortrows(F, 'matdate_start');
if ~ismember('pro_measured', F.Properties.VariableNames), F.pro_measured = zeros(height(F), 1); end

pick = [];
for pm = [0 1]
    ii = find(F.pro_measured == pm);
    n = nPerCruise / (1 + pm);
    if ~isempty(ii)
        pick = [pick; ii(unique(round(linspace(1, numel(ii), min(n, numel(ii))))))]; %#ok<AGROW>
    end
end

beads = load_beads(bp);
uw = load_uw(bp);

rows = {};
for k = pick'
    fcsname = F.fcslist{k};
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fcsname));
    c = load(fullfile(bp, 'class', regexprep(fcsname, '\.fcs$', '.mat')), 'class');
    cls = c.class(:, 1);
    pn = {fcshdr.par.name};

    r = struct();
    r.cruise = string(cruise);
    r.file = string(fcsname);
    r.time = F.matdate_start(k);
    r.pro_measured = F.pro_measured(k);
    r.trigger = string(fcshdr.tr1_par);
    r.trigger_level = fcshdr.tr1_level;
    r.nevents = size(fcsdat, 1);
    for ch = chans
        j = strcmp(pn, [ch{1} '-H']);
        r.(['hv_' ch{1}]) = nan_if_empty([fcshdr.par(j).hv]);
    end
    for cc = classes
        idx = cls == cc.code;
        r.(['n_' cc.name]) = sum(idx);
        for ch = chans
            j = strcmp(pn, [ch{1} '-H']);
            v = NaN;
            if any(j) && sum(idx) >= 20
                x = fcsdat(idx, j); x = x(x > 0);
                if numel(x) >= 20, v = median(log10(x)); end
            end
            r.([cc.name '_' ch{1}]) = v;
        end
    end

    % nearest bead run in time
    r.bead_dt_days = NaN;
    for b = 1:3
        for ch = chans, r.(sprintf('bead%d_%s', b, ch{1})) = NaN; end
    end
    for ch = chans, r.(['beadhv_' ch{1}]) = NaN; end
    if ~isempty(beads)
        % nearest bead run in time with that population QC'd (load_beadstat layout)
        P = {'05' '1' '6'};
        [dt, bi] = min(abs(datenum(beads.time) - r.time));
        r.bead_dt_days = dt;
        for ch = chans
            for b = 1:3
                q = find(beads.(['QC' P{b}]));
                if isempty(q), continue, end
                [~, k] = min(abs(datenum(beads.time(q)) - r.time));
                val = beads.(['peak' P{b}]).([ch{1} '_H'])(q(k));
                if val > 0, r.(sprintf('bead%d_%s', b, ch{1})) = log10(val); end
            end
            r.(['beadhv_' ch{1}]) = beads.hv.(ch{1})(bi);
        end
    end

    % location and temperature
    r.lat = NaN; r.lon = NaN; r.temp = NaN;
    if ~isempty(uw)
        u = find(strcmp(uw.file, fcsname), 1);
        if ~isempty(u), r.lat = uw.lat(u); r.lon = uw.lon(u); r.temp = uw.temp(u); end
    end
    rows{end+1, 1} = r; %#ok<AGROW>
end
end

function beads = load_beads(bp)
beads = [];
try
    beads = load_beadstat(bp);
catch
end
end
