function [M, E] = merge_discrete_tables(cruises, outpath, opt)
% function [M, E] = merge_discrete_tables(cruises, outpath, Name=Value)
% Merge the per-Niskin preserved (discrete) sample tables of several cruises
% (DiscreteTable.mat from process_cruise(cruise, Sample="preserved")) into one
% table, and make the EDI table of the discrete product from it.
%
%   cruises   cruise names ("EN727"), paths to DiscreteTable.mat files, or "all"
%             (default): every <Root>\*\preserved\bead_calibrated\DiscreteTable.mat
%   outpath   folder for the outputs; "" (default) writes nothing
%
% Outputs
%   M   every column of every DiscreteTable (the union; missing ones NaN / "" / NaT),
%       one row per cruise, cast, niskin (underway discretes by sample_label) and sizefrac
%   E   EDI table. The columns of the published package (knb-lter-nes.29) come first,
%       in the same order and units, then the extra columns chosen with Extra.
%   Written to outpath (when given): Attune_Discrete_Table_v2.mat (M, E,
%   merge_metadata) and the csv named by CsvName.
%
% Mapping to the published columns (v2 class codes, see preserved_table)
%   syn_*                  Syn, class 2 (PE run, else CHL)
%   redeuk_leq_Xum_*       Euk without PE, class 1, ESD <= X um (CHL run, else PE)
%   redeuk_volume_analyzed_ml, redeuk_filename   the run used for euk
%   hetprok_cells_per_ml   SYBR run: Hetprok="incl_pro" class 7 + 9 (as the published
%                          numbers, whose bacteria gate included Pro) or
%                          "heterotrophs" class 9 only
%   hetprok_carbon_concentration   hetprok_cells_per_ml x 20 fg C per cell (Lee & Fuhrman 1987)
%   potential_temperature_c = potemp090c; date_processed = date sampled + storage days
%   of the euk run (else Syn, Pro, SYBR run), as a date
%   samplelog_QC           flag from the SampleLogQC csv (cruisename, cast, niskin, flag);
%                          NaN without it
%   station_distance       NaN (not in DiscreteTable yet)
%
% Options (Name=Value)
%   Root        folder of cruise folders (default \\sosiknas1\Lab_data\Attune\cruise_data\)
%   Hetprok     "incl_pro" (default, as published) or "heterotrophs"
%   Extra       groups of columns added after the published ones (default: all):
%               "pro"             pro_cells_per_ml, _biovolume_concentration, _carbon_concentration,
%                                 _volume_analyzed_ml, _filename, pro_lower_bound_cells_per_ml,
%                                 pro_lower_bound_filename
%               "crypto"          crypto_dim_* and crypto_bright_* (classes 3, 4; ESD <= 20 um)
%               "euk_with_syn"    euk_with_syn_* (class 5; ESD <= 20 um)
%               "euk_pe_shift"    euk_pe_shift_* (class 6: Euks whose PE channels are raised by
%                                 preservation; all sizes, already in the redeuk columns)
%               "redeuk_all"      redeuk_cells_per_ml etc.: class 1, all sizes
%               "syn_on_particle" syn_on_particle_cells_per_ml (class 8)
%               "syn_clump"       syn_clump_cells_per_ml (class 11: Syn clumps, not in syn)
%               "syn_types"       syn_typeA/B/C/D/ambiguous/unusable_cells_per_ml
%               "sybr"            pro_sybr_cells_per_ml and the other hetprok choice
%                                 (hetprok_excl_pro_ or hetprok_incl_pro_cells_per_ml)
%               "hna"             hetprok_lna_cells_per_ml, hetprok_hna_cells_per_ml (class 9
%                                 heterotrophs split by bact_type; sum = hetprok excluding Pro)
%                                 and hetprok_hna_method ("valley": cut at the run's SYBR
%                                 density valley; "fixed": no valley, HNA and LNA merged and
%                                 the split nominal)
%               "storage_days"    syn/redeuk/pro/hetprok_storage_days
%               "runs"            n_CHL_runs, n_PE_runs, n_pro_runs, n_hbac_runs
%               "sample"          sample_label, sizefrac
%               "processing"      processing: "v2 classes" or "hand gates" (Legacy)
%   Legacy      cruises taken from the old pipeline (<cruise>\preserved\outputs\EDI_table.mat)
%               instead of a DiscreteTable; their extra columns are NaN (default none)
%   SampleLogQC csv of sample-log flags (default "": samplelog_QC NaN)
%   CsvName     default "attune-transect-discrete-samples.csv"
%   AllowNAS    false (default): refuse to write into the EDI package folder or the
%               cruise_data root, where the published and merged tables live
%
% Nov 2026, process_cruise preserved mode

arguments
    cruises string = "all"
    outpath (1,1) string = ""
    opt.Root (1,1) string = string([char([92 92]) 'sosiknas1\Lab_data\Attune\cruise_data\'])
    opt.Hetprok (1,1) string {mustBeMember(opt.Hetprok, ["incl_pro" "heterotrophs"])} = "incl_pro"
    opt.Extra string = ["pro" "crypto" "euk_with_syn" "euk_pe_shift" "redeuk_all" "syn_on_particle" "syn_clump" "syn_types" ...
        "sybr" "hna" "storage_days" "runs" "sample" "processing"]
    opt.Legacy string = strings(0)
    opt.SampleLogQC (1,1) string = ""
    opt.CsvName (1,1) string = "attune-transect-discrete-samples.csv"
    opt.AllowNAS (1,1) logical = false
end
known = ["pro" "crypto" "euk_with_syn" "euk_pe_shift" "redeuk_all" "syn_on_particle" "syn_clump" "syn_types" "sybr" "hna" ...
    "storage_days" "runs" "sample" "processing"];
bad = setdiff(opt.Extra, known);
if ~isempty(bad), error('merge_discrete_tables:extra', 'unknown Extra "%s"; known: %s', bad(1), strjoin(known, ', ')), end
root = char(opt.Root); if ~endsWith(root, filesep), root = [root filesep]; end

% ---- DiscreteTables
files = find_tables(cruises, root, opt.Legacy);
parts = cell(numel(files), 1); notes = cell(numel(files), 1);
for k = 1:numel(files)
    S = load(files(k), 'T', 'table_metadata');
    parts{k} = S.T; notes{k} = [files(k); string(S.table_metadata(:))];
    fprintf('merge: %s %d rows (%s)\n', S.T.cruise(1), height(S.T), files(k))
end
M = stack(parts);
if ~isempty(M)
    key = M.cruise + "|" + M.cast + "|" + M.niskin + "|" + M.sample_label + "|" + M.sizefrac;
    [u, ~, j] = unique(key);
    dup = u(accumarray(j, 1) > 1);
    if ~isempty(dup), error('merge_discrete_tables:duplicate', 'duplicate samples, e.g. %s', dup(1)), end
    M = sortrows(M, {'cruise' 'cast' 'niskin' 'sample_label' 'sizefrac'});
end

% ---- EDI table
E = edi_table(M, opt);
if ~isempty(opt.Legacy)
    E = stack({E; legacy_table(opt.Legacy, root, E)});
end
if opt.SampleLogQC ~= "" && ~isempty(E)
    Q = readtable(opt.SampleLogQC, 'TextType', 'string');
    [ok, j] = ismember(upper(E.cruise) + "|" + E.cast + "|" + E.niskin, upper(Q.cruisename) + "|" + Q.cast + "|" + Q.niskin);
    E.samplelog_QC(ok) = Q.flag(j(ok));
end
if ~isempty(E), E = sortrows(E, {'date_sampled' 'cruise' 'cast' 'niskin'}); end

merge_metadata = {"merge_discrete_tables"; string(datetime('now')); ...
    "Hetprok = " + opt.Hetprok; "Extra = " + strjoin(opt.Extra, ', '); "Legacy = " + strjoin(opt.Legacy, ', '); ...
    "SampleLogQC = " + opt.SampleLogQC; "sources, each followed by its table_metadata:"};
merge_metadata = [merge_metadata; vertcat(notes{:})];

% ---- write
if outpath == "", return, end
out = char(outpath); if ~endsWith(out, filesep), out = [out filesep]; end
lo = lower(strrep(out, '/', '\'));
if ~opt.AllowNAS && (contains(lo, 'edi_data_packages') || strcmp(lo, lower(strrep(root, '/', '\'))))
    error('merge_discrete_tables:protected', ['%s holds the published or merged discrete table; ' ...
        'pass AllowNAS=true to write there'], out)
end
if ~exist(out, 'dir'), mkdir(out), end
save([out 'Attune_Discrete_Table_v2.mat'], 'M', 'E', 'merge_metadata')
writetable(csv_ready(E), [out char(opt.CsvName)])
fprintf('merge: %d samples in %d cruises -> %s\n', height(E), numel(unique(E.cruise)), out)
end


function files = find_tables(cruises, root, legacy)
sub = fullfile('preserved', 'bead_calibrated', 'DiscreteTable.mat');
if isscalar(cruises) && lower(cruises) == "all"
    d = dir(fullfile(root, '*', sub));
    files = string(fullfile({d.folder}, {d.name}))';
    if ~isempty(legacy) && ~isempty(files) % a Legacy cruise comes from its old table only
        c = regexp(files, '_([^_\\/]+)[\\/]preserved', 'tokens', 'once');
        c = upper(cellfun(@(x) string(x{1}), c));
        files = files(~ismember(c, upper(legacy)));
    end
    return
end
files = strings(0, 1);
for c = cruises(:)'
    if endsWith(c, ".mat", 'IgnoreCase', true)
        f = c;
    else
        f = string(fullfile(cruise_folder(c, root), sub));
    end
    if ~isfile(f), error('merge_discrete_tables:missing', 'no DiscreteTable for %s (%s)', c, f), end
    files(end+1, 1) = f; %#ok<AGROW>
end
end


function f = cruise_folder(c, root)
d = dir(fullfile(root, ['*_' char(c)]));
d = d([d.isdir]);
if numel(d) ~= 1, error('merge_discrete_tables:folder', 'found %d folders for cruise %s in %s', numel(d), c, root), end
f = fullfile(d.folder, d.name);
end


function E = edi_table(M, opt)
% published columns in published order, then extras
n = height(M);
E = table();
E.cruise = col(M, 'cruise', n);
E.cast = col(M, 'cast', n);
E.niskin = col(M, 'niskin', n);
E.latitude = col(M, 'latitude', n);
E.longitude = col(M, 'longitude', n);
E.depth_m = col(M, 'depth_m', n);
E.salinity = col(M, 'salinity', n);
E.potential_temperature_c = col(M, 'potemp090c', n);
E.date_sampled = col(M, 'date_sampled', n, NaT(n, 1));
E.date_processed = processed_date(M, E.date_sampled);
for v = ["cells_per_ml" "biovolume_concentration" "carbon_concentration"]
    E.("syn_" + v) = col(M, "syn_" + v, n);
end
E.syn_volume_analyzed_ml = col(M, 'syn_volume_analyzed_ml', n);
E.syn_filename = col(M, 'syn_filename', n, strings(n, 1));
for x = [2 3 5 10 20]
    for v = ["cells_per_ml" "biovolume_concentration" "carbon_concentration"]
        E.(sprintf('redeuk_leq_%dum_%s', x, v)) = col(M, sprintf('euk_leq%dum_%s', x, v), n);
    end
end
E.redeuk_volume_analyzed_ml = col(M, 'euk_volume_analyzed_ml', n);
E.redeuk_filename = col(M, 'euk_filename', n, strings(n, 1));
incl = col(M, 'bacteria_incl_pro_cells_per_ml', n); het = col(M, 'bacteria_cells_per_ml', n);
if opt.Hetprok == "incl_pro", E.hetprok_cells_per_ml = incl; else, E.hetprok_cells_per_ml = het; end
E.hetprok_carbon_concentration = E.hetprok_cells_per_ml * 20e-6; % 20 fg C per cell -> ug/L
E.hetprok_volume_analyzed_ml = col(M, 'hbac_volume_analyzed_ml', n);
E.hetprok_filename = col(M, 'hbac_filename', n, strings(n, 1));
E.samplelog_QC = NaN(n, 1);
E.nearest_station = col(M, 'nearest_station', n, strings(n, 1));
E.station_distance = col(M, 'station_distance', n);

ex = opt.Extra;
if ismember("pro", ex)
    for v = ["cells_per_ml" "biovolume_concentration" "carbon_concentration" "volume_analyzed_ml"]
        E.("pro_" + v) = col(M, "pro_" + v, n);
    end
    E.pro_filename = col(M, 'pro_filename', n, strings(n, 1));
    E.pro_lower_bound_cells_per_ml = col(M, 'pro_lower_bound_cells_per_ml', n);
    E.pro_lower_bound_filename = col(M, 'pro_lower_bound_filename', n, strings(n, 1));
end
grp = {"crypto", ["crypto_dim" "crypto_bright"]; "euk_with_syn", "euk_with_syn"; "euk_pe_shift", "euk_pe_shift"; "redeuk_all", "euk"};
for g = 1:size(grp, 1)
    if ~ismember(grp{g, 1}, ex), continue, end
    for nm = grp{g, 2}
        for v = ["cells_per_ml" "biovolume_concentration" "carbon_concentration"]
            e = nm; if nm == "euk", e = "redeuk"; end
            E.(e + "_" + v) = col(M, nm + "_" + v, n);
        end
    end
end
if ismember("syn_on_particle", ex), E.syn_on_particle_cells_per_ml = col(M, 'syn_on_particle_cells_per_ml', n); end
if ismember("syn_clump", ex), E.syn_clump_cells_per_ml = col(M, 'syn_clump_cells_per_ml', n); end
if ismember("syn_types", ex)
    for t = ["A" "B" "C" "D" "_ambiguous" "_unusable"]
        E.("syn_type" + t + "_cells_per_ml") = col(M, "syn_type" + t + "_cells_per_ml", n);
    end
end
if ismember("sybr", ex)
    E.pro_sybr_cells_per_ml = col(M, 'pro_sybr_cells_per_ml', n);
    if opt.Hetprok == "incl_pro", E.hetprok_excl_pro_cells_per_ml = het; else, E.hetprok_incl_pro_cells_per_ml = incl; end
end
if ismember("hna", ex)
    E.hetprok_lna_cells_per_ml = col(M, 'bacteria_LNA_cells_per_ml', n);
    E.hetprok_hna_cells_per_ml = col(M, 'bacteria_HNA_cells_per_ml', n);
    E.hetprok_hna_method = col(M, 'bacteria_hna_method', n, strings(n, 1));
end
if ismember("storage_days", ex)
    for g = ["syn" "euk" "pro" "hbac"; "syn" "redeuk" "pro" "hetprok"]
        E.(g(2) + "_storage_days") = col(M, g(1) + "_storage_days", n);
    end
end
if ismember("runs", ex)
    for r = ["CHL" "PE" "pro" "hbac"], E.("n_" + r + "_runs") = col(M, "n_" + r + "_runs", n); end
end
if ismember("sample", ex)
    E.sample_label = col(M, 'sample_label', n, strings(n, 1));
    E.sizefrac = col(M, 'sizefrac', n, strings(n, 1));
end
if ismember("processing", ex), E.processing = repmat("v2 classes", n, 1); end
end


function x = col(M, name, n, empty)
% column of M, or the empty value (NaN by default) when M lacks it
if nargin < 4, empty = NaN(n, 1); end
if ismember(char(name), M.Properties.VariableNames), x = M.(char(name)); else, x = empty; end
if iscellstr(x) || ischar(x), x = string(x); end %#ok<ISCLSTR>
end


function d = processed_date(M, ts)
% run date from date_sampled + storage days of the euk run, else syn, pro, SYBR
n = height(M); sd = NaN(n, 1);
for g = ["euk" "syn" "pro" "hbac"]
    s = col(M, g + "_storage_days", n);
    sd(isnan(sd)) = s(isnan(sd));
end
d = dateshift(ts + days(sd), 'start', 'day');
d.Format = 'yyyy-MM-dd';
end


function L = legacy_table(cruises, root, E)
% old pipeline EDI_table.mat (preserved\outputs, or preserved2 for EN655), as published
parts = cell(numel(cruises), 1);
for k = 1:numel(cruises)
    base = cruise_folder(cruises(k), root);
    f = fullfile(base, 'preserved', 'outputs', 'EDI_table.mat');
    if ~isfile(f), f = fullfile(base, 'preserved2', 'outputs', 'EDI_table.mat'); end
    S = load(f, 'EDI_table'); T = S.EDI_table;
    T.cruise = upper(string(T.cruise));
    for v = string(T.Properties.VariableNames)
        if iscell(T.(v)), T.(v) = string(T.(v)); T.(v)(T.(v) == "NaN") = ""; end
    end
    if any(T.syn_volume_analyzed_ml > 10 | T.redeuk_volume_analyzed_ml > 10)
        warning('merge_discrete_tables:legacyvol', ['%s: volume analyzed > 10 ml in the old table ' ...
            '(microlitres?); concentrations are suspect'], cruises(k))
    end
    T.processing = repmat("hand gates", height(T), 1);
    parts{k} = T;
    fprintf('merge: legacy %s %d rows (%s)\n', cruises(k), height(T), f)
end
L = stack(parts);
L = L(:, intersect(L.Properties.VariableNames, E.Properties.VariableNames, 'stable'));
end


function T = stack(parts)
% vertcat tables with different columns: missing ones NaN, "", NaT, or 0 for logicals
parts = parts(~cellfun(@isempty, parts));
if isempty(parts), T = table(); return, end
names = {};
for k = 1:numel(parts), names = [names setdiff(parts{k}.Properties.VariableNames, names, 'stable')]; end %#ok<AGROW>
for k = 1:numel(parts)
    P = parts{k}; n = height(P);
    for v = setdiff(names, P.Properties.VariableNames)
        ref = [];
        for q = 1:numel(parts), if ismember(v{1}, parts{q}.Properties.VariableNames), ref = parts{q}.(v{1}); break, end, end
        if isdatetime(ref), P.(v{1}) = NaT(n, 1);
        elseif isstring(ref) || iscell(ref), P.(v{1}) = strings(n, 1);
        else, P.(v{1}) = NaN(n, 1);
        end
    end
    for v = names, if iscell(P.(v{1})), P.(v{1}) = string(P.(v{1})); end, end
    if ismember('date_sampled', names) && isdatetime(P.date_sampled), P.date_sampled.TimeZone = ''; end
    parts{k} = P(:, names);
end
T = vertcat(parts{:});
end


function E = csv_ready(E)
% formats for the csv: lat/lon 4 decimals, other numbers 6 significant digits,
% dates as in the published package, missing text "NaN" (its missing value code)
for v = string(E.Properties.VariableNames)
    x = E.(v);
    if isstring(x)
        x(ismissing(x) | x == "") = "NaN"; E.(v) = x; % published missing value code
    elseif isdatetime(x)
        if v == "date_processed", x.Format = 'yyyy-MM-dd'; else, x.Format = 'yyyy-MM-dd HH:mm:ss'; end
        E.(v) = x;
    elseif isfloat(x)
        if ismember(v, ["latitude" "longitude"]), E.(v) = round(x, 4);
        else, E.(v) = round(x, 6, 'significant');
        end
    end
end
end
