function out = perf_beadstats(cruises, V, R, outdir)
% function out = perf_beadstats(cruises, V, R, outdir)
% FCB-equivalent bead tables for cruises without FCB bead runs, from the
% performance tests (pt_perf_table) run during each cruise:
%   log10 FCB peak (channel, population) = PT perf peak + offset
% at the performance test's own voltages (stored in hv, so bead_reference
% scales to each sample's voltages with cfg.k as usual). offset = median of
% (FCB peak - PT perf peak - k log10(hv_run / hv_perf)) over FCB bead runs
% paired with performance tests in the same optical state
% (validate_perf_transfer: layout/OD2 placement).
%
% Output tables have the beadstat format (process_beads_PT_adjust_3 /
% load_beadstat), one row per good performance test from 3 days before the
% cruise's first sample file to its last; meta.source = "PT performance test
% transfer". They are written to outdir, NOT to the cruise's bead_calibrated
% folder.
%
% Caveat (meta.notes, row field filterset): Jan - early Apr 2018 and early
% Jul 2018 (BL3 perf voltage ~388 V, physical GL4 ~368 V) had BL3 695/40 and
% green 695/40 filters (settings history: EN608, AR28B). The offsets were
% measured with 680/40, so BL3 and GL4 assume the FCB/PT ratio is the same
% through both filters.
%
% cruises - cruise folder names (default: EN608, AR28B, AR29, EN617, AR31A, EN627)
% V       - validate_perf_transfer output (default: load ../output/perf_transfer_validation.mat)
% R       - pt_perf_table output (default: load ../output/pt_perf.mat)
% outdir  - default ../output/pt_transfer
%
% Nov 2026, PT transfer for early cruises

here = fileparts(mfilename('fullpath'));
dev = fileparts(here);
addpath(dev, fileparts(dev))
if nargin < 1 || isempty(cruises)
    cruises = {'20180131_EN608' '20180404_AR28B' '20180414_AR29' '20180720_EN617' '20181020_AR31A' '20190201_EN627'};
end
if nargin < 2 || isempty(V), s = load(fullfile(dev, 'output', 'perf_transfer_validation.mat'), 'V'); V = s.V; end
if nargin < 3 || isempty(R), s = load(fullfile(dev, 'output', 'pt_perf.mat'), 'R'); R = s.R; end
if nargin < 4 || isempty(outdir), outdir = fullfile(dev, 'output', 'pt_transfer'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end
S = V.Properties.CustomProperties.summary;
R = R(R.ok, :);
R.state = R.layout + "/" + R.od2;
root = ['\' '\sosiknas1\lab_data\Attune\cruise_data\'];
chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
names = [strcat(chans, '_H') {'SSC_A' 'GL1_A'}];
P = {'05' '1' '6'};
for ic = 1:numel(cruises)
    cdir = cruises{ic}; cruise = extractAfter(cdir, '_');
    [t0, t1] = cruise_window(fullfile(root, cdir, 'FCS'));
    t0 = t0 - days(3);
    conf = attune_config(cruise);
    % only tests in the cruise's own optical state (e.g. EN627: not the
    % tests from before the OD2 filter came off SSC on 1 Feb 2019)
    Q = R(R.time >= t0 & R.time <= t1 & R.od2 == string(conf.OD2), :);
    beadstat = table;
    for r = 1:height(Q)
        x = table;
        x.filename = Q.file(r);
        x.time = Q.time(r);
        hv = struct();
        for c = chans, hv.(c{1}) = Q.(['hv_' c{1}])(r); end
        x.hv = struct2table(hv);
        x.state = Q.state(r);
        x.filterset = "680/40";
        % 695/40 needs a higher BL3 voltage relative to BL1 (BL3/BL1 ~0.96 vs ~0.91;
        % BL1 drifts with the blue laser, so use the ratio, not BL3 alone)
        if Q.layout(r) == "v2019" && Q.hv_BL3(r) / Q.hv_BL1(r) > 0.935, x.filterset = "695/40"; end
        for ip = 1:3
            pk = nan(1, numel(chans));
            for k = 1:numel(chans)
                o = S.offset(S.chan == chans{k} & S.pop == P{ip} & S.state == Q.state(r));
                if ~isempty(o), pk(k) = Q.(['pk_' chans{k}])(r) + o; end
            end
            x.(['n' P{ip}]) = Q.nsing(r);
            x.(['QC' P{ip}]) = isfinite(pk(2)); % SSC offset available for this state
            x.(['why' P{ip}]) = "";
            if ~x.(['QC' P{ip}]), x.(['why' P{ip}]) = "no offset for state " + Q.state(r); end
            x.(['peak' P{ip}]) = array2table(10 .^ [pk NaN NaN], 'VariableNames', names);
            x.(['med' P{ip}]) = array2table(nan(1, numel(names)), 'VariableNames', names);
            x.(['rsd' P{ip}]) = array2table(nan(1, numel(names)), 'VariableNames', names);
        end
        beadstat = [beadstat; x]; %#ok<AGROW>
    end
    meta = struct();
    meta.source = "PT performance test transfer";
    meta.config = attune_config(cruise);
    meta.offsets = S;
    meta.cruise = string(cruise);
    meta.processed = datetime('now');
    meta.code = "perf_beadstats (Nov 2026)";
    meta.notes = ["FCB-equivalent peaks = PT performance-test peak + (FCB - PT perf) offset per channel, population and optical state, from FCB bead runs paired with performance tests (AR34B, RB1904: no OD2; TN368: OD2 on SSC)";
        "hv = the performance test's own voltages; bead_reference scales to sample voltages with cfg.k";
        "filterset 695/40: BL3 and green-chl (GL4) filters differed from the 680/40 used for the offsets; FCB/PT ratio assumed unchanged"];
    out.(matlab.lang.makeValidName(cruise)) = beadstat;
    fprintf('%s: %d performance tests (%s), filtersets %s\n', cruise, height(beadstat), ...
        strjoin(unique(beadstat.state), '/'), strjoin(unique(beadstat.filterset), '/'))
    meta.window = [t0 + days(3) t1];
    if height(beadstat) > 0
        save(fullfile(outdir, [cruise '_beadstat.mat']), 'beadstat', 'meta')
    end
end
end

function [t0, t1] = cruise_window(fcsdir)
% acquisition time range from the $DATE/$BTIM of up to 60 sample files spread
% through the folder (file modification dates are not reliable: some cruise
% folders were copied later)
F = dir(fullfile(fcsdir, '*.fcs'));
F = F(~contains({F.name}, {'bead' 'Bead' 'test' 'Test'}));
pick = unique(round(linspace(1, numel(F), min(60, numel(F)))));
t = NaT(numel(pick), 1);
for i = 1:numel(pick)
    fid = fopen(fullfile(fcsdir, F(pick(i)).name));
    txt = char(fread(fid, 20000, 'uint8')'); fclose(fid);
    d = regexp(txt, '\$DATE/([^/]*)/', 'tokens', 'once');
    b = regexp(txt, '\$BTIM/([^/]*)/', 'tokens', 'once');
    if isempty(d) || isempty(b), continue, end
    try, t(i) = datetime([d{1} ' ' b{1}], 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'); catch, end
end
t0 = min(t); t1 = max(t);
end
