function R = pt_perf_table(folder, outfile)
% function R = pt_perf_table(folder, outfile)
% Table of Attune performance tests (automated PT-bead setup runs exported as
% FCS, "GetFinalStats_FlowRate100Iteration1*.fcs"). In each test the Attune
% sets every detector voltage so that the PT bead reaches a fixed target, so
% the voltage records the channel's sensitivity at that time.
%
% For each test: time (#PTDATE + $BTIM; these files have no $DATE), laser
% delays and area scaling factors, number of singlet beads, and per physical
% channel (current-layout names, as in beadstat; attune_config gives the
% filter layout for the date): voltage hv, log10 PT peak (FSC/SSC: the single
% scatter peak; fluorescence: the brightest PT intensity level), robust CV.
% od2 = 'SSC' when the SSC voltage shows the OD2 filter was on SSC (> 330 V
% before the Aug 2019 layout change), else as attune_config.
%
% folder  - default \\sosiknas1\Lab_data\Attune\performance_test_fcs
% outfile - optional .mat to save R in
%
% Nov 2026, PT transfer for early cruises

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
if nargin < 1 || isempty(folder), folder = ['\' '\sosiknas1\Lab_data\Attune\performance_test_fcs']; end
F = dir(fullfile(folder, '*.fcs'));
chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
rows = cell(numel(F), 1);
for k = 1:numel(F)
    fn = fullfile(F(k).folder, F(k).name);
    try
        [d, h] = fca_readfcs(fn);
    catch
        continue
    end
    fid = fopen(fn); txt = char(fread(fid, 8192, 'uint8')'); fclose(fid);
    kw = @(n) regexp(txt, [n '/([^/]*)/'], 'tokens', 'once');
    d0 = kw('#PTDATE');
    if isempty(d0), continue, end
    r = struct();
    r.file = string(F(k).name);
    r.time = datetime([d0{1} ' ' h.starttime], 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    r.delay1 = str2double(kw('#LASER1DELAY')); r.delay2 = str2double(kw('#LASER2DELAY'));
    r.asf1 = str2double(kw('#LASER1ASF')); r.asf2 = str2double(kw('#LASER2ASF'));
    names = {h.par.name};
    conf = attune_config('', datenum(r.time));
    r.layout = string(conf.layout);
    % singlets: central cluster on FSC-H x SSC-H
    fs = log10(max(d(:, strcmp(names, 'FSC-H')), 1));
    ss = log10(max(d(:, strcmp(names, 'SSC-H')), 1));
    g = abs(fs - median(fs)) < 0.15 & abs(ss - median(ss)) < 0.15;
    r.n = size(d, 1); r.nsing = sum(g);
    for j = 1:numel(chans)
        ch = chans{j}; p = conf.chan.(ch);
        c = find(strcmp(names, p), 1);
        if isempty(p) || isempty(c)
            r.(['hv_' ch]) = NaN; r.(['pk_' ch]) = NaN; r.(['rcv_' ch]) = NaN;
            continue
        end
        r.(['hv_' ch]) = h.par(c).hv;
        x = d(g, c);
        if ismember(ch, {'FSC' 'SSC'})
            v = x;
        else % brightest PT intensity level
            lx = log10(max(x, 1));
            v = x(lx > prctile(lx, 70));
        end
        r.(['pk_' ch]) = log10(median(v));
        r.(['rcv_' ch]) = iqr(v) / 1.349 / median(v);
    end
    if conf.layout == "v2019" && r.hv_SSC > 330
        r.od2 = "SSC";
    elseif conf.layout == "v2019"
        r.od2 = "None";
    else
        r.od2 = "GL1";
    end
    rows{k} = r;
end
R = struct2table([rows{~cellfun(@isempty, rows)}]);
R = sortrows(R, 'time');
% failed or aborted tests: too few singlets, broad SSC, or channels far off target
R.ok = R.nsing > 2500 & R.rcv_SSC < 0.07 & abs(R.pk_SSC - median(R.pk_SSC, 'omitnan')) < 0.05;
if nargin > 1 && ~isempty(outfile), save(outfile, 'R'); end
end
