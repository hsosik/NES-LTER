function B = collect_bead_runs(outfile)
% function B = collect_bead_runs(outfile)
% Gather every FCB bead run from all processed cruises into one table:
% cruise, time, detector voltage and log10 peak position (NaN unless that
% population passed QC) for each of the 3 bead populations in each height
% channel. Reads beadstat.mat via load_beadstat.
% Used to fit detector response vs voltage (beads are a constant target).
%
% Nov 2026, classifier redesign exploration

root = '\\sosiknas1\lab_data\Attune\cruise_data\';
chans = {'FSC' 'SSC' 'BL1' 'BL2' 'BL3' 'GL1' 'GL2' 'GL3' 'GL4'};
D = dir([root '20*_*']);
D = D([D.isdir]);
rows = {};
for i = 1:numel(D)
    bp = fullfile(root, D(i).name, 'bead_calibrated');
    try
        bs = load_beadstat(bp);
    catch
        continue
    end
    P = {'05' '1' '6'};
    for r = 1:height(bs)
        x = struct();
        x.cruise = string(extractAfter(D(i).name, '_'));
        x.time = datenum(bs.time(r));
        for ch = chans
            x.(['hv_' ch{1}]) = bs.hv.(ch{1})(r);
            for b = 1:3
                v = NaN;
                m = bs.(['peak' P{b}]).([ch{1} '_H'])(r);
                if bs.(['QC' P{b}])(r) && m > 0, v = log10(m); end
                x.(sprintf('bead%d_%s', b, ch{1})) = v;
            end
        end
        rows{end+1, 1} = x; %#ok<AGROW>
    end
end
B = struct2table([rows{:}]');
if exist('outfile', 'var') && ~isempty(outfile)
    save(outfile, 'B', 'chans')
end
end
