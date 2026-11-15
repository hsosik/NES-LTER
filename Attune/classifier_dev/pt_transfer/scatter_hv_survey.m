function H = scatter_hv_survey(nfiles)
% function H = scatter_hv_survey(nfiles)
% SSC and GL1 voltages used on each cruise with the current filter layout
% (folders dated 2019-08-19 on), from the FCS headers ($PnN/$PnV) of up to
% nfiles sample files and all bead files per cruise. Reads the header text
% only (fast). Needed because the SSC/GL1 merge window (cfg.gl1_merge) is set
% in raw GL1 counts, so it moves in bead units when the GL1 voltage changes.
%
% H - one row per cruise x file type (sample/bead) x (SSC hv, GL1 hv) with
%     the number of files
%
% Nov 2026, PT transfer / volume calibration

if nargin < 1 || isempty(nfiles), nfiles = 60; end
root = ['\' '\sosiknas1\lab_data\Attune\cruise_data\'];
D = dir(root);
D = D([D.isdir] & ~cellfun(@isempty, regexp({D.name}, '^\d{8}_', 'once')));
D = D(cellfun(@(n) str2double(n(1:8)), {D.name}) >= 20190819);
rows = {};
for i = 1:numel(D)
    F = dir(fullfile(root, D(i).name, 'FCS', '*.fcs'));
    if isempty(F), continue, end
    isb = contains({F.name}, {'bead' 'Bead'});
    S = F(~isb); B = F(isb);
    pick = S(unique(round(linspace(1, numel(S), min(nfiles, numel(S))))));
    if isempty(S), pick = S; end
    for typ = ["sample" "bead"]
        if typ == "sample", L = pick; else, L = B; end
        for k = 1:numel(L)
            fid = fopen(fullfile(L(k).folder, L(k).name));
            txt = char(fread(fid, 20000, 'uint8')'); fclose(fid);
            hv = struct('SSC', NaN, 'GL1', NaN);
            for c = ["SSC" "GL1"]
                n = regexp(txt, ['\$P(\d+)N/' char(c) '-H/'], 'tokens', 'once');
                if isempty(n), continue, end
                v = regexp(txt, ['\$P' n{1} 'V/([^/]*)/'], 'tokens', 'once');
                if ~isempty(v), hv.(c) = str2double(v{1}); end
            end
            rows{end+1, 1} = struct('cruise', string(extractAfter(D(i).name, '_')), 'folder', string(D(i).name), ...
                'type', typ, 'hv_SSC', hv.SSC, 'hv_GL1', hv.GL1); %#ok<AGROW>
        end
    end
end
R = struct2table([rows{:}]');
H = groupcounts(R, {'folder', 'cruise', 'type', 'hv_SSC', 'hv_GL1'});
H = H(:, {'folder', 'cruise', 'type', 'hv_SSC', 'hv_GL1', 'GroupCount'});
end
