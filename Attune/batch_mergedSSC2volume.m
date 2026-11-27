function batch_mergedSSC2volume(p)
% function batch_mergedSSC2volume(p)
% process_cruise volume step (underway): cell volumes for every class file in
% p.classpath, added to the class files, and mergeTable.mat (per-file scatter
% voltages and SSC/GL1 overlap statistics) in p.outpath.
% Settings and bead positions: volume_setup; per-file calculation: volume_file
% (also used by the classify step, which computes volumes as it classifies when
% the classify and volume steps are run together).
% Files without a bead position (e.g. no QC'd bead run at the file's SSC
% voltage) are skipped with a warning.

load([p.outpath 'FCSfileinfo'], 'FCSfileinfo')
fcslist = dir([p.classpath '*.mat']); %just get the ones in the class folder (not the beads, etc.)
fcslist = regexprep({fcslist.name}', '.mat', '.fcs');
[~, ia] = ismember(fcslist, FCSfileinfo.fcslist);
filetime = datetime(FCSfileinfo.matdate_start(ia), 'ConvertFrom', 'datenum');

V = volume_setup(p, min(filetime));
n = numel(fcslist);
rows = cell(n, 1); nskip = 0;
for filecount = 1:n
    if ~rem(filecount, 50)
        disp(['volume calculating ' num2str(filecount) ' of ' num2str(n)])
    end
    cf = [p.classpath regexprep(fcslist{filecount}, '.fcs', '.mat')];
    try
        [fcsdat, fcshdr] = fca_readfcs([p.fpath, fcslist{filecount}]);
        c = load(cf);
        [c, rows{filecount}] = volume_file(c, fcsdat, fcshdr, V, filetime(filecount));
        save(cf, '-struct', "c")
    catch ME
        nskip = nskip + 1;
        warning('batch_mergedSSC2volume:file', '%s: volume not calculated: %s', fcslist{filecount}, ME.message)
    end
end
volume_merge_table(p, V, fcslist, filetime, rows, false);
if nskip, fprintf('volume: %d of %d files skipped\n', nskip, n); end
end
