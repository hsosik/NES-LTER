function mergeT = volume_merge_table(p, V, files, filetime, rows, keep_others)
% function mergeT = volume_merge_table(p, V, files, filetime, rows, keep_others)
% Writes mergeTable.mat in p.outpath: one row per file with the fields
% volume_file returned (scatter voltages, SSC/GL1 overlap statistics or the
% decompression flag count; preserved samples: voltages only, so cruise_status can
% see which files have volumes). Files whose volume failed (empty row) get NaN.
% keep_others true (classify on some files): rows of other files already in
% mergeTable.mat are kept.

files = cellstr(files(:)); filetime = filetime(:);
ok = ~cellfun(@(r) isempty(r) || isempty(fieldnames(r)), rows);
if any(ok)
    R = struct2table(vertcat(rows{ok}), 'AsArray', true);
    T = array2table(NaN(numel(files), width(R)), 'VariableNames', R.Properties.VariableNames);
    T(ok, :) = R;
else
    T = table;
end
mergeT = [table(files, filetime, 'VariableNames', {'filename' 'filetime'}) T];
fn = [p.outpath 'mergeTable.mat'];
if keep_others && exist(fn, 'file')
    try
        o = load(fn, 'mergeT');
        oT = o.mergeT(~ismember(o.mergeT.filename, mergeT.filename), :);
        if height(oT)
            for v = setdiff(mergeT.Properties.VariableNames, oT.Properties.VariableNames), oT.(v{1}) = NaN(height(oT), 1); end
            for v = setdiff(oT.Properties.VariableNames, mergeT.Properties.VariableNames), mergeT.(v{1}) = NaN(height(mergeT), 1); end
            mergeT = sortrows([oT(:, mergeT.Properties.VariableNames); mergeT], 'filetime');
        end
    catch ME
        warning('volume_merge_table:merge', 'earlier mergeTable not merged: %s', ME.message)
    end
end
if ismember(V.mode, {'merge' 'preserved'})
    GL1_merge_bead = V.GL1_merge_bead;
    save(fn, "mergeT", "GL1_merge_bead")
else
    OD2setting = V.OD2setting;
    save(fn, "mergeT", "OD2setting")
end
end
