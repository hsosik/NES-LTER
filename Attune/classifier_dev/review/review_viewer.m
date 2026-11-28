function review_viewer(outdir, Q, R, cruise, decfile, build)
% function review_viewer(outdir, Q, R, cruise, decfile, build)
% Write the static cruise review viewer into outdir: index.html (copied from
% review_viewer_template.html) and review_data.js (settings and per-file QC).
% Frames are expected in outdir/frames/. Q is a cell array of structs from
% review_qc, each with added fields file, time (datenum), frame, nrow.
% Open index.html in a browser; it works from the network share (no server).
% decfile: the review decisions file (default review_decisions.csv in the folder above
% outdir, i.e. bead_calibrated). It is never written here: its current contents are
% copied into review_data.js so earlier decisions show after a rebuild, and the viewer
% saves new decisions straight to it (Chrome or Edge; see review_decisions).
% build: label of the classification run shown (default: now); decisions record it.
% R.keep_order true: files stay in the order of Q (default: sorted by time), e.g. the
% held-out class views (holdout_class_view), which follow their A/B page's order.
% R.link (optional): {href, text} link shown next to the title.
%
% Nov 2026, classification review prototype

if nargin < 3 || isempty(R), R = review_config(); end
outdir = regexprep(char(outdir), '[\/]+$', '');
if nargin < 5 || isempty(decfile), decfile = fullfile(fileparts(outdir), 'review_decisions.csv'); end
if nargin < 6 || isempty(build), build = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')); end
if ~(isfield(R, 'keep_order') && R.keep_order)
    t = cellfun(@(q) q.time, Q);
    [~, o] = sort(t); Q = Q(o);
end
for i = 1:numel(Q)
    Q{i}.time = char(datetime(Q{i}.time, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
end
V = rmfield(R, intersect({'views' 'view_lim' 'type_views'}, fieldnames(R)));
V.views = cell(size(R.views, 1), 1);
for v = 1:size(R.views, 1)
    lx = R.lim.(R.views{v, 1}); ly = R.lim.(R.views{v, 2});
    if ~isempty(R.view_lim{v}), lx = R.view_lim{v}(1, :); ly = R.view_lim{v}(2, :); end
    tp = ''; if isfield(R, 'type_views'), tp = R.type_views{v}; end % Syn type view overlay: '', 'diag' or 'vbox'
    V.views{v} = struct('x', R.views{v, 1}, 'y', R.views{v, 2}, 'xlab', R.views{v, 3}, 'ylab', R.views{v, 4}, 'xlim', lx, 'ylim', ly, 'type', tp);
end
D = struct('cruise', cruise, 'made', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), 'R', V);
D.build = build;
[dd, dn, de] = fileparts(decfile); [~, ddn] = fileparts(dd);
D.decisions_path = decfile; D.decisions_name = [dn de]; D.decisions_dir = dd; D.decisions_dir_name = ddn;
D.decisions = num2cell(table2struct(review_decisions(decfile, true))); % cell: always a JSON array
D.files = Q;
fid = fopen(fullfile(outdir, 'review_data.js'), 'w', 'n', 'UTF-8');
fprintf(fid, 'const REVIEW = %s;\n', jsonencode(D));
fclose(fid);
% index.html: the template with a version stamp on review_data.js, so browsers don't
% show a cached copy after a rerun
html = fileread(fullfile(fileparts(mfilename('fullpath')), 'review_viewer_template.html'));
html = strrep(html, '__STAMP__', char(datetime('now', 'Format', 'yyyyMMddHHmmss')));
fid = fopen(fullfile(outdir, 'index.html'), 'w', 'n', 'UTF-8'); fwrite(fid, unicode2native(html, 'UTF-8')); fclose(fid);
end
