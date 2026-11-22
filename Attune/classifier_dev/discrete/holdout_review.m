function holdout_review(sets, outdir, refresh)
% function holdout_review(sets, outdir, refresh)
% Build the static review viewer for the held-out check (step 4,
% holdout_review_design.md in the shared notes). For every scored held-out run
% outside the pass range (failures_<set>.csv from holdout_report), plus a seeded
% ~10% of the passing runs as controls, draws three cytogram rows: hand gates,
% frozen classifier, and the events where they disagree. Writes outdir/index.html
% (from holdout_review_template.html), review_data.js and frames/*.png.
% In the viewer the two labellings are shown as A/B (blind mode, default) in an
% order fixed per run here; verdicts save to outdir/holdout_verdicts.csv
% (read back with holdout_verdicts).
% sets   - default {'2a' '2a_prime' '2b'}
% outdir - default output/discrete/holdout/review
% refresh - true: remake the frames of the runs already in outdir/review_data.js with
%           the current classifier and hand classes, keeping each run's id and A/B order
%           so saved verdicts still apply (old frames copied to frames_<old build date>);
%           the run list is not redrawn from the current failures (default false)
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fileparts(fileparts(here)))
od = fullfile(fileparts(here), 'output', 'discrete');
if nargin < 1 || isempty(sets), sets = {'2a' '2a_prime' '2b'}; end
if nargin < 2 || isempty(outdir), outdir = fullfile(od, 'holdout', 'review'); end
if nargin < 3, refresh = false; end
prev = [];
if refresh
    js = fileread(fullfile(outdir, 'review_data.js'));
    prev = jsondecode(regexprep(strtrim(js), '^window\.HOLDOUT\s*=\s*|;$', ''));
    if iscell(prev.runs), prev.runs = [prev.runs{:}]; end
    bk = fullfile(outdir, ['frames_' regexprep(char(prev.made), '[^0-9]', '')]);
    if exist(fullfile(outdir, 'frames'), 'dir') && ~exist(bk, 'dir'), copyfile(fullfile(outdir, 'frames'), bk); end
end
fd = fullfile(outdir, 'frames'); if ~exist(fd, 'dir'), mkdir(fd), end
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
L = readtable(fullfile(od, 'holdout', 'holdout_lists.csv'), 'TextType', 'string');

% runs to review: failures, plus controls (refresh: the runs already in the viewer)
sel = table;
if refresh
    pr = prev.runs(:);
    sel = table(string({pr.set})', string({pr.file})', string({pr.why})', [pr.control]', ...
        string({pr.id})', string({pr.auto_is})', 'VariableNames', {'set' 'file' 'why' 'control' 'id' 'auto_is'});
    sets = {};
end
for s = 1:numel(sets)
    H = load(fullfile(od, 'holdout', ['score_' sets{s} '.mat'])).H;
    Fl = readtable(fullfile(od, 'holdout', ['failures_' sets{s} '.csv']), 'TextType', 'string');
    bad = unique(Fl.file);
    why = strings(size(bad));
    for b = 1:numel(bad)
        k = find(Fl.file == bad(b));
        why(b) = strjoin(arrayfun(@(j) sprintf('%s %.2f (classifier/hand)', Fl.class(j), Fl.ratio(j)), k(:)', 'UniformOutput', false), '; ');
    end
    ok = setdiff(H.file, bad);
    rng(2029 + s, 'twister');
    ctl = ok(randperm(numel(ok), max(1, round(0.1 * numel(H.file)))));
    sel = [sel; table(repmat(string(sets{s}), numel(bad), 1), bad, why, false(numel(bad), 1), 'VariableNames', {'set' 'file' 'why' 'control'}); ...
        table(repmat(string(sets{s}), numel(ctl), 1), ctl(:), repmat("", numel(ctl), 1), true(numel(ctl), 1), 'VariableNames', {'set' 'file' 'why' 'control'})]; %#ok<AGROW>
end
if ~refresh, rng(2040, 'twister'); sel = sel(randperm(height(sel)), :); end %review order: failures and controls mixed
fprintf('%d runs to review (%d controls)\n', height(sel), sum(sel.control));

V.phyto = {'size' 'BL3'; 'BL3' 'GL2'; 'size' 'PEchl'; 'PEchl' 'GL2GL3'; 'size' 'GL2'};
V.hbac = {'SSC' 'BL1'; 'BL1' 'BL3'; 'BL1' 'GL2'; 'SSC' 'BL3'; 'SSC' 'GL2'};
lim = struct('size', [-2.6 2.8], 'SSC', [-3.2 2.6], 'BL1', [-1.6 3.2], 'BL3', [-1.3 3.2], 'GL2', [-3.5 1.7], ...
    'PEchl', [-4.4 1.8], 'GL2GL3', [-2.5 1.5]);
names = struct('size', 'size', 'SSC', 'SSC', 'BL1', 'SYBR (BL1)', 'BL3', 'chl (BL3)', 'GL2', 'PE (GL2)', 'PEchl', 'PE:chl', 'GL2GL3', 'GL2:GL3');
runs = {};
Gc = struct();
rng(2050, 'twister');
for i = 1:height(sel)
    j = find(L.file == sel.file(i), 1);
    cdir = char(L.cruisedir(j)); fn = char(L.file(j)); r = L.runtype(j);
    cr = extractAfter(cdir, '_'); if isempty(cr), cr = cdir; end
    key = matlab.lang.makeValidName(cdir);
    if ~isfield(Gc, key)
        G = load(fullfile(root, cdir, 'preserved', 'outputs', 'Gated_Table.mat')).gated_table;
        G.Properties.VariableNames = lower(G.Properties.VariableNames); Gc.(key) = G;
    end
    G = Gc.(key); gi = find(strcmp(G.fcslist, fn), 1);
    try
        [dat, hdr] = fca_readfcs(fullfile(root, cdir, 'preserved', 'FCS', fn));
        cfg = discrete_config(); cfg.cruise = cr; cfg.gl1_is_od2_scatter = true;
        ref = discrete_bead_reference(lab, datenum([hdr.date ' ' hdr.starttime]), cfg);
        if r == "hbac", [auto, info] = assign_class_hbac(dat, hdr, ref, cfg);
        else, [auto, info] = assign_class_discrete(dat, hdr, ref, cfg, char(r)); end
    catch err
        fprintf('%s: %s\n', fn, err.message); continue
    end
    hand = double(G.class{gi}(:));
    % hand classes rebuilt with the wider gate-name map, as in holdout_score
    if all(ismember({'gate_names' 'gate_assignments'}, G.Properties.VariableNames)) && ~isempty(G.gate_assignments{gi})
        hr = hand_class_rebuild(G.gate_names{gi}, G.gate_assignments{gi}, true);
        if numel(hr) == numel(hand), hand = hr(:); end
    end
    if r == "hbac"
        h = hand; h(ismember(h, [3 4])) = 97; h(~ismember(h, [1 2 97])) = 0;
        a = double(auto); a(ismember(a, [7 9])) = 97; a(~ismember(a, [1 2 97])) = 0;
        views = V.hbac;
    else
        h = hand; h(h == 4) = 7; h(ismember(h, [5 6])) = 34; h(h == 3) = 0;
        if contains(fn, 'EN727') && contains(fn, ["_C19N13" "_C19N02" "_C09N01" "_C06N06"]), h(h == 7) = 0; end
        a = double(auto); a(ismember(a, [3 4 5])) = 34;
        views = V.phyto;
    end
    late = true(size(h)); if any(h ~= 0), late = dat(:, 1) >= min(dat(h ~= 0, 1)); end
    F = info.features; F = F(late, :); h = h(late); a = a(late);
    if ~ismember('size', F.Properties.VariableNames) && ismember('SSC', F.Properties.VariableNames), F.size = F.SSC; end
    if ~ismember('SSC', F.Properties.VariableNames) && ismember('size', F.Properties.VariableNames), F.SSC = F.size; end
    if ~ismember('PEchl', F.Properties.VariableNames) && all(ismember({'GL2' 'BL3'}, F.Properties.VariableNames)), F.PEchl = F.GL2 - F.BL3; end
    if ~ismember('GL2GL3', F.Properties.VariableNames) && all(ismember({'GL2' 'GL3'}, F.Properties.VariableNames)), F.GL2GL3 = F.GL2 - F.GL3; end
    id = sprintf('r%03d', i); if refresh, id = char(sel.id(i)); end
    frame(F, h, views, lim, fullfile(fd, [id '_h.png']));
    frame(F, a, views, lim, fullfile(fd, [id '_a.png']));
    framediff(F, h, a, views, lim, fullfile(fd, [id '_d.png']));
    cls = [1 2 7 34 97 8];
    cn = struct(); for c = cls, cn.(sprintf('c%d', c)) = [sum(h == c) sum(a == c)]; end
    dep = L.depth_m(j);
    runs{end + 1} = struct('id', id, 'set', char(sel.set(i)), 'file', fn, 'cruise', cr, 'runtype', char(r), ...
        'depth', dep, 'control', sel.control(i), 'why', char(sel.why(i)), 'auto_is', 'A', ...
        'counts', cn, 'n', sum(late), 'agree', sum(h == a & h ~= 0) / max(1, sum(h ~= 0 | a ~= 0)), ...
        'views', {arrayfun(@(k) sprintf('%s vs %s', names.(views{k, 2}), names.(views{k, 1})), 1:size(views, 1), 'UniformOutput', false)}); %#ok<AGROW>
    if rand < 0.5, runs{end}.auto_is = 'B'; end %blind order: which row (A/B) is the classifier
    if refresh, runs{end}.auto_is = char(sel.auto_is(i)); end %keep the order the verdicts were made in
    fprintf('%d/%d %s\n', i, height(sel), fn);
end
bname = ['v2-holdout-test ' datestr(now, 'yyyy-mm-dd HH:MM')];
if refresh, bname = [bname ' (frames remade, runs from ' char(prev.made) ')']; end
D = struct('build', bname, 'made', datestr(now, 'yyyy-mm-dd HH:MM'), ...
    'decisions_name', 'holdout_verdicts.csv', 'decisions_path', fullfile(outdir, 'holdout_verdicts.csv'), ...
    'decisions_dir', outdir, 'runs', {runs});
[~, D.decisions_dir_name] = fileparts(outdir);
vf = fullfile(outdir, 'holdout_verdicts.csv');
D.decisions = {};
if exist(vf, 'file'), T = holdout_verdicts(vf); D.decisions = table2struct(T)'; end
fid = fopen(fullfile(outdir, 'review_data.js'), 'w', 'n', 'UTF-8');
fprintf(fid, 'window.HOLDOUT = %s;\n', jsonencode(D)); fclose(fid);
tpl = fileread(fullfile(here, 'holdout_review_template.html'));
tpl = strrep(tpl, 'review_data.js"', ['review_data.js?v=' datestr(now, 'yyyymmddHHMMSS') '"']);
fid = fopen(fullfile(outdir, 'index.html'), 'w', 'n', 'UTF-8'); fwrite(fid, tpl, 'char'); fclose(fid);
fprintf('Viewer: %s\n', fullfile(outdir, 'index.html'));
end

function frame(F, c, views, lim, outfile)
% one row of panels, classes on the shared scheme
W = 220; g = 4; nv = size(views, 1);
img = ones(W, nv * W + (nv - 1) * g, 3);
cols = containers.Map('KeyType', 'double', 'ValueType', 'any');
cols(1) = [0 .6 0]; cols(2) = [1 .55 0]; cols(7) = [0 .35 1]; cols(97) = [0 .6 .6];
for v = 1:nv
    [ix, ok] = binxy(F.(views{v, 1}), F.(views{v, 2}), lim.(views{v, 1}), lim.(views{v, 2}), W);
    P = ones(W, W, 3);
    P = dens(P, ix, ok & c == 0, [.45 .45 .45], .55);
    ab = [1 2 7 97]; [~, o] = sort(arrayfun(@(k) sum(c == k), ab), 'descend');
    for k = ab(o), P = dens(P, ix, ok & c == k, cols(k), .9); end
    P = pts(P, ix, ok & c == 34, [.9 0 .75]);
    P = pts(P, ix, ok & c == 8, [.85 0 0]);
    img(:, (v - 1) * (W + g) + (1:W), :) = P;
end
writepng(img, outfile)
end

function framediff(F, h, a, views, lim, outfile)
% events labelled the same (grey density) and differently (red points), without
% saying which labelling is which
W = 220; g = 4; nv = size(views, 1);
img = ones(W, nv * W + (nv - 1) * g, 3);
same = h == a; d = ~same;
for v = 1:nv
    [ix, ok] = binxy(F.(views{v, 1}), F.(views{v, 2}), lim.(views{v, 1}), lim.(views{v, 2}), W);
    P = ones(W, W, 3);
    P = dens(P, ix, ok & same & h == 0, [.7 .7 .7], .4);
    P = dens(P, ix, ok & same & h ~= 0, [.35 .35 .35], .6);
    P = pts(P, ix, ok & d, [.85 0 0]);
    img(:, (v - 1) * (W + g) + (1:W), :) = P;
end
writepng(img, outfile)
end

function writepng(img, outfile)
img = uint8(round(img * 255));
[X, map] = rgb2ind(img, 128, 'nodither');
imwrite(X, map, outfile)
end

function [ix, ok] = binxy(x, y, lx, ly, W)
ok = isfinite(x) & isfinite(y);
c = 1 + floor((x - lx(1)) / diff(lx) * W); c = min(max(c, 1), W);
r = W - floor((y - ly(1)) / diff(ly) * W); r = min(max(r, 1), W);
ix = r + (c - 1) * W; ix(~ok) = 1;
end

function P = dens(P, ix, m, col, amax)
if ~any(m), return, end
W = size(P, 1); amin = 0.3;
n = accumarray(ix(m), 1, [W * W 1]);
a = zeros(W * W, 1); k = n > 0;
a(k) = amax * (amin + (1 - amin) * log1p(n(k)) / log1p(max(n)));
a = reshape(a, W, W);
for ch = 1:3, P(:, :, ch) = P(:, :, ch) .* (1 - a) + col(ch) * a; end
end

function P = pts(P, ix, m, col)
if ~any(m), return, end
W = size(P, 1);
[r, c] = ind2sub([W W], unique(ix(m)));
[dr, dc] = meshgrid(-1:1);
rr = min(max(r + dr(:)', 1), W); cc = min(max(c + dc(:)', 1), W);
k = rr(:) + (cc(:) - 1) * W;
for ch = 1:3, P(k + (ch - 1) * W * W) = col(ch); end
end
