function live_holdout_review(sets, outdir, refresh)
% function live_holdout_review(sets, outdir, refresh)
% Static review viewer for the live held-out check, as the discrete one
% (discrete/holdout_review): every scored file outside the pass range
% (failures_<set>.csv from live_holdout_score) plus a seeded ~10% of the passing
% files as controls, with three cytogram rows per file: old classes, v2, and the
% events where they differ. Shown blind as A/B in an order fixed per file here.
% Two pages, each with its own verdict file:
%   outdir/review        files failing on anything other than PE euk, plus the controls
%   outdir/review_peeuk  files whose only failure is the PE euk ratio (mostly the
%                        intended move of Euk+Syn, class 5, into the PE euk group;
%                        the tuning files fail it too)
% Verdicts (v2 wrong / old wrong / ambiguous / OK) save to <page>/live_verdicts.csv;
% read them with
%   T = holdout_verdicts(fullfile(outdir, 'review', 'live_verdicts.csv'))
% (columns as the discrete viewer's; verdict hand_wrong there reads old_wrong here).
% sets    - default {'L1' 'L2'} (T, the tuning files, can be added as a check)
% outdir  - default output/live_holdout
% refresh - true: remake the frames of the files already in the two pages with the
%           current classifier, keeping ids, pages and A/B order (default false)
%
% Nov 2026, classifier redesign (held-out check, live samples)

here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(here))
if nargin < 1 || isempty(sets), sets = {'L1' 'L2'}; end
if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'output', 'live_holdout'); end
if nargin < 3, refresh = false; end
rdir = {fullfile(outdir, 'review') fullfile(outdir, 'review_peeuk')};
ptitle = {'Live held-out review' 'Live held-out review: PE euk only'};
prev = cell(1, 2); made0 = '';
sel = table;
for p = 1:2 * refresh
    f = fullfile(rdir{p}, 'review_data.js');
    if ~exist(f, 'file'), continue, end
    q = jsondecode(regexprep(strtrim(fileread(f)), '^window\.HOLDOUT\s*=\s*|;$', ''));
    if iscell(q.runs), q.runs = [q.runs{:}]; end
    prev{p} = q; made0 = char(q.made);
    bk = fullfile(rdir{p}, ['frames_' regexprep(made0, '[^0-9]', '')]);
    if exist(fullfile(rdir{p}, 'frames'), 'dir') && ~exist(bk, 'dir'), copyfile(fullfile(rdir{p}, 'frames'), bk); end
    pr = q.runs(:);
    if isempty(pr), continue, end
    sel = [sel; table(string({pr.set})', string({pr.file})', string({pr.why})', [pr.control]', ...
        string({pr.id})', string({pr.auto_is})', 'VariableNames', {'set' 'file' 'why' 'control' 'id' 'auto_is'})]; %#ok<AGROW>
end
if refresh, sets = {}; end
for p = 1:2, if ~exist(fullfile(rdir{p}, 'frames'), 'dir'), mkdir(fullfile(rdir{p}, 'frames')), end, end
L = readtable(fullfile(outdir, 'holdout_lists.csv'), 'TextType', 'string');
for s = 1:numel(sets)
    H = load(fullfile(outdir, ['score_' sets{s} '.mat'])).H;
    Fl = readtable(fullfile(outdir, ['failures_' sets{s} '.csv']), 'TextType', 'string');
    bad = unique(Fl.file);
    why = strings(size(bad));
    for b = 1:numel(bad)
        k = find(Fl.file == bad(b));
        why(b) = strjoin(arrayfun(@(j) sprintf('%s %.2f (v2/old)', Fl.class(j), Fl.ratio(j)), k(:)', 'UniformOutput', false), '; ');
        why(b) = regexprep(why(b), 'agreement ([\d.]+) \(v2/old\)', 'agreement $1');
    end
    ok = setdiff(H.file, bad);
    rng(3029 + s, 'twister');
    ctl = ok(randperm(numel(ok), min(numel(ok), max(1, round(0.1 * height(H))))));
    sel = [sel; table(repmat(string(sets{s}), numel(bad), 1), bad, why, false(numel(bad), 1), 'VariableNames', {'set' 'file' 'why' 'control'}); ...
        table(repmat(string(sets{s}), numel(ctl), 1), ctl(:), repmat("", numel(ctl), 1), true(numel(ctl), 1), 'VariableNames', {'set' 'file' 'why' 'control'})]; %#ok<AGROW>
end
if ~refresh, rng(3040, 'twister'); sel = sel(randperm(height(sel)), :); end
fprintf('%d files to review (%d controls)\n', height(sel), sum(sel.control));

views = {'size' 'BL3'; 'BL3' 'GL2'; 'size' 'PEchl'; 'PEchl' 'GL2GL3'; 'size' 'GL2'; 'chlSSC' 'GL2GL3'; 'SSC' 'BL3'};
lim = struct('size', [-2.6 2.8], 'SSC', [-3.2 2.6], 'BL3', [-1.3 3.2], 'GL2', [-3.5 1.7], ...
    'PEchl', [-4.4 1.8], 'GL2GL3', [-2.5 1.5], 'chlSSC', [-1 3]);
names = struct('size', 'size', 'SSC', 'SSC', 'BL3', 'chl (BL3)', 'GL2', 'PE (GL2)', 'PEchl', 'PE:chl', 'GL2GL3', 'GL2:GL3', 'chlSSC', 'chl:size');
runs = {{} {}};
% page 2: files whose only failure is PE euk
pe_only = arrayfun(@(w) w ~= "" && all(startsWith(strtrim(split(w, ';')), "peeuk")), sel.why);
rng(3050, 'twister');
for i = 1:height(sel)
    pg = 1 + pe_only(i); fd = fullfile(rdir{pg}, 'frames');
    j = find(L.file == sel.file(i) & L.set == sel.set(i), 1);
    cd0 = char(L.cruisedir(j)); fn = char(L.file(j));
    try
        P = cruise_paths(cd0);
        cfg = classifier_config(); cfg.cruise = extractAfter(cd0, '_');
        cfg.gl1_is_od2_scatter = strcmpi(attune_config(cfg.cruise).OD2, 'GL1');
        [dat, hdr] = fca_readfcs(fullfile(P.basepath, 'FCS', fn));
        ref = bead_reference(load_beadstat(P.bead_dir), L.time(j), cfg);
        [a, info] = assign_class_v2(dat, hdr, ref, cfg, L.pro_measured(j));
        h = load_old_class(P, fn);
    catch err
        fprintf('%s: %s\n', fn, err.message); continue
    end
    h = double(h(:)); h(~ismember(h, [0 1 2 3 4 7])) = 0; h(ismember(h, [3 4])) = 34;
    a = double(a(:)); a(ismember(a, [3 4 5])) = 34;
    F = info.features;
    id = sprintf('r%03d', i); if refresh, id = char(sel.id(i)); end
    frame(F, h, views, lim, fullfile(fd, [id '_h.png']));
    frame(F, a, views, lim, fullfile(fd, [id '_a.png']));
    framediff(F, h, a, views, lim, fullfile(fd, [id '_d.png']));
    cn = struct(); for c = [1 2 7 34 8], cn.(sprintf('c%d', c)) = [sum(h == c) sum(a == c)]; end
    cr = extractAfter(cd0, '_');
    r = struct('id', id, 'set', char(sel.set(i)), 'file', fn, 'cruise', cr, ...
        'runtype', ternary(L.pro_measured(j) == 1, 'Pro-capable', 'standard'), 'when', datestr(L.time(j), 'yyyy-mm-dd HH:MM'), ...
        'control', sel.control(i), 'why', char(sel.why(i)), 'auto_is', 'A', ...
        'counts', cn, 'n', numel(a), 'agree', sum(h == a & h ~= 0) / max(1, sum(h ~= 0 | a ~= 0)), ...
        'flags', sprintf('v2 flags: Syn-pigment particle %d, below Euk chl gap %d, trigger-edge pile %d', sum(info.flag == 1), sum(info.flag == 2), sum(info.flag == 3)), ...
        'views', {arrayfun(@(k) sprintf('%s vs %s', names.(views{k, 2}), names.(views{k, 1})), 1:size(views, 1), 'UniformOutput', false)});
    if rand < 0.5, r.auto_is = 'B'; end
    if refresh, r.auto_is = char(sel.auto_is(i)); end
    runs{pg}{end + 1} = r;
    fprintf('%d/%d %s\n', i, height(sel), fn);
end
[~, commit] = system(['git -C "' here '" rev-parse --short HEAD']);
bname = ['v2 live held-out ' strtrim(commit) ' ' datestr(now, 'yyyy-mm-dd HH:MM')];
if refresh, bname = [bname ' (frames remade, files from ' made0 ')']; end
write_pages(here, rdir, ptitle, runs, bname, datestr(now, 'yyyy-mm-dd HH:MM'));
end

function write_pages(here, rdir, ptitle, runs, bname, made)
% review_data.js and index.html for each page, each linking to the other
lnk = {'../review_peeuk/index.html' '../review/index.html'};
lname = {'PE euk only page' 'main page'};
for p = 1:2
    rd = rdir{p};
    if ~exist(rd, 'dir'), mkdir(rd), end
    D = struct('build', bname, 'made', made, 'title', ptitle{p}, 'other_href', lnk{p}, 'other_name', lname{p}, ...
        'decisions_name', 'live_verdicts.csv', 'decisions_path', fullfile(rd, 'live_verdicts.csv'), ...
        'decisions_dir', rd, 'runs', {runs{p}});
    [~, D.decisions_dir_name] = fileparts(rd);
    vf = fullfile(rd, 'live_verdicts.csv');
    D.decisions = {};
    if exist(vf, 'file'), T = holdout_verdicts(vf); D.decisions = table2struct(T)'; end
    fid = fopen(fullfile(rd, 'review_data.js'), 'w', 'n', 'UTF-8');
    fprintf(fid, 'window.HOLDOUT = %s;\n', jsonencode(D)); fclose(fid);
    tpl = fileread(fullfile(here, 'live_holdout_review_template.html'));
    tpl = strrep(tpl, 'review_data.js"', ['review_data.js?v=' datestr(now, 'yyyymmddHHMMSS') '"']);
    fid = fopen(fullfile(rd, 'index.html'), 'w', 'n', 'UTF-8'); fwrite(fid, tpl, 'char'); fclose(fid);
    fprintf('%s: %d files, %s\n', ptitle{p}, numel(runs{p}), fullfile(rd, 'index.html'));
end
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end

function frame(F, c, views, lim, outfile)
% one row of panels, classes on the discrete viewer's colours
W = 220; g = 4; nv = size(views, 1);
img = ones(W, nv * W + (nv - 1) * g, 3);
cols = containers.Map('KeyType', 'double', 'ValueType', 'any');
cols(1) = [0 .6 0]; cols(2) = [1 .55 0]; cols(7) = [0 .35 1];
for v = 1:nv
    [ix, ok] = binxy(F.(views{v, 1}), F.(views{v, 2}), lim.(views{v, 1}), lim.(views{v, 2}), W);
    P = ones(W, W, 3);
    P = dens(P, ix, ok & c == 0, [.45 .45 .45], .55);
    ab = [1 2 7]; [~, o] = sort(arrayfun(@(k) sum(c == k), ab), 'descend');
    for k = ab(o), P = dens(P, ix, ok & c == k, cols(k), .9); end
    P = pts(P, ix, ok & c == 34, [.9 0 .75]);
    P = pts(P, ix, ok & c == 8, [.85 0 0]);
    img(:, (v - 1) * (W + g) + (1:W), :) = P;
end
writepng(img, outfile)
end

function framediff(F, h, a, views, lim, outfile)
% same label (grey density) and different label (red points), without saying which is which
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
