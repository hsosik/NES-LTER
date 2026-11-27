function S = bead_review_status(bp)
% function S = bead_review_status(bp)
% Whether a cruise's bead review (bead_review page) has been marked complete for
% the bead table that is there now. bp: the folder with beadstat.mat and
% bead_review_decisions.csv (live: <cruise>\bead_calibrated; preserved:
% <cruise>\preserved\bead_calibrated).
%
% The page's "mark review complete" button saves one row in
% bead_review_decisions.csv with file = '__review_complete__' and status 'complete'
% ('' after "reopen"), plus reviewer, saved (UTC) and build (the bead stats the
% page showed: meta.processed of beadstat.mat). Same rules as the page: the review
% counts as complete only when
%   - that row says complete,
%   - its build matches the beadstat.mat now in bp (not rebuilt since), and
%   - no run decision was changed after it.
%
% S.complete   true/false
% S.state      'complete', 'no bead table', 'not reviewed' (no decisions file),
%              'not marked complete', 'reopened', 'bead stats rebuilt since',
%              'decisions changed since'
% S.reviewer, S.saved (datetime, UTC), S.build: from the completion row ('' / NaT)
% S.message    one line for warnings, e.g. from process_cruise
%
% Nov 2026, bead review

S = struct('complete', false, 'state', '', 'reviewer', '', 'saved', NaT, 'build', '', 'message', '');
bp = char(bp);
fn = fullfile(bp, 'beadstat.mat');
if ~isfile(fn)
    S.state = 'no bead table';
    S.message = sprintf('no beadstat.mat in %s', bp);
    return
end
m = load(fn, 'meta');
build = char(datetime(m.meta.processed, 'Format', 'yyyy-MM-dd HH:mm')); % as bead_review's D.build
dec = fullfile(bp, 'bead_review_decisions.csv');
page = mapped(fullfile(bp, 'bead_review', 'index.html')); % the page saves decisions only when opened from a mapped drive
if ~isfile(dec)
    S.state = 'not reviewed';
    S.message = sprintf('bead review not done (no %s); open from NAS mapped drive, e.g., %s', 'bead_review_decisions.csv', page);
    return
end
T = review_decisions(dec, true);
k = find(strcmp(T.file, '__review_complete__'), 1, 'last');
if isempty(k)
    S.state = 'not marked complete';
    S.message = sprintf('bead review not marked complete: open from NAS mapped drive, e.g., %s, and use "mark review complete" when done', page);
    return
end
S.reviewer = T.reviewer{k}; S.build = T.build{k};
S.saved = iso(T.saved{k});
who = S.reviewer; if isempty(who), who = 'someone'; end
if ~strcmp(T.status{k}, 'complete')
    S.state = 'reopened';
    S.message = sprintf('bead review reopened by %s on %s and not marked complete again: %s', who, char(S.saved, 'yyyy-MM-dd'), page);
    return
end
if ~strcmp(S.build, build)
    S.state = 'bead stats rebuilt since';
    S.message = sprintf('bead review was completed for older bead stats (%s; beadstat.mat now %s): review the new table at %s', S.build, build, page);
    return
end
runs = ~strcmp(T.file, '__review_complete__');
later = runs & cellfun(@(s) iso(s) > S.saved, T.saved);
if any(later)
    S.state = 'decisions changed since';
    S.message = sprintf('%d bead run decision(s) changed after the review was marked complete (%s, %s): mark it complete again at %s', ...
        sum(later), who, char(S.saved, 'yyyy-MM-dd'), page);
    return
end
S.complete = true;
S.state = 'complete';
S.message = sprintf('bead review complete (%s, %s)', who, char(S.saved, 'yyyy-MM-dd HH:mm'));
end

function f = mapped(f)
% \\sosiknas1\Lab_data\... -> <drive>:\... (nas_drive: this computer's mapping of
% the NAS share, default Z)
if exist('nas_drive', 'file')
    f = nas_drive('path', f);
elseif startsWith(f, '\\sosiknas1\Lab_data\', 'IgnoreCase', true)
    f = ['Z:\' f(numel('\\sosiknas1\Lab_data\')+1:end)];
end
end

function t = iso(s)
% ISO 8601 UTC from the page (toISOString) -> datetime; NaT if blank or unreadable
t = NaT;
if isempty(s), return, end
try
    t = datetime(s, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', 'TimeZone', 'UTC');
catch
    try, t = datetime(s, 'TimeZone', 'UTC'); catch, end
end
if ~isnat(t), t.TimeZone = ''; end
end
