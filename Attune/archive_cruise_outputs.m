function L = archive_cruise_outputs(cruises, tag, opt)
% function L = archive_cruise_outputs(cruises, tag, Name=Value)
% Archive each cruise's processing output before reprocessing with process_cruise:
% rename <cruise>\bead_calibrated to <cruise>\bead_calibrated_archive_<tag>, make
% a new empty bead_calibrated, and copy into it the inputs that are not rebuilt
% from the FCS files:
%   bead_review_decisions.csv  bead runs excluded in the bead review (hand decisions)
%   beadstat.mat               cruises without FCB bead runs only (attune_config
%                              fcb_beads false: EN608 to EN627), whose bead table is a
%                              PT transfer or, for EN627, a hand-made PT conversion
%   (Options.CarryReview) review_decisions.csv, classification review decisions
% Everything else (FCSfileinfo, class files, volumes, tables, review pages, movies)
% is rebuilt by process_cruise. The archive gets ARCHIVE_README.txt: date, tag,
% code version used for archiving, what made the output (process_cruise log if
% present), what was carried forward, and the archived file list.
%
% Nothing is changed unless DryRun=false.
%
% cruises - cruise names, or "all" (every row of attune_cruise_settings.csv)
% tag     - archive suffix, e.g. "2026-11" -> bead_calibrated_archive_2026-11
% Options
%   DryRun       true (default): list the actions only
%   Sample       "underway" (default), "preserved" (<cruise>\preserved\bead_calibrated)
%                or "both"
%   CarryReview  false (default): also carry review_decisions.csv forward
%   Root         folder of cruise folders (default as attune_cruise_settings)
%
% L - table, one row per cruise and sample type: cruise, sample, folder, action,
%     status ("dry run", "done", "skipped", "failed"), message
%
% e.g. L = archive_cruise_outputs("all", "2026-11");                  % look first
%      L = archive_cruise_outputs("EN727", "2026-11", DryRun=false);   % then do it
%
% Nov 2026, processing wrapper redesign (reprocessing all cruises with v2)

arguments
    cruises string
    tag (1,1) string {mustBeNonzeroLengthText}
    opt.DryRun (1,1) logical = true
    opt.Sample (1,1) string {mustBeMember(opt.Sample, ["underway" "preserved" "both"])} = "underway"
    opt.CarryReview (1,1) logical = false
    opt.Root (1,1) string = ""
end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'classifier_dev'))
if isscalar(cruises) && lower(cruises) == "all"
    T = readtable(fullfile(here, 'attune_cruise_settings.csv'), 'TextType', 'string', 'Delimiter', ',');
    cruises = T.cruise';
end
samples = opt.Sample; if samples == "both", samples = ["underway" "preserved"]; end
[~, h] = system(['git -C "' here '" rev-parse --short HEAD']); version = strtrim(h);
arch = "bead_calibrated_archive_" + tag;

L = table;
for c = cruises
    try
        s = attune_cruise_settings(c, char(opt.Root));
    catch ME
        L = [L; row(c, "", "", "", "failed", ME.message)]; %#ok<AGROW>
        continue
    end
    for smp = samples
        if smp == "preserved", base = s.preserved_path; else, base = s.basepath; end
        src = fullfile(base, 'bead_calibrated'); dst = fullfile(base, char(arch));
        if ~exist(src, 'dir')
            L = [L; row(c, smp, src, "none", "skipped", "no bead_calibrated folder")]; %#ok<AGROW>
            continue
        end
        if exist(dst, 'dir')
            L = [L; row(c, smp, src, "none", "skipped", "archive already exists: " + dst)]; %#ok<AGROW>
            continue
        end
        carry = "bead_review_decisions.csv";
        if smp == "underway" && ~s.config.fcb_beads, carry = [carry "beadstat.mat"]; end %#ok<AGROW>
        if opt.CarryReview, carry = [carry "review_decisions.csv"]; end %#ok<AGROW>
        carry = carry(arrayfun(@(f) exist(fullfile(src, f), 'file') == 2, carry));
        d = dir(src); d = d(~ismember({d.name}, {'.' '..'}));
        act = sprintf('rename to %s; new bead_calibrated; carry forward: %s', arch, listtxt(carry));
        if opt.DryRun
            L = [L; row(c, smp, src, act, "dry run", sprintf('%d items in bead_calibrated', numel(d)))]; %#ok<AGROW>
            continue
        end
        try
            [ok, msg] = movefile(src, dst);
            if ~ok, error('archive:rename', '%s', msg); end
            mkdir(src);
            for f = carry, copyfile(fullfile(dst, f), fullfile(src, f)); end
            write_readme(dst, s, smp, tag, version, carry, d);
            L = [L; row(c, smp, src, act, "done", "")]; %#ok<AGROW>
        catch ME
            % leave things as found where possible: a renamed folder without a new one is
            % reported, not undone, so nothing is lost
            L = [L; row(c, smp, src, act, "failed", ME.message)]; %#ok<AGROW>
        end
    end
end
disp(L(:, {'cruise' 'sample' 'status' 'message'}))
fprintf('%d done, %d dry run, %d skipped, %d failed\n', sum(L.status == "done"), sum(L.status == "dry run"), ...
    sum(L.status == "skipped"), sum(L.status == "failed"))
end


function r = row(c, smp, folder, act, status, msg)
r = table(string(c), string(smp), string(folder), string(act), string(status), string(msg), ...
    'VariableNames', {'cruise' 'sample' 'folder' 'action' 'status' 'message'});
end


function write_readme(dst, s, smp, tag, version, carry, d)
made = "unknown (no process_cruise_log.mat: process_wrapper era)";
lf = fullfile(dst, 'process_cruise_log.mat');
if exist(lf, 'file')
    try
        g = load(lf, 'runs'); r = g.runs(end);
        made = sprintf('process_cruise, last run %s, steps %s', string(r.date), strjoin(string(r.p.steps), ' '));
    catch
        made = "process_cruise (log unreadable)";
    end
end
fid = fopen(fullfile(dst, 'ARCHIVE_README.txt'), 'w');
c = onCleanup(@() fclose(fid));
fprintf(fid, 'Archived processing output for %s (%s samples)\n', s.cruise, smp);
fprintf(fid, 'Archived: %s, tag %s, by archive_cruise_outputs (code %s)\n', string(datetime('now')), tag, version);
fprintf(fid, 'Made by: %s\n', made);
fprintf(fid, 'Replaced by a new bead_calibrated folder for reprocessing with process_cruise (assign_class_v2).\n');
fprintf(fid, 'Carried forward into the new folder: %s\n\n', listtxt(carry));
fprintf(fid, 'Contents at archiving (%d items):\n', numel(d));
for k = 1:numel(d)
    if d(k).isdir, fprintf(fid, '  %s\\\n', d(k).name);
    else, fprintf(fid, '  %s  %d bytes  %s\n', d(k).name, d(k).bytes, d(k).date);
    end
end
end


function t = listtxt(x)
if isempty(x), t = "(nothing)"; else, t = strjoin(x, ', '); end
end
