function P = cruise_paths(cruisedir, opt)
% function P = cruise_paths(cruisedir, Name=Value)
% Folders for scoring v2 against the old classifier once a cruise's output has
% been archived (archive_cruise_outputs renames <cruise>\bead_calibrated to
% bead_calibrated_archive_<tag> and process_cruise rebuilds bead_calibrated).
%
% cruisedir - cruise folder name, e.g. '20250124_EN727'
% Options
%   Tag     archive tag to use for the old classes (default: the only archive; the
%           latest by name, with a warning, if there are several)
%   Sample  "underway" (default) or "preserved"
%   Root    folder of cruise folders (default \\sosiknas1\lab_data\Attune\cruise_data\)
%
% P.basepath   cruise folder (FCS files in P.basepath\FCS)
% P.current    <cruise>\bead_calibrated (new output)
% P.archived   true if an archive folder was found
% P.tag        archive tag used ('' if none)
% P.old        folder holding the old (pre-v2) output: the archive if there is one,
%              otherwise bead_calibrated itself (cruise not yet archived)
% P.old_class  P.old\class (old class files, one per FCS file)
% P.info_dir   folder to read FCSfileinfo.mat from: P.current once process_cruise
%              has rebuilt it, otherwise P.old
% P.bead_dir   folder to read beadstat.mat from (same rule)
% P.baseline_class  P.current\class when the cruise is archived and v2 class files
%              exist there (the reprocessed v2 baseline), otherwise ''
%
% Old class files are plain class matrices; v2 class files from process_cruise
% carry class_meta. Use cruise_paths_check_old(file) before treating a file as old.
%
% Nov 2026, classifier redesign (scoring after archiving)

arguments
    cruisedir (1,:) char
    opt.Tag (1,1) string = ""
    opt.Sample (1,1) string {mustBeMember(opt.Sample, ["underway" "preserved"])} = "underway"
    opt.Root (1,:) char = '\\sosiknas1\lab_data\Attune\cruise_data\'
end
P.basepath = fullfile(opt.Root, cruisedir);
if opt.Sample == "preserved", P.basepath = fullfile(P.basepath, 'preserved'); end
P.current = fullfile(P.basepath, 'bead_calibrated');
d = dir(fullfile(P.basepath, 'bead_calibrated_archive_*'));
d = d([d.isdir]);
tags = erase(string({d.name}), "bead_calibrated_archive_");
if opt.Tag ~= ""
    if ~ismember(opt.Tag, tags)
        error('cruise_paths:tag', 'no archive with tag %s in %s', opt.Tag, P.basepath)
    end
    P.tag = char(opt.Tag);
elseif isempty(tags)
    P.tag = '';
else
    tags = sort(tags);
    if numel(tags) > 1
        warning('cruise_paths:tags', '%s has %d archives (%s); using %s', cruisedir, numel(tags), strjoin(tags, ', '), tags(end))
    end
    P.tag = char(tags(end));
end
P.archived = ~isempty(P.tag);
if P.archived, P.old = fullfile(P.basepath, ['bead_calibrated_archive_' P.tag]);
else, P.old = P.current;
end
P.old_class = fullfile(P.old, 'class');
P.info_dir = P.old; if exist(fullfile(P.current, 'FCSfileinfo.mat'), 'file'), P.info_dir = P.current; end
P.bead_dir = P.old; if exist(fullfile(P.current, 'beadstat.mat'), 'file'), P.bead_dir = P.current; end
P.baseline_class = '';
if P.archived && isfolder(fullfile(P.current, 'class')), P.baseline_class = fullfile(P.current, 'class'); end
end
