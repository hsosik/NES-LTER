function tag = classifier_version(classifier)
% function tag = classifier_version(classifier)
% Version tag of a classifier's output. Bump the tag (date + letter) whenever a
% change alters the class files a classifier writes: class codes, boundaries,
% syn_type, flags, or class_meta fields read downstream. Do not bump for review,
% plotting or documentation changes.
%
% process_cruise saves the tag in class_meta.classifier_tag and in class_log;
% cruise_status marks a cruise's classify step stale when its files were made
% with another tag (or before tags existed).
%
% classifier  'assign_class_v2' (underway), 'assign_class_discrete' (preserved
%             CHL / PE / pro runs), 'assign_class_hbac' (preserved SYBR runs);
%             no input: struct of all tags
%
% Nov 2026, cruise console (processing status)

T = struct( ...
    'assign_class_v2',       "2026-11-28a", ...
    'assign_class_discrete', "2026-11-28a", ...
    'assign_class_hbac',     "2026-11-27a");
if nargin == 0
    tag = T;
    return
end
classifier = char(classifier);
if isfield(T, classifier)
    tag = T.(classifier);
else
    tag = "";
end
end
