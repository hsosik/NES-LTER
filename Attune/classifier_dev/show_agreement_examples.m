function E = show_agreement_examples(tags, nexamples, outdir)
% function E = show_agreement_examples(tags, nexamples, outdir)
% Pool score_agreement output for the cruises in tags (e.g. {'AR92' 'EN727'}),
% pick nexamples files evenly spaced in rank from highest to lowest
% agreement, and draw plot_class_cytograms for each into outdir/examples.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
A = table;
for t = tags
    a = load(fullfile(outdir, ['agreement_' t{1} '.mat']), 'A'); A = [A; a.A]; %#ok<AGROW>
end
A = sortrows(A, 'agreement', 'descend');
fprintf('%d files: agreement median %.3f, range %.3f-%.3f\n', height(A), median(A.agreement), min(A.agreement), max(A.agreement))
E = A(unique(round(linspace(1, height(A), nexamples))), :);
disp(E(:, {'cruise' 'file' 'agreement' 'n_phyto' 'n_added' 'n_dropped' 'n_relabelled'}))
exdir = fullfile(outdir, 'examples');
if ~exist(exdir, 'dir'), mkdir(exdir), end
for i = 1:height(E)
    cpath = fullfile(root, char(E.cruise(i)));
    P = cruise_paths(char(E.cruise(i)), Root=root); % archive-aware
    bs = load_beadstat(P.bead_dir);
    fn = char(E.file(i));
    [fcsdat, fcshdr] = fca_readfcs(fullfile(cpath, 'FCS', fn));
    ref = bead_reference(bs, E.time(i), cfg);
    [new, info] = assign_class_v2(fcsdat, fcshdr, ref, cfg, E.pro_measured(i));
    old = load_old_class(P, fn);
    ttl = sprintf('example %d of %d: %s   agreement %.3f (rank %d of %d)   v2 adds %d, drops %d, relabels %d of %d phyto events   trigger %s %g', ...
        i, height(E), fn, E.agreement(i), find(A.file == E.file(i)), height(A), E.n_added(i), E.n_dropped(i), E.n_relabelled(i), E.n_phyto(i), ...
        fcshdr.tr1_par, fcshdr.tr1_level);
    plot_class_cytograms(fcsdat, fcshdr, old, new, info, ttl, fullfile(exdir, sprintf('example%d_%s.png', i, extractAfter(char(E.cruise(i)), '_'))));
end
end
