function Q = examine_discrete_events(tag)
% function Q = examine_discrete_events(tag)
% Feature distributions of hand-gated populations in preserved runs (event
% subsample E from score_discrete), and where the automatic classifier
% disagrees. Prints quantiles; figure output/discrete/events_<tag>.png
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
s = load(fullfile(od, ['score_' tag '.mat']));
E = s.E;
E.BL2GL2 = E.BL2 - E.GL2;
feat = {'size' 'BL3' 'GL2' 'PEchl' 'GL2GL3' 'BL2GL2' 'BL1' 'chlSSC'};
grp = {'hand Euk' E.hand == 1; 'hand Syn' E.hand == 2; 'hand Pro' E.hand == 7; 'hand PE euk' E.hand == 34; ...
    'hand 0' E.hand == 0; 'auto Syn, hand 0' E.auto == 2 & E.hand == 0; 'auto Syn, hand Euk' E.auto == 2 & E.hand == 1; ...
    'auto PEeuk, hand Euk' E.auto == 34 & E.hand == 1; 'auto PEeuk, hand Syn' E.auto == 34 & E.hand == 2; 'auto 0, hand Syn' E.auto == 0 & E.hand == 2};
Q = table;
for g = 1:size(grp, 1)
    x = E(grp{g, 2} & E.runtype ~= "PE", :); %CHL and pro runs (PE runs trigger on PE)
    row = table(string(grp{g, 1}), height(x), 'VariableNames', {'group' 'n'});
    for f = feat
        q = prctile(x.(f{1}), [5 50 95]);
        row.(f{1}) = sprintf('%5.2f %5.2f %5.2f', q);
    end
    Q = [Q; row]; %#ok<AGROW>
end
disp('CHL and pro runs: 5th / 50th / 95th percentiles')
disp(Q)
cols = [0 .6 0; 1 0 0; .5 .5 .5; 1 0 1; .8 .8 .8];
fig = figure('visible', 'off', 'position', [50 50 1500 800]);
pairs = {'size' 'BL3'; 'BL3' 'GL2'; 'size' 'PEchl'; 'PEchl' 'GL2GL3'; 'GL2' 'BL2GL2'; 'BL3' 'BL1'};
x = E(E.runtype ~= "PE", :);
for row = 1:2
    for k = 1:size(pairs, 1)
        subplot(2, size(pairs, 1), (row - 1) * size(pairs, 1) + k), hold on
        if row == 1, lab = x.hand; ttl = 'hand gates'; else, lab = x.auto; ttl = 'auto'; end
        cs = [0 1 2 7 34];
        for c = 1:5
            i = lab == cs(c);
            if cs(c) == 0, i = i & rand(size(i)) < 0.3; end
            plot(x.(pairs{k, 1})(i), x.(pairs{k, 2})(i), '.', 'color', cols(mod(c - 2, 5) + 1, :), 'markersize', 1)
        end
        xlabel(pairs{k, 1}), ylabel(pairs{k, 2}), title(ttl, 'fontsize', 8), grid on
    end
end
annotation('textbox', [0 0.95 1 0.05], 'string', 'green Euk, red Syn, gray Pro, magenta PE euk, light gray 0 (CHL and pro runs, pooled)', 'edgecolor', 'none', 'horizontalalignment', 'center')
exportgraphics(fig, fullfile(od, ['events_' tag '.png']), 'resolution', 80), close(fig)
end
