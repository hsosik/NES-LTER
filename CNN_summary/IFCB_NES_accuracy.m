load("\\vortex\omics\sosik\training-output\20201022_NES\results.mat")

figure
subplot(1,2,2)
imagesc(confusion_matrix)
axis square
%for ii = 1:numel(confusion_matrix), if confusion_matrix(ii), [a,b] = ind2sub(size(confusion_matrix),ii); text(b,a,num2str(confusion_matrix(ii)), 'HorizontalAlignment', 'center', 'fontsize', 8), end; end
set(gca, 'ytick', 1:1:numel(ind), 'yticklabel', regexprep(class_labels(classes_by_recall), '_', ' '), 'fontsize', 8)
set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels(classes_by_recall), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8, 'XAxisLocation', 'top')

subplot(2,2,1)
[~,ind] = sort(f1_perclass);
bar(1:numel(ind),f1_perclass(ind))
set(gca, 'xtick', 1:1:numel(ind),'xticklabel', []) %, 'xticklabel', regexprep(extractAfter(class_labels(ind), '_'), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8)
ylabel('F1-score', 'fontsize', 12)

subplot(2,2,3)
[~,ind] = sort(f1_perclass);
bar(1:numel(ind),counts_perclass(ind))
set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels(ind), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 10)
ylabel('Count', 'fontsize', 12)
p = get(gca,'position');
p(2) = p(2)*2;
set(gca,'position', p)

set(gcf, 'Position', get(0, 'Screensize'))
%%
figure
group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
[~,ia,ib] = intersect(group_table.CNN_classlist, class_labels);

[~,ind] = sort(f1_perclass);
bar(1:numel(ind),f1_perclass(ind), 'edgecolor', 'none')
ylabel('F1-score', 'fontsize', 16)
set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels(ind), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8)
set(gca, 'position', [.13 .58 .78 .34], 'fontweight', 'bold')
set(gcf, 'Position', get(0, 'Screensize'))
ax = gca;
ax.YAxis.FontSize = 16;

pause

class_ind = ib(find(group_table.Diatom(ia(ind))));
for ii = 1:length(class_ind), ax.XTickLabel{class_ind(ii)} = ['\color{red}' ax.XTickLabel{class_ind(ii)}]; end

class_ind = ib(find(group_table.Dinoflagellate(ia(ind))));
for ii = 1:length(class_ind), ax.XTickLabel{class_ind(ii)} = ['\color{blue}' ax.XTickLabel{class_ind(ii)}]; end

class_ind = ib(find(group_table.Nano(ia(ind))));
for ii = 1:length(class_ind), ax.XTickLabel{class_ind(ii)} = ['\color{green}' ax.XTickLabel{class_ind(ii)}]; end
