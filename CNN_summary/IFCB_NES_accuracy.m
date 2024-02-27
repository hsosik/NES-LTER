%load("\\vortex\omics\sosik\training-output\20211229_Jan2022_NES_2.2\results.mat")
%load("\\vortex\omics\sosik\training-output\20220101_Gdel_sub_1_0\results.mat")
%load("\\vortex\omics\sosik\training-output\20220105_Jan2022_NES_2.1_iv3plus\results.mat")
%load("\\vortex\omics\sosik\training-output\20220113_EXPORTS_pacific_Dec2021_3\results");
%f1load("C:\work\CNN_models\20201102_NES_norm_seed__THIS IS really 20210429_Dec2020_NES\results.mat")
%load("C:\work\CNN_models\20210430_Dec2020_NES_norm_1.7\results.mat")
%load("C:\work\CNN_models\202106_Dec2020_NES_1.6\results.mat")
%load('\\sosiknas1\training_sets\IFCB\training-output\20220416_Delmar_NES_1\results')
load('\\sosiknas1\training_sets\IFCB\training-output\20220209_Jan2022_NES_2.4\results')
%ad('\\vast\proj\omics\sosik\training-output\20231215_NBP2022_NES_1.0\results')
%load('\\vast\proj\omics\sosik\training-output\20231218_NBP2022_NES_1.0\results')

figure
subplot(1,2,2)
imagesc(confusion_matrix)
axis square
%for ii = 1:numel(confusion_matrix), if confusion_matrix(ii), [a,b] = ind2sub(size(confusion_matrix),ii); text(b,a,num2str(confusion_matrix(ii)), 'HorizontalAlignment', 'center', 'fontsize', 8), end; end
[~,ind] = sort(f1_perclass);
%set(gca, 'ytick', 1:1:numel(ind), 'yticklabel', regexprep(class_labels(classes_by_recall), '_', ' '), 'fontsize', 8)
%set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels(classes_by_recall), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8, 'XAxisLocation', 'top')
set(gca, 'ytick', 1:1:numel(ind), 'yticklabel', regexprep(class_labels, '_', ' '), 'fontsize', 8)
set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels, '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8, 'XAxisLocation', 'top')

subplot(2,2,1)
[~,ind] = sort(f1_perclass);
bar(1:numel(ind),f1_perclass(ind))
set(gca, 'xtick', 1:1:numel(ind)-1,'xticklabel', []) %, 'xticklabel', regexprep(extractAfter(class_labels(ind), '_'), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8)
ylabel('F1-score', 'fontsize', 12)

subplot(2,2,3)
[~,ind] = sort(f1_perclass);
bar(1:numel(ind),counts_perclass(ind))
set(gca, 'xtick', .5:1:numel(ind), 'xticklabel', regexprep(class_labels(ind), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 10)
ylabel('Count', 'fontsize', 12)
p = get(gca,'position');
p(2) = p(2)*2;
set(gca,'position', p)

set(gcf, 'Position', get(0, 'Screensize'))
%%
group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
[~,ia,ib] = intersect(group_table.CNN_classlist, class_labels);
%%
figure
[~,ind] = sort(f1_perclass);
bar(1:numel(ind),f1_perclass(ind), 'edgecolor', 'none')
ylabel('F1-score', 'fontsize', 16)
set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels(ind), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8)
%set(gca, 'position', [.13 .58 .78 .34], 'fontweight', 'bold')
set(gca,'position', [.065 .6 .93 .34], 'fontweight', 'bold')
set(gcf, 'Position', get(0, 'Screensize'))
ax = gca;
ax.YAxis.FontSize = 14;
grid on

figure
[~,ind] = sort(f1_perclass);
bar(1:numel(ind),counts_perclass(ind))
set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels(ind), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8)
ylabel('Count', 'fontsize', 12)
set(gca,'position', [.065 .6 .93 .34], 'fontweight', 'bold')
set(gcf, 'Position', get(0, 'Screensize'))
ax = gca;
ax.YAxis.FontSize = 14;

%%
figure
[~,ind] = sort(f1_perclass,'descend');
bar(1:numel(ind),f1_perclass(ind), 'edgecolor', 'none')
ylabel('F1-score', 'fontsize', 16)
set(gcf, 'position', [325 425 700 200])
xlabel('Class number', 'fontsize', 16)
set(gca,'position', [.088 .24 .9 .7], 'fontweight', 'bold', 'xtick', 0:20:140)
grid on
%%
figure
[~,ind] = sort(f1_perclass, 'descend');
bar(1:numel(ind),f1_perclass(ind), 'edgecolor', 'none')
%set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels(ind), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8)
set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(class_labels(ind), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 6, 'XTickLabelRotationMode', 'auto')
ylabel({'ML classifier performance';'(F1-score)'}, 'fontsize', 12, 'fontweight', 'bold')
set(gcf, 'position', [325 425 700 200])
%xlabel('Class number', 'fontsize', 16)
set(gca,'position', [.065 .4 .93 .5])
set(gcf, 'position', [125 125 900 300])


%%

return

class_ind = ib(find(group_table.Diatom(ia(ind))));
for ii = 1:length(class_ind), ax.XTickLabel{class_ind(ii)} = ['\color{red}' ax.XTickLabel{class_ind(ii)}]; end

class_ind = ib(find(group_table.Dinoflagellate(ia(ind))));
for ii = 1:length(class_ind), ax.XTickLabel{class_ind(ii)} = ['\color{blue}' ax.XTickLabel{class_ind(ii)}]; end

class_ind = ib(find(group_table.Nano(ia(ind))));
for ii = 1:length(class_ind), ax.XTickLabel{class_ind(ii)} = ['\color{green}' ax.XTickLabel{class_ind(ii)}]; end

class_ind = ib(find(group_table.Ciliate(ia(ind))));
for ii = 1:length(class_ind), ax.XTickLabel{class_ind(ii)} = ['\color{cyan}' ax.XTickLabel{class_ind(ii)}]; end

class_ind = ib(find(group_table.Coccolithophore(ia(ind))));
for ii = 1:length(class_ind), ax.XTickLabel{class_ind(ii)} = ['\color{magenta}' ax.XTickLabel{class_ind(ii)}]; end
