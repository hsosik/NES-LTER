%load('training-output\20201120_OTZzooscan\results.mat')
%load('\\vortex\share\otz-data\zooscan\training-output\20201207_OTZzooscan\results')
%load('\\vortex\share\otz-data\zooscan\training-output\20201214_OTZzooscan_jpg\results')

load('\\vortex\share\otz-data\zooscan\training-output\OTZ_zooscan_photic_20220318_1.4\results')
%load('\\vortex\share\otz-data\zooscan\training-output\20210302_OTZzooscan_1.1jpg\results')
%load('\\vortex\share\otz-data\zooscan\training-output\20210302_OTZzooscan_1.2jpg\results')

figure
subplot(1,2,2)
confusion_matrix = confusion_matrix';
imagesc(confusion_matrix)
axis square
for ii = 1:numel(confusion_matrix)
    if confusion_matrix(ii), [a,b] = ind2sub(size(confusion_matrix),ii); 
        text(b,a,num2str(confusion_matrix(ii)), 'HorizontalAlignment', 'center', 'fontsize', 10, 'color', 'w', 'fontweight', 'bold')
    end
end
%set(gca, 'ytick', 1:1:numel(class_labels), 'yticklabel', regexprep(extractAfter(class_labels(classes_by_recall), '_'), '_', ' '), 'yTickLabelRotation', -20, 'fontsize', 10)
%set(gca, 'xtick', 1:1:numel(class_labels), 'xticklabel', regexprep(extractAfter(class_labels(classes_by_recall), '_'), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 10, 'XAxisLocation', 'top')
set(gca, 'ytick', 1:1:numel(class_labels), 'yticklabel', regexprep(class_labels, '_', ''), 'yTickLabelRotation', -20, 'fontsize', 10)
set(gca, 'xtick', 1:1:numel(class_labels), 'xticklabel', regexprep(class_labels, '_', ''), 'XTickLabelRotation', 80, 'fontsize', 10, 'XAxisLocation', 'top')
xlabel('Manual class', 'fontsize', 20)
ylabel('Predicted class', 'fontsize', 20)

subplot(2,2,1)
[~,ind] = sort(f1_perclass);
bar(1:numel(ind),f1_perclass(ind))
%set(gca, 'xtick', 1:1:numel(ind),'xticklabel', []) %, 'xticklabel', regexprep(extractAfter(class_labels(ind), '_'), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 8)
set(gca, 'xtick', 1:1:numel(ind),'xticklabel', []) %, 'xticklabel', regexprep(class_labels(ind), '_', ''), 'XTickLabelRotation', 80, 'fontsize', 8)
ylabel('F1-score', 'fontsize', 12)
set(gca, 'ygrid', 'on')

subplot(2,2,3)
[~,ind] = sort(f1_perclass);
bar(1:numel(ind),counts_perclass(ind))
%set(gca, 'xtick', 1:1:numel(ind), 'xticklabel', regexprep(extractAfter(class_labels(ind), '_'), '_', ' '), 'XTickLabelRotation', 80, 'fontsize', 10)
set(gca, 'xtick', .5:1:numel(ind), 'xticklabel', regexprep(class_labels(ind), '_', ''), 'XTickLabelRotation', 80, 'fontsize', 10)
ylabel('Count', 'fontsize', 12)
p = get(gca,'position');
p(2) = p(2)*2;
set(gca,'position', p)
set(gca, 'ygrid', 'on')

set(gcf, 'Position', get(0, 'Screensize'))