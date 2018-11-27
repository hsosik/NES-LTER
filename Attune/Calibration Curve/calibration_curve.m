%% calibration curve
% filename = '\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\calibration_curve_data_attune2.csv';
filename = 'E:\Attune_Data\Size_calibration_July2018\calibration_curve_data_attuneFINAL';
sheet = 'Sheet1';
range = 'C2:D14' ;
curve=xlsread(filename,sheet,range)

name_range = 'B2:B14';
[~,names]  = xlsread(filename,sheet, name_range)
pnames = {'\itPavlova Lutheri';...
    '\itPyrmnesium Parvum';'\itAlexandrium Minutum';...
    '\itPeridinium';'\itScrippsiella Trochoidea';'\itGymnodium galatheanum';...
    '\itHeterosigman akashiwo';'\itHeterocapsa Triquetra';'\itIsochrysis Galbana';...
    '\itHeterosigman Akashiwo';'\itHeterocapsa Triquetra';'\itCryptophora';'\itSynechecoccus'}
ssc = curve(:,1)
cc = curve(:,2)

vol = (4/3)*pi.*(cc./2).^3;

ssc_log = log10(ssc);
vol_log = log10(vol);
figure
plot(ssc_log,vol_log,'o','MarkerSize',17,'MarkerFaceColor','g')
xlabel('Log10 of SSC-A','FontSize',20)
ylabel('Log10 of CC Vol (um^3)','FontSize',20)
set(gca, 'FontSize', 20)
set(gcf, 'Position', get(0, 'Screensize'));
% legend(pnames)
pbaspect([1 1 1])
volume = (1.3.*ssc)-2.9;

%%
% Open a new figure
figure;
% Turn hold on so all the points can 
% be plotted individually
hold on;
% Calculation for the amount by which the
% label should be displaced in the 'y'
% direction

map = brewermap(14,'Set1'); 
% plot and label the individual points
for i = 1:length(ssc_log)
     plot(ssc_log(i),vol_log(i),'o','MarkerSize',20,'MarkerFaceColor',map(i+1,:),'MarkerEdgeColor',map(i,:));
     pbaspect([1 1 1])
     % Label the points with the index
%      text(ssc_log(i)+0.04,vol_log(i)+0.06,pnames(i),'FontSize',16);
end
legend(pnames,'Location','eastoutside')
set(gca,'FontSize',24)
set(gcf, 'Position', get(0, 'Screensize'));
hold on
x = 1:1:6
y = 1.3.*x-2.9
plot(y,'LineWidth',3,'Color','k','HandleVisibility','off')
xlabel('Log10 of Side Scattering','FontSize', 25)
ylabel('Log10 of CC Vol (um^3)','FontSize', 25)
dim = [0.2,0.75,0.1,0.1];
annotation('textbox',dim,'String','y = 1.3*x - 2.9','FontSize',36,'LineStyle','none');
