% function [fcsplot] = fcs_plot_signals(singlefcsfilepath)

%basepath = '\\sosiknas1\Lab_data\Attune\EN608';
%basepath = '\\sosiknas1\Backup\SPIROPA\20180414_AR29\Attune\';
basepath= '\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS';
index = 1;
fcsfil

%% Loading in data
if ~exist([basepath '\Summary\compiled_stats.mat'],'file')
    open compile_attune
else
load([basepath '\Summary\compiled_stats.mat'])
end

fcs_path =  [basepath '\ExportedFCS\']

%%
filename = 'E:\Attune_Data\EN608\ExportedFCS\NESLTER_EN608_02Feb2018C_Group_day0_Sample(10).fcs'

%% Synchecoccus Graph
% [~,fcshdr,fcsdatscaled] =fca_readfcs(filename);
filename = 'E:\Attune_Data\EN608\ExportedFCS\NESLTER_EN608_31Jan2018B_Group_day0_Sample(1).fcs';
[~,fcshdr,fcsdatscaled] =fca_readfcs(filename);
%this defines the edges of the rectange for synechecoccus
min = 10^2
max =  10^6
%for syn count box
xmin= 200
xmax= 10^4
ymin= 10^3
ymax= 10^5
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

figure
loglog(fcsdatscaled(:,11),fcsdatscaled(:,19),'k.','HandleVisibility','off')
xlim([min max])
ylim([min max])
hold on
title('\itSynchecoccus')
in_syn = inpolygon(fcsdatscaled(:,11),fcsdatscaled(:,19),x_rect,y_rect);
syn_count = length(find(in_syn==1));
fsc_signal = fcsdatscaled(:,11);
txt1 = ['Syn: ',num2str(syn_count)];
text(xmin+100,10^5.15,txt1,'Color','r')
xx = fcsdatscaled(:,11);
yy = fcsdatscaled(:,19);
hold on
loglog(xx(in_syn),yy(in_syn),'r.')
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itSynechococcus');
xlabel('Forward Scattering')
ylabel('Phycoerythrin')

%Plotting Histogram of the scattering signal
figure
[n, xout] = hist(fsc_signal(find(in_syn==1)),syn_count);
bar(xout, n, 'barwidth', 1, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Forward Scattering');
ylabel('Count');
title('Histogram of Forward Scattering');
%set(gca,'yscale','log')

%% Small Eukaryotes
[~,fcshdr,fcsdatscaled] =fca_readfcs(filename);
ssc_signal = fcsdatscaled(:,12);
y_signal =fcsdatscaled(:,15);

plot_xmin = 10^1;
plot_xmax = 10 ^6;
plot_ymin = 10^2;
plot_ymax = 10^6;

x_polygon = [25  50 10^4 10^6 10^6 10^5 10^4 25];
y_polygon = [300 3500 10^6 10^6 10^5 10^4 10^3 300];

figure
loglog(fcsdatscaled(:,12),fcsdatscaled(:,15),'k.','HandleVisibility','off')
xlim([plot_xmin plot_xmax]);
ylim([plot_ymin plot_ymax]);

in_euk = inpolygon(ssc_signal,y_signal,x_polygon,y_polygon);

hold on
loglog(ssc_signal(in_euk),y_signal(in_euk),'.','Color',[0 0.75 0])
xlim([plot_xmin plot_xmax]);
ylim([plot_ymin plot_ymax]);


euk_count = length(find(in_euk==1));
txt1 = ['small Euk: ',num2str(euk_count)];
text(10^2,10^5.15,txt1,'Color',[0 0.75 0]);

hold on
loglog(ssc_signal(in_syn),y_signal(in_syn),'r.')
xlim([plot_xmin plot_xmax]);
ylim([plot_ymin plot_ymax]);

hold on
loglog(x_polygon,y_polygon,'LineWidth',1,'LineStyle','--','Color',[0 0.75 0]);
xlim([plot_xmin plot_xmax]);
ylim([plot_ymin plot_ymax]);
xlabel('Side Scattering');
ylabel('GL-1');
title('Chlorophyll Signal for Small Eukaryotes')
lh = legend( 'Small eukaryotes','\itSynechococcus');

%histogram
figure
[n, xout] = hist(ssc_signal(find(in_euk==1)),euk_count);
bar(xout, n, 'barwidth', 0.5, 'basevalue', 1,'FaceColor',[0 0.75 0]);
set(gca,'yscale','log','xscale','log');
xlabel('Side Scattering');
ylabel('Count');
title('Histogram of Side Scattering for Small Eukaryotes');


