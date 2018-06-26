%% Counts
if ~exist('\\sosiknas1\Lab_data\Attune\EN608\Summary\compiled_stats.mat','file')
    open compile_attune
else
load \\sosiknas1\Lab_data\Attune\EN608\Summary\compiled_stats.mat;
end
%% Synchecoccus Graph
fcs_path = '\\sosiknas1\Lab_data\Attune\EN608\ExportedFCS\'
[~,fcshdr,fcsdatscaled] =fca_readfcs(char(fullfile(fcs_path,fcsfile_syn(2076))));

%this defines the edges of the rectange for synechecoccus
min = 10^2
max =  10 ^6
%for syn count box
xmin= 200
xmax= 10^4
ymin= 200
ymax= 10^5
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];
figure
loglog(fcsdatscaled(:,11),fcsdatscaled(:,19),'k.','HandleVisibility','off')
xlim([min max])
ylim([min max])
hold on
title('\itSynchecoccus')
%title(['PE Signal for ',char(fullfile(fcs_path,fcsfile_syn(2076)))])
in_syn = inpolygon(fcsdatscaled(:,11),fcsdatscaled(:,19),x_rect,y_rect);
syn_count = length(find(in_syn==1));
fsc_signal = fcsdatscaled(:,11);
txt1 = ['Syn: ',num2str(syn_count)];
text(xmin+100,10^5.15,txt1,'Color','r')
xx = fcsdatscaled(:,11);
yy = fcsdatscaled(:,19);
hold on
loglog(xx(in_syn),yy(in_syn),'g.')
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r')
lh = legend('\itSynechococcus');

%Plotting Histogram of the scattering signal
figure
[n, xout] = hist(fsc_signal(find(in_syn==1)),syn_count);
bar(xout, n, 'barwidth', 1, 'basevalue', 0,'FaceColor','g');
xlabel('FSC');
ylabel('Count');
title('FSC vs. Count');
%set(gca,'yscale','log')

%% Small Eukaryotes
[~,fcshdr,fcsdatscaled] =fca_readfcs(char(fullfile(fcs_path,fcsfile_syn(2076))));
ssc_signal = fcsdatscaled(:,12);
y_signal =fcsdatscaled(:,15);
plot_xmin = 10^1;
plot_xmax =  10 ^6;
plot_ymin = 10^2;
plot_ymax = 10^6;
x_polygon = [25  50 10^4 10^6 10^4 900 25];
y_polygon = [300 3500 10^6 10^6 10^4 2000 300];
figure
loglog(fcsdatscaled(:,12),fcsdatscaled(:,15),'k.','HandleVisibility','off')
xlim([plot_xmin plot_xmax]);
ylim([plot_ymin plot_ymax]);
in = inpolygon(fcsdatscaled(:,12),fcsdatscaled(:,15),x_polygon,y_polygon);
hold on
loglog(ssc_signal(in),y_signal(in),'g.')
euk_count = length(find(in==1));
txt1 = ['small Euk: ',num2str(euk_count)];
text(10^2,10^5.15,txt1,'Color','g');
hold on
loglog(ssc_signal(in_syn),y_signal(in_syn),'r.')
hold on
loglog(x_polygon,y_polygon,'LineWidth',1);
title('Chlorophyll Signal for Small Eukaryotes')
lh = legend( 'Small eukaryotes','\itSynechococcus');
%histogram
figure
[n, xout] = hist(ssc_signal(find(in==1)),euk_count);
bar(xout, n, 'barwidth', 0.5, 'basevalue', 1);
set(gca,'yscale','log','xscale','log');
xlabel('SSC');
ylabel('Count');
title('Small Eukaryotes SSC vs. Count');

