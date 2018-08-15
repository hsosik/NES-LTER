%plot_attune

%check gating for syn
%check gating for small euks
%attune cell conc
%attune biovolume
%ifcb cell conc
%ifcb biovolume
%biovolume from both instruments (no overlap)

basepath = '\\sosiknas1\Lab_data\Attune\EN608';
load([basepath '\Summary\Attune'])


%% check gating for Syn
% [~,fcshdr,fcsdatscaled] =fca_readfcs(filename);
% filename = 'E:\Attune_Data\EN608\ExportedFCS\NESLTER_EN608_31Jan2018B_Group_day0_Sample(1).fcs';
filename = '\\sosiknas1\Lab_data\Attune\EN608\ExportedFCS\NESLTER_EN608_31Jan2018B_Group_day0_Sample(1).fcs';

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

%% Check gates for Small Eukaryotes
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

%% Attune Cell Concentration

attunefig = figure;
[~,ii] = sort(Attune.FCSfileinfo.matdate_start)
att_conc.lesstwo = (Attune.Count.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_conc.twoten = (Attune.Count.twoten)./(Attune.vol_analyzed.*10^(-6));
att_conc.tentwen = (Attune.Count.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_conc.twen = (Attune.Count.twen)./(Attune.vol_analyzed.*10^(-6));
att_conc.syn = (Attune.Count.SynTotal)./(Attune.vol_analyzed.*10^(-6));


plot(Attune.fcsmatch.lat,att_conc.lesstwo,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.lat,att_conc.twoten,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
plot(Attune.fcsmatch.lat,att_conc.tentwen,'.','MarkerSize',25,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
plot(Attune.fcsmatch.lat,att_conc.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
plot(Attune.fcsmatch.lat,att_conc.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
set(gca, 'YScale', 'log','FontSize',24)

lh = legend('<2 \mum', '2-10 \mum','10-20 \mum','>20 \mum','\itSynechecoccus')
title(lh,'Size Classes')
title(['Concentration of Cells by Size Class during Cruise ' Attune.cruiseName],'FontSize',24)
xlabel('Latitude (Degrees North)')
ylabel('[Cell Concentration] (mL^{-1})')
set(gcf, 'Position', get(0, 'Screensize'));

savefig(attunefig,[basepath '\Figures\AttuneCellConc.fig'])
%% Attune Biovolume

attunefig = figure;
[~,ii] = sort(Attune.FCSfileinfo.matdate_start)
att_bv.lesstwo = (Attune.Biovol.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_bv.twoten = (Attune.Biovol.twoten)./(Attune.vol_analyzed.*10^(-6));
att_bv.tentwen = (Attune.Biovol.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_bv.twen = (Attune.Biovol.twen)./(Attune.vol_analyzed.*10^(-6));
att_bv.syn = (Attune.Biovol.Syn)./(Attune.vol_analyzed.*10^(-6));


plot(Attune.fcsmatch.lat,att_bv.lesstwo,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.lat,att_bv.twoten,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
plot(Attune.fcsmatch.lat,att_bv.tentwen,'.','MarkerSize',25,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
plot(Attune.fcsmatch.lat,att_bv.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
plot(Attune.fcsmatch.lat,att_bv.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
set(gca, 'YScale', 'log','FontSize',24)

lh = legend('<2 \mum', '2-10 \mum','10-20 \mum','>20 \mum','\itSynechecoccus')
title(lh,'Size Classes')
title(['Biovolume of Cells by Size Class during Cruise ' Attune.cruiseName],'FontSize',24)
xlabel('Latitude (Degrees North)')
ylabel('[Biovolume (\mum^{3}mL^{-1})')
set(gcf, 'Position', get(0, 'Screensize'));

savefig(attunefig,[basepath '\Figures\AttuneCellConc.fig'])

%% IFCB
basepath = '\\sosiknas1\IFCB_products\NESLTER_transect\summary\'
load([basepath 'IFCB_biovolume_size_classes_manual_14Aug2018'])
load([basepath 'count_biovol_size_manual_14Aug2018'])

%% IFCB Cell Concentration

N0_10_Conc = N0_10_phyto./ml_analyzed
N10_20_Conc = N10_20_phyto./ml_analyzed
N20_inf_Conc = N20_inf_phyto./ml_analyzed

figure
plot(ifcblat,N0_10_Conc,'.','MarkerSize',10,'MarkerEdgeColor',[0.2,0,0.2],'MarkerFaceColor',[0.2,0,0.2])
hold on
plot(ifcblat,N10_20_Conc,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.3,0.5],'MarkerFaceColor',[0,0.3,0.5])
plot(ifcblat,N20_inf_Conc,'.','MarkerSize',30,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])
ylabel('Cell Concentrat0ion (ml^{-1})')
xlabel('Latitude(Degrees North)')
lh = legend('0-10 \mum','10-20 \mum', '>= 20 \mum')
set(gca,'YScale','log')


%% IFCB Biovolume

biovol0_10_Conc = biovol0_10_phyto./ml_analyzed
biovol10_20_Conc = biovol10_20_phyto./ml_analyzed
biovol20_inf_Conc = biovol20_inf_phyto./ml_analyzed

figure
plot(ifcblat,biovol0_10_Conc,'.','MarkerSize',10,'MarkerEdgeColor',[0.2,0,0.2],'MarkerFaceColor',[0.2,0,0.2])
hold on
plot(ifcblat,biovol10_20_Conc,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.3,0.5],'MarkerFaceColor',[0,0.3,0.5])
plot(ifcblat,biovol20_inf_Conc,'.','MarkerSize',30,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])
lh = legend('0-10 \mum','10-20 \mum', '>= 20 \mum')
ylabel('Biovolume (\mum^{3}ml^{-1})')
xlabel('Latitude(Degrees North)')
set(gca,'YScale','log')

%% Combined biovol w/out overlap vs. Latitude
figure('units','normalized','outerposition',[0 0 1 1])
att_bv.lesstwo = (Attune.Biovol.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_bv.twoten = (Attune.Biovol.twoten)./(Attune.vol_analyzed.*10^(-6));
att_bv.tentwen = (Attune.Biovol.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_bv.twen = (Attune.Biovol.twen)./(Attune.vol_analyzed.*10^(-6));
att_bv.syn = (Attune.Biovol.Syn)./(Attune.vol_analyzed.*10^(-6));

plot(Attune.fcsmatch.lat,att_bv.lesstwo,'.','MarkerSize',5,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.lat,att_bv.twoten,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
% plot(Attune.fcsmatch.lat,att_bv.tentwen,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
% plot(Attune.fcsmatch.lat,att_bv.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
plot(Attune.fcsmatch.lat,att_bv.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])

biovol0_10_Conc = biovol0_10_phyto./ml_analyzed
biovol10_20_Conc = biovol10_20_phyto./ml_analyzed
biovol20_inf_Conc = biovol20_inf_phyto./ml_analyzed


% plot(ifcblat,biovol0_10_Conc,'.','MarkerSize',10,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

plot(ifcblat,biovol10_20_Conc,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.3,0.5],'MarkerFaceColor',[0,0.3,0.5])
plot(ifcblat,biovol20_inf_Conc,'.','MarkerSize',30,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

ylabel('Biovolume (\mum^{3}ml^{-1})')
xlabel('Latitude(Degrees North)')
set(gca, 'YScale', 'log','FontSize',24)

lh = legend('<2 \mum (attune)', '2-10 \mum (attune)','\itSynechecoccus (attune)','10-20 \mum (IFCB)','>20 \mum (IFCB)','Location','eastoutside')
title(lh,'Size Classes')


%% Combined biovol w/ overlap vs. Latitude
figure('units','normalized','outerposition',[0 0 1 1])
att_bv.lesstwo = (Attune.Biovol.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_bv.twoten = (Attune.Biovol.twoten)./(Attune.vol_analyzed.*10^(-6));
att_bv.tentwen = (Attune.Biovol.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_bv.twen = (Attune.Biovol.twen)./(Attune.vol_analyzed.*10^(-6));
att_bv.syn = (Attune.Biovol.Syn)./(Attune.vol_analyzed.*10^(-6));

plot(Attune.fcsmatch.lat,att_bv.lesstwo,'.','MarkerSize',5,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.lat,att_bv.twoten,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
plot(Attune.fcsmatch.lat,att_bv.tentwen,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
% plot(Attune.fcsmatch.lat,att_bv.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
plot(Attune.fcsmatch.lat,att_bv.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])

biovol0_10_Conc = biovol0_10_phyto./ml_analyzed
biovol10_20_Conc = biovol10_20_phyto./ml_analyzed
biovol20_inf_Conc = biovol20_inf_phyto./ml_analyzed


plot(ifcblat,biovol0_10_Conc,'.','MarkerSize',10,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

plot(ifcblat,biovol10_20_Conc,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.3,0.5],'MarkerFaceColor',[0,0.3,0.5])
plot(ifcblat,biovol20_inf_Conc,'.','MarkerSize',30,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

ylabel('Biovolume (\mum^{3}ml^{-1})')
xlabel('Latitude(Degrees North)')
set(gca, 'YScale', 'log','FontSize',24)

lh = legend('<2 \mum', '2-10 \mum','10-20 \mum','\itSynechecoccus','0-10 \mum','10-20 \mum','>20 \mum','Location','eastoutside')
title(lh,'Size Classes')


%% Combined biovol w/ overlap vs. time
figure('units','normalized','outerposition',[0 0 1 1])
att_bv.lesstwo = (Attune.Biovol.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_bv.twoten = (Attune.Biovol.twoten)./(Attune.vol_analyzed.*10^(-6));
att_bv.tentwen = (Attune.Biovol.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_bv.twen = (Attune.Biovol.twen)./(Attune.vol_analyzed.*10^(-6));
att_bv.syn = (Attune.Biovol.Syn)./(Attune.vol_analyzed.*10^(-6));

plot(Attune.fcsmatch.mdate_start,att_bv.lesstwo,'.','MarkerSize',5,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.mdate_start,att_bv.twoten,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
% plot(Attune.fcsmatch.mdate_start,att_bv.tentwen,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
% plot(Attune.fcsmatch.lat,att_bv.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
plot(Attune.fcsmatch.mdate_start,att_bv.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])

biovol0_10_Conc = biovol0_10_phyto./ml_analyzed
biovol10_20_Conc = biovol10_20_phyto./ml_analyzed
biovol20_inf_Conc = biovol20_inf_phyto./ml_analyzed


% plot(matdate,biovol0_10_Conc,'.','MarkerSize',10,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

plot(matdate,biovol10_20_Conc,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.3,0.5],'MarkerFaceColor',[0,0.3,0.5])
plot(matdate,biovol20_inf_Conc,'.','MarkerSize',30,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

ylabel('Biovolume (\mum^{3}ml^{-1})')
xlabel('Latitude(Degrees North)')
set(gca, 'YScale', 'log','FontSize',24)

lh = legend('<2 \mum (attune)', '2-10 \mum (attune)','\itSynechecoccus (attune)','10-20 \mum (IFCB)','>20 \mum (IFCB)','Location','eastoutside')
title(lh,'Size Classes')

%% Combined biovol w/ overlap vs. time
figure('units','normalized','outerposition',[0 0 1 1])
att_bv.lesstwo = (Attune.Biovol.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_bv.twoten = (Attune.Biovol.twoten)./(Attune.vol_analyzed.*10^(-6));
att_bv.tentwen = (Attune.Biovol.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_bv.twen = (Attune.Biovol.twen)./(Attune.vol_analyzed.*10^(-6));
att_bv.syn = (Attune.Biovol.Syn)./(Attune.vol_analyzed.*10^(-6));

plot(Attune.fcsmatch.mdate_start,att_bv.lesstwo,'.','MarkerSize',5,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.mdate_start,att_bv.twoten,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
plot(Attune.fcsmatch.mdate_start,att_bv.tentwen,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
% plot(Attune.fcsmatch.lat,att_bv.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
plot(Attune.fcsmatch.mdate_start,att_bv.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])

biovol0_10_Conc = biovol0_10_phyto./ml_analyzed
biovol10_20_Conc = biovol10_20_phyto./ml_analyzed
biovol20_inf_Conc = biovol20_inf_phyto./ml_analyzed


plot(matdate,biovol0_10_Conc,'.','MarkerSize',10,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

plot(matdate,biovol10_20_Conc,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.3,0.5],'MarkerFaceColor',[0,0.3,0.5])
plot(matdate,biovol20_inf_Conc,'.','MarkerSize',30,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

ylabel('Biovolume (\mum^{3}ml^{-1})')
xlabel('Latitude(Degrees North)')
set(gca, 'YScale', 'log','FontSize',24)

lh = legend('<2 \mum', '2-10 \mum','10-20 \mum','\itSynechecoccus','0-10 \mum','10-20 \mum','>20 \mum','Location','eastoutside')
title(lh,'Size Classes')

%% Combined Concentration w/out overlap vs. lat

attunefig = figure;
[~,ii] = sort(Attune.FCSfileinfo.matdate_start)
att_conc.lesstwo = (Attune.Count.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_conc.twoten = (Attune.Count.twoten)./(Attune.vol_analyzed.*10^(-6));
att_conc.tentwen = (Attune.Count.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_conc.twen = (Attune.Count.twen)./(Attune.vol_analyzed.*10^(-6));
att_conc.syn = (Attune.Count.SynTotal)./(Attune.vol_analyzed.*10^(-6));


plot(Attune.fcsmatch.lat,att_conc.lesstwo,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.lat,att_conc.twoten,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
% plot(Attune.fcsmatch.lat,att_conc.tentwen,'.','MarkerSize',25,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
% plot(Attune.fcsmatch.lat,att_conc.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
% plot(Attune.fcsmatch.lat,att_conc.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
set(gca, 'YScale', 'log','FontSize',24)

title(['Concentration of Cells by Size Class during Cruise ' Attune.cruiseName],'FontSize',24)
xlabel('Latitude (Degrees North)')
ylabel('[Cell Concentration] (mL^{-1})')
set(gcf, 'Position', get(0, 'Screensize'));


N0_10_Conc = N0_10_phyto./ml_analyzed
N10_20_Conc = N10_20_phyto./ml_analyzed
N20_inf_Conc = N20_inf_phyto./ml_analyzed

% plot(ifcblat,N0_10_Conc,'.','MarkerSize',10,'MarkerEdgeColor',[0.2,0,0.2],'MarkerFaceColor',[0.2,0,0.2])
hold on
plot(ifcblat,N10_20_Conc,'.','MarkerSize',30,'MarkerEdgeColor',[0,0,0.5],'MarkerFaceColor',[0,0,0.5])
plot(ifcblat,N20_inf_Conc,'.','MarkerSize',40,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])
ylabel('Cell Concentrat0ion (ml^{-1})')
xlabel('Latitude(Degrees North)')
lh = legend('<2\mum (Attune)','2-10 \mum (Attune)','10-20 \mum', '>= 20 \mum','Location','eastoutside')
set(gca,'YScale','log')


%% Combined Concentration w overlap vs. lat

attunefig = figure;
[~,ii] = sort(Attune.FCSfileinfo.matdate_start)
att_conc.lesstwo = (Attune.Count.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_conc.twoten = (Attune.Count.twoten)./(Attune.vol_analyzed.*10^(-6));
att_conc.tentwen = (Attune.Count.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_conc.twen = (Attune.Count.twen)./(Attune.vol_analyzed.*10^(-6));
att_conc.syn = (Attune.Count.SynTotal)./(Attune.vol_analyzed.*10^(-6));


plot(Attune.fcsmatch.lat,att_conc.lesstwo,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.lat,att_conc.twoten,'*','MarkerSize',5,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
plot(Attune.fcsmatch.lat,att_conc.tentwen,'.','MarkerSize',5,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
% plot(Attune.fcsmatch.lat,att_conc.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
plot(Attune.fcsmatch.lat,att_conc.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
set(gca, 'YScale', 'log','FontSize',24)

title(['Concentration of Cells by Size Class during Cruise ' Attune.cruiseName],'FontSize',24)
xlabel('Latitude (Degrees North)')
ylabel('[Cell Concentration] (mL^{-1})')
set(gcf, 'Position', get(0, 'Screensize'));


N0_10_Conc = N0_10_phyto./ml_analyzed
N10_20_Conc = N10_20_phyto./ml_analyzed
N20_inf_Conc = N20_inf_phyto./ml_analyzed

plot(ifcblat,N0_10_Conc,'*','MarkerSize',5,'MarkerEdgeColor',[0.2,0,0.2],'MarkerFaceColor',[0,0,0.2])
hold on
plot(ifcblat,N10_20_Conc,'.','MarkerSize',5,'MarkerEdgeColor',[0,0,0.5],'MarkerFaceColor',[0,0,0.5])
plot(ifcblat,N20_inf_Conc,'.','MarkerSize',40,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])
ylabel('Cell Concentrat0ion (ml^{-1})')
xlabel('Latitude(Degrees North)')
lh = legend('<2\mum (Attune)','2-10 \mum (Attune)','10-20 \mum (Attune)','\itSyn (attune)','0-10\mum (ifcb)','10-20 \mum (IFCB)', '>= 20 \mum (IFCB)','Location','eastoutside')
set(gca,'YScale','log')


%% combined conce w/overlap vs. time

attunefig = figure;
[~,ii] = sort(Attune.FCSfileinfo.matdate_start)
att_conc.lesstwo = (Attune.Count.lesstwo)./(Attune.vol_analyzed.*10^(-6));
att_conc.twoten = (Attune.Count.twoten)./(Attune.vol_analyzed.*10^(-6));
att_conc.tentwen = (Attune.Count.tentwen)./(Attune.vol_analyzed.*10^(-6));
att_conc.twen = (Attune.Count.twen)./(Attune.vol_analyzed.*10^(-6));
att_conc.syn = (Attune.Count.SynTotal)./(Attune.vol_analyzed.*10^(-6));


plot(Attune.fcsmatch.mdate_start,att_conc.lesstwo,'.','MarkerSize',10,'MarkerEdgeColor',[0,0.2,0],'MarkerFaceColor',[0,0.2,0])
hold on
plot(Attune.fcsmatch.mdate_start,att_conc.twoten,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.5,0],'MarkerFaceColor',[0,0.5,0])
plot(Attune.fcsmatch.mdate_start,att_conc.tentwen,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.7,0],'MarkerFaceColor',[0,0.7,0])
% plot(Attune.fcsmatch.mdate_start,att_conc.twen,'.','MarkerSize',35,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0])
% plot(Attune.fcsmatch.mdate_start,att_conc.syn,'.','MarkerSize',10,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
set(gca, 'YScale', 'log','FontSize',24)

title(['Concentration of Cells by Size Class during Cruise ' Attune.cruiseName],'FontSize',24)
xlabel('Latitude (Degrees North)')
ylabel('[Cell Concentration] (mL^{-1})')
set(gcf, 'Position', get(0, 'Screensize'));


N0_10_Conc = N0_10_phyto./ml_analyzed
N10_20_Conc = N10_20_phyto./ml_analyzed
N20_inf_Conc = N20_inf_phyto./ml_analyzed

plot(matdate,N0_10_Conc,'.','MarkerSize',20,'MarkerEdgeColor',[0.2,0,0.2],'MarkerFaceColor',[0.2,0,0.2])
hold on
plot(matdate,N10_20_Conc,'.','MarkerSize',30,'MarkerEdgeColor',[0,0,0.5],'MarkerFaceColor',[0,0,0.5])
plot(matdate,N20_inf_Conc,'.','MarkerSize',40,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])
ylabel('Cell Concentrat0ion (ml^{-1})')
xlabel('Time')
lh = legend('<2\mum (Attune)','2-10 \mum (Attune)','10-20 \mum (Attune)','0-10\mum (IFCB)','10-20 \mum (IFCB)', '>= 20 \mum (IFCB)','Location','eastoutside')
set(gca,'YScale','log')
