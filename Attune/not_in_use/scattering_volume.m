%get cell concentratiosn for 2-10 microns and then 10 -20 microns and then
%2 microns
% cell concentration of 2-10, 10 -20 and < 10 and greater than 20
%biovolume across the cruise
%biovolume
%refer to compile_attune as a function
%fcsdat vs. fcsdatscaled?

%basepath ='\\sosiknas1\Lab_data\Attune\EN608\';%at the lab
basepath='E:\Attune_Data\EN608\';

%outpath= '\\sosiknas1\Lab_data\Attune\EN608\Summary'
% basepath = '\\sosiknas1\Backup\SPIROPA\20180414_AR29\Attune\';
% basepath = '';%hard drive
% basepath = '';%
% basepath = '';%
compile_attune(basepath)

%% Converting Scattering to Effective Spherical Diameter for Small Eukaryotes
basepath = 'E:\Attune_Data\EN608\';
% basepath = '\\sosiknas1\Backup\LTER\20180404_AR28\'
fpath = [basepath '\ExportedFCS\'];
outpath = [basepath '\Summary\'];

% Extracting files out of the directory sorts NES out from SFD
%first it will populate with NES titled files but if empty will go for SFD
%file string
filelist = dir([fpath 'NES*']);
filelist = {filelist.name}';
flistchar = char(filelist);
dstr = flistchar(:,15:end-28);
mdate = datenum(dstr);
[~,s] = sort(mdate);
sortedlist = filelist(s);

Count.lesstwo = [];
Count.twoten = [];
Count.tentwen =[];
Count.twen = [];
Count.vol_analyzed = [];
Biovol.lesstwo = [];
Biovol.twoten = [];
Biovol.tentwen =[];
Biovol.twen = [];
for count = 1:length(sortedlist)
    disp(sortedlist{count});
    filename = [fpath sortedlist{count}];
    [~,fcshdr,fcsdatscaled] =fca_readfcs(filename);
    Count.vol_analyzed = [Count.vol_analyzed; fcshdr.VOL];
% filelist = filelist(s);
   ssc_signal = fcsdatscaled(:,3);
   y_signal =fcsdatscaled(:,15);
%    plot_xmin = 10^1;
%    plot_xmax = 10 ^6;
%    plot_ymin = 10^2;
%    plot_ymax = 10^6;
   x_polygon = [25  50 10^4 10^6 10^6 10^5 10^4 25];
   y_polygon = [300 3500 10^6 10^6 10^5 10^4 10^3 300];
%    figure
%    loglog(ssc_signal,y_signal,'k.','HandleVisibility','off')
%    xlim([plot_xmin plot_xmax]);
%    ylim([plot_ymin plot_ymax]);
   in_euk = inpolygon(ssc_signal,y_signal,x_polygon,y_polygon);
   ssc = log10(ssc_signal(find(in_euk==1)));
   y = y_signal(find(in_euk==1));
%    hold on
%    loglog(ssc,y,'.','Color',[0 0.75 0])
%    xlim([plot_xmin plot_xmax]);
%    ylim([plot_ymin plot_ymax]);
%    euk_count = length();
%    txt1 = ['small Euk: ',num2str(euk_count)];
%    text(10^2,10^5.15,txt1,'Color',[0 0.75 0]);
%    hold on
%    loglog(x_polygon,y_polygon,'LineWidth',1,'LineStyle','--','Color',[0 0.75 0]);
%    xlim([plot_xmin plot_xmax]);
%    ylim([plot_ymin plot_ymax]);
%    xlabel('Side Scattering');
%    ylabel('BL3-H');
   volume = 1.3.*ssc - 2.9;
   lin_vol = 10.^(volume);
   diameter = ((lin_vol/pi).*(3/4)).^(1/3);
   
 size2 = [];
 size2_10 =[];
 size10_20 = [];
 size20 = [];
 for ii = 1:length(diameter)
    if diameter(ii) <= 2
        size2 = [size2;diameter(ii)];
    elseif diameter(ii)>= 2 & diameter(ii) <= 10
        size2_10 = [size2_10;diameter(ii)];
    elseif  diameter(ii)>= 10 & diameter(ii) <= 20
        size10_20 =[size10_20; diameter(ii)];
    elseif diameter(ii) >= 20 
        size20 =[ size20;diameter(ii)];
    end
 end
   Biovol.lesstwo = [Biovol.lesstwo; sum((4/3).*pi.*((size2./2).^3))];
   Biovol.twoten = [Biovol.twoten; sum((4/3).*pi.*(size2_10./2).^3)];
   Biovol.tentwen =[Biovol.tentwen; sum((4/3).*pi.*(size10_20./2).^3)];
   Biovol.twen = [Biovol.twen; sum((4/3).*pi.*(size20./2).^3)];
   Count.lesstwo = [Count.lesstwo ; length(size2)];
   Count.twoten = [Count.twoten ; length(size2_10)];
   Count.tentwen =[Count.tentwen ; length(size10_20)];
   Count.twen = [Count.twen ; length(size20)];
end

%% Synechecoccus
basepath = 'E:\Attune_Data\EN608\';
%basepath = 'E:\Attune_data\AR28\';
%basepath = 'E:\Attune_data\AR29\';
%basepath = 'E:\Attune_data\EN617\';

fpath = [basepath '\ExportedFCS\'];
outpath = [basepath '\Summary\'];

% Extracting files out of the directory sorts NES out from SFD
%first it will populate with NES titled files but if empty will go for SFD
%file string
filelist = dir([fpath 'NES*']);
filelist = {filelist.name}'
flistchar = char(filelist);
dstr = flistchar(:,15:end-28);
mdate = datenum(dstr);
[~,s] = sort(mdate);
sortedlist = filelist(s);

Count.Syn = [];
Biovol.Syn = [];
Count.vol_analyzed = [];
for count = 1:length(sortedlist)
    disp(sortedlist{count});
    filename = [fpath sortedlist{count}];
    [~,fcshdr,fcsdatscaled] =fca_readfcs(filename);
    Count.vol_analyzed = [Count.vol_analyzed; fcshdr.VOL];
    %this sets up the 4 vertices of the gate for Synechecoccus
    xmin= 200;
    xmax= 10^4;
    ymin= 10^3;
    ymax= 10^5;
    x_rect = [xmin xmin xmax xmax xmin];
    y_rect = [ymin ymax ymax ymin ymin];
    
   ssc_signal = fcsdatscaled(:,11);
   y_signal =fcsdatscaled(:,19);
   in_syn = inpolygon(ssc_signal,y_signal,x_rect,y_rect);
   ssc = log10(ssc_signal(find(in_syn ==1)));
   volume = 1.3.*ssc - 2.9;
   lin_vol = 10.^(volume);
   diameter = ((lin_vol/pi).*(3/4)).^(1/3);
   Count.Syn = [Count.Syn; length(ssc)];
   Biovol.Syn = [Biovol.Syn; sum((4/3).*pi.*(diameter./2).^3)];

%    save(outpath,Count)
end 

%%
filename = 'E:\Attune_Data\EN608\ExportedFCS\NESLTER_EN608_31Jan2018C_Group_day0_Sample(60).fcs';
[~,fcshdr,fcsdatscaled] =fca_readfcs(filename);

%this defines the edges of the rectange for synechecoccus
min = 10^2;
max =  10^6;
%for syn count box
xmin= 200;
xmax= 10^4;
ymin= 10^3;
ymax= 10^5;
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

figure
loglog(fcsdatscaled(:,11),fcsdatscaled(:,19),'k.','HandleVisibility','off')
xlim([min max])
ylim([min max])
hold on
title('\itSynchecoccus')
in_syn = inpolygon(fcsdatscaled(:,3),fcsdatscaled(:,19),x_rect,y_rect);
syn_count = length(find(in_syn==1));
fsc_signal = fcsdatscaled(:,3);
txt1 = ['Syn: ',num2str(syn_count)];
text(xmin+100,10^5.15,txt1,'Color','r')
xx = fcsdatscaled(:,3);
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
%%

figure
loglog(fcsdatscaled(:,11),fcsdatscaled(:,19),'k.','HandleVisibility','off')
xlim([min max])
ylim([min max])
hold on

syn_count = length(find(in_syn==1));
fsc_signal = fcsdatscaled(:,11)
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
%%
%Plotting Histogram of the scattering signal
figure
[n, xout] = hist(fsc_signal(find(in_syn==1)),syn_count);
bar(xout, n, 'barwidth', 1, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Forward Scattering');
ylabel('Count');
title('Histogram of Forward Scattering');

%% ATTUNE CONCENTRATION
figure
markersize = 20;
% att_conc.lesstwo = (Count.lesstwo)./(Count.vol_analyzed.*10^(-6));
% att_conc.twoten = (Count.twoten)./(Count.vol_analyzed.*10^(-6));
att_conc.syn =(Count.Syn)./(Count.vol_analyzed.*10^(-6));
% plot(attunelat,att_conc.lesstwo,'.','MarkerSize',15,'MarkerEdgeColor',[0,0.2,0])
% hold on
% plot(attunelat,att_conc.twoten,'.','MarkerSize',20,'MarkerEdgeColor',[0,0.6,0],'MarkerFaceColor',[0,0.2,0])
% plot(attunelat,(Count.tentwen./(Count.vol_analyzed.*10^(-6))),'.','MarkerSize',30,'MarkerEdgeColor',[0,1,0])
% plot(attunelat,(Count.twen./(Count.vol_analyzed.*10^(-6))),'.','MarkerSize',30,'MarkerEdgeColor',[0,1,0])
plot(attunelat,att_conc.syn,'.','MarkerSize',markersize,'MarkerEdgeColor','k')
set(gca, 'YScale', 'log','FontSize',19)

lh = legend('Synechecoccus')
% title(lh,'Size Classes')
% title('Cncentration of Cells by Size Class on the North-South LTER Transect in Winter','FontSize',24)
xlabel('Latitude (Degrees)','FontSize',24)
ylabel('[Cell Concentration] (mL^{-1})','FontSize',19)
set(gcf, 'Position', get(0, 'Screensize'));


%% ATTUNE BIOVOLUME
figure
index = find(attunelon >= -70.9 & attunelon <= -70.87)'
markersize = 20;
% plot(attunelat(index),Biovol.lesstwo(index)./(Count.vol_analyzed(index).*10^(-6)),'.','MarkerSize',15,'MarkerEdgeColor',[0,0.2,0])
hold on
% plot(attunelat(index),Biovol.twoten(index)./(Count.vol_analyzed(index).*10^(-6)),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0])
% plot(attunelat(index),Biovol.tentwen(index)./(Count.vol_analyzed(index).*10^(-6)),'.','MarkerSize',30,'MarkerEdgeColor',[0,1,0])
plot(attunelat,Biovol.Syn/Count.vol_analyzed,'.','MarkerEdgeColor','k','MarkerSize',20)
set(gca, 'YScale', 'log')
% xlim([10^(-3) 10^(-1)])

lh = legend('Synechecoccus')

xlabel('Latitude (Degrees)','FontSize',24)
ylabel('Biovolume (\mum^{3}mL^{-1})','FontSize',24)
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',24)


%% Example Code
%filename = [basepath filelist(2).name];
filename = 'E:\Attune_Data\EN608\ExportedFCS\NESLTER_EN608_31Jan2018C_Group_day0_Sample(63).fcs'
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

% hold on
% loglog(ssc_signal(in_syn),y_signal(in_syn),'r.')
% xlim([plot_xmin plot_xmax]);
% ylim([plot_ymin plot_ymax]);

hold on
loglog(x_polygon,y_polygon,'LineWidth',1,'LineStyle','--','Color',[0 0.75 0]);
xlim([plot_xmin plot_xmax]);
ylim([plot_ymin plot_ymax]);
xlabel('Side Scattering');
ylabel('GL-1');
title('Chlorophyll Signal for Small Eukaryotes')
lh = legend( 'Small eukaryotes');

%% figure 4 temperature plot

load([outpath 'compiled_stats'], 'fcsfile*', 'SynConc', 'EukConc','SynCount','SynYcv', 'EukCount','EukYcv');
%return
figure
plot(time,SynConc*1000, 'b.','MarkerSize',12)
hold on
plot(time,EukConc*1000, 'g.','MarkerSize',12)
%xlim([0 2563])
ylabel('Cell concentration (ml^{-1})','FontSize',20)
xlabel('2-minute sample resolution, 31-Jan to 5-Feb 2018','FontSize',20)
%xlabel('2-minute sample resolution, 16-Apr to 29-Apr 2018')
lh = legend('\itSynechococcus', 'Small eukaryotes', 'location', 'northwest');
%title('onshore                      \leftarrow                     offshore            \rightarrow                             onshore')
set(lh, 'fontsize', 20)



%%
figure
plot(fcsmatch.mdate_start,fcsmatch.lat)
datetick