%%
basepath = 'E:\Attune_Data\Size_calibration_July2018\ExportedFCS\';
filelist = dir([basepath 'cruiseMethod*']); 

char(fullfile(basepath,filelist(1).name))

%% Gymgal Data

[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_GYMGAL_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_GYMGAL_(1).fcs');

disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^4;
xmax= 10^6;
ymin= 10^4;
ymax= 10^6;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Gymgal: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itGymgal');
xlabel('SSC-A');
ylabel('BL3-H');

gym = log10(ssca(find(in_gate==1)));
gym2 = ssca(find(in_gate==1))
%scattering is converted to volume with the this function from the
%calibration curve
gymVol = (1.3.*gym)-2.9;

lin_Vol = 10.^(gymVol)

%this converts volume(that has been translated to linear space to diamter
GYMGAL = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.GYMGAL.attune = GYMGAL

%Plotting Histogram of the scattering signal
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
% [gymN,gymxout] = hist(gymD,gate_count);
% bar(gymxout, gymN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
% xlabel('');
% ylabel('Count');
% title('Histogram of Side Scattering for Gymgal');

%% Htriq
%reading in the fcs file
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_Heterocapsa_HTriq_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_Heterocapsa_HTriq_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^4;
xmax= 10^6;
ymin= 10^5;
ymax= 10^6;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Htriq: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itHtriq');
xlabel('SSC-A');
ylabel('BL3-H');

htriq = log10(ssca(find(in_gate==1)));
htriq1 = ssca(find(in_gate==1))
%scattering is converted to volume with the this function from the
%calibration curve
htriqVol = (1.3.*htriq)-2.9;
lin_Vol = 10.^(htriqVol)
HTriq = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.HTriq.attune = HTriq


%Plotting Histogram of the scattering signal
% figure
% histogram(alexVol,2,'BinWidth',0.2 ,'EdgeColor','b','FaceColor','b')
% %bar(xout, n, 'barwidth', 3, 'basevalue', 0,'FaceColor',[1 0 0]);
% xlabel('Log of volume (um^3)');
% ylabel('Count');
% title('Histogram of Side Scattering for Alexandrium');

figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[alexN,alexout] = hist(alexVol,gate_count);
bar(alexout, alexN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Alexandrium');


%% Alexandrium Data

% Alexandrium 3
%reading in the fcs file
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_Alexandrium_minutum_AMBOPO14_(3).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_Alexandrium_minutum_AMBOPO14_(3).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^4.5;
xmax= 10^5.8;
ymin= 10^5;
ymax= 10^6;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Alexandrium: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itAlexandrium');
xlabel('SSC-A');
ylabel('BL3-H');

alex = log10(ssca(find(in_gate==1)));
alex2 = ssca(find(in_gate==1))
mode(alex2)

%scattering is converted to volume with the this function from the
%calibration curve
alexVol = (1.3.*alex)-2.9;
lin_Vol = 10.^(alexVol)
AMBOPO = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.AMBOPO.attune = AMBOPO


%Plotting Histogram of the scattering signal
% figure
% histogram(alexVol,2,'BinWidth',0.2 ,'EdgeColor','b','FaceColor','b')
% %bar(xout, n, 'barwidth', 3, 'basevalue', 0,'FaceColor',[1 0 0]);
% xlabel('Log of volume (um^3)');
% ylabel('Count');
% title('Histogram of Side Scattering for Alexandrium');

figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[alexN,alexout] = hist(alexVol,gate_count);
bar(alexout, alexN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Alexandrium');


%% Crypto

%reading in the fcs file
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_crypto_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_crypto_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^3.5;
xmax= 10^6;
ymin= 10^5;
ymax= 10^6;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Crypto: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r')

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itCrypto');
xlabel('SSC-A')
ylabel('BL3-H')

crypto = log10(ssca(find(in_gate==1)));
crypto1 = ssca(find(in_gate==1))

%scattering is converted to volume with the this function from the
%calibration curve
cryptoVol = (1.3.*crypto)-2.9;

lin_Vol = 10.^(cryptoVol)
CRYPTO_SPMC = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.CRYPTO_SPMC.attune = CRYPTO_SPMC


figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[crypN,crypout] = hist(cryptoVol,gate_count);
bar(crypout, crypN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Crypto');

%% H. Akashiwo
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_HAkashiwo_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_HAkashiwo_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;
%the bounds for the rectangular gate
xmin= 10^4;
xmax= 10^6;
ymin= 10^4;
ymax= 10^6;
% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['HAkashiwo ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itHAkashiwo');
xlabel('SSC-A');
ylabel('BL3-H');

hak = log10(ssca(find(in_gate==1)));
hak1 = ssca(find(in_gate==1))
%scattering is converted to volume with the this function from the
%calibration curve
hakVol = (1.3.*hak)-2.9;

lin_Vol = 10.^(hakVol)
H_Akashiwo_CCMP3374 = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.H_Akashiwo_CCMP3374.attune = H_Akashiwo_CCMP3374

figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[hakN,hakout] = hist(hakVol,gate_count);
bar(hakout, hakN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for HK');


%% 'cruiseMethod_HK_(1).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_HK_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_HK_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^4;
xmax= 10^5.5;
ymin= 10^4.8;
ymax= 10^6;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['HK: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itHK');
xlabel('SSC-A');
ylabel('BL3-H');

hk = log10(ssca(find(in_gate==1)));
hk1 = ssca(find(in_gate==1))
mode(hk1)
%scattering is converted to volume with the this function from the
%calibration curve
hkVol = (1.3.*hk)-2.9;

lin_Vol = 10.^(hkVol)
HK = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.HK.attune = HK


figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[hkN,hkout] = hist(hkVol,gate_count);
bar(hkout, hkN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for HK');

%% 'cruiseMethod_Heterocapsa_HET_(1).fcs'

[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_Heterocapsa_HET_(1).fcs');
% [~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_Heterocapsa_HET_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^4.5;
xmax= 10^6;
ymin= 10^4.5;
ymax= 10^6;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Heterocapsa: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itHeterocapsa');
xlabel('SSC-A');
ylabel('BL3-H');

het1 = log10(ssca(find(in_gate==1)));
het2=ssca(find(in_gate==1))
mode(het2)
%scattering is converted to volume with the this function from the
%calibration curve
het1Vol = (1.3.*het1)-2.9;


lin_Vol = 10.^(het1Vol)
HET = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.HET.attune = HET


figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[het1N,het1out] = hist(het1Vol,gate_count);
bar(het1out, het1N, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Heterocapsa');

%% 'cruiseMethod_Igalbana_(1).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_Igalbana_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_Igalbana_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^2.8;
xmax= 10^4.5;
ymin= 10^4.5;
ymax= 10^6;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Ilgalbana: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itIlgalbana');
xlabel('SSC-A');
ylabel('BL3-H');

Ilgal= log10(ssca(find(in_gate==1)));
Ilgal1 = ssca(find(in_gate==1))
mode(Ilgal1)

%scattering is converted to volume with the this function from the
%calibration curve
galVol = (1.3.*Ilgal)-2.9;

lin_Vol = 10.^(galVol);
I_galbana = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.I_galbana.attune = I_galbana

figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[galN,galout] = hist(galVol,gate_count);
bar(galout, galN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Il Galbana');


%% 'cruiseMethod_PRYMN_(2).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_PRYMN_(2).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_PRYMN_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^3.5;
xmax= 10^5.2;
ymin= 10^4.2;
ymax= 10^5.8;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Pyrmn: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itPyrmn');
xlabel('SSC-A');
ylabel('BL3-H');

prymn = log10(ssca(find(in_gate==1)));
prymn1 = ssca(find(in_gate==1))
mode(prymn1)

%scattering is converted to volume with the this function from the
%calibration curve
prymnVol = (1.3.*prymn)-2.9;

lin_Vol = 10.^(prymnVol)
PRYMN_PAV = (((lin_Vol/pi).*(3/4)).^(1/3)).*2
sizestr.PRYMN_PAV.attune = PRYMN_PAV

figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[pN,pout] = hist(prymnVol,gate_count);
bar(pout, pN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Il Pyrmn');


%% 'cruiseMethod_Pavlova_lutheri_(1).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_Pavlova_lutheri_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_Pavlova_lutheri_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^3.2;
xmax= 10^5.2;
ymin= 10^3.8;
ymax= 10^5.5;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Pavlova: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itPavlova');
xlabel('SSC-A');
ylabel('BL3-H');

pav = log10(ssca(find(in_gate==1)));
pav1 = ssca(find(in_gate==1))
mode(pav1)

%scattering is converted to volume with the this function from the
%calibration curve
pavVol = (1.3.*pav)-2.9;

lin_Vol = 10.^(pavVol)
Pavlova = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.Pavlova.attune = Pavlova


figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[pN,pout] = hist(pavVol,gate_count);
bar(pout, pN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Pavlova');


%% cruiseMethod_Peridinium_PTPP02_(1).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_Peridinium_PTPP02_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_Peridinium_PTPP02_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^4.7;
xmax= 10^6;
ymin= 10^5;
ymax= 10^6.2;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Peridinium: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itPeridinium');
xlabel('SSC-A');
ylabel('BL3-H');

per = log10(ssca(find(in_gate==1)));
per1 = ssca(find(in_gate==1))
mode(per1)

%scattering is converted to volume with the this function from the
%calibration curve
perVol = (1.3.*per)-2.9;

lin_Vol = 10.^(perVol)
PPTP02 = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.PPTP02.attune = PPTP02


figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[pN,pout] = hist(perVol,gate_count);
bar(pout, pN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Peridinium');


%% 'cruiseMethod_Scrippsiella_SA2_(1).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_Scrippsiella_SA2_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_Scrippsiella_SA2_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10^2;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^4.7;
xmax= 10^6;
ymin= 10^5;
ymax= 10^6.2;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
xlim([min max]);
ylim([min max]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Scrippsiella: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itScrippsiella');
xlabel('SSC-A');
ylabel('BL3-H');

scripp = log10(ssca(find(in_gate==1)));
scripp1 = ssca(find(in_gate==1))
mode(scripp1)

%scattering is converted to volume with the this function from the
%calibration curve
scripVol = (1.3.*scripp)-2.9;

lin_Vol = 10.^(scripVol)
SA2 = (((lin_Vol/pi).*(3/4)).^(1/3)).*2
sizestr.SA2.attune = SA2


figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[sN,sout] = hist(SA2,gate_count);
bar(sout, sN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel(' Diameter (um)');
ylabel('Count');
title('Histogram of Side Scattering for Scrippsiella');


%% 'cruiseMethod_Syn_10.1H_(1).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_Syn_10.1H_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_Syn_10.1H_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^1.2;
xmax= 10^3.2;
ymin= 10^3.5;
ymax= 10^5;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,19);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
% xlim([axismin aximax]);
% ylim([axismin aximax]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Syn: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itSyn');
xlabel('SSC-A');
ylabel('GL3-H');

syn = log10(ssca(find(in_gate==1)));
syn1 = ssca(find(in_gate==1))
mode(syn1)

%scattering is converted to volume with the this function from the
%calibration curve
synVol = (1.3.*syn)-2.9;
lin_Vol = 10.^(synVol)
Syn = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.Syn.attune = Syn

figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[sN,sout] = hist(synVol,gate_count);
bar(sout, sN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Syn');

%% 'cruiseMethod_dun_(1).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_dun_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_dun_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^3.7;
xmax= 10^5.2;
ymin= 10^5;
ymax= 10^6.2;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
% xlim([axismin aximax]);
% ylim([axismin aximax]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Duniella: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itDuniella');
xlabel('SSC-A');
ylabel('BL3-H');

dun = log10(ssca(find(in_gate==1)));

%scattering is converted to volume with the this function from the
%calibration curve
dunVol = (1.3.*dun)-2.9;

lin_Vol = 10.^(dunVol)

Dun7_12 = (((lin_Vol/pi).*(3/4)).^(1/3)).*2


sizestr.Dun7_12.attune = Dun7_12

figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[sN,sout] = hist(dunVol,gate_count);
bar(sout, sN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Duniella');

%% Micromonas 'cruiseMethod_micro_(1).fcs'
[~,fcshdr,fcsdatscaled] =fca_readfcs('\\sosiknas1\Lab_data\Attune\Size_calibration_July2018\ExportedFCS\cruiseMethod_micro_(1).fcs');
%[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_micro_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^1.8;
xmax= 10^3.8;
ymin= 10^3.5;
ymax= 10^4.5;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
% xlim([axismin aximax]);
% ylim([axismin aximax]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Micromonas: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itMicromonas');
xlabel('SSC-A');
ylabel('BL3-H');

micro = log10(ssca(find(in_gate==1)));

%scattering is converted to volume with the this function from the
%calibration curve
microVol = (1.3.*micro)-2.9;
lin_Vol = 10.^(microVol)
Micromonas_cold = (((lin_Vol/pi).*(3/4)).^(1/3)).*2

sizestr.Micromonas_cold.attune = Micromonas_cold


figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[sN,sout] = hist(microVol,gate_count);
bar(sout, sN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Duniella');

%% 'cruiseMethod_nano_(1).fcs'

[~,fcshdr,fcsdatscaled] =fca_readfcs('E:\Attune_Data\Size_calibration_July2018\ExportedFCS\cruiseMethod_nano_(1).fcs');
disp(fcshdr.filename)
%the bounds for the axis for the scatter plot of the signal
axismin = 10;
aximax =  10^8;

%the bounds for the rectangular gate
xmin= 10^1.8;
xmax= 10^3.8;
ymin= 10^3.5;
ymax= 10^4.5;

% rectangular gate
x_rect = [xmin xmin xmax xmax xmin];
y_rect = [ymin ymax ymax ymin ymin];

ssca = fcsdatscaled(:,3);
bl3a = fcsdatscaled(:,15);
%plotting the scattering signals from different lasers
figure
loglog(ssca,bl3a,'k.','HandleVisibility','off')
% xlim([axismin aximax]);
% ylim([axismin aximax]);
hold on
in_gate = inpolygon(ssca,bl3a,x_rect,y_rect);
gate_count = length(find(in_gate==1));
txt1 = ['Nanochloris: ',num2str(gate_count)];
text(xmin+100,10^4.15,txt1,'Color','r');

%plotting the points in the gate in red
hold on
loglog(ssca(in_gate),bl3a(in_gate),'r.')
%plotting the gate itself in red 
hold on
loglog(x_rect,y_rect,'LineWidth',2,'Color','r','LineStyle','--')
lh = legend('\itNanochloris');
xlabel('SSC-A');
ylabel('BL3-H');

nano = ssca(find(in_gate==1));

%scattering is converted to volume with the this function from the
%calibration curve
nanoVol = (1.3.*nano)-3;

lin_Vol = 10.^(nanoVol)
Nano = (((lin_Vol/pi).*(3/4)).^(1/3)).*2
sizestr.Nano.attune = NaN

figure
%histogram(gymVol,'BinWidth',0.2 ,'EdgeColor','r','FaceColor','r')
[sN,sout] = hist(nanoVol,gate_count);
bar(sout, sN, 'barwidth', 2, 'basevalue', 0,'FaceColor',[1 0 0]);
xlabel('Log of volume (um^3)');
ylabel('Count');
title('Histogram of Side Scattering for Nanochloris');


%%

scattering = ssca(find(in_gate == 1));
%calibration equation
volume = (1.3.*scattering) - 2.9
lin_vol = 10.^(volume)
diameter = (((lin_Vol/pi).*(3/4)).^(1/3)).*2
sizestr.Nano.attune = diameter




%% Super Histogram
map = brewermap(20,'Set1'); 

binwid = 0.05;
binedge = 1:0.1:50
figure
histogram(log(((gymVol/pi)*(3/4)^(1/3))*2) ,'EdgeColor',map(1,:),'FaceColor',map(1,:),'FaceAlpha',0.5,'BinEdges',binedge);
hold on
histogram(log(((alexVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid,'EdgeColor',map(2,:),'FaceColor',map(2,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((hakVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor',map(4,:),'FaceColor',map(4,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((het1Vol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor',map(5,:),'FaceColor',map(5,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((htriqVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor',map(6,:),'FaceColor',map(6,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((galVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor',map(7,:),'FaceColor',map(7,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((prymnVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor',map(8,:),'FaceColor',map(8,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((pavVol/pi)*(3/4)^(1/3))*2),'BinWidth', binwid,'EdgeColor',map(9,:),'FaceColor',map(9,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((perVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor',map(10,:),'FaceColor',map(10,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((scripVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid,'EdgeColor',map(11,:),'FaceColor',map(11,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((dunVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor','r','FaceColor',map(12,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((microVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor',map(13,:),'FaceColor',map(13,:),'FaceAlpha',0.5,'BinEdges',binedge);
histogram(log(((cryptoVol/pi)*(3/4)^(1/3))*2),'BinWidth',binwid ,'EdgeColor',map(20,:),'FaceColor',map(20,:),'FaceAlpha',0.5,'BinEdges',binedge);
title = {'Gymgal','Alex','H.Akashiwo','Heterocapsa','HTriquetra','IlGalbana','Prymnesium','Pavlova','Peridinium', 'Scripsiella','Duniella','Micromonas','Crypto'}
legend(title')
xlabel('log of Diameter (um)');
barstrings = num2str(n');
ylabel('Count');
