
%function [ccmean, ccmode] = ccreadraw(datapath, CCfile2analyze)
%  -- reads either old (missing last 8 channels) or new CC files (acquired with Coulter2.bas or .exe)
%  -- corrects for data acquired with Full set on Log Linear Diameter 

% datapath = 'C:\Cabrations2013\Coulter_vs_Acquri_Jun2013\'
% CCfile2analyze = 'TEST50BD.002'
%[CCfile2analyze, datapath] = uigetfile('C:\Cabrations2013\Coulter_vs_Acquri_Jun2013\TEST*');   % for Coulter TEST exp  25-26 Jun 2013
%[CCfile2analyze, datapath] = uigetfile('C:\Cabrations2013\Coulter_sep13\acc*');   % for Coulter exp 12 Sep 2013
[CCfile2analyze, datapath] = uigetfile('C:\cal2018\*');   % for Coulter exp 12 Sep 2013
CCfile2analyze
% if CCfile2analyze == 'n',
%     ccmean = -1; ccmode = -1;
%     return
% end
smooth_points = 10, % 16 or 10

fid = fopen([datapath CCfile2analyze]);

% better ways to read data... (Heidi)
% t = textscan(fid, '%c'); % or t = textscan(fid, '%s') % to read as strings, for header
% s = t{1}(125:end-1); % for header
% u = str2num(reshape(s, 5,256)'); % for data
 
%read header (6 lines)
for j = 1:6
    a{j} = fgetl(fid);
    %a{j} = fgets(fid);
end
%'Get info from Header
dte = [num2str(a{1}(1:2)) '/' num2str(a{1}(3:4)) '/' a{1}(5:6)];
orifice = str2num(a{1}(7:9));
kd = str2num(a{1}(20:25));
size = str2num(a{2}(1:5));
units = str2num(a{2}(7));

%'Find Current and Gain Setting
indic = 0; row = 2; first = 9; last = 11;
indic = zerone(a,row,first,last);
if indic == 9, cg = 'Auto'; end
if indic == 10, cg = 'Manual';end
if indic == 11, cg = 'Off';end

%'Find Value of the Current
current = str2num(a{2}(12:15));

%'Find Value of the Gain
row = 2; first = 19; last = 23;
indic = zerone(a,row,first,last);
gain = 2 ^ (indic - 19);

%'Find Polarity
row = 3; first = 2; last = 4;
indic = zerone(a,row,first,last);
if indic == 2, polarity = '+'; end
if indic == 3, polarity = '-';end
if indic == 4, polarity = 'Alternating';end

%'Find Control i.e. Sampling Mode
row = 3; first = 5; last = 9;
indic = zerone(a,row,first,last);
if indic == 5, control = 'manual'; end
if indic == 6, control = 'time'; end
if indic == 7, control = 'siphon'; end
if indic == 8, control = 'channel'; end
if indic == 9, control = 'total'; end

%'Find Time Specified under sampling time mode
T = a{3}(10:13);

%'Find Channel Count Setting
chcount = a{3}(17:20);

%'Find Total Count Seting
totcount = str2num(a{3}(22:24)); exponent = str2num(a{3}(25));
tot = totcount * 10 ^ (exponent);

%'Find Number of Channels Used
row = 4; first = 1; last = 5;
indic = zerone(a,row,first,last);
channel = 2 ^ (indic + 3);


%'Determine if Data was Corrected for Coincidence
row = 4; first = 10; last = 11;
indic = zerone(a,row,first,last);
if indic == 10, coicor = 'Coincidence Correction ON'; end
if indic == 11, coicor = 'Coincidence Correction OFF'; end

%'Find Volume Siphoned
vol = a{4}(12:15);

%'Determine if Sampling was done in Full, Narrow or Window mode
menu = str2num(a{5}(1));
if menu == 4, menu = 'Full'; end
if menu == 5, menu = 'Narrow'; end
if menu == 6, menu = 'Window'; end

%'Find Accumulation Time
actime = str2num(a{5}(2:7));

%Find Raw Total Count
rawct = str2num(a{5}(8:14));

%'Find Coicidence Corrected Count
coinct = str2num(a{5}(15:21));

%'Find Sample Number
sample = str2num(a{5}(22:24));

%'Find Channel Number Position of Left Cursor in Full Mode
lhf = str2num(a{5}(25:27));

%'Find Channel Number Position of Right Cursor in Full Mode
rhf =str2num(a{5}(28:30));

%'Find the Accumulation Law used in Full Mode
aclawf = str2num(a{5}(34));
if aclawf == 0, aclawf = 'Off';end
if aclawf == 1, aclawf = 'Linear Diameter';end
if aclawf == 2, aclawf = 'Volume';end
if aclawf == 3, aclawf = 'Area';end
if aclawf == 4, aclawf = 'Log Linear Diameter';end

%'Find Accumulation Law for Narrow
aclawn = str2num(a{5}(38));
if aclawn == 0, aclawn = 'Off';end
if aclawn == 1, aclawn = 'Linear Diameter';end
if aclawn == 2, aclawn = 'Volume';end
if aclawn == 3, aclawn = 'Area';end
if aclawn == 4, aclawn = 'Log Linear Diameter';end

%'Find Channel Number for Left Cursor in Narrow
lhnar = str2num(a{5}(39:41));

%'Find Channel Number for Right Hand Cursor in Narrow
rhnar = str2num(a{5}(42:44));



%'Convert strings to values (frequency)
%'Gets five digit numbers from string of numbers
dilute = 1;

t = textscan(fid, '%c');  %this way can read different length files (see comment below about coulter.exe and coulter2.exe)
tchar = t{1}(1:end-1); %remove the final '1', which was dilution factor added at end for some reason.  
t = reshape(tchar,5,length(tchar)/5);
n = str2num(t');
if length(n)<256
    n(length(n)+1:256) = 0;
end
if length(n)>256   
    n = n(1:256);
end

% n=zeros(256,1);
% for t = 1:248  %coulter.exe (coulter.bas) didn't put last 8 channels into files (only 25 lines read instead of 26).  Rob made coulter2.exe (28 Oct 09) to fix this.
% % for t = 1:256  % For use with files made by coulter2.exe
%     n(t) = dilute * str2num(fscanf(fid,'%5s',1));
% end


%  'Mark's CONVERT subroutine
%
if strcmp(menu,'Full')
    if strcmp(aclawf,'Linear Diameter')
        dleft = (kd * .5) / (channel * ((gain * current) ^ (1 / 3)));
        dright = kd * (channel - .5) / (channel * ((gain * current) ^ (1 / 3)));
    end
end
if strcmp(menu,'Narrow')
    if strcmp(aclawn,'Linear Diameter')
        dlnar = kd * (lhf - .5) / (channel * ((gain * current) ^ (1 / 3)));
        drnar = kd * (rhf - .5) / (channel * ((gain * current) ^ (1 / 3)));
    end
end

%'find diameter of left and right cursor positions that the narrow was done on

if strcmp(menu,'Narrow')
    if strcmp(aclawn,'Volume')
        dlnar = kd * (lhf - .5) / (channel * ((gain * current) ^ (1 / 3)));
        dlnar = pi * dlnar ^ 3 / 6;
    end
end
if strcmp(menu,'Narrow')
    if strcmp(aclawn,'Volume')
        drnar = kd * (rhf - .5) / (channel * ((gain * current) ^ (1 / 3)));
        drnar = pi * drnar ^ 3 / 6;
    end
end

%'find volume of left and right cursor positions that the narrow was done on

%'find the diameter all channels
if strcmp(menu,'Narrow')
    dl = dlnar;
    dr = drnar;
end
if strcmp(menu,'Full')
    dl = dleft;
    dr = dright;
end
for j = 1 : channel
    dn(j) = ((j - 1) * (dr - dl) / (channel - 1)) + dl;
end

%         'Coincidence corrected cell density
ncell = 0;
for j = 1 : channel
    ncell = ncell + n(j) / dilute;
end
correction = coinct / ncell;
if strcmp(coicor,'Coincidence Correction OFF')
    correction = rawct / ncell;
end

% 'column 1 = size (um); column 2 = count corrected for deadtime
% '       and coinc corr, if applicable (assumes dilution =1)
%disp([dn' n * correction])
fclose(fid);

% % to correct for when CC was set to Log (by mistake) on Full instead of Linear
diamcorr = zeros(256,1);
if findstr(aclawf,'Log')
    for chan = [1 256]
        diamcorr(chan) = (kd/(current*gain)^(1/3))*(2^(1/50)^(256*((dn(chan)*256*(current*gain)^(1/3)/kd)/256-1)));
    end
    for chan = 2:255
        diamcorr(chan) = ((chan-1)*(diamcorr(256)-diamcorr(1))/(256-1))+diamcorr(1);
    end
    dn = diamcorr';
    disp('corrected diameters: ')
    %disp([dn' n * correction])
end


figure(10)
clf
tempL = []; tempU = [];
screen = get(0, 'ScreenSize');
width = screen(3);
height = screen(4); %turtlescreen = [1 1 1400 1050]; pos on turtle = [190 -7 1330 976];
set(gcf,'position',[width*0 height*0.05 width*.8 height*.85])
%plot(dn,n, 'g.-')
plot(dn,n, 'b.-')
hold on
%nsm = smooth(n,10);  %original Rob's
%nsm = smooth(n,16);  % Alexi for 6.2 beads doublets 
nsm = smooth(n,smooth_points);  % Alexi for 6.2 beads doublets 
plot(dn,nsm, 'r.-')
line([2.139 2.139],[0 max(nsm)])  %mark expected location of 2.14 um beads (internal standard)
xlabel(aclawn)
ylabel('frequency')
%axis([0.6 max(dn) 0 max(nsm(dn>0.6))]) % don't plot all the junk near the baseline
%axis([0.6 max(dn) 0 max(nsm(dn>0.65))]) % plot even less junk for 2010
%axis([0.6 max(dn) 0 max(nsm(dn>0.9))]) % for big cells only

axis([0.6 max(dn) 0 1.1*max(n(dn>0.9))]) % for big cells only

%find mode, mean inside window
title([CCfile2analyze '    Choose LOWER bound for stats......'])
while isempty(tempL)  % make the user pick 2 channels before going on...
    tempL = ginput(1); %tempL=tempL(1); %pick lower and upper channels to bracket peak
end
lower = (max(find(dn<tempL(1))));
line([dn(lower) dn(lower)],[0 max(nsm)],'color','r')
title([CCfile2analyze '    Choose UPPER bound for stats......'])
while isempty(tempU)  % make the user pick 2 channels before going on...
    tempU = ginput(1); %tempU = tempU(1); %pick upper channel to bracket peak
end
upper = (min(find(dn>tempU(1))));
line([dn(upper) dn(upper)],[0 max(nsm)],'color','b')
title(CCfile2analyze)
dnul = dn(lower:upper)';
nul = nsm(lower:upper);
ccmode = dnul(find(nul==(max(nul))));
Smooth_mode = ccmode

nulraw = n(lower:upper);
ccmoderaw = dnul(find(nulraw==(max(nulraw))));
Raw_mode = ccmoderaw

%ccmode = dnul(find(nul==(max(nul(lower:upper)))));
if length(ccmode)>1, ccmode = mean(ccmode);end
line([ccmode  ccmode],[0 max(nsm)],'color','g','linestyle','--','linewidth',3)
ccmean = sum(dn(lower:upper)' .* nsm(lower:upper)) / (sum(nsm(lower:upper)))
if length(ccmean)>1, ccmean = mean(ccmean);end
line([ccmean  ccmean],[0 max(nsm)],'color','k','linestyle','--')
    

% SUB zerone (a$(), row, first, last, indic)
%
%         indic = 0
%         FOR j = first TO last
%         rj$ = MID$(a$(row), j, 1)
%         IF rj$ = "1" THEN indic = j
%         NEXT
%
% END SUB





