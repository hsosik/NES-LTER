% basepath = 'E:\Attune_Data\EN608\';
basepath = '\\sosiknas1\Lab_data\Attune\EN608'
% basepath = '\\sosiknas1\Backup\LTER\20180404_AR28\'

fpath = [basepath '\ExportedFCS\'];
outpath = [basepath '\Summary\'];

% Extracting files out of the directory sorts NES out from SFD
%first it will populate with NES titled files but if empty will go for SFD
%file string
filelist = dir([fpath 'NES*']);
filelist = {filelist.name}';

[Attune.FCSfileinfo] = FCS_DateTimeList(fpath)

Attune.Count.lesstwo = [];
Attune.Count.twoten = [];
Attune.Count.tentwen =[];
Attune.Count.twen = [];
Attune.Count.Syn = [];
Attune.vol_analyzed = [];
Attune.Biovol.lesstwo = [];
Attune.Biovol.twoten = [];
Attune.Biovol.tentwen =[];
Attune.Biovol.twen = [];

for count = 1:length(sortedlist)
    disp(sortedlist{count});
    filename = [fpath sortedlist{count}];
    [~,fcshdr,fcsdatscaled] =fca_readfcs(filename);
    
    %
    synSignal = fcsdatscaled(:,11);
    %
    fscSignal =fcsdatscaled(:,19);
   
   %Defining the rectangular gate for the Synechecoccus Signal
    SynXmin= 200;
    SynXmax= 10^4;
    SynYmin= 10^3;
    SynYmax= 10^5;
    x_rect = [SynXmin SynXmin SynXmax SynXmax SynXmin];
    y_rect = [SynYmin SynYmax SynYmax SynYmin SynYmin];
    
   in_syn = inpolygon(synSignal,fscSignal,x_rect,y_rect);
   SynSsc = log10(synSignal(find(in_syn ==1)));
   SynY = log10(fscSignal(find(in_syn ==1)));
   %volume from scattering conversion
   SynVolume = 1.3.*SynSsc - 2.9;
   lin_SynVol = 10.^(SynVolume);
   SynDiameter = ((lin_SynVol/pi).*(3/4)).^(1/3);
   
   %S
   ssc_signal = fcsdatscaled(:,3);
   %SSC
   y_signal =fcsdatscaled(:,15);
   
   %defing the polygon gate for the Small Eukaryote Signal
   x_polygon = [25  50 10^4 10^6 10^6 10^5 10^4 25];
   y_polygon = [300 3500 10^6 10^6 10^5 10^4 10^3 300];

   in_euk = inpolygon(ssc_signal,y_signal,x_polygon,y_polygon);
   ssc = log10(ssc_signal(find(in_euk==1)));
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
 
   Attune.Biovol.lesstwo = [Attune.Biovol.lesstwo; sum((4/3).*pi.*((size2./2).^3))];
   Attune.Biovol.twoten = [Attune.Biovol.twoten; sum((4/3).*pi.*(size2_10./2).^3)];
   Attune.Biovol.tentwen =[Attune.Biovol.tentwen; sum((4/3).*pi.*(size10_20./2).^3)];
   Attune.Biovol.twen = [Attune.Biovol.twen; sum((4/3).*pi.*(size20./2).^3)];
   Attune.Biovol.Syn = [Attune.Biovol.Syn; sum((4/3).*pi.*(SynDiameter./2).^3)];
   Attune.Count.lesstwo = [Attune.Count.lesstwo ; length(size2)];
   Attune.Count.twoten = [Attune.Count.twoten ; length(size2_10)];
   Attune.Count.tentwen =[Attune.Count.tentwen ; length(size10_20)];
   Attune.Count.twen = [Attune.Count.twen ; length(size20)];
   Attune.Count.SynTotal = [Attune.Count.Syn; length(SynSsc)];
   Attune.Count.EukTotal = [Attune.Count.total; length(diameter)]
   Attune.Count.vol_analyzed = [Attune.Count.vol_analyzed; fcshdr.VOL];
 
end

save([basepath '\Summary'],Attune)