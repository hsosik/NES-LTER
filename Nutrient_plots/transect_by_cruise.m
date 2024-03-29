% Heidi M. Sosik, NES-LTER transect plotting of discrete nutrient concentrations
% June 2021

% Missing nuts (so far) : 'AR16' 'AR34B' 'AR39A' 'EN661' 'AR44' 'AR48A' 'AR48B' 

%cruises = {'AR22' 'AR24A' 'AR24B' 'AR24C' 'EN608' 'AR28A' 'AR28B' 'EN617'...
%    'AR31A' 'AR31B' 'AR31C' 'AR32' 'EN627' 'AR34A' 'AR34B' 'AR38' 'EN644'...
%    'EN655' 'EN657' 'AR39B'};  

%cruises = {'EN657'}; %use this option to run one cruise OR above case for saving and plottting multiple cruises
%cruises = {'en627' 'ar34b' 'en644' 'ar39b'};

%cruises = readtable("C:\Users\heidi\Downloads\cruise_metadata.csv");
opt = weboptions('Timeout', 60);
cruises = webread('https://nes-lter-data.whoi.edu/api/cruises/metadata.csv', opt);
cruises = cruises.cruise;
%%
if 1 %1 to read from the APIs, 0 to load the stored (multi-cruise) file
    nut = table;
    for count1 = 1:length(cruises)
        disp(cruises{count1})
        opt = weboptions('Timeout', 120);
        try
            %opt = detectImportOptions(['https://nes-lter-data.whoi.edu/api/nut/' cruises{count1} '.csv']);
            %opt.VariableTypes(strcmp(opt.VariableNames, 'alternate_sample_id')) = {'char'};
            %opt.Timeout = 120;
            n = webread(['https://nes-lter-data.whoi.edu/api/nut/' cruises{count1} '.csv'], opt);
            n.alternate_sample_id = []; %move this column since the type doesn't match between all cruises
            %n.nearest_station = [];
            %n.distance_km = [];
            v = n.Properties.VariableNames;
    %        if ~ismember('alternate_sample_id',v)
    %            n.alternate_sample_id(:) = {''};
    %        end
            if ~ismember('nearest_station',v)
                n.nearest_station(:) = {''};
            end
            if ~ismember('distance_km',v)
                n.distance_km(:) = NaN;
            end
            if iscell(n.distance_km)
                n.distance_km_temp = n.distance_km;
                n = removevars(n,"distance_km");
                n.distance_km(strcmp(n.distance_km_temp,'NA')) = NaN;
                if sum(~strcmp(n.distance_km_temp,'NA'))
                    n.distance_km(~strcmp(n.distance_km_temp,'NA')) = str2num(char(n.distance_km_temp(~strcmp(n.distance_km_temp,'NA'))));
                end
                n = removevars(n,"distance_km_temp");
            end
            nut = [nut; n];
        catch
            disp(cruises{count1})
            disp('no data found')
        end
    end
    %nut.mdate = datenum(nut.date, 'yyyy-mm-dd hh:MM:ss+00:00');
    nut.datetime(~strcmp(nut.date,'NA')) = datetime(nut.date(~strcmp(nut.date,'NA')),'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');

    save('c:\work\LTER\nut_all', 'nut', 'cruises')
else
    load('c:\work\LTER\nut_all')
end

%%
ilat = 39.5:.05:41.5;
ilon = ones(size(ilat)).*-70.8855;
idpth = 0:1:200;
bdata = load('ngdc2.mat'); %from Gordon
ibdpth = griddata(bdata.lon',bdata.lat,bdata.h,ilon(:),ilat(:));

%cstr = 'EN617';
tstr = {'NO_3+NO_2' 'NH_4' 'PO_4' 'SIO_4'};
cmax = [10 5 2 15];
latfactor = 100;
    
for count = 1:length(cruises)
    cstr = cruises{count};
    nind = find(ismember(lower(nut.cruise), lower(cstr)));
if ~isempty(nind)
    inut(:,:,1) = griddata(nut.latitude(nind)*latfactor,nut.depth(nind),nut.nitrate_nitrite(nind),ilat*latfactor,idpth');
    inut(:,:,2) = griddata(nut.latitude(nind)*latfactor,nut.depth(nind),nut.ammonium(nind),ilat*latfactor,idpth');
    inut(:,:,3) = griddata(nut.latitude(nind)*latfactor,nut.depth(nind),nut.phosphate(nind),ilat*latfactor,idpth');
    inut(:,:,4) = griddata(nut.latitude(nind)*latfactor,nut.depth(nind),nut.silicate(nind),ilat*latfactor,idpth');

    for ii = 1:length(ilat)
        inut(idpth>ibdpth(ii),ii,:) = NaN; 
    end

    figure
    for cc = 1:4
        subplot(4,1,cc)
        pcolor(ilat, idpth', squeeze(inut(:,:,cc)))
        hold on
        plot(nut.latitude(nind), nut.depth(nind), '+')
        plot(ilat,ibdpth,'k','linewidth',2);
        shading interp, set(gca, 'ydir', 'rev', 'xdir', 'rev')
        %if count == 1, title(cstr), end
        if cc == 1
            title([cstr ' ' datestr(min(floor(nut.mdate(nind))),'dd mmm') '-' datestr(max(floor(nut.mdate(nind))),'dd mmm yyyy')])
        end
        text(41.4,175,tstr(cc)) 
        colorbar
        caxis([0 cmax(cc)])
    end
    set(gcf, 'position', [488 41.8 560 740.8])
    print(['c:\work\lter\nutrient_sections\' cruises{count}], '-dpng')
end
end


return
%%
figure
set(gcf', 'position', [0 350 1200 400])
%tstr = [{'Winter 2019'; 'NO_3+NO_2'} {'Spring 2019'; 'NO_3+NO_2'} {'Summer 2019'; 'NO_3+NO_2'} {'Fall 2019'; 'NO_3+NO_2'}];
tstr = [{'Winter 2019'} {'Spring 2019'} {'Summer 2019'} {'Fall 2019'}];
tiledlayout(2,2,'TileSpacing', 'compact')
for count = 1:length(cruises)
    cstr = cruises{count};
    nind = find(ismember(lower(nut.cruise), lower(cstr)));

    inut(:,:,1) = griddata(nut.latitude(nind)*latfactor,nut.depth(nind),nut.nitrate_nitrite(nind),ilat*latfactor,idpth');

    for ii = 1:length(ilat)
        inut(idpth>ibdpth(ii),ii,:) = NaN; 
    end

        nexttile
        pcolor(ilat, idpth', squeeze(inut))
        hold on
        plot(nut.latitude(nind), nut.depth(nind), '+')
        plot(ilat,ibdpth,'k','linewidth',2);
        shading interp, set(gca, 'ydir', 'rev', 'xdir', 'rev')
        text(41.4,150,tstr(:,count), 'fontsize', 16) 
        %colorbar
        caxis([0 12])
end
cbh = colorbar;
title(cbh, 'NO_3+NO_2')