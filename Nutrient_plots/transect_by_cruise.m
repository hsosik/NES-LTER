% Heidi M. Sosik, NES-LTER transect plotting of discrete nutrient concentrations
% June 2021

% Missing nuts (so far) : 'AR16' 'AR34B' 'AR39A' 'EN661' 'AR44' 'AR48A' 'AR48B' 

%cruises = {'AR22' 'AR24A' 'AR24B' 'AR24C' 'EN608' 'AR28A' 'AR28B' 'EN617'...
%    'AR31A' 'AR31B' 'AR31C' 'AR32' 'EN627' 'AR34A' 'AR34B' 'AR38' 'EN644'...
%    'EN655' 'EN657' 'AR39B'};  

cruises = {'EN657'}; %use this option to run one cruise OR above case for saving and plottting multiple cruises

if 1 %1 to read from the APIs, 0 to load the stored (multi-cruise) file
    nut = table;
    for count1 = 1:length(cruises)
        disp(cruises{count1})
        n = webread(['https://nes-lter-data.whoi.edu/api/nut/' cruises{count1} '.csv']);
        n.alternate_sample_id = []; %move this column since the type doesn't match between all cruises
        nut = [nut; n];
    end
    nut.mdate = datenum(nut.date, 'yyyy-mm-dd hh:MM:ss+00:00');
    save('c:\work\LTER\nut_all', 'nut', 'cruises')
else
    load('c:\work\LTER\nut_all')
end

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
    nind = find(ismember(nut.cruise, cstr));

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
