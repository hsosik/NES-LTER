% Heidi M. Sosik, NES-LTER transect plotting of discrete POC, PON, and C:N
% Sept 2022
p = 'C:\work\LTER\POC\';
load([p 'NESLTER_CHN_table'])

cruises = unique(CHNtable.Cruise);
%%
ilat = 39.5:.05:41.5;
ilon = ones(size(ilat)).*-70.8855;
idpth = 0:1:200;
bdata = load('ngdc2.mat'); %from Gordon
ibdpth = griddata(bdata.lon',bdata.lat,bdata.h,ilon(:),ilat(:));

tstr = {'NO_3+NO_2' 'NH_4' 'PO_4' 'SIO_4'};
tstr = {'POC (\mug l^{-1})' 'PON (\mug l^{-1})' 'C:N (molar)' };
cmax = [600 100 12];
latfactor = 100;

%%
for count = 1:length(cruises)
    cstr = cruises{count};
    nind = find(ismember(lower(CHNtable.Cruise), lower(cstr)));
if ~isempty(nind)
    inut(:,:,1) = griddata(CHNtable.latitude(nind)*latfactor,CHNtable.depth(nind),CHNtable.POC_ugperL(nind),ilat*latfactor,idpth');
    inut(:,:,2) = griddata(CHNtable.latitude(nind)*latfactor,CHNtable.depth(nind),CHNtable.PON_ugperL(nind),ilat*latfactor,idpth');
    inut(:,:,3) = griddata(CHNtable.latitude(nind)*latfactor,CHNtable.depth(nind),CHNtable.C_to_N_molar_ratio(nind),ilat*latfactor,idpth');

    for ii = 1:length(ilat)
        iCHNtable(idpth>ibdpth(ii),ii,:) = NaN; 
    end

    figure
    for cc = 1:3
        subplot(3,1,cc)
        pcolor(ilat, idpth', squeeze(inut(:,:,cc)))
        hold on
        plot(CHNtable.latitude(nind), CHNtable.depth(nind), '+')
        plot(ilat,ibdpth,'k','linewidth',2);
        shading interp, set(gca, 'ydir', 'rev', 'xdir', 'rev')
        %if count == 1, title(cstr), end
        if cc == 1
            title([cstr ' ' datestr(min(dateshift(CHNtable.datetime(nind),'start', 'day')),'dd mmm') '-' datestr(max(dateshift(CHNtable.datetime(nind),'start', 'day')),'dd mmm yyyy')])
        end
        text(41.4,175,tstr(cc)) 
        colorbar
        caxis([0 cmax(cc)])
    end
    set(gcf, 'position', [488 41.8 560 740.8])
    print([p 'plots' filesep datestr(min(dateshift(CHNtable.datetime(nind),'start', 'day')),'yyyymmdd') '_' cruises{count}], '-dpng')
end
    pause
end

