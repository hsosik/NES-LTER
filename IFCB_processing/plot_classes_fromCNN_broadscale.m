yr = 2013:2019;
for count = 1:length(yr)
    %s(count) = load(['\\sosiknas1\IFCB_products\NESLTER_broadscale\summary\summary_biovol_allHDF_min20_' num2str(yr(count))]);
    s(count) = load(['c:\work\IFCB_products\NESLTER_broadscale\summary\summary_biovol_allHDF_min20_' num2str(yr(count))]);
end

IFCBsum = table;
slist = {'filelist' 'classcount' 'ml_analyzed'};
for count = 1:length(slist)
    t = slist{count}; IFCBsum.(t) = cat(1,s.(t)); 
end
class2use = s(1).class2use;
%meta = readtable('\\sosiknas1\IFCB_products\NESLTER_broadscale\exported_metadata.csv');
meta = readtable('c:\work\IFCB_products\NESLTER_broadscale\exported_metadata.csv');
[~,a,b] = intersect(IFCBsum.filelist, meta.pid); 
IFCBsum.lat(a) = meta.latitude(b);
IFCBsum.lon(a) = meta.longitude(b);
mdate = IFCB_file2date(IFCBsum.filelist);
[~,a,b] = intersect(IFCBsum.filelist, meta.pid);
IFCBsum.ml_analyzed(a) = meta.ml_analyzed(b);

clear t count yr 

%%
cc = strmatch('Hemiaulus', class2use);
%cc = strmatch('Guinardia_delicatula', class2use, 'exact');
cc = strmatch('Guinardia_delicatula', class2use);
ind = find(mdate>datenum('1-15-2017') & mdate<datenum('3-1-2017'));
figure
[ASIT(1), ASIT(2)] = m_ll2xy(-70.5667,41.325);
%m_gshhs_h('save','gumby3')
m_proj('Mollweide','long',[-76.5 -65],'lat',[36 44]);
m_usercoast('gumby3','patch',[.8 .8 .8],'edgecolor','k');
m_grid('box','fancy','tickdir','out','fontsize', 10, 'fontname','Times New Roman');
m_elev('contour',[-200 -200],'edgecolor','b');
[X,Y] = m_ll2xy(IFCBsum.lon(ind), IFCBsum.lat(ind)); %convert positions of IFCB manual points to map coordinates
Z = log10(sum(IFCBsum.classcount(ind,cc),2)./IFCBsum.ml_analyzed(ind));
ii = find(isinf(Z));
plot(X(ii),Y(ii), 'ok', 'markersize',4)
hold on
Z(ii) = NaN;
scatter(X,Y,20,Z,'filled') %microg per liter
cbh = colorbar;
caxis([log10(.1) log10(100)])
set(cbh, 'ytick', log10([.1 1 10 100]))
set(cbh, 'yticklabel', 10.^(get(cbh, 'ytick')))
set(cbh, 'position', [.88 .12 .03 .6])
title(cbh, {'chains' ; 'ml^{-1}'}, 'fontsize', 14)
plot(ASIT(1), ASIT(2), 'rpentagram', 'markersize', 20, 'markerfacecolor', 'r')
%text(ASIT(1), ASIT(2), '  MVCO', 'fontsize', 14, 'color', 'r')
title('HB1701; 11-23 Feb 2017')

return

latlim = [35 45];
lonlim = [-76 -65];
figure
ax = usamap(latlim, lonlim)
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
scatterm(IFCBsum.lat(ind), IFCBsum.lon(ind),10,log10(sum(IFCBsum.classcount(ind,cc),2)./IFCBsum.ml_analyzed(ind)), 'filled') %microg per liter
ASIT = [41.325 -70.5667];
cbh = colorbar;
caxis([log10(.1) log10(100)])
set(cbh, 'ytick', log10([.5 5 50]))
set(cbh, 'yticklabel', 10.^(get(cbh, 'ytick')))
plotm(ASIT(1), ASIT(2), 'rpentagram', 'markersize', 20, 'markerfacecolor', 'r')
textm(ASIT(1), ASIT(2), '  MVCO', 'fontsize', 14, 'color', 'r')
title('HB1701; 11-23 Feb 2017')

   