%%Example script to plot coastline and transect data for NES-LTER
% uses m_map toolbox:
%   Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", 
%   version 1.4m, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html
%
% Heidi M. Sosik, Woods Hole Oceanographic Institution, May 2020
%%
% load some NCP data
%load C:\work\LTER\Stanley_data\ncplterEn617 ; titlestr = 'EN617 20-25 Jul 2018';
load C:\work\LTER\Stanley_data\ncplterEn644 ; titlestr = 'EN644 19-25 Aug 2019';


figure
m_proj('Mollweide','long',[-72.5 -69.5],'lat',[39.5 43]);
m_usercoast('gumby3','patch',[.8 .8 .8],'edgecolor','k');
m_grid('box','fancy','tickdir','out','fontsize', 10, 'fontname','Times New Roman');
m_tbase('contour',[-200 -200],'edgecolor','k');

[ASIT(1), ASIT(2)] = m_ll2xy(-70.5667,41.325); %MVCO tower lon/lat as map coordinates

[X,Y] = m_ll2xy(ncplter(:,6), ncplter(:,5)); %convert lon/lat positions of sample points to map coordinates
Z = ncplter(:,9); %NCP
scatter(X,Y,20,Z,'filled') %microg per liter
cbh = colorbar;
set(cbh, 'position', [.76 .15 .03 .6])
th = title(cbh, {'NCP' ; 'mmol O_2 m^{-2} d^{-1}'}, 'fontsize', 14);
set(th, 'position', [6.3 -49 0])
plot(ASIT(1), ASIT(2), 'rpentagram', 'markersize', 20, 'markerfacecolor', 'r')
text(ASIT(1), ASIT(2)*1.003, '  MVCO', 'fontsize', 14, 'color', 'r')
title(titlestr)
caxis([0 50])
