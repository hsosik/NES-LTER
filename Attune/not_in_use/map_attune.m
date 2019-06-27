function [] = mapper()
%
%
% addpath 
basepath = 'C:\Users\mps48\Desktop\Summary\Attune'
load(basepath)
addpath 'C:\Users\mps48\Desktop\m_map'
addpath 'C:\Program Files\MATLAB\R2018a\toolbox\add_ons'
figure
m_proj('UTM','longitude',[-72.5 -69.5],'latitude',[39.5 42.5])
m_gshhs_i('patch', [0 0.5 0]);
% m_gshhs_i('speckle','color','k');
m_grid('backgroundcolor',[0.5843 0.8157 .98],'box','fancy','tickdir','in');
% m_track(fcsmatch.lon,fcsmatch.lat,'linewdith',1.5);
% m_pcolor(long,lati,c)
hold on
c = cmocean(color);
cmap = colormap(c);

m_scatter(Attune.fcsmatch.lon,Attune.fcsmatch.lat,20,Attune.Count.EukTotal,'filled');
colorbar('eastoutside');
title('EukTotal');
m_proj get;

