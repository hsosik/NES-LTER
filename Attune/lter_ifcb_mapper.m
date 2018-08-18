function  lter_ifcb_mapper(basepath,varagin,color,plotName,ylabel)
% basepath is the filepath to the structure
%
%
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

m_scatter(IFCB.match.lon,IFCB.match.lat,20,varagin,'filled');
colorbar('eastoutside');
title(plotName);
ylabel(ylabel)
y = ylabel(a ,y_label,'Rotation', -90)
set(y, 'position', get(y,'position')-[-1,0,0])
m_proj get;

end
