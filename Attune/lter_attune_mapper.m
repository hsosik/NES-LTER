function  lter_attune_mapper(basepath,varagin,color,plotName,y_label)
%basepath: path to cruise folder with summary inside
%varagin is a variable to be plotted against lat and lon (e.g. Attune.fcsmatch.temperature or Attune.Count.EukTotal)
%color options (for colorbar)  :'thermal','haline','solar','ice','gray','deep','dense','algae','oxy','matter','turbid'
%(see https://matplotlib.org/cmocean/ for the full range of premade color maps)
%plotName should be a string e.g. 'Sea Surface Temperature for EN608 Cruise'
%ylabel (string) e.g. 'Temperature (Degrees C)'

load([basepath '\Summary\Attune']);
addpath 'C:\Users\mps48\Desktop\m_map';
addpath 'C:\Program Files\MATLAB\R2018a\toolbox\add_ons';

figure
m_proj('UTM','longitude',[-72.5 -69.5],'latitude',[39.5 42.5]);
m_gshhs_i('patch', [0 0.5 0]);
m_grid('backgroundcolor',[0.5843 0.8157 .98],'box','fancy','tickdir','in');

hold on
c = cmocean(color);
cmap = colormap(c);

m_scatter(Attune.fcsmatch.lon,Attune.fcsmatch.lat,20,varagin,'filled');
a = colorbar('eastoutside');

title(plotName);
y = ylabel(a ,y_label,'Rotation', -90);
set(y, 'position', get(y,'position')-[-1,0,0]);

m_proj get;

end