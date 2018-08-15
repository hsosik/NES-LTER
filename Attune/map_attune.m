%function [map] = map_attune()
%% Map of Sea Surface Temperature
%needs fcsmatch.lat and fcsmatch.lon

lonlim = [-72.5 -69.5];
latlim = [39.5 42.5];

states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
stateName = 'Massachusetts';

ma = states(strcmp(states.Name, stateName));
figure('units','normalized','outerposition',[0 0 1 1])
% ax = usamap('ma',latlim,lonlim);
ax = usamap(latlim,lonlim);
oceanColor = [.5 .7 .9];
setm(ax, 'FFaceColor', oceanColor,'FontSize', 20,'FontName','Calibri')
geoshow(states)

%Colormap is defined as a 3 column matrix, each row being an RGB triplet
% map = zeros(numel(fcsmatch.temperature),3);
% map(:,1)=1;
% map(:,2)=0;
% map(:,3)=fcsmatch.temperature./max(fcsmatch.temperature);
% %Set the current Colormap
% cmap = colormap(map)
% colorbar
%Display Colorbar
[~,I] = sort(Attune.fcsmatch.mdate_start);
geoshow(Attune.fcsmatch.lat(I),Attune.fcsmatch.lon(I),'LineWidth',3)
colorbar

ylabel('Sea Surface Temperature (Celsius)','FontName','Calibri')
%Annotations

% Coordinates for first arrow
x = [0.57 0.57] ;
y = [0.2 0.55] ;

%Coordinates for second arrow
x2 = [0.475 0.475 ] ;
y2 = [0.55 0.2] ;

annotation('arrow',x,y, 'LineStyle', '--','LineWidth',4, 'HeadStyle','plain','Color','white')
annotation('arrow',x2,y2','LineWidth',4, 'HeadStyle','plain','Color','white')
% text(0.5,0.5,'Northbound')
annotation('textbox',...
    [0.375000000 0.320247933884298 0.0620914826498422 0.0826446280991737],...
    'Color',[1 1 1],...
    'String',{'Southbound','Jan 31- Feb 3'},...
    'LineStyle','none',...
    'FontName','Calibri',...
    'FitBoxToText','off','FontSize',19);
annotation('textbox',...
    [0.580 0.34297520661157 0.0623543638275507 0.0562024793388431],...
    'Color',[1 1 1],...
    'String',{'Northbound','Feb 3 - Feb 5'},...
    'LineStyle','none',...
    'FontName','Calibri',...
    'FitBoxToText','off','FontSize',19);

hold on

map = zeros(numel(Attune.fcsmatch.temperature),3);
map(:,1)=1;
map(:,2)=0;
map(:,3)=Attune.fcsmatch.temperature./max(Attune.fcsmatch.temperature);
%Set the current Colormap
c= colormap(map);
%Display Colorbar
scatter(Attune.fcsmatch.lon,Attune.fcsmatch.lat,1,c)
colorbar
