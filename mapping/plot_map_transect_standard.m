figure(1)
clf
set(gcf,'PaperPosition',[0 0 8 11]);
hold on,set(gca,'LineWidth',1.5,'FontSize',18,'Layer','top')

xlim1 = [-71.5 -69.5]; ylim1 = [39.5 42];
axis([xlim1 ylim1])
% plot bathymetry; pick your contour levels
load bathy lonb_big latb_big zb_big
[cb,hb] = contour(lonb_big,latb_big,zb_big,[20 40 60 80 100],'k-');
[cb2,hb2] = contour(lonb_big,latb_big,zb_big,[200 600 1000 1400 1800 2200 2600],'b-');
% load coastline
load small_MVCO_coastline_f.mat

%% set color and thickness of contour lines
%set(hb,'Color',[.4 .4 .4],'LineWidth',1)
% pick contour label positions manually, if you want
% hl = clabel(cb,hb,'manual','FontSize',18,'Color',.4*[1 1 1]);
% set(hb,'Color',[.4 .4 .4],'LineWidth',1.5)
% otherwise, auto-label
clabel(cb,hb,'FontSize',12,'Color',.4*[1 1 1],'LineWidth',1,'LabelSpacing',400);
clabel(cb2,hb2,'FontSize',12,'Color','b','LineWidth',1,'LabelSpacing',400);

%% plot islands as separate patches of grey (color .7*[1 1 1])
ind = find(isnan(ncst(:,1)));
for k = 1:length(ind)-1
    fill(ncst(ind(k)+1:ind(k+1)-1,1),ncst(ind(k)+1:ind(k+1)-1,2),.7*[1 1 1])
end
% outline the islands in black
h = plot(ncst(:,1),ncst(:,2),'Color','k','LineWidth',1);

% fix map to correct for deformation (so 1 km on y axis is the same size as 1 km on x axis)
% this routine is by Carlos Moffat.  It calls "sw_dist.m"
c_basescale(mean(ylim1));

%%
N12 = [41.3366 -70.5564];  %12-m node MVCO
ASIT = [41 19.5 70 34.0];  %tower MVCO
pioneer = [39.9371 -70.887; 39.9365, -70.8802; 40.36719, -70.8818; 40.3619, -70.87829; 39.9394, -70.7708;...
    40.3649, -70.76999; 40.0963, -70.87889; 40.2267, -70.87819; 40.13409, -70.7701; 40.1334, -70.7785];
GSO = [41 29.485 71 25.432];

LTER = [41 11.8 70 53; 41 1.8 70 53; 40 51.8 70 53; 40 41.8 70 53; 40 30.8 70 53; 40 21.8 70 53; 40 13.6 70 53; ...
    40 08.2 70 46.5; 40 5.9 70 53; 39 56.4 70 53; 39 46.4 70 53; 39 56.4 70 46.5; 40 21.9 70 46.5];
LTER_labels = {'L1'; 'L2'; 'L3'; 'L4'; 'L5'; 'L6'; 'L7'; 'L8'; 'L9'; 'L10'; 'L11'; 'L12'; 'L13'};

LTERup = LTER([2,2,4,4,6,6,9,9,11,11],:); LTERdown = LTERup;
LTERup(:,4) =   [39.7 46.35 39.8 46.4 39.9 46.45 39.9 46.45 40 46.5]; %nudge up and down stream stations to ~10 nm by longitude 
LTERdown(:,4) = [66.3 59.65 66.2 59.6 66.1 59.55 66.1 59.55 66 59.5]; 
extra = [LTERup LTERdown]; clear LTERup LTERdown
extra = reshape(extra',4,length(extra(:))/4)';
extra_labels = {'u2b'; 'd2b'; 'u2a'; 'd2a'; 'u4b'; 'd4b'; 'u4a'; 'd4a'; 'u6b'; 'd6b'; 'u6a'; 'd6a'; 'u9b'; 'd9b'; 'u9a'; 'd9a'; 'u11b'; 'd11b'; 'u11a'; 'd11a'};
%extra_labels = {''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''; ''};

extra = [extra; ASIT];
extra_labels = [extra_labels; 'MVCO'];
 
temp = [GSO; LTER; extra];
station_table = table(['GSO'; LTER_labels; extra_labels], temp(:,1), temp(:,2), temp(:,3), temp(:,4), 'VariableNames', {'StnLabel' 'lat_deg' 'lat_min' 'lon_deg' 'lon_min'});

stn_order = {'GSO'    'L1'     'L2'     'u2b'    'u2a'    'd2a'   'd2b'    'L2'     'L3'      'L7'     'L8'     'L9'     'u9b'    'u9a'    'd9a'    'd9b'    'L9'     'L10'    'L12'    'L11'    'u11b'   'u11a'   'd11a'   'd11b'   'L11'     'L6'     'u6b'    'L13'   'd6a'    'd6b'     'L6'     'L5'      'L4'     'u4b'    'u4a'    'd4a'    'd4b'    'L4'     'L1'     'MVCO'   'L1'     'GSO'};
[~,stn_ind] = ismember(stn_order,station_table.StnLabel);

%decimal degrees
station_table.lat_dec_deg = station_table.lat_deg+station_table.lat_min/60;
station_table.lon_dec_deg = -station_table.lon_deg-station_table.lon_min/60;

H(1) = plot(station_table.lon_dec_deg, station_table.lat_dec_deg, '.r', 'markersize', 20);
H(2) = plot(station_table.lon_dec_deg(stn_ind), station_table.lat_dec_deg(stn_ind), '--r', 'linewidth', 2);

% label = station_table.StnLabel(stn_ind); 
% %next two lines to skip labeling up/down stns
% label(strmatch('u', label)) = {''};
% label(strmatch('d', label)) = {''};
% 
% text_offset = -2/60;
% %T = text(station_table.lon_dec_deg(stn_ind)+text_offset, station_table.lat_dec_deg(stn_ind), label, 'horizontalalignment', 'right', 'color', 'm', 'fontsize', 12, 'fontweight', 'bold');
% T = text(station_table.lon_dec_deg(stn_ind)+text_offset, station_table.lat_dec_deg(stn_ind), label, 'horizontalalignment', 'right', 'verticalalignment', 'top', 'color', 'm', 'fontsize', 12, 'fontweight', 'bold');

stn2label = {'GSO' 'L1' 'L2' 'L3' 'L4' 'L5' 'L6' 'L7' 'L9' 'L10' 'L11'};
[~,text_ind] = ismember(stn2label,station_table.StnLabel);
text_offset = -2/60;
T = text(station_table.lon_dec_deg(text_ind)+text_offset, station_table.lat_dec_deg(text_ind), station_table.StnLabel(text_ind), 'horizontalalignment', 'right', 'verticalalignment', 'top', 'color', 'm', 'fontsize', 12, 'fontweight', 'bold');

stn2label = {'MVCO' 'L12' 'L13'};
[~,text_ind] = ismember(stn2label,station_table.StnLabel);
text_offset = 2/60;
T = text(station_table.lon_dec_deg(text_ind)+text_offset, station_table.lat_dec_deg(text_ind), station_table.StnLabel(text_ind), 'horizontalalignment', 'left', 'verticalalignment', 'top', 'color', 'm', 'fontsize', 12, 'fontweight', 'bold');

stn2label = {'L8'};
[~,text_ind] = ismember(stn2label,station_table.StnLabel);
text_offset = 2/60;
T = text(station_table.lon_dec_deg(text_ind)+text_offset, station_table.lat_dec_deg(text_ind), station_table.StnLabel(text_ind), 'horizontalalignment', 'left', 'verticalalignment', 'bottom', 'color', 'm', 'fontsize', 12, 'fontweight', 'bold');

set(gca, 'box', 'on')
grid

H(11) = plot(pioneer(:,2), pioneer(:,1), 'b^', 'markersize', 10);
H(12) = plot(-GSO(:,3)-GSO(:,4)/60, GSO(:,1)+GSO(:,2)/60, 'gs', 'markersize', 10, 'markerfacecolor', 'g');

legend(H([1 11 12]), 'Station', 'Pioneer', 'GSO', 'location', 'southeast')
set(gca, 'fontsize', 10)
set(gcf, 'position', [400 40 550 750])
