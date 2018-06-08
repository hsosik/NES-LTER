load \\10.100.100.109\work\SFD\code\FCM_match_summary
uw = load('\\multiproc\science_share\Underway_data\AR29_underway');

%Stations A1-A16
Stn_lat = [40.725 NaN NaN NaN NaN 40.400 40.335 40.270 40.205 40.140 40.075 40.010 39.945 39.880 39.815 39.750 39.685 39.620];
Stn_lon = [-70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 -70.83 ];

%1. 2018/04/17 13:15:04, Cast 1
%2. 2018/04/18 03:10:04, EK80
%3. 2018/04/18 10:30:04, Cast 15
%4. 2018/04/19 02:01:34, Cast 22
%5. 2018/04/19 10:40:04, Cast 29
%6. 2018/04/19 23:01:49, EK80
%7. 2018/04/20 04:18:00, EK80

transect_start = datenum({...
%'2018/04/17 13:15:04'... %use this if you want to go all the way back to shore
'2018/04/17 05:00:00'...
'2018/04/18 03:10:04'...
'2018/04/18 10:30:04'...
'2018/04/18 23:30:00'...
'2018/04/19 7:30:00'...
'2018/04/19 13:30:00'...
'2018/04/19 23:01:49'...
'2018/04/20 04:18:00'...
'2018/04/20 10:00:00'...
'2018/04/20 17:20:00'...
'2018/04/21 06:35:00'...
'2018/04/21 19:00:00'...
'2018/04/22 07:00:00'...
'2018/04/22 19:00:00'...
'2018/04/22 23:40:00'...
'2018/04/23 09:00:00'...
'2018/04/24 01:00:00'...
'2018/04/24 09:00:00'...
'2018/04/24 16:12:00'...
'2018/04/25 01:20:00'...
'2018/04/25 07:45:00'...
'2018/04/25 16:25:00'...
'2018/04/26 10:00:00'...
'2018/04/26 19:20:00'...
'2018/04/27 00:00:00'...
'2018/04/27 11:30:00'...
'2018/04/28 08:00:00'...
'2018/04/29 02:00:00'...
'2018/04/29 06:30:00'...
'2018/04/29 11:00:00'...

    });

for ii = 1:length(transect_start)-1
    aa{ii} = find(fcsmatch.mdate_start >=transect_start(ii) & fcsmatch.mdate_start < transect_start(ii+1));
    bb{ii} = find(uw.mdate >=transect_start(ii) & uw.mdate < transect_start(ii+1));
end

figure
plot(fcsmatch.mdate_start(aa{1}), fcsmatch.lat(aa{1}), 'linewidth', 2), hold on
for ii = 2:length(transect_start)-1
    plot(fcsmatch.mdate_start(aa{ii}), fcsmatch.lat(aa{ii}), 'linewidth', 2)
end;

km_per_degree = sw_dist([39 40], [-70.83 -70.83], 'km');
for ii = 1:length(transect_start)-2
    distrange = [15 120];
    if ii <= 11
        flrrange = [50 100];
        salrange = [32.5 36];
        temprange = [5 18];
        picorange = [0 3e4];
        nanorange = [0 8e3];
        synrange = [0 4e3];
    else
        flrrange = [50 100];
        salrange = [32.5 36.5];
        temprange = [5 20];
        picorange = [0 3e4];
        nanorange = [0 8e3];
        synrange = [0 1e4];
    end
    if ii>25
    figure, set(gcf, 'position', [488.2 41.8 560 740.8])
    subplot(6,1,1)
    d = (Stn_lat(1)-fcsmatch.lat(aa{ii}))*km_per_degree;
    d2 = (Stn_lat(1)-uw.lat(bb{ii}))*km_per_degree;
    labeld = (Stn_lat(1)-Stn_lat(6:16))*km_per_degree;
    plot(d2, uw.sbe48T(bb{ii}), 'linewidth', 2), xlim(distrange), ylim([temprange]), ylabel('Temperature \circC'), set(gca, 'xgrid', 'on')
    th = title(['Transect ' num2str(ii) '  ' datestr(fcsmatch.mdate_start(aa{ii}(1))) ' - ' datestr(fcsmatch.mdate_start(aa{ii}(end)))], 'fontsize', 8);
    tp = get(th, 'position'); tp(2) = tp(2)*1.2;
    set(th, 'position', tp)
    text(labeld,ones(11,1)*temprange(2)+1.5, {'A6' 'A7' 'A8' 'A9' 'A10' 'A11' 'A12' 'A13' 'A14' 'A15' 'A16'}, 'fontsize', 10)
    subplot(6,1,2)
    plot(d2, uw.sbe45S(bb{ii}), 'linewidth', 2), xlim(distrange), ylim(salrange), ylabel('Salinity'), set(gca, 'xgrid', 'on')
    subplot(6,1,3)
    plot(d2, uw.flr(bb{ii}), 'linewidth', 2), xlim(distrange), ylim(flrrange), ylabel('Chl fluor'),set(gca, 'xgrid', 'on')
    subplot(6,1,4)
    plot(d, SynConc(aa{ii})*1000,  '.', 'markersize', 10), xlim(distrange), ylim(synrange), ylabel('Syn ml^{-1}'),set(gca, 'xgrid', 'on')
    subplot(6,1,5)
    plot(d, (PicoeukConc(aa{ii})+UltraeukConc(aa{ii}))*1000,  '.', 'markersize', 10), xlim(distrange), ylim(picorange), ylabel('PicoEuk ml^{-1}'),set(gca, 'xgrid', 'on')
    subplot(6,1,6)
    plot(d, NanoeukConc(aa{ii})*1000, '.', 'markersize', 10), xlim(distrange), ylim(nanorange), ylabel('NanoEuk ml^{-1}'),set(gca, 'xgrid', 'on')
    xlabel('Distance offshore from A1 (km)')
    
%     subplot(5,2,2)
%     plot(fcsmatch.lat(aa{ii}), PicoeukConc(aa{ii})*1000, 'linewidth', 2), xlim(latrange), ylim([0 4e4]), ylabel('Picoeuk ml^{-1}')
%     subplot(5,2,4)
%     plot(fcsmatch.lat(aa{ii}), UltraeukConc(aa{ii})*1000, 'linewidth', 2), xlim(latrange), ylim([0 6e3]), ylabel('Ultraeuk ml^{-1}')
%     subplot(5,2,6)
%     plot(fcsmatch.lat(aa{ii}), NanoeukConc(aa{ii})*1000, 'linewidth', 2), xlim(latrange), ylim([0 1e4]), ylabel('Nanoeuk ml^{-1}')
    orient tall
    end
end
