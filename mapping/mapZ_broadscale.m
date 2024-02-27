function [cbh] = mapZ_broadscale(input2plot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
persistent broad_bathy_Z broad_bathy_lat broad_bathy_lon
    lon = [-76.5 -65];
    lat = [36 44];
    m_proj('Mollweide','long',lon,'lat',lat);
    m_usercoast('gumby3','patch',[.8 .8 .8],'edgecolor','k');
    m_grid('box','fancy','tickdir','out','fontsize', 10, 'fontname','Times New Roman'); hold on
    if isempty(broad_bathy_Z)
        [broad_bathy_Z,broad_bathy_lat,broad_bathy_lon]=mygrid_sand2([lon lat ]); broad_bathy_lon = broad_bathy_lon-360;
    end
    m_contour(broad_bathy_lon,broad_bathy_lat,broad_bathy_Z,[-100 -200],'k-')
       
    [X,Y] = m_ll2xy(input2plot.lon, input2plot.lat); %convert positions of IFCB manual points to map coordinates
    Z = input2plot.Z;
    ii = find(isinf(Z) | Z==0);
    plot(X(ii),Y(ii), 'o', 'markersize',4, 'color', [.5 .5 .5])
    hold on
    Z(ii) = NaN;
    scatter(X,Y,20,Z,'filled'); 
    cbh = colorbar(gca);
    caxis(input2plot.caxis_range);
    if input2plot.caxis_log
       %caxis([log10(.1) log10(100)])
       %set(cbh, 'ytick', log10([.1 1 10 100]))
       t = get(cbh, 'ylim');
       set(cbh, 'ytick', [floor(t(1)):ceil(t(2))])
       set(cbh, 'yticklabel', 10.^(get(cbh, 'ytick')))
    end
    %set(cbh, 'position', [.88 .12 .03 .6])
    title(cbh, input2plot.colorbar_title, 'fontsize', 10)
    title(input2plot.title)
end

