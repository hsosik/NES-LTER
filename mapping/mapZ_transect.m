function [] = mapZ_broadscale(input2plot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    m_proj('Mollweide','long',[-72.5 -69.5],'lat',[39.5 43]);
  %  m_proj('Mollweide','long',[-73 -69],'lat',[38.5 43]); %spiropa
    m_usercoast('gumby3','patch',[.8 .8 .8],'edgecolor','k');
    m_grid('box','fancy','tickdir','out','fontsize', 10, 'fontname','Times New Roman');
    m_elev('contour',[-200 -100],'edgecolor','b');
    [X,Y] = m_ll2xy(input2plot.lon, input2plot.lat); %convert positions of IFCB manual points to map coordinates
    Z = input2plot.Z;
    ii = find(isinf(Z));
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

