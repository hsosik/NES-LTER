%run summarize_broadscale.m first
ind = (IFCBsum_meta.sample_type=='underway' & ~IFCBsum_meta.skip & IFCBsum_meta.ml_analyzed<5);

lat_smooth = round(IFCBsum_meta.latitude(ind),1);
ilat = unique(lat_smooth);
ilat = ilat(~isnan(ilat));
mn = month(IFCBsum_meta.datetime(ind));
yr = year(IFCBsum_meta.datetime(ind));

%Z2 = totalCconc/1000; titlestr = 'Plankton C (\mug l^{-1})'; dl = [0 200];
zstr = 'Ditylum_brightwellii'; dl = [0 2]; 
zstr = 'Guinardia_delicatula'; dl = [0 15];
zstr = 'Cerataulina_pelagica'; dl = [0 1];

%Z2 = IFCBcount.(zstr)(ind)./IFCBsum_meta.ml_analyzed(ind); 
%ystr = 'Concentration (ml^{-1})'; dl = [0 10];
zstr = 'Pseudo-nitzschia'; zstr2 = 'pennate_Pseudo-nitzschia'; dl = [0 1];
zstr = 'Hemiaulus'; zstr2 = 'Hemiaulus_membranaceus'; dl = [0 20];
%Z2 = (IFCBcount.(zstr)(ind)+IFCBcount.(zstr2)(ind))./IFCBsum_meta.ml_analyzed(ind);
Z2 = (IFCBCarbon.(zstr)(ind)+IFCBCarbon.(zstr2)(ind))./IFCBsum_meta.ml_analyzed(ind)/1000;
Z2 = IFCBCarbon.(zstr)(ind)./IFCBsum_meta.ml_analyzed(ind)/1000; 
ystr = 'Carbon (\mug ml^{-1})';

%%
tstr = {'Jan-Mar' 'Apr-Jun' 'Jul-Sep' 'Oct-Dec'};
%tstr = {'Feb-Mar' 'Apr-May' 'Jun-Aug' 'Sep-Nov'};
figure 
set(gcf,'position', [360 40 500 600])
tl = tiledlayout(4,1, 'TileSpacing', 'compact')
ccmonth = {[1:3] [4:6] [7:9] [10:12]};
%ccmonth = {[2:3] [4:5] [6:8] [9:11]};
for cc = 1:length(ccmonth) %2:2:9
    nexttile
    ii = find((mn>=ccmonth{cc}(1) & mn<=ccmonth{cc}(end)));
%    boxplot([Z2(ind(ii)); NaN(size(ilat))],[lat_smooth(ind(ii)); ilat], 'whisker', 3, 'datalim', [0 1.5e6], 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    boxplot([Z2(ii); NaN(size(ilat))],[lat_smooth(ii); ilat], 'whisker', 3, 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    lines = findobj(gca, 'type', 'line', 'Tag', 'Median');
    hold on
    xMed = mean(vertcat(lines.XData),2); % n x 1
    yMed = vertcat(lines.YData); % n x 2 (duplicate columns)
    plot(xMed, yMed(:,1), 'r-')
    set(gca, 'xdir', 'rev')
    if cc < length(ccmonth)
        set(gca,'xticklabel', [])
    end
    ylim(dl)
    text(5, dl(2)*.8, tstr{cc}, 'fontsize', 14)
end
set(gca, 'XTickLabelRotation',90)
set(gcf, 'paperposition', [.25 .25 6 10.5])
ylabel(tl, ystr )
xlabel(tl, 'Latitude')
title(tl, '2018-2022', 'fontsize', 14)

%%

tstr = {'Jan-Mar' 'Apr-Jun' 'Jul-Sep' 'Oct-Dec'};
%tstr = {'Feb-Mar' 'Apr-May' 'Jun-Aug' 'Sep-Nov'};
titlestr = '2007-2022';
for yy = 2017:2022
figure 
set(gcf,'position', [360 40 500 600])
tl = tiledlayout(4,1, 'TileSpacing', 'compact')
ccmonth = {[1:3] [4:6] [7:9] [10:12]};
%ccmonth = {[2:3] [4:5] [6:8] [9:11]};
for cc = 1:length(ccmonth) %2:2:9
    nexttile
    ii = find((mn>=ccmonth{cc}(1) & mn<=ccmonth{cc}(end)) & yr == yy);
    boxplot([Z2(ii); NaN(size(ilat))],[lat_smooth(ii); ilat], 'whisker', 3, 'datalim', dl, 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
    lines = findobj(gca, 'type', 'line', 'Tag', 'Median');
    hold on
    xMed = mean(vertcat(lines.XData),2); % n x 1
    yMed = vertcat(lines.YData); % n x 2 (duplicate columns)
    plot(xMed, yMed(:,1), 'r-')
    set(gca, 'xdir', 'rev')
    if cc < length(ccmonth)
        set(gca,'xticklabel', [])
    end
    ylim(dl)
    text(6, dl(2)*.7, tstr{cc}, 'fontsize', 14)
end
set(gca, 'XTickLabelRotation',90)
set(gcf, 'paperposition', [.25 .25 6 10.5])
ylabel(tl, ystr)
xlabel(tl, 'Latitude')
title(tl, yy, 'fontsize', 14)
end

%%

%%
lat_edges = [41.3 40.980 40.327 39.923 39.5];
group = zeros(size(Z2));
for g = 1:4
    group(IFCBsum_meta.latitude(ind)<=lat_edges(g) & IFCBsum_meta.latitude(ind)>lat_edges(g+1)) = g;
end
gstr = {'inner shelf', 'mid shelf', 'outer shelf', 'slope'};
%%
figure
tl = tiledlayout(1,4);
set(gcf, 'position', [100 250 1100 225])
for g = 1:4
    nexttile
    boxplot([Z2(group==g); NaN(12,1)],[mn(group==g); (1:12)'], 'whisker', 3, 'datalim', dl, 'extrememode', 'compress', 'notch', 'on');
    title(gstr(g))
    set(gca, 'XTickLabel', ['JFMAMJJASOND']')
    ylim(dl*1.2075)
    lines = findobj(gca, 'type', 'line', 'Tag', 'Median');
    hold on
    xMed = mean(vertcat(lines.XData),2); % n x 1
    yMed = vertcat(lines.YData); % n x 2 (duplicate columns)
    plot(xMed, yMed(:,1), 'r-', 'linewidth', 1.5)
end
title(tl, ['\it{' regexprep(zstr, '_', ' ') '} \rmbroadscale'])

ylabel(tl, ystr)