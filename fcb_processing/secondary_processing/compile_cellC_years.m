summarypath = '\\sosiknas1\Lab_data\MVCO\FCB\summary\';
filelist = dir([summarypath 'cellC*.mat']);

matdate = [];
ml_analyzed = [];
sumC = [];
sumnum = [];
sumvol = [];

for count = 1:length(filelist)
    fname = filelist(count).name;
    temp = load([summarypath fname]);
    matdate = [matdate; temp.cellresults_all(:,1)];
    ml_analyzed = [ml_analyzed; temp.cellresults_all(:,3)];
    sumC = [sumC; temp.sumC_all];
    sumnum = [sumnum; temp.sumnum_all];
    sumvol = [sumvol; temp.sumvol_all];
end
clear count

[matdate, sind ]= sort(matdate);
ml_analyzed = ml_analyzed(sind);
sumC = sumC(sind,:);
sumnum = sumnum(sind,:);
sumvol = sumvol(sind,:);
sumtitles = temp.sumtitles;
dt = datetime(matdate, 'ConvertFrom', 'datenum');
TTsumC = table2timetable(array2table(sumC,'VariableNames',sumtitles), 'RowTimes', dt);
TTsumC.ml = ml_analyzed;
TTsumnum = table2timetable(array2table(sumnum,'VariableNames',sumtitles), 'RowTimes', dt);
TTsumnum.ml = ml_analyzed;
TTsumC(isnat(TTsumC.Time),:) = [];
TTsumnum(isnat(TTsumnum.Time),:) = [];

notes = {'Carbon units, picograms'; 'Created with compile_cellC_years.m after cellC_save3.m'; 'Heidi M. Sosik, Woods Hole Oceanographic Institution'; datestr(now)}
save([summarypath 'FCB_compiledC_tables'], 'TTsumC', 'TTsumnum', 'notes')

ind = find(~isnan(ml_analyzed));
[ mdate_mat, SynC_mat, yearlist, yd ] = timeseries2ydmat( matdate(ind), sumC(ind,1));
[ mdate_mat, ml_mat, yearlist, yd ] = timeseries2ydmat( matdate(ind), ml_analyzed(ind));
figure
plot(mdate_mat(:), SynC_mat(:)./ml_mat(:)/1000, '-')

%%
[ mdate_mat, Synnum_mat, yearlist, yd ] = timeseries2ydmat( matdate(ind), sumnum(ind,1));
[ mdate_mat, ml_mat, yearlist, yd ] = timeseries2ydmat( matdate(ind), ml_analyzed(ind));
figure
semilogy(mdate_mat(:), Synnum_mat(:)./ml_mat(:), 'b-','linewidth',1)
set(gcf, 'position', [50 400 1200 220])
datetick
ylim([10 1e6])
set(gca, 'ytick', [1e2 1e4 1e6], 'xgrid', 'on', 'fontsize', 12)
set(gca, 'position', [0.07    0.1335    0.9    0.8083])
ylabel('Cells ml^{-1}')
text(datenum(2003,9,1), 100, '\itSynechococcus', 'fontsize', 16, 'fontweight', 'bold')
grid on

%%
ind = find(~isnan(ml_analyzed));
[ mdate_mat, Euknum_mat, yearlist, yd ] = timeseries2ydmat( matdate(ind), sumnum(ind,2));
[ mdate_mat, ml_mat, yearlist, yd ] = timeseries2ydmat( matdate(ind), ml_analyzed(ind));
figure
semilogy(mdate_mat(:), Euknum_mat(:)./ml_mat(:), 'b-','linewidth',1)
set(gcf, 'position', [50 400 1200 220])
datetick
ylim([1e3 2e5])
set(gca, 'ytick', [1e2 1e4 1e6], 'xgrid', 'on', 'fontsize', 12)
set(gca, 'position', [0.07    0.1335    0.9    0.8083])
ylabel('Cells ml^{-1}')
text(datenum(2003,9,1), 2000, 'Eukaryotes', 'fontsize', 16, 'fontweight', 'bold')
grid on
%%

return
%semilogy(time_syn(:), daily_syn(:), 'b', 'linewidth', 1)
semilogy(matdate, sumnum(:,1)./ml_analyzed, 'b', 'linewidth', 1)
ylim(([10 1e6]))
set(gca, 'ytick', [1e2 1e4 1e6], 'xgrid', 'on', 'fontsize', 14)
ylabel('\itSynechococcus\rm, ml^{-1}')
set(gcf, 'position', [288.2000  400  940.8000  238.4000])
set(gca, 'position', [0.0884    0.1335    0.8806    0.8083])
set(gca, 'xlim', [datenum('1-0-2003') datenum('1-1-2018')])
datetick

figure

load("\\sosiknas1\Lab_data\MVCO\EnvironmentalData\CTD_tables.mat")
load('\\sosiknas1\Lab_data\MVCO\EnvironmentalData\Tall_day.mat')

ii = find(MVCO_ctd_table.dates> datetime('1-Jan-2019'));
unqday = unique(MVCO_ctd_table.dates(ii));
for iii = 1:length(unqday)
    Tday_new(iii) = mean(MVCO_ctd_table.temperature(ii(MVCO_ctd_table.dates(ii) == unqday(iii))));
end


T = [Tday2; Tday_new'];
T_mdate = [mdate2; datenum(unqday)];
Tinterp = interp1(T_mdate, T, matdate);

ind = find(~isnan(ml_analyzed));
[ mdate_mat, Syn_mat, yearlist, yd ] = timeseries2ydmat( matdate(ind), sumnum(ind,1));
[ mdate_mat, ml_mat, yearlist, yd ] = timeseries2ydmat( matdate(ind), ml_analyzed(ind));
Tinterp = interp1(T_mdate, T, mdate_mat(:));

dv = datevec(mdate_mat(:));

%%
figure
scatter(Tinterp, Syn_mat(:)./ml_mat(:),10,dv(:,2), 'filled')
cbh = colorbar
set(cbh, 'xdir', 'rev')
caxis([1 12])
set(cbh, 'ticklabels', {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'}')
set(gca, 'yscale', 'log')
ylabel('\itSynechococcus\rm concentration (ml^{-1})', 'fontsize', 14)
xlabel('Temperature (\circC)', 'fontsize', 14)
set(gca, 'box', 'on')
title('Inner shelf, MVCO 2003-2021', 'fontsize', 14)
axis([-2 28 20 1e6])

