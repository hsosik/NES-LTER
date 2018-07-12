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




return
%semilogy(time_syn(:), daily_syn(:), 'b', 'linewidth', 1)
semilogy(matdate, sumnum(:,1)./ml_analyzed, 'b', 'linewidth', 1)
ylim(([10 1e6]))
set(gca, 'ytick', [1e2 1e4 1e6], 'xgrid', 'on', 'fontsize', 14)
ylabel('\itSynechococcus\rm, ml^{-1}')
set(gcf, 'position', [288.2000  447.4000  940.8000  238.4000])
set(gca, 'position', [0.0884    0.1335    0.8806    0.8083])
set(gca, 'xlim', [datenum('1-0-2003') datenum('1-1-2018')])
datetick


