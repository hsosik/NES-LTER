t1 = load('\\queenrose\mvco\mvco_jan2012\data\processed\beads\beadresults.mat')
g1 = load('\\queenrose\mvco\mvco_jan2012\data\processed\grouped\groupsum')
figure
plot(g1.cellresultsall(:,1), g1.cellNUMall(:,1)./g1.cellresultsall(:,3), 'g.-')
hold on
plot(g1.cellresultsall(:,1), g1.cellNUMall(:,4)./g1.cellresultsall(:,3), '.-')
line([t1.beadresults(:,2)'; t1.beadresults(:,2)'], ylim, 'color', 'r')