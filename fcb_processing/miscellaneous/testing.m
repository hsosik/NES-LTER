loglog(partialdatmerged(datbins2,5), partialdatmerged(datbins2,2), '.')
ii = find(partialdatmerged(datbins2,4)./partialdatmerged(datbins2,2)>=5);
hold on
loglog(partialdatmerged(datbins2(ii),5), partialdatmerged(datbins2(ii),2), '.g')

maxvalue = 1e6;
bins = 10.^(0:log10(maxvalue)/63:log10(maxvalue));
for c = 1:length(bins)-1
    iii = find(partialdatmerged(datbins2(ii),5) > bins(c) & partialdatmerged(datbins2(ii),5) < bins(c+1)); 
    nm(c) = length(iii); 
    m2(c) = median(partialdatmerged(datbins2(ii(iii)),2));
end
plot(bins(1:end-1),100+smooth(m2',5))


maxvalue = 5e6;
bins = 10.^(0:log10(maxvalue)/31:log10(maxvalue));
ii2 = find(partialdatmerged(datbins2,4)./partialdatmerged(datbins2,2)<=4 & partialdatmerged(datbins2,2)>100 & partialdatmerged(datbins2,4)./partialdatmerged(datbins2,2)>.1);
h = histmulti5([partialdatmerged(datbins2(ii2),4),partialdatmerged(datbins2(ii2),2)],[bins; bins]');
[n,p] = max(h);
t = find(n>5);
fitobject = fit(log10(bins(t))', log10(bins(p(t)))', 'poly1');
ci = confint(fitobject);
c = coeffvalues(fitobject);
figure
plot(log10(bins(p(t))),log10(bins(t)), '.')
hold on
yest = bins(t).^(c(1)).*10^(c(2));
plot(log10(yest), log10(bins(t)), 'r-')
yest = (bins(t)).^(c(1)).*10^(ci(2,2));
plot(log10(yest), log10(bins(t)), 'g-')


ii = find(partialdatmerged(datbins2,4)>= partialdatmerged(datbins2,2).^c(1).*10^(c(2)+.5));
h = histmulti5([partialdatmerged(datbins2(ii),5),partialdatmerged(datbins2(ii),2)],[bins; bins]');
[nt,pt] = max(h');
mask = ones(size(h));
for idx = 1:length(bins), mask(idx,1:pt(idx)) = NaN; end;
h2 = h.*mask;
h3 = zeros(size(h2));
minv = max([repmat(4,length(bins),1)'; min(h2')*1.2])';
v = repmat(minv,1, 32);
h3(h2<v) = 1;
h3(h2>=v) = 0;
[~,ph] = max(h3')
s = find(nt>5);
figure
loglog((partialdatmerged(datbins2(ii),5)), (partialdatmerged(datbins2(ii),2)), '.r')
hold on
plot(bins(s), bins(ph(s)), '*g')

yint = interp1(bins(s), bins(ph(s)), bins, 'pchip')
fitobject = fit(log10(bins)', log10(yint)', 'smoothingspline', 'SmoothingParam', 0.99);
%fitobject = fit(log10(bins(s))', log10(bins(ph(s)))', 'smoothingspline', 'SmoothingParam', 0.99);
y = feval(fitobject, log10(bins));
plot(bins, 10.^y, 'r-')
