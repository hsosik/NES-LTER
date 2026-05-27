p = '\\sosiknas1\Lab_data\Attune\cruise_data\IFCB_Attune_merge\summary_files\';
pout = 'C:\work\IFCB_products\Attune_IFCB_merge\';
flist = dir([p '*_IFCB_Attune_merge_sum.mat']);
flist = {flist.name}';
flist(contains(flist', 'AR78')) = []; %skip for now (missing IFCB metadata)

outname1 = regexprep(flist, '_IFCB_Attune_merge_sum.mat', '');
%outname = strcat(outname1, '_IFCB_Attune_merge_byGroup');
%%
uwlist = dir(['\\sosiknas1\IFCB_data\\*\match_up\*uw_match.mat']);
uwlist(contains({uwlist.name}', 'STAINING')) = [];

uwfile = cell(size(flist));
for count = 1:length(flist)
    tempStr = outname1{count}(10:end);
    tempind = contains({uwlist.name}', tempStr);
    tempind = find(tempind);
    temp = fullfile(uwlist(tempind(1)).folder,uwlist(tempind(1)).name);
    for ii = 2:length(tempind)
        temp = [temp; fullfile(uwlist(tempind(ii)).folder,uwlist(tempind(ii)).name)];
    end
    uwfile{count} = temp; clear temp 
end
%%
for count = 1:length(flist)
    %cruiseStr = fList{count};
    %f = dir([p flist{count}]);
    T{count} = load([p flist{count}]);
    T{count}.f = flist(count);
    if contains(flist(count), '_EN627')
        varstr = {'tsg2_temperature' 'tsg2_salinity'};
    elseif contains(flist(count), '_AE2426')
        varstr = {'sbe45pri_temp_c' 'sbe45pri_sal_psu'};    
    elseif contains(flist(count), '_EN')
        varstr = {'tsg1_temperature' 'tsg1_salinity'};
    elseif contains(flist(count), '_AR43')
        varstr = {'sbe45t' 'sbe45s'};
    elseif contains(flist(count),'_AR') | contains(flist(count), '_AT46')
        varstr = {'sbe48t' 'sbe45s'};
    elseif contains(flist(count), '_GU') | contains(flist(count), '_PC') | contains(flist(count), '_TN')...
            | contains(flist(count), '_HB') | contains(flist(count), '_RB')
        varstr = {'ts' 'ssps'};
    elseif contains(flist(count), 'HRS')
        varstr = {'water_temperature_degree_c' 'salinity_psu'};
    end
    uw = load(uwfile{count}(1,:));
    for ii = 2:size(uwfile{count},1)
        uw2 = load(uwfile{count}(ii,:));
        uw.IFCB_match_uw_results = [uw.IFCB_match_uw_results; uw2.IFCB_match_uw_results];
        clear uw2
    end
    T{count}.ifcbHmeta.dt = datetime(T{count}.ifcbHmeta.sample_time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');
    T{count}.ifcbHmeta.mdate = datenum(T{count}.ifcbHmeta.dt);
    if isfield(T{count},'attunemlSyn') %case for alternating run (e.g., EN688)
        T{count}.attuneml = T{count}.attunemlEuk;
        %adjust the Syn sums to the equivalent ml run for the Euks so everything can be combined with IFCB in a consistent way
        T{count}.attuneHcarbon.Syn = T{count}.attuneHcarbon.Syn./T{count}.attunemlSyn'.*T{count}.attunemlEuk';
        T{count}.attuneHcount.Syn = T{count}.attuneHcount.Syn./T{count}.attunemlSyn'.*T{count}.attunemlEuk';
    end

% trnum = 1:length(tstime);
% T{count}.ifcbHmeta.transect = floor(interp1(tstime,trnum, T{count}.ifcbHmeta.mdate));
% clear trnum tstime
[aa,bb] = ismember(T{count}.ifcbHmeta.pid,uw.IFCB_match_uw_results.pid);
    T{count}.uw_match = uw.IFCB_match_uw_results(bb,:);
    T{count}.uw_match.Properties.VariableNames = lower(T{count}.uw_match.Properties.VariableNames);
    T{count}.ifcbHmeta.temp = T{count}.uw_match.(varstr{1});
    T{count}.ifcbHmeta.sal = T{count}.uw_match.(varstr{2});
 end
clear f count uw
%%
binedges = T{1}.binedges;
binrange = [binedges(1:end-1); binedges(2:end); 1:length(binedges)-1];
load('\\sosiknas1\Lab_data\Attune\cruise_data\IFCB_Attune_merge\summary_files\cutoffTable_bycruise')
%%
conc = [];
Cconc = [];
meta = table;
Cconc_taxa = table; 
conc_taxa = table;
sumT = table;
for count = 1:length(flist)
    A = T{count};
    cutind = strmatch(outname1(count), cutoffTable.cruise);
    bin2use = zeros(2,size(binrange,2)); 
    bin2use(1,1:find(binrange(2,:)==cutoffTable.Attunemax(cutind))) = 1;
    bin2use(2,find(binrange(1,:)==cutoffTable.IFCBmin(cutind)):end) = 1;
    attuneml_bybin = bin2use(1,:).*A.attuneml';
    ifcbml_bybin = bin2use(2,:).*A.ifcbHmeta.ml_analyzed;
    totalml_bybin = attuneml_bybin+ifcbml_bybin;
    totalml_bybin(attuneml_bybin(:,1)==0,:) = NaN; %skip cases with no Attune match up
    
   % ifcbHcount = [A.ifcbHcount.diatom+A.ifcbHcount.dino+A.ifcbHcount.nano+A.ifcbHcount.ciliate+A.ifcbHcount.otherPhyto];
   % ifcbHcarbon = [A.ifcbHcarbon.diatom+A.ifcbHcarbon.dino+A.ifcbHcarbon.nano+A.ifcbHcarbon.ciliate+A.ifcbHcarbon.otherPhyto];
    ifcbHcount = A.ifcbHcount.protist_tricho;
    ifcbHcarbon = A.ifcbHcarbon.protist_tricho;
    attuneHcount = A.attuneHcount.Syn+A.attuneHcount.Euk;
    attuneHcarbon = A.attuneHcarbon.Syn+A.attuneHcarbon.Euk;
    if ~isreal(attuneHcarbon)
        attuneHcarbon = real(attuneHcarbon);
        disp(['imaginary values in ' flist(count)]);
    end
    totalHcount = attuneHcount.*bin2use(1,:)+ifcbHcount.*bin2use(2,:);
    totalHcarbon = attuneHcarbon.*bin2use(1,:)+ifcbHcarbon.*bin2use(2,:);
    conc = [conc; totalHcount./totalml_bybin];
    Cconc = [Cconc; totalHcarbon./totalml_bybin];
    %%conc = totalHcount./totalml_bybin;
    %%Cconc = totalHcarbon./totalml_bybin;
    % TT = array2table([sum(A.ifcbHcarbon.diatom,2)./A.ifcbHmeta.ml_analyzed sum(A.ifcbHcarbon.dino,2)./A.ifcbHmeta.ml_analyzed...
    % sum(A.ifcbHcarbon.ciliate,2)./A.ifcbHmeta.ml_analyzed A.attuneHcarbon.Syn(:,1)./A.attuneml']/1000, 'VariableNames',{'diatom' 'dino' 'ciliate' 'Syn'});    
    % Cconc_taxa = TT;
    % TT = array2table([sum(A.ifcbHcount.diatom,2)./A.ifcbHmeta.ml_analyzed sum(A.ifcbHcount.dino,2)./A.ifcbHmeta.ml_analyzed...
    % sum(A.ifcbHcount.ciliate,2)./A.ifcbHmeta.ml_analyzed A.attuneHcount.Syn(:,1)./A.attuneml'], 'VariableNames',{'diatom' 'dino' 'ciliate' 'Syn'});    
    % conc_taxa = TT;
    sumT = [sumT; A.ifcbHmeta(:,["cast" "cruise" "depth" "dt" "ifcb" "latitude" "longitude" "mdate" "ml_analyzed" ...
        "niskin" "pid" "sal" "sample_type" "skip" "temp"])];
end
labels = regexprep(cellstr(strcat('phytoC_', num2str(binrange(1,:)'), '_', num2str(binrange(2,:)'))), ' ', '');
sumT = [sumT array2table(Cconc/1000, 'VariableNames', labels)];
sumT.phytoC_total = sum(sumT{:,find(strncmp('phyto', cellstr(sumT.Properties.VariableNames),5))},2);
itemp = startsWith(sumT.Properties.VariableNames, {'phytoC_2_3' 'phytoC_3_5' 'phytoC_5_7' 'phytoC_7_10' 'phytoC_10_15' 'phytoC_15_20'});
sumT.phytoC_2_20 = sum(sumT{:,itemp},2);
itemp = startsWith(sumT.Properties.VariableNames, {'phytoC_20_30' 'phytoC_30_40' 'phytoC_40_50' 'phytoC_50_100'});
sumT.phytoC_20_100 = sum(sumT{:,itemp},2);
itemp = startsWith(sumT.Properties.VariableNames, {'phytoC_0_2' 'phytoC_2_3'});
sumT.phytoC_0_3 = sum(sumT{:,itemp},2);
itemp = startsWith(sumT.Properties.VariableNames, {'phytoC_3_5' 'phytoC_5_7' 'phytoC_7_10'});
sumT.phytoC_3_10 = sum(sumT{:,itemp},2);
itemp = startsWith(sumT.Properties.VariableNames, {'phytoC_10_15' 'phytoC_15_20' 'phytoC_20_30' 'phytoC_30_40' 'phytoC_40_50' 'phytoC_50_100'});
sumT.phytoC_10_100 = sum(sumT{:,itemp},2);
sumT.depth = [];
sumT.cast = [];
sumT.niskin = [];
sumT.skip = [];

phytoC_all = sumT;
notes = {'Carbon concentration in micrograms per liter'; "Cell size ranges in micrometers equivalent spherical diameter"};
save([p 'Attune_IFCB_merge_class_summary_uw_allcruises.mat'], 'phytoC_all', 'notes')
%save("c:\work\IFCB_products\NESLTER_transect\summary\Attune_IFCB_merge_class_summary_uw_allcruises.mat", 'phytoC_all', 'notes')

phytoC_spiropa = sumT(ismember(sumT.cruise, {'AR29' 'RB1904' 'TN368'}),:);
notes = {'Carbon concentration in micrograms per liter'; "Cell size ranges in micrometers equivalent spherical diameter"};
save([p 'Attune_IFCB_merge_class_summary_uw_spiropa'], 'phytoC_spiropa', 'notes')
%save("c:\work\IFCB_products\NESLTER_transect\summary\Attune_IFCB_merge_class_summary_uw_spiropa", 'phytoC_spiropa', 'notes')


return
%%
cruiseStr = 'AR29'; % 'AR29' 'RB1904' 'TN368' 
ii = find(meta.cruise==cruiseStr & meta.longitude>-70.84 & meta.longitude<-70.82 & meta.latitude<40.8);
w = 0.2;
cstrList = {'pico' 'nano' 'micro' 'Syn' 'diatom' 'dino' };
Z2_label = 'Carbon (\mug l^{-1})'; Z2str = 'Cconc';
%Z2_label = 'Concentration (ml^{-1})'; Z2str = 'conc';
for count = 1:length(cstrList)
    cstr = cstrList{count};
switch cstr
    case {'pico' 'nano' 'micro'}
        Z2 = eval([cstr Z2str]); dl = [0 inf];
    case 'nano'
         dl = [0 200];
    case 'micro'
         dl = [0 20];
    case {'Syn' 'diatom' 'dino'}
        Z2 = eval([Z2str '_taxa.(cstr)']); dl = [0 inf];   
    end

figure('position', [100 200 1100 430]), tlh = tiledlayout(1,3);

tl(1) = nexttile;
%scatter(uw1_match.SSPS(ii),uw1_match.lat(ii), 10,uw1_match.matdate(ii), 'filled'), xlim([32 36])
scatter(meta.sal(ii),meta.latitude(ii), 10,meta.mdate(ii), 'filled'), xlim([32 36])
ylim([39.5 40.8])
line([34.5 34.5], ylim, 'color', 'k')
ylabel('Latitude'), xlabel('Salinity')
tl(2) = nexttile;
scatter(Z2(ii),meta.latitude(ii),10,meta.mdate(ii),'filled')
ylim([39.5 40.8])
ylabel('Latitude'), xlabel(Z2_label)
tl(3) = nexttile;
scatter(meta.sal(ii), Z2(ii),10,meta.mdate(ii),'filled'), xlim([32 36])
hold on
X = meta.sal(ii); Y = Z2(ii);
Xi = min(X):.05:max(X); Yi = NaN(size(Xi));
for cc =1:length(Xi), tt = find(X>Xi(cc)-w & X<Xi(cc)+w); Yi(cc) = median(Y(tt),'omitnan');end
plot(Xi,Yi,'-k', 'linewidth', 2)
line([34.5 34.5], ylim, 'color', 'k')
ylabel(Z2_label), xlabel('Salinity')

cb = colorbar(tl(end)); datetick(cb, 'y')

set(tl, 'Colormap', winter, 'CLim', cb.Limits)
cb.Layout.Tile = 'east'; 
title(tlh, [cruiseStr ' ' cstr])
print([pout cruiseStr '_' cstr '_' Z2str '_scatter'], '-dpng')

figure('Position', [350 180 400 400])
tl = tiledlayout(2,1);
nexttile
wl = 0.1;
ilat = (39.5:wl:41)';
[u,v] = discretize(meta.latitude, ilat-wl/2);
lat_smooth = NaN(size(meta.latitude));
lat_smooth((~isnan(u))) = ilat(u(~isnan(u)));

boxplot([Z2(ii); NaN(size(ilat))],[lat_smooth(ii); ilat], 'datalim', dl, 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
lines = findobj(gca, 'type', 'line', 'Tag', 'Median');
hold on
xMed = mean(vertcat(lines.XData),2); % n x 1
yMed = vertcat(lines.YData); % n x 2 (duplicate columns)
plot(xMed, yMed(:,1), 'r-')
set(gca, 'xdir', 'rev')
xlabel('Latitude')
ylim([0 inf])

nexttile
ws = 0.3;
isal = (32:ws:36.2)';
[u,v] = discretize(meta.sal, isal-ws/2);
sal_smooth = NaN(size(meta.sal));
sal_smooth((~isnan(u))) = isal(u(~isnan(u)));

boxplot([Z2(ii); NaN(size(isal))],[sal_smooth(ii); isal], 'datalim', dl, 'extrememode', 'compress', 'notch', 'on', 'LabelOrientation', 'horizontal');
lines = findobj(gca, 'type', 'line', 'Tag', 'Median');
hold on
xMed = mean(vertcat(lines.XData),2); % n x 1
yMed = vertcat(lines.YData); % n x 2 (duplicate columns)
plot(xMed, yMed(:,1), 'r-', 'linewidth',2)
xlabel('Salinity')
xint = interp1(isal,get(gca, 'xtick'), 34.5);
yl = ylim; yl = [0 yl(2)];
ylim(yl)
hold on, line(xint*[1 1], ylim, 'color', 'k')
title(tl, [cruiseStr ' ' cstr])
ylabel(tl, Z2_label)
print([pout cruiseStr '_' cstr '_' Z2str '_box'], '-dpng')

end

%%
figure('position', [200 200 900 400]), tlh = tiledlayout(1,3)
for count = 1:3
    nexttile, cruiseStr = cruiseList{count};
    ind = find(meta.cruise==cruiseStr & meta.longitude>-70.84 & meta.longitude<-70.82);
    plot(meta.longitude(meta.cruise==cruiseStr), meta.latitude(meta.cruise==cruiseStr), '.')
    hold on
    plot(meta.longitude(ind), meta.latitude(ind), '.')
    title(cruiseStr)
    axis([-71.5 -70 38.5 41.5])
end

return
%%
%RB tstime indices for main transect
%2 has a jog to the east, then a slight offline diagonal
%3 and 4 are short segments on the north of the line
%27 is on main line but only two IFCB samples
%33,34,38,39,40 on main line but short and few samples
%seems like 43 and 44 could be merged for underway
%seems like 49 and 50 could be merged for underway
[1:6 12:14 32 43 47:51 64]
ii = ismember(ifcbHmeta.transect, [1:6 12:14 32 43 47:51 64]);

%%
%AR29
% transects 1-9, 12, 13
% 14 maybe, 20, 21, 28, 29, 34, 35, 36, 37, 38
% 39
% 40 (good transect plus a crossing??)
% 51
% [1:9 12:14 20 21 28 29 34:40 51]
% 
% ii = ismember(ifcbHmeta.transect, [1:9 12:14 20 21 28 29 34:40 51]);