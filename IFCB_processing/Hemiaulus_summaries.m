%ystr = '2018'; cruise = 'EN617';
ystr = '2019'; cruise = 'EN644';
%ystr = '2020'; cruise = 'EN655';
%ystr = '2021'; cruise = 'EN668';
%ystr = '2022'; cruise = 'EN687';

plot_flag = {'transect_uw'};
%plot_flag = {'cast'};
base = 'C:\work\IFCB_products\NESLTER_transect\summary\'; %cruise = 'EN644';%cruise = 'EN644'; 'EN617'
%yrstr = '2019';base = 'C:\work\IFCB_products\NESLTER_broadscale\summary\'; plot_flag = {'broadscale_uw'}; cruise = 'HB1902'; %'GU1905' 'GU1902' 'HB1902'
%yrstr = '2019'; base = 'C:\work\IFCB_products\SPIROPA\summary\'; cruise = 'TN368'; 
%yrstr = '2019';base = 'C:\work\IFCB_products\OTZ\summary\'; cruise = 'HB1907'; yrstr = '2019';
%
L = load([base 'summary_biovol_allHDF_min20_' ystr 'lists.mat']);
load([base 'summary_biovol_allHDF_min20_' ystr], 'meta_data', 'mdate');
ml_new = load([regexprep(base, 'C:\\work\\', '\\\\sosiknas1\\') 'ml_new_' ystr]);
[ia,ib] = ismember(meta_data.pid, ml_new.filelist);
if ~isequal(sum(ia),length(mdate))
    disp('missing some?')
    keyboard
end
meta_data.ml_analyzed(ia) = ml_new.ml_analyzed(ib);
L.mdate = mdate;
L.meta_data = meta_data;
clear mdate meta_data

base2 = regexprep(base, 'C:\\work\\IFCB_products', '\\\\sosiknas1\\IFCB_data');
base2 = regexprep(base2, 'summary', 'match_up');
if strcmp(plot_flag, {'cast'})
     f = dir([base2 '*' cruise '_cast_match.mat']);
else
    f = dir([base2 '*' cruise '*uw_match.mat']);
end
load([base2 f.name])


if strcmp(plot_flag, {'cast'})
    castind = find(L.meta_data.cruise==cruise & L.meta_data.sample_type=='cast' & ~L.meta_data.skip);
    typeind = castind; outfile = [cruise '_cast_Hemiaulus'];
else
    if strcmp(cruise, 'HB1907')
        uwind = find(strcmp(L.meta_data.cruise, cruise) & strcmp(L.meta_data.sample_type,'underway') & ~L.meta_data.skip); %OTZ
    else
        uwind = find(L.meta_data.cruise==cruise & L.meta_data.sample_type=='underway' & ~L.meta_data.skip);
    end
    typeind = uwind; outfile = [cruise '_underway_Hemiaulus'];
end

if strcmp(plot_flag, {'cast'})
    [ii,ia] = ismember(L.filelist(castind), IFCB_match_btl_results.pid);
    match = IFCB_match_btl_results(ia,:);
    match.Temp = match.t090c;
    match.Sal = match.sal00;
else
    [ii,ia] = ismember(L.filelist(uwind), IFCB_match_uw_results.pid);
    match = IFCB_match_uw_results(ia,:);
    switch cruise
    case {'EN617' 'EN644' 'EN655' 'EN668' 'EN687'}
        match.Temp = match.tsg1_temperature;
        match.Sal = match.tsg1_salinity;
    case {'GU1902' 'TN368' 'HB1907' 'GU1905' 'HB1902'}
        match.Temp = match.TS;
        match.Sal = match.SSPS;
    end
end

%if ~isequal(sum(ii), length(uwind))
%    disp('missing data?')
%    keyboard
%end
%match = IFCB_match_uw_results(ia,:);

%%
group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
[~,ia,ib] = intersect(group_table.CNN_classlist, L.class2use);
diatom_ind = ib(find(group_table.Diatom(ia)));
diatom_ind = setdiff(diatom_ind, strmatch('Thalassiosira_TAG_external_detritus', L.class2use));
dino_ind = ib(find(group_table.Dinoflagellate(ia)));
nano_ind = ib(find(group_table.Nano(ia) | group_table.flagellate(ia) | group_table.Coccolithophore(ia)));
otherPhyto_ind = [strmatch('Trichodesmium', L.class2use); strmatch('Pseudochattonella_farcimen', L.class2use); strmatch('Phaeocystis', L.class2use, 'exact')];
ciliate_ind = ib(find(group_table.Ciliate(ia)));
artifact_ind = ib(find(group_table.IFCBArtifact(ia)));
particle_ind = setdiff(1:length(L.class2use), artifact_ind);

check_ind = setdiff(1:length(L.class2use), [diatom_ind; dino_ind; nano_ind; otherPhyto_ind; ciliate_ind]);

%cc = strmatch('Hemiaulus', L.class2use);
%temp = cat(1,L.classFeaList{uwind,cc(1)}); 

%%
T = L.meta_data(typeind,:);
var2keep = {'sample_time' 'latitude' 'longitude' 'depth' 'pid' 'ml_analyzed' 'sample_type' 'cast' 'niskin'};
var2del = setdiff(T.Properties.VariableNames, var2keep);
T = removevars(T,var2del);
T.matdate = L.mdate(typeind);
T = movevars(T, 'matdate', 'after', 'sample_time');
T.Temperature = match.Temp;
T.Salinity = match.Sal;
ESDmin = 5;
scoremin_totalPhyto = 0;
scoremin = 0.5;
scoremin_Hemiaulus = 0.9; %maximizes F1 on broadscale evaluation set
for ibin = 1:length(typeind)
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),[diatom_ind; dino_ind; nano_ind; otherPhyto_ind]}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin_totalPhyto);
    T.totalPhytoC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter 
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),dino_ind}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin);
    T.dinoflagellateC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter;
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),diatom_ind}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin);
    T.diatomC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter;
    cc = strmatch('Hemiaulus', L.class2use);
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),cc}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin_Hemiaulus);
    T.HemiaulusC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter;
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),ciliate_ind}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin);
    T.ciliateC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter;
end

if exist('castind', 'var')
    [G, Gcast, Gniskin] = findgroups(T.cast, T.niskin);
    var2pool = {'totalPhytoC' 'dinoflagellateC' 'diatomC' 'HemiaulusC' 'ciliateC'};
    var2keep = {'sample_time' 'latitude' 'longitude' 'depth' 'pid' 'ml_analyzed' 'sample_type' 'cast' 'niskin' 'Temperature' 'Salinity'};
    var2del = setdiff(T.Properties.VariableNames, var2keep);
    [~,ia] = unique(G);
    Tpooled = T(ia,:);
    Tpooled = removevars(Tpooled,var2del);
    Tpooled.ml_analyzed = splitapply(@sum,T.ml_analyzed,G);
    for count = 1:length(var2pool)
        Tpooled.(var2pool{count}) = splitapply(@sum,T.(var2pool{count}).*T.ml_analyzed,G)./Tpooled.ml_analyzed;
    end
end

%%
notes = {['Only includes cells with ESD >= ' num2str(ESDmin) ' microns']; 'Phytoplankton carbon values in micrograms per liter'; 'diatom includes Hemiaulus'; 'total phyto includes diatoms and dinoflagellates'; 'created with Hemiaulus_summaries.m'};
IFCB_carbon_concentration = T;
classifier = L.classpath_generic;
save([base outfile], 'notes', 'IFCB_carbon_concentration', 'scoremin', 'scoremin_Hemiaulus', 'scoremin_totalPhyto', 'classifier')

%%
%broadscale
if strcmp(plot_flag, {'broadscale_uw'})
input2plot.caxis_range = [-1.5 2];
input2plot.caxis_log = 1;
input2plot.colorbar_title = {'Carbon'; '(mg L^{-1})'};
input2plot.lat = T.latitude;
input2plot.lon = T.longitude;

figure
tiledlayout(2,3)
nexttile
input2plot.title = { [cruise ' total phytoplankton'] };
input2plot.Z = log10(T.totalPhytoC);
mapZ_broadscale(input2plot)
nexttile
input2plot.title = { [cruise ' total diatom'] };
input2plot.Z = log10(T.diatomC);
mapZ_broadscale(input2plot)
nexttile
input2plot.title = { [cruise ' \itHemiaulus'] };
input2plot.Z = log10(T.HemiaulusC);
mapZ_broadscale(input2plot)
nexttile
input2plot.title = { [cruise ' dinoflagellate'] };
input2plot.Z = log10(T.dinoflagellateC);
mapZ_broadscale(input2plot)
nexttile
input2plot.title = { [cruise ' ciliate'] };
input2plot.Z = log10(T.ciliateC);
%input2plot.caxis_range = [-1.5 1];
mapZ_broadscale(input2plot)
nexttile
input2plot.title = { [cruise ' salinity'] };
input2plot.Z = match.SSPS;
input2plot.caxis_log = 0;
input2plot.caxis_range = [28 37];
input2plot.colorbar_title = {'Salinity'};
mapZ_broadscale(input2plot)

set(gcf,'WindowState','fullscreen')
print([base outfile], '-dpng', '-r300')
end

%%
%transect
if strcmp(plot_flag, {'transect_uw'})
input2plot.caxis_range = [-1.5 2];
input2plot.caxis_log = 1;
input2plot.colorbar_title = {'Carbon'; '(mg L^{-1})'};
input2plot.lat = T.latitude;
input2plot.lon = T.longitude;

figure
tiledlayout(2,3)
nexttile
input2plot.title = { [cruise ' total phytoplankton'] };
input2plot.Z = log10(T.totalPhytoC);
mapZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' total diatom'] };
input2plot.Z = log10(T.diatomC);
mapZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' \itHemiaulus'] };
input2plot.Z = log10(T.HemiaulusC);
mapZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' dinoflagellate'] };
input2plot.Z = log10(T.dinoflagellateC);
mapZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' ciliate'] };
input2plot.Z = log10(T.ciliateC);
%input2plot.caxis_range = [-1 1];
mapZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' salinity'] };
input2plot.Z = match.Sal;
input2plot.caxis_log = 0;
input2plot.caxis_range = [28 37];
input2plot.colorbar_title = {'Salinity'};
mapZ_transect(input2plot)

set(gcf,'WindowState','fullscreen')
print([base outfile], '-dpng', '-r300')
end

%%
%Transect casts
if strcmp(plot_flag, {'cast'})
input2plot.caxis_range = [0 2];
input2plot.caxis_log = 1;
input2plot.colorbar_title = {'Carbon'; '(mg L^{-1})'};
input2plot.lat = Tpooled.latitude;
input2plot.lon = Tpooled.longitude;
input2plot.depth = Tpooled.depth;

tl = tiledlayout(2,3);
nexttile
input2plot.title = { [cruise ' total phytoplankton'] };
input2plot.Z = log10(Tpooled.totalPhytoC);
sectionZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' total diatom'] };
input2plot.Z = log10(Tpooled.diatomC);
sectionZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' \itHemiaulus'] };
input2plot.Z = log10(Tpooled.HemiaulusC);
sectionZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' dinoflagellate'] };
input2plot.Z = log10(Tpooled.dinoflagellateC);
sectionZ_transect(input2plot)
nexttile
input2plot.title = { [cruise ' ciliate'] };
input2plot.Z = log10(Tpooled.ciliateC);
%input2plot.caxis_range = [-1 1];
sectionZ_transect(input2plot)

ylabel(tl, 'Depth (m)')
set(gcf,'WindowState','fullscreen')
print([base outfile], '-dpng', '-r300')
%%
cruise = 'EN644'; TT = load(['C:\work\IFCB_products\NESLTER_transect\summary\' cruise '_underway_Hemiaulus.mat']);
%cruise = 'EN687'; TT = load(['C:\work\IFCB_products\NESLTER_transect\summary\' cruise '_underway_Hemiaulus.mat']);
%cruise = 'HB1907'; TT = load(['C:\work\IFCB_products\OTZ\summary\' cruise '_underway_Hemiaulus.mat']);
%cruise = 'GU1902'; TT = load(['C:\work\IFCB_products\NESLTER_broadscale\summary\' cruise '_underway_Hemiaulus.mat']); %'GU1905' 'GU1902' 'HB1902'
%cruise = 'TN368'; TT = load(['C:\work\IFCB_products\SPIROPA\summary\' cruise '_underway_Hemiaulus.mat']);
T = TT.IFCB_carbon_concentration;
figure
input2plot.title = { [cruise ' \itHemiaulus'] };
input2plot.Z = log10(T.HemiaulusC);
input2plot.caxis_range = [-1.5 1.8];
input2plot.caxis_log = 1;
input2plot.colorbar_title = {'Carbon'; '(mg L^{-1})'};
input2plot.xlim = [28 37];
input2plot.ylim = [10 30];
input2plot.x = T.Salinity;
input2plot.y = T.Temperature;
ts_Z(input2plot)
print(['C:\work\IFCB_products\NESLTER_transect\summary\' cruise '_Hemiaulus_TS'], '-dpng', '-r300')

end
%%
return

%%
%Hemiaulus boxplot
cruise = 'EN644'; TT = load(['C:\work\IFCB_products\NESLTER_transect\summary\' cruise '_underway_Hemiaulus.mat']);
T = TT.IFCB_carbon_concentration;
figure('Position',[350 300 550 300])
wl = 0.1;
ilat = (39.8:wl:41.5)';
[u,v] = discretize(T.latitude, ilat-wl/2);
lat_smooth = NaN(size(T.latitude));
lat_smooth((~isnan(u))) = ilat(u(~isnan(u)));
boxplot([T.HemiaulusC; NaN(size(ilat))],[lat_smooth; ilat], 'notch', 'on', 'Whisker', 3, 'LabelOrientation','horizontal')
lines = findobj(gca, 'type', 'line', 'Tag', 'Median');
hold on
xMed = mean(vertcat(lines.XData),2); % n x 1
yMed = vertcat(lines.YData); % n x 2 (duplicate columns)
plot(xMed, yMed(:,1), 'r-', 'linewidth', 1.5)
set(gca, 'xdir', 'rev')
xlabel('Latitude')
ylim([0 75])
ylabel('\itHemiaulus\rm carbon (\mug l^{-1})')
%print([base 'EN644_HemiaulusC_boxplot'], '-dpng', '-r300')
print([base 'EN644_HemiaulusC_boxplot'], '-depsc', '-r600')

%%
ilat = 39.5:.05:41.5;
ilon = ones(size(ilat)).*-70.8855;
idpth = 0:1:200;
bdata = load('ngdc2.mat'); %from Gordon
ibdpth = griddata(bdata.lon',bdata.lat,bdata.h,ilon(:),ilat(:));
%%
figure('Position',[350 300 550 300])
input2plot.caxis_range = [0 70];
input2plot.caxis_log = 0;
input2plot.colorbar_title = {'Carbon'; '(mg L^{-1})'};
input2plot.lat = Tpooled.latitude;
input2plot.lon = Tpooled.longitude;
input2plot.depth = Tpooled.depth;
input2plot.title = { [cruise ' \itHemiaulus'] };
input2plot.Z = (Tpooled.HemiaulusC);
sectionZ_transect(input2plot)
plot(ilat,ibdpth,'k','linewidth',2);
cb = colorbar('Location', 'south');
set(cb, 'position', [.15 .25 .4 .07])
title(cb, '\itHemiaulus\rm carbon (\mug L^{-1})', 'fontsize', 10)
xlabel('Latitude')
ylabel('Depth (m)')
xlim([39.75 41.55])
ylim([0 150])
title('')
disp('delete the default colorbar from sectionZ_transect; then hit enter to save png')
pause
%print([base 'EN644_HemiaulusC_section'], '-dpng', '-r300')
print([base 'EN644_HemiaulusC_section'], '-depsc', '-r600')
