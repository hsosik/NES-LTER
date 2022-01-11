ystr = '2019';
%base = 'C:\work\IFCB_products\NESLTER_transect\summary\'; cruise = 'EN644';
base = 'C:\work\IFCB_products\NESLTER_broadscale\summary\'; cruise = 'GU1905';
%base = 'C:\work\IFCB_products\SPIROPA\summary\'; cruise = 'TN368';

L = load([base 'summary_biovol_allHDF_min20_' ystr 'lists.mat']);
load([base 'summary_biovol_allHDF_min20_' ystr], 'meta_data', 'mdate');
L.mdate = mdate;
L.meta_data = meta_data;
clear mdate meta_data

uwind = find(L.meta_data.cruise==cruise & L.meta_data.sample_type=='underway' & ~L.meta_data.skip);
%castind = find(L.meta_data.cruise=='TN368' & L.meta_data.sample_type=='cast' & ~L.meta_data.skip);
typeind = uwind; outfile = [cruise '_underway_Hemiaulus'];
%typeind = castind;

%%
group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
[~,ia,ib] = intersect(group_table.CNN_classlist, L.class2use);
diatom_ind = ib(find(group_table.Diatom(ia)));
dino_ind = ib(find(group_table.Dinoflagellate(ia)));
nano_ind = ib(find(group_table.Nano(ia) | group_table.flagellate(ia) | group_table.Coccolithophore(ia)));
otherPhyto_ind = [strmatch('Pseudochattonella_farcimen', L.class2use); strmatch('Phaeocystis', L.class2use, 'exact')];
ciliate_ind = ib(find(group_table.Ciliate(ia)));
artifact_ind = ib(find(group_table.IFCBArtifact(ia)));
particle_ind = setdiff(1:length(L.class2use), artifact_ind);

check_ind = setdiff(1:length(L.class2use), [diatom_ind; dino_ind; nano_ind; otherPhyto_ind; ciliate_ind]);

%cc = strmatch('Hemiaulus', L.class2use);
%temp = cat(1,L.classFeaList{uwind,cc(1)}); 

%%
T = L.meta_data(typeind,:);
var2keep = {'sample_time' 'latitude' 'longitude' 'depth' 'pid' 'ml_analyzed' 'sample_type'};
var2del = setdiff(T.Properties.VariableNames, var2keep);
T = removevars(T,var2del);
T.matdate = L.mdate(typeind);
T = movevars(T, 'matdate', 'after', 'sample_time');
ESDmin = 5;
scoremin = 0;
for ibin = 1:length(typeind)
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),[diatom_ind; dino_ind; nano_ind; otherPhyto_ind]}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin);
    T.totalPhytoC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter 
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),dino_ind}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin);
    T.dinoflagellateC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter;
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),diatom_ind}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin);
    T.diatomC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter;
    cc = strmatch('Hemiaulus', L.class2use);
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),cc}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin);
    T.HemiaulusC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter;
    temp = array2table(cat(1,L.classFeaList{typeind(ibin),ciliate_ind}), 'VariableNames', L.classFeaList_variables);
    ind = (temp.ESD>=ESDmin & temp.score>=scoremin);
    T.ciliateC(ibin) = sum(temp.cellC(ind))./T.ml_analyzed(ibin)/1000; %micrograms per liter;
end
notes = {['Only includes cells with ESD >= ' num2str(ESDmin) ' microns']; 'Phytoplankton carbon values in micrograms per liter'; 'diatom includes Hemiaulus'; 'total phyto includes diatoms and dinoflagellates'; 'created with Hemiaulus_summaries.m'};
IFCB_carbon_concentration = T;
save([base outfile], 'notes', 'IFCB_carbon_concentration')

%%
%broadscale
if 1
input2plot.caxis_range = [0 2];
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
input2plot.caxis_range = [-1 1];
mapZ_broadscale(input2plot)

set(gcf,'WindowState','fullscreen')
print([base outfile], '-dpng', '-r300')
end

%%
%transect
if 1
input2plot.caxis_range = [0 2];
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

set(gcf,'WindowState','fullscreen')
print([base outfile], '-dpng', '-r300')
end
