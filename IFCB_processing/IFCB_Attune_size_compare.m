%%2021:'EN661' 'EN668'  2020:'EN649' 'EN655' 'EN657' 2019:'EN627' 'EN644' 2018:'EN608' 'EN617'
%SPIROPA 'TN368'; 'RB1904' 'AR29'
cruiseStr = 'AR39';
%ifcbbase = '\\sosiknas1\IFCB_products\NESLTER_transect\summary\';
ifcbbase = 'C:\work\IFCB_products\NESLTER_transect\summary\';
%ifcbbase = 'C:\work\IFCB_products\SPIROPA\summary\';
attunebase = '\\sosiknas1\Lab_data\Attune\cruise_data\';
alist = dir([attunebase '*' cruiseStr]);
yr = alist.name(1:4);
%load([attunebase alist.name '\bead_calibrated\AttuneTable_synmod'])
load([attunebase alist.name '\bead_calibrated\AttuneTable'])
%attunebase2 = [attunebase alist.name '\bead_calibrated\class_newVolEst\'];
%attunebase_out = [attunebase2 'sizeDistPlots_newVolEst\'];
attunebase2 = [attunebase alist.name '\bead_calibrated\class\'];
attunebase_out = [attunebase2 'sizeDistPlots\'];
if ~exist(attunebase_out, 'dir')
    mkdir(attunebase_out)
end
%%
ifcbL = load([ifcbbase 'summary_biovol_allHDF_min20_' yr 'lists']);
ifcb_yr = load([ifcbbase 'summary_biovol_allHDF_min20_' yr]);

%%

%ifcb_uwind = find(ifcb_yr.meta_data.cruise==cruiseStr & ifcb_yr.meta_data.sample_type == 'underway' & ~ifcb_yr.meta_data.skip);
% handle case for example AR39 (attune all AR39, IFCB AR39A, AR39B
ifcb_uwind = find(startsWith(ifcb_yr.meta_data.cruise,cruiseStr) & ifcb_yr.meta_data.sample_type == 'underway' & ~ifcb_yr.meta_data.skip);

%%
group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
[~,ia,ib] = intersect(group_table.CNN_classlist, ifcb_yr.class2use);
diatom_ind = ib(find(group_table.Diatom(ia)));
dino_ind = ib(find(group_table.Dinoflagellate(ia)));
nano_ind = ib(find(group_table.Nano(ia) | group_table.flagellate(ia) | group_table.Coccolithophore(ia)));
otherPhyto_ind = [strmatch('Pseudochattonella_farcimen', ifcb_yr.class2use); strmatch('Phaeocystis', ifcb_yr.class2use, 'exact')];
phyto_ind = [diatom_ind; dino_ind; nano_ind; otherPhyto_ind];
ciliate_ind = ib(find(group_table.Ciliate(ia)));
artifact_ind = ib(find(group_table.IFCBArtifact(ia)));
particle_ind = setdiff(1:length(ifcbL.class2use),artifact_ind); %all except artifacts
%%
Twin = 12/60/24; %12 minutes as days
diambins = logspace(-.5,log10(50),100);
for count = 467:length(ifcb_uwind)
    temp = array2table(cat(1,ifcbL.classFeaList{ifcb_uwind(count),[phyto_ind; ciliate_ind]}), 'VariableNames', ifcbL.classFeaList_variables);
    %temp = array2table(cat(1,ifcbL.classFeaList{ifcb_uwind(count),:}), 'VariableNames', ifcbL.classFeaList_variables);
    ifcbH = histcounts(temp.ESD,[diambins inf]);
    figure(1), clf
    semilogy(diambins,ifcbH./ifcb_yr.meta_data.ml_analyzed(ifcb_uwind(count)), 'b.-')
    hold on
    temp = array2table(cat(1,ifcbL.classFeaList{ifcb_uwind(count),particle_ind}), 'VariableNames', ifcbL.classFeaList_variables);
    ifcbH = histcounts(temp.ESD,[diambins inf]);
    semilogy(diambins,ifcbH./ifcb_yr.meta_data.ml_analyzed(ifcb_uwind(count)), 'b-.')
    t = datenum(AttuneTable.StartDate);
    attune_ind = find(t>ifcb_yr.mdate(ifcb_uwind(count))-Twin & t<ifcb_yr.mdate(ifcb_uwind(count))+Twin);  
    attuneH = zeros(size(ifcbH));
    attuneml = 0;
    h = NaN(length(attune_ind),length(diambins));
    h2 = h;
    for count2 = 1:length(attune_ind)
        f = char(regexprep(AttuneTable.Filename(attune_ind(count2)), '.fcs', '.mat'));
        c = load([attunebase2 f]);
        v = real(c.volume(c.class>=1&c.class<=4)); %leave out the large coincident class 5,6
        %v = real(c.rel_volume(c.class>=1&c.class<=4)); %leave out the large coincident class 5,6
        d = real(biovol2esd(v));
        h(count2,:) = histcounts(d,[diambins inf]);
        attuneml = attuneml + AttuneTable.VolAnalyzed_ml(attune_ind(count2));
    end
    attuneH = sum(h,1);
    semilogy(diambins,attuneH./attuneml, 'g.-')
    if ~isempty(h)
        [~,ii] = min(abs(t(attune_ind)-ifcb_yr.mdate(ifcb_uwind(count))));
        semilogy(diambins,h(ii,:)./.4, 'c-.') %closest Attune sample to IFCB sample
        %semilogy(diambins,max(h)./.4, 'g-.')
        %semilogy(diambins,min(h)./.4, 'g-.')
        %semilogy(diambins,sum(h2,1)./attuneml, 'r.-', 'linewidth', 1)
    end
    title(ifcb_yr.filelist(ifcb_uwind(count)), 'interpreter', 'none')
    legend(['IFCB ' num2str(round(ifcb_yr.meta_data.ml_analyzed(ifcb_uwind(count)),2)) ' ml'], 'IFCB w/detritus', ['Attune +/-12min ' num2str(attuneml) ' ml'], 'Attune nearest 0.4 ml')
    ylabel('Cells ml^{-1} per bin')
    xlabel('ESD (\mum)')
    grid
    set(gca, 'xtick',0:5:50, 'xlim',[0 50])
   print([attunebase_out ifcb_yr.filelist{ifcb_uwind(count)} '.png'], '-dpng')
end