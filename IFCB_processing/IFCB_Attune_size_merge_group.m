function IFCB_Attune_size_merge_group(clist)
%% 2022:'AT46' 2021:'EN661' 'EN668'  2020:'EN649' 'EN655' 'EN657' 2019:'EN627' 'EN644' 2018:'EN608' 'EN617'
%% 2019: AR34 AR39  2021: AR61B AR62 AR63
%SPIROPA 'TN368'; 'RB1904' 'AR29'
% OTZ 'AR43';
%clist = {'EN608' 'EN617' 'EN627' 'EN644' 'EN649' 'EN655' 'EN657' 'EN661' 'EN668' 'AT46' 'AR66B' 'EN687' 'EN685' 'EN695' 'HRS2303' 'EN706'};
%clist = {'AR78' 'AR79' 'EN712' 'AR82A' 'EN715' 'EN720' 'AE2426' 'EN727'};  %need to check/redo EN695 HRS2303 EN706 %need Attune AR70B AR77 AR78 AR79
%clist = {'AR99'};
if ischar(clist), clist = cellstr(clist); end

ifcbbase = '\\sosiknas1\IFCB_products\NESLTER_transect\summary\';
%clist = {'EN608' 'EN617' 'EN627' 'EN644' 'EN649' 'EN655' 'EN657' 'EN661' 'EN668' 'AT46'};
%clist = {'AR28' 'AR31' 'AR34' 'AR39' 'AR61B' 'AR62' 'AR63' 'AR66B'};

%clist = {'AR29' 'RB1904' 'TN368'}; ifcbbase = '\\sosiknas1\IFCB_products\SPIROPA\summary\'; %'AR29' 'RB1904' 'TN368'
%clist = {'AR43' 'EN688'}; ifcbbase = '\\sosiknas1\IFCB_products\OTZ\summary\';
attunebase = '\\sosiknas1\Lab_data\Attune\cruise_data\';
pout = '\\sosiknas1\Lab_data\Attune\cruise_data\IFCB_Attune_merge\summary_files\';

for clistn = 1:length(clist)

    cruiseStr = clist{clistn};
    disp(cruiseStr)
    alist = dir([attunebase '*' cruiseStr]);
    yr = alist.name(1:4);
    %load([attunebase alist.name '\bead_calibrated\AttuneTable_synmod'])
    if ismember(cruiseStr, {'EN608' 'EN617' 'AR28' 'AR29' 'AR31'})
        load([attunebase alist.name '\ifcb_calibrated\AttuneTable'])
        attunebase2 = [attunebase alist.name '\ifcb_calibrated\class\'];
    else
        load([attunebase alist.name '\bead_calibrated\AttuneTable'])
        attunebase2 = [attunebase alist.name '\bead_calibrated\class\'];
    end
    AttuneTable = AttuneTable(AttuneTable.QC_flag==1,:); %omit the rows with bad QC flag
    t_datenum = datenum(AttuneTable.StartDate);

    attunebase_out = [attunebase2 'sizeDistPlots\'];
    if ~exist(attunebase_out, 'dir')
        mkdir(attunebase_out)
    end
    %%
    ifcbL = load([ifcbbase 'summary_biovol_allHDF_min20_' yr 'lists_group']);
    ifcb_yr = load([ifcbbase 'summary_biovol_allHDF_min20_' yr]);
    %%

    %ifcb_uwind = find(ifcb_yr.meta_data.cruise==cruiseStr & ifcb_yr.meta_data.sample_type == 'underway' & ~ifcb_yr.meta_data.skip);
    % handle case for example AR39 (attune all AR39, IFCB AR39A, AR39B
    %ifcb_uwind = find(startsWith(ifcb_yr.meta_data.cruise,cruiseStr) & ifcb_yr.meta_data.sample_type == 'underway' & ~ifcb_yr.meta_data.skip);
    if strcmp(cruiseStr, 'HB1907') | strcmp(cruiseStr, 'AR43') | strcmp(cruiseStr, 'EN688') %OTZ
        ifcb_uwind = find(strcmp(ifcb_yr.meta_data.cruise, cruiseStr) & strcmp(ifcb_yr.meta_data.sample_type,'underway') & ~ifcb_yr.meta_data.skip); %OTZ
    else
        ifcb_uwind = find(startsWith(ifcb_yr.meta_data.cruise, cruiseStr) & strcmp(ifcb_yr.meta_data.sample_type,'underway') & ~ifcb_yr.meta_data.skip);
    end
    %if isequal(cruiseStr, 'AR62')
    %    dt = datetime(ifcb_yr.meta_data.sample_time, 'InputFormat', 'yyyy-MM-dd hh:mm:ss+00:00');
    %    ifcb_uwind = find(dt>datetime(2021,11,17) & dt<datetime(2021,11,24));
    %    disp('FIX THIS LATER--special case for missing metadata for AR62')
    %end
    if isequal(cruiseStr, 'AR29')
        disp('Using only IFCB127 on AR29')
        ifcb_uwind = find(startsWith(ifcb_yr.meta_data.cruise,cruiseStr) & ifcb_yr.meta_data.ifcb==127 & ifcb_yr.meta_data.sample_type == 'underway' & ~ifcb_yr.meta_data.skip);
    end
    if isequal(cruiseStr, 'EN661')
        disp('Using only IFCB102 on EN661 for now')
        ifcb_uwind = find(startsWith(ifcb_yr.meta_data.cruise,cruiseStr) & ifcb_yr.meta_data.ifcb==102 & ifcb_yr.meta_data.sample_type == 'underway' & ~ifcb_yr.meta_data.skip);
    end
    %%
    %   group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
    %     [~,ia,ib] = intersect(group_table.CNN_classlist, ifcb_yr.class2use);
    %     diatom_ind = ib(find(group_table.Diatom(ia)));
    %     diatom_ind = setdiff(diatom_ind, strmatch('Thalassiosira_TAG_external_detritus', ifcb_yr.class2use));
    %     dino_ind = ib(find(group_table.Dinoflagellate(ia)));
    %     nano_ind = ib(find(group_table.Nano(ia) | group_table.flagellate(ia) | group_table.Coccolithophore(ia)));
    %     otherPhyto_ind = [strmatch('Pseudochattonella_farcimen', ifcb_yr.class2use); strmatch('Phaeocystis', ifcb_yr.class2use, 'exact'); ; strmatch('Trichodesmium', ifcb_yr.class2use, 'exact') ];
    %     phyto_ind = [diatom_ind; dino_ind; nano_ind; otherPhyto_ind];
    %     ciliate_ind = ib(find(group_table.Ciliate(ia)));
    %     artifact_ind = ib(find(group_table.IFCBArtifact(ia)));
    %     particle_ind = setdiff(1:length(ifcbL.class2use),artifact_ind)'; %all except artifacts
    %
    %%
    Twin = 12/60/24; %12 minutes as days
    %diambins = logspace(-.5,log10(50),100);
    binedges = [0 2 3 5 7 10 15 20 30 40 50 100 inf];
    nbins = length(binedges)-1;
    %cases2do = {'diatom' 'dino' 'nano' 'otherPhyto' 'ciliate' 'particle'};
    cases2do = {'protist_tricho' 'Detritus' 'Diatom_noDetritus' 'Dinoflagellate' 'Ciliate' 'NanoFlagCocco'};
    ifcbHmeta = ifcb_yr.meta_data(ifcb_uwind,:);

    blankH = zeros(length(ifcb_uwind),nbins);
    for casei = 1:length(cases2do)
        ifcbHcount.(cases2do{casei}) = blankH;
    end
    ifcbHcarbon = ifcbHcount;
    ifcbHvol = ifcbHcount;
    attuneHcount.Syn = blankH;
    attuneHcount.Euk = blankH;
    attuneHcount.Pro = blankH;
    attuneHcarbon = attuneHcount;
    attuneHvol = attuneHcount;
    attuneHsa = attuneHcount;
    attuneml = nan(1,length(ifcb_uwind));
    attunemlPro = nan(1,length(ifcb_uwind));
    for count = 1:length(ifcb_uwind)
        if ~rem(count,20)
            disp([num2str(count) ' of ' num2str(length(ifcb_uwind))])
        end
        for casei = 1:length(cases2do)
            %            ind = eval([cases2do{casei} '_ind']);
            %tempFea = ifcbL.classFeaList{ifcb_uwind(count),ind};
            %temp = array2table(cat(1,tempFea{:}), 'VariableNames', ifcbL.classFeaList_variables);
            tempFea = ifcbL.groupFeaList.(cases2do{casei})(ifcb_uwind(count));
            % temp = array2table(cat(1,tempFea{:}), 'VariableNames', ifcbL.groupFeaList_variables);
            temp = array2table(tempFea{:}, 'VariableNames', ifcbL.groupFeaList_variables);
            temp2 = discretize(temp.ESD,binedges);
            if ~isempty(temp2)
                ifcbHcount.(cases2do{casei})(count,:) = histcounts(temp2, 1:nbins+1);
                ifcbHcarbon.(cases2do{casei})(count,unique(temp2)) = splitapply(@sum , temp.cellC, findgroups(temp2))';
                ifcbHvol.(cases2do{casei})(count,unique(temp2)) = splitapply(@sum , temp.summedBiovolume, findgroups(temp2))';
                ifcbHsa.(cases2do{casei})(count,unique(temp2)) = splitapply(@sum , temp.summedSurfaceArea, findgroups(temp2))';                
            end
        end
       % t = datenum(AttuneTable.StartDate);
        attune_ind = find(t_datenum>ifcb_yr.mdate(ifcb_uwind(count))-Twin & t_datenum<ifcb_yr.mdate(ifcb_uwind(count))+Twin);
        ml = 0;
        mlPro = 0;
        hSyn = zeros(length(attune_ind),nbins);
        hSynC = hSyn;
        hSynV = hSyn;
        hSynSA = hSyn;
        hEuk = hSyn;
        hEukC = hSyn;
        hEukV = hSyn;
        hEukSA = hSyn;
        hPro = hSyn;
        hProC = hSyn;
        hProV = hSyn;
        hProSA = hSyn;
        
        for count2 = 1:length(attune_ind)  %2:length(attune_ind) why was this starting at 2??
            f = char(regexprep(AttuneTable.Filename(attune_ind(count2)), '.fcs', '.mat'));
            if exist([attunebase2 f], 'file')
                c = load([attunebase2 f]);
                %if ~isnan(c.bead_value)
                if isfield(c, 'beadSSCmean')
                %if ~isnan(c.beadSSCmean.mean_SSCA_1micron)
                    v = real(c.volume_cubic_microns(c.class==2)); %Syn
                    carbon = biovol2carbon(v,0);
                    d = real(biovol2esd(v));
                    sa = pi*(d/2).^2;
                    temp2 = discretize(d,binedges);
                    if ~isempty(temp2)
                        hSyn(count2,:) = histcounts(temp2,1:nbins+1); %histcounts(d,binedges);
                        hSynC(count2,unique(temp2)) = splitapply(@sum , carbon  , findgroups(temp2))';
                        hSynV(count2,unique(temp2)) = splitapply(@sum , v  , findgroups(temp2))';
                        hSynSA(count2,unique(temp2)) = splitapply(@sum , sa  , findgroups(temp2))';
                    end
                    %v = real(c.rel_volume(c.class>=1&c.class<=4)); %leave out the large coincident class 5,6
                    v = c.volume_cubic_microns(ismember(c.class, [1,3,4])); %Euk
                    %v = c.volume_cubic_microns(ismember(c.class, [1,4])); %Euk
                    carbon = biovol2carbon(v,0);
                    d = real(biovol2esd(v));
                    sa = pi*(d/2).^2;
                    temp2 = discretize(d,binedges);
                    if ~isempty(temp2)
                        hEuk(count2,:) = histcounts(temp2,1:nbins+1); %histcounts(d,binedges);
                        hEukC(count2,unique(temp2)) = splitapply(@sum , carbon  , findgroups(temp2))';
                        hEukV(count2,unique(temp2)) = splitapply(@sum , v  , findgroups(temp2))';
                        hEukSA(count2,unique(temp2)) = splitapply(@sum , sa  , findgroups(temp2))';
                    end
                    if ~isnan(AttuneTable.Pro_count(attune_ind(count2)))
                    v = c.volume_cubic_microns(c.class==7); %Pro
                    d = real(biovol2esd(v));
                    sa = pi*(d/2).^2;
                    hPro(count2,1) = length(v);
                    hProC(count2,1) = hPro(count2,1)*50/1000; %round from Bertilsson et al. 2003 
                    hProV(count2,1) = sum(v);
                    hProSA(count2,1) = sum(sa);
                    mlPro = mlPro + AttuneTable.VolAnalyzed_ml(attune_ind(count2));
                    end
                    ml = ml + AttuneTable.VolAnalyzed_ml(attune_ind(count2));
                else
                    disp([AttuneTable.Filename(attune_ind(count2)) ' NaN bead_value'])
                end
            end
        end
        attuneml(1,count) = ml;
        attunemlPro(1,count) = mlPro;
        attuneHcount.Syn(count,:) = sum(hSyn,1);
        attuneHcarbon.Syn(count,:) = sum(hSynC,1);
        attuneHvol.Syn(count,:) = sum(hSynV,1);
        attuneHsa.Syn(count,:) = sum(hSynSA,1);
        attuneHcount.Euk(count,:) = sum(hEuk,1);
        attuneHcarbon.Euk(count,:) = sum(hEukC,1);
        attuneHvol.Euk(count,:) = sum(hEukV,1); 
        attuneHsa.Euk(count,:) = sum(hEukSA,1);
        attuneHcount.Pro(count,1) = sum(hPro(:,1),1,'omitmissing');
        attuneHcarbon.Pro(count,:) = sum(hProC,1,'omitmissing');
        attuneHvol.Pro(count,:) = sum(hProV,1,'omitmissing');
        attuneHsa.Pro(count,:) = sum(hProSA,1,'omitmissing');
    end

    clear ind ifcb_uwind attune_ind

    %class2use = ifcbL.class2use;

    save([pout char(alist.name) '_IFCB_Attune_merge_sum'], 'attuneH*', 'attuneml*', 'ifcbH*', 'binedges')
end
end
