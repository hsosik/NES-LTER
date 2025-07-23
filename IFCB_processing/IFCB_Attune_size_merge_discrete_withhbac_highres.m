%clist = {'EN608' 'EN617' 'EN627' 'EN644' 'EN649' 'EN655' 'EN657' 'EN661' 'EN668' 'AT46' 'AR66B' 'EN687' 'EN685' 'EN695' 'HRS2303' 'EN706'};
%clist = {'AT46'};  %need to check/redo EN695 HRS2303 EN706 %need Attune AR70B AR77 AR78 AR79
clist = {'EN706'};
attunebase = '\\sosiknas1\Lab_data\Attune\cruise_data\';
t = dir([attunebase '\\*\preserved*']);
t = t(cat(1,t.isdir));
clist_all = {t.folder}';

%SKIP SOME FOR NOW since output not available from attune proc
clist_all(contains(clist_all, 'AR28B')) = [];
clist_all(contains(clist_all, 'AR32')) = [];
clist_all(contains(clist_all, 'SR2018')) = [];
clist_all(contains(clist_all, 'AR66B')) = [];
clist_all(contains(clist_all, 'AR77')) = [];
clist_all(contains(clist_all, 'AE2426')) = [];

%SKIP SOME MORE That need different IFCB files
clist_all(contains(clist_all, 'SG2105')) = [];
clist_all(contains(clist_all, 'HB1907')) = [];
clist_all(contains(clist_all, 'TN368')) = [];
clist_all(contains(clist_all, 'EN688')) = [];
clist_all(contains(clist_all, 'MVCO')) = [];
clist_all(contains(clist_all, '20190303_AL')) = [];


cruiseStr_all = split(clist_all,'_');
cruiseStr_all = cruiseStr_all(:,end);
yr = split(clist_all,filesep);
yr = yr(:,end);
yr = extractBefore(yr,5);

ifcbbase = '\\sosiknas1\IFCB_products\NESLTER_transect\summary\';

%clist = {'AR29' 'RB1904' 'TN368'}; ifcbbase = '\\sosiknas1\IFCB_products\SPIROPA\summary\'; %'AR29' 'RB1904' 'TN368'
%clist = {'AR43' 'EN688'}; ifcbbase = '\\sosiknas1\IFCB_products\OTZ\summary\';
%just do this cruise
clist_all = clist_all(contains(clist_all, 'EN727')); yr = {'2025'}; cruiseStr_all = {'EN727'};
yrlist = unique(yr);

pout = '\\sosiknas1\Lab_data\Attune\cruise_data\IFCB_Attune_merge\summary_files_discrete\';
attune_euk_class = [1,5,6];
%
%Lee and Fuhrman 1987, on hbac: "Surprisingly, in six cultures with average per-cell biovolumes ranging from 0.036 to 0.073 μm3, 
% the average per-cell carbon biomass was relatively constant at 20 ± 0.08 fg of C (mean ± standard error of the mean)."
hba_cell_cubic_microns = mean([0.036 0.073]);
hba_cell_esd_microns = 2*(hba_cell_cubic_microns*3/4/pi).^(1/3);
hba_cell_carbon_picograms = 20/1000;
%
%%
for ylistn = 1:length(yrlist) %6 %just 2023 for now %3:length(yrlist)
    t = strcmp(yrlist(ylistn),yr);
    clist = clist_all(t);
    cruiseStr_now = cruiseStr_all(t);
    clear t
    %
    disp('loading IFCB year summary...')
    ifcbL = load([ifcbbase 'summary_biovol_allHDF_min20_' yrlist{ylistn} 'lists_group']);
    ifcb_yr = load([ifcbbase 'summary_biovol_allHDF_min20_' yrlist{ylistn}]);
%%
for clistn = 1:length(clist)
    cruiseStr = cruiseStr_now(clistn);
    disp(cruiseStr)
   load([clist{clistn} '\preserved\outputs\SummaryTable']);
   T = load([clist{clistn} '\preserved\outputs\EDI_table']);
   CNTable = renamevars(CNTable, {'Cast' 'Niskin'}, {'cast' 'niskin'});
   if strcmp('TN368', cruiseStr)
       CNTable(CNTable.cast==0,:) = [];
   end
   load([clist{clistn} '\preserved\outputs\EDI_table']); %load this to get the vol_analyzed
   TT = T.EDI_table(:,{'cast' 'niskin' 'depth_m' 'syn_volume_analyzed_ml' 'redeuk_volume_analyzed_ml' 'hetprok_volume_analyzed_ml'});
   CNTable = join(CNTable, TT);
   clear T TT 
   attunebase2 = [clist{clistn} '\preserved\outputs\class\'];
    %    AttuneTable = AttuneTable(AttuneTable.QC_flag==1,:); %omit the rows with bad QC flag

    attunebase_out = [attunebase2 'sizeDistPlots\'];
    if ~exist(attunebase_out, 'dir')
        mkdir(attunebase_out)
    end
   
    binedges = [0 2 3 5 7 10 15 20 30 40 50 100 inf]; saveStr = 'carbon2'; IFCBvar = 'ESD'; Attunevar = 'd'; hba_bin = discretize(hba_cell_esd_microns, binedges);
    %binedges = [0:1:100 inf]; saveStr = 'highres';
    %keep following option for "NBSS_vertint_FCM_IFCB..." to combine with bongo for LTER synthesis WG
    %binedges = [0 2.^(0:1:8)]; saveStr = 'base2'; IFCBvar = 'ESD'; Attunevar = 'd'; hba_bin = discretize(hba_cell_esd_microns, binedges);
   
%    binedges = [0 2.^(-5:1:18)]; saveStr = 'carbon'; IFCBvar = 'cellC'; Attunevar = 'carbon'; hba_bin = discretize(hba_cell_carbon_picograms, binedges);

    nbins = length(binedges)-1;
    cases2do = {'protist_tricho' 'Detritus'};
    %ifcbHmeta = ifcb_yr.meta_data(ifcb_uwind,:);
    attuneHmeta = CNTable(:, ["cruise" "cast" "niskin" "depth_m" "date_sampled" "latitude" "longitude" "nearest_station" "salinity" "potemp090c"]);
    Nrows = size(CNTable,1);
    blankH = nan(Nrows,nbins);
    for casei = 1:length(cases2do)
        ifcbHcount.(cases2do{casei}) = blankH;
    end
    ifcbHcarbon = ifcbHcount;
    attuneHcount.Syn = blankH;
    attuneHcount.Euk = blankH;
    attuneHcount.Hba = blankH;
    attuneHcarbon = attuneHcount;
    attuneml_syn = CNTable.syn_volume_analyzed_ml; %nan(Nrows,1);
    attuneml_euk = CNTable.redeuk_volume_analyzed_ml; %nan(Nrows,1);
    attuneml_hba = CNTable.hetprok_volume_analyzed_ml;  %nan(Nrows,1);
    %attuneml(:) = .32; %FIX LATER
    ifcbml = nan(size(attuneml_syn));
    for count = 1:Nrows
        disp([num2str(count) ' of ' num2str(Nrows)])
        ifcb_ind = find(startsWith(ifcb_yr.meta_data.cruise, cruiseStr) & ifcb_yr.meta_data.cast==num2str(CNTable.cast(count)) & ifcb_yr.meta_data.niskin==CNTable.niskin(count) & ~ifcb_yr.meta_data.skip);
        if ~isempty(ifcb_ind)
            ifcbml(count) = sum(ifcb_yr.meta_data.ml_analyzed(ifcb_ind));
            for casei = 1:length(cases2do)
                tempFea = ifcbL.groupFeaList.(cases2do{casei})(ifcb_ind);
                temp = array2table(cat(1,tempFea{:}), 'VariableNames', ifcbL.groupFeaList_variables);
                temp2 = discretize(temp.(IFCBvar),binedges);
                if ~isempty(temp2)
                    ifcbHcount.(cases2do{casei})(count,:) = histcounts(temp2, 1:nbins+1);
                    ifcbHcarbon.(cases2do{casei})(count,unique(temp2)) = splitapply(@sum , temp.cellC, findgroups(temp2))';
                end
            end
        end
        hSyn = nan(1,nbins);
        hSynC = hSyn;
        hEuk = hSyn;
        hEukC = hSyn;
        hHba = hSyn;
        hHbaC = hSyn;
        if ~isempty(CNTable.BacteriaFile{count})
            hHba(:) = 0;
            hHbaC(:) = 0;
            hHba(hba_bin) = CNTable.bac_per_ml(count)*CNTable.hetprok_volume_analyzed_ml(count); %counts
            hHbaC(hba_bin) = hHba(hba_bin)*hba_cell_carbon_picograms; % C in picograms
            attuneHcount.Hba(count,:) = hHba;
            attuneHcarbon.Hba(count,:) = hHbaC;
        end
        if ~isempty(CNTable.Synfile{count})
            f = char(regexprep(CNTable.Synfile(count), '.fcs', '.mat'));
            c = load([attunebase2 f]);
            if ~isnan(c.bead_value)
                v = real(c.volume(c.class==2)); %Syn
                carbon = biovol2carbon(v,0);
                x = real(biovol2esd(v));
                if strcmp(saveStr, 'carbon')
                    x = carbon;
                end
                temp2 = discretize(x,binedges);
                if ~isempty(temp2)
                    hSyn(1,:) = histcounts(temp2,1:nbins+1); %histcounts(d,binedges);
                    hSynC(1,unique(temp2)) = splitapply(@sum , carbon  , findgroups(temp2))';
                end
            else
                disp([f ' NaN bead_value'])
            end
            attuneHcount.Syn(count,:) = hSyn;
            attuneHcarbon.Syn(count,:) = hSynC;
        end
        if ~isempty(CNTable.Eukfile{count})
            f = char(regexprep(CNTable.Eukfile(count), '.fcs', '.mat'));
            c = load([attunebase2 f]);
            if ~isnan(c.bead_value)
                v = c.volume(ismember(c.class, attune_euk_class)); %Euk for discretes, including with hi/low PE
                carbon = biovol2carbon(v,0);
                x = real(biovol2esd(v));
                if strcmp(saveStr, 'carbon')
                    x = carbon;
                end
                temp2 = discretize(x,binedges);
                if ~isempty(temp2)
                    hEuk(1,:) = histcounts(temp2,1:nbins+1); %histcounts(d,binedges);
                    hEukC(1,unique(temp2)) = splitapply(@sum , carbon  , findgroups(temp2))';
                end
            else
                disp([f ' NaN bead_value'])
            end
            attuneHcount.Euk(count,:) = hEuk;
            attuneHcarbon.Euk(count,:) = hEukC;
        end

    end

    clear ind ifcb_uwind

    save([pout char(cruiseStr) '_IFCB_Attune_merge_sum_' saveStr], 'attuneH*', 'attuneml*', 'ifcbml', 'ifcbH*', 'binedges', 'attune_euk_class')
    disp('results saved:')
    disp([pout char(cruiseStr) '_IFCB_Attune_merge_sum_' saveStr])
end
%%
end
