 %Generates the attune table from Class files and FCSfileinfo.mat only
%slightly different from EDI table
%namely, we have counts of all particles, and piceouks_with_PE < 3

%Saves a table of stats about the fcs files as well as a variable 
%table_metadata which saves date generated and assign class function used

%classpath should end in \

function [] = generate_attune_table(classpath, FCSfileinfopath)

Attune = load([FCSfileinfopath]);

PE_euk_celltypes = [3 4];%input('which class numbers are PE containing Euks for this cruise?');

% identify bead run files and determine mean bead SSC-H for cruise
% may need to adjust epsilon and/or minpts depending on cruise

% read in and process Attune data files
%classlist = dir([classpath, '*.mat']); 
filelist = regexprep(Attune.FCSfileinfo.fcslist,'fcs', 'mat'); 
AttuneTable = table(Attune.FCSfileinfo.fcslist, datetime(Attune.FCSfileinfo.matdate_start, 'ConvertFrom', 'datenum'), datetime(Attune.FCSfileinfo.matdate_stop, 'ConvertFrom', 'datenum'), Attune.FCSfileinfo.vol_analyzed/1e6, 'VariableNames', {'Filename' 'StartDate' 'StopDate' 'VolAnalyzed_ml'});

%go ahead and add these to AttuneTable for all cruises. Their purpose is
%information for masking alternating settings. So after that only run
%masking on cruises that sampled with alternating settings.
AttuneTable.trigger1_parameter = Attune.FCSfileinfo.trigger1_parameter;
AttuneTable.trigger1_threshhold = Attune.FCSfileinfo.trigger1_threshhold;
AttuneTable.trigger2_parameter = Attune.FCSfileinfo.trigger2_parameter;
AttuneTable.trigger2_threshhold = Attune.FCSfileinfo.trigger2_threshhold;

% Creating the variables
%we want Syn, Euk<=2, Euk<=3, Euk<=5, Euk<=10, Euk<=20, PEeuk<=2, PEeuk<=3,
%PEeuk<5, PEeuk<10, PEeuk<20

EukSizes = [0 2 3 5 10 20];
%numBins = 1+2*(length(EukSizes)-1);
numBins = 1+2*(length(EukSizes))-1; %this is brittle, added a column for Pro
Count = NaN(length(filelist),numBins-1);
Biovol = Count;
Carbon = Count;
rem_ind = [];
Scatter_hv = NaN(length(filelist), 1);
% classes and per-event codes added with assign_class_v2 (Nov 2026):
% count, biovolume, carbon of single classes; Syn-type and flag counts
% (NaN for class files without syn_type / flag, i.e. older classifiers)
newclass = [3 4 5 10];
newclassnames = {'Crypto_dimX' 'Crypto_brightX' 'Euk_with_SynX' 'Euk_like_low_chlX'}; %10: low-chl Euk-like cells, identity uncertain (not in the Euk columns)
NewCount = NaN(length(filelist), numel(newclass));
NewBiovol = NewCount; NewCarbon = NewCount;
SynPart_count = NaN(length(filelist), 1); %class 8: Syn on particle, count only (not Syn size or carbon)
syntype_classes = {2 'Syn'; 5 'Euk_with_Syn'; 8 'Syn_on_particle'};
syntype_codes = {1 'typeA'; 2 'typeB'; 3 'typeC'; 4 'typeD'; 9 'type_unusable'}; %4 types since 26 Nov 2026 (syn_type_current)
SynTypeCount = NaN(length(filelist), size(syntype_classes, 1)*size(syntype_codes, 1));
flag_codes = {1 'Flag_lowchl_Syn_pigment'; 2 'Flag_below_Euk_chl_gap'; 3 'Flag_chl_edge_noise'; 4 'Flag_Pro_like'; ...
    5 'Flag_Euk_chl_size_cut'; 6 'Flag_split_Syn_merged_green_part'; 7 'Flag_split_Syn_unpaired_green_part'; ...
    8 'Flag_split_Syn_unpaired_blue_part'; 9 'Flag_split_other_green_part'};
%flags 6-9: laser dropout, a Syn recorded as two events. Pairs are merged into one event (flag 6 = the
%merged green part). Flag 7 (unpaired green part) is a separate Syn count, not in Syn_count (Syn total =
%Syn_count + Flag_split_Syn_unpaired_green_part_count); flags 6, 8 and 9 are not counted
Pro_near_trigger_fraction = NaN(length(filelist), 1); %class_meta.pro_near_trigger_frac (Pro-setting files)
ProMeasured = NaN(length(filelist), 1); %class_meta.pro_measured: the setting the file was classified with
FlagCount = NaN(length(filelist), size(flag_codes, 1));
Num_particles= NaN(length(filelist), 1); 

for count = 1:length(filelist) %go through each of the files in the FCSfileinfo
    if ~rem(count,10)
        disp([num2str(count) ' of ' num2str(length(filelist))])
    end
    filename = [classpath filelist{count}]; %load class file 
    if exist(filename)
    
    load(filename)
    if ~exist('volume_cubic_microns', 'var') %% if we haven't done size calibration yet, function should still run 
        eval('volume_cubic_microns = NaN.*class;'); 
    end
    
    Num_particles(count) = length(volume_cubic_microns); 

    %volume = real(volume); 
    volume = volume_cubic_microns;
    carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
    carbon = real(carbon); %having issues with formatting, keeps having valus with + 0i. 

    if exist('file_hv', 'var')
        Scatter_hv(count) = file_hv.SSC;
    end
    
    eval(['rename_class = class;']) %issues with class as a variable name since its a function in matlab 
    clear class 
    
    
    %SynFirst
        Count(count,1) = sum(rename_class==2);
        Biovol(count,1) = nansum(volume(rename_class==2));
        Carbon(count,1) = nansum(carbon(rename_class==2));

    

    %now Euks
    celltype = 1; 
    diam = (volume*3/4/pi).^(1/3)*2; %equivalent spherical diam, micrometers
    % pre-Aug 2019 no-OD2 cruises: cells beyond the correctable SSC range
    % (vol_flag, batch_mergedSSC2volume) carry the largest correctable volume
    % as a lower bound (~9.6 um ESD at 160 V); count them as >10 um for the
    % size bins (biovolume/carbon keep the lower-bound volume). Nov 2026
    if exist('vol_flag', 'var') && numel(vol_flag) == numel(diam)
        diam(vol_flag) = max(diam(vol_flag), 10.001);
    end
    for ii = 1:length(EukSizes)-1
        ind = find(diam<=EukSizes(ii+1) & rename_class==celltype);
        Count(count,ii+1) = size(ind,1);
        Biovol(count,ii+1) = nansum(volume(ind));
        Carbon(count,ii+1) = nansum(carbon(ind));
    end

    %now PE containing Euks
    celltype = PE_euk_celltypes;
    for ii = 2:length(EukSizes)-1
        ind = find(diam<=EukSizes(ii+1) & ismember(rename_class, celltype));
        Count(count,ii+length(EukSizes)-1) = size(ind,1);
        Biovol(count,ii+length(EukSizes)-1) = nansum(volume(ind));
        Carbon(count,ii+length(EukSizes)-1) = nansum(carbon(ind));
    end

    %Pro next (added April 2024)
        Count(count,11) = sum(rename_class==7);
        Biovol(count,11) = nansum(volume(rename_class==7));
        Carbon(count,11) = nansum(carbon(rename_class==7));

    %single classes split out of the groups above, and Syn on particle (Nov 2026)
    for ii = 1:numel(newclass)
        ind = rename_class==newclass(ii);
        NewCount(count,ii) = sum(ind);
        NewBiovol(count,ii) = nansum(volume(ind));
        NewCarbon(count,ii) = nansum(carbon(ind));
    end
    SynPart_count(count) = sum(rename_class==8);
    %syn_type and flag read by name: a bare load does not shadow functions
    %(flag is a colormap function), so read them into a struct
    ev = load(filename, 'syn_type', 'flag', 'class_meta');
    if isfield(ev, 'class_meta') && isfield(ev.class_meta, 'pro_near_trigger_frac')
        Pro_near_trigger_fraction(count) = ev.class_meta.pro_near_trigger_frac;
    end
    if isfield(ev, 'class_meta') && isfield(ev.class_meta, 'pro_measured')
        ProMeasured(count) = double(ev.class_meta.pro_measured);
    end
    if isfield(ev, 'syn_type') && numel(ev.syn_type) == numel(rename_class)
        if isfield(ev, 'class_meta'), ev.syn_type = syn_type_current(ev.syn_type, ev.class_meta); end
        col = 0;
        for ii = 1:size(syntype_classes, 1)
            for jj = 1:size(syntype_codes, 1)
                col = col + 1;
                SynTypeCount(count,col) = sum(rename_class==syntype_classes{ii,1} & ev.syn_type==syntype_codes{jj,1});
            end
        end
    end
    if isfield(ev, 'flag') && numel(ev.flag) == numel(rename_class)
        for ii = 1:size(flag_codes, 1)
            FlagCount(count,ii) = sum(ev.flag==flag_codes{ii,1});
        end
    end

    else
        rem_ind = [rem_ind count]; %indeces to remove from final table
        disp(['skipped ', filelist{count}])
    end
    
    clear volume volume_cubic_microns vol_flag syn_type ev
end



%define column names
classnames = {'SynX'};
for ii = 2:length(EukSizes)
    classnames = [classnames; ['Euk_without_PE_leq' num2str(EukSizes(ii)) 'umX']]; 
end
for ii = 3:length(EukSizes)
    classnames = [classnames; ['Euk_w_PE_leq' num2str(EukSizes(ii)) 'umX']] ;
end

classnames = [classnames; {'ProX'}];

AttuneTable = [AttuneTable array2table(Count, 'VariableNames', regexprep(classnames, 'X', '_count'))];
AttuneTable = [AttuneTable array2table(Biovol, 'VariableNames', regexprep(classnames, 'X', '_biovolume'))];
AttuneTable = [AttuneTable array2table(Carbon, 'VariableNames', regexprep(classnames, 'X', '_carbon'))];
% pro_measured: as the file was classified (class_meta, process_cruise), else from
% FCSfileinfo (older class files; older FCSfileinfo files lack it)
pm = ProMeasured;
if ismember('pro_measured', Attune.FCSfileinfo.Properties.VariableNames)
    fpm = double(Attune.FCSfileinfo.pro_measured);
    pm(isnan(pm)) = fpm(isnan(pm));
end
if any(~isnan(pm))
    % Pro not measured: Pro columns empty, and the lower bound is the Pro above the chl
    % trigger that v2 calls class 7 there (assign_class_v2, pro_in_euk_runs) plus the
    % remaining Pro-like events (flag 4); in Pro-setting files the lower bound is empty
    % (Nov 2026)
    nopro = pm == 0;
    pro_eukrun = AttuneTable.Pro_count; pro_eukrun(isnan(pro_eukrun)) = 0;
    AttuneTable.Pro_count(nopro) = NaN;
    AttuneTable.Pro_biovolume(nopro) = NaN;
    AttuneTable.Pro_carbon(nopro) = NaN;
    AttuneTable.Pro_lower_bound_count = FlagCount(:, 4) + pro_eukrun;
    AttuneTable.Pro_lower_bound_count(~nopro) = NaN;
    AttuneTable.Pro_measured = pm;
else
    AttuneTable.Pro_lower_bound_count = NaN(height(AttuneTable), 1);
    AttuneTable.Pro_measured = pm;
end
AttuneTable.Pro_near_trigger_fraction = Pro_near_trigger_fraction;

AttuneTable = [AttuneTable array2table(NewCount, 'VariableNames', regexprep(newclassnames, 'X', '_count'))];
AttuneTable = [AttuneTable array2table(NewBiovol, 'VariableNames', regexprep(newclassnames, 'X', '_biovolume'))];
AttuneTable = [AttuneTable array2table(NewCarbon, 'VariableNames', regexprep(newclassnames, 'X', '_carbon'))];
AttuneTable.Syn_on_particle_count = SynPart_count;
stnames = {};
for ii = 1:size(syntype_classes, 1)
    for jj = 1:size(syntype_codes, 1)
        stnames{end+1} = [syntype_classes{ii,2} '_' syntype_codes{jj,2} '_count'];
    end
end
AttuneTable = [AttuneTable array2table(SynTypeCount, 'VariableNames', stnames)];
AttuneTable = [AttuneTable array2table(FlagCount, 'VariableNames', strcat(flag_codes(:,2)', '_count'))];

AttuneTable.Num_particles = Num_particles; 
AttuneTable.QC_flag = Attune.FCSfileinfo.QC_flag;
AttuneTable.QC_flowrates = Attune.FCSfileinfo.QC_flowrates; 
AttuneTable.QC_dataintegrity = Attune.FCSfileinfo.QC_dataintegrity;
% final post-review flag (Nov 2026): 1 = marked exclude in the class review
% (bead_calibrated\review_decisions.csv). Every row is kept here; the EDI tables leave
% these out. After a later review, update with apply_review_excludes(cruise)
AttuneTable.Review_exclude = review_exclude_column(AttuneTable.Filename, fullfile(classpath, '..', 'review_decisions.csv'));

AttuneTable.Scatter_hv = Scatter_hv; 

AttuneTable(rem_ind,:) = []; 
AttuneTable = sortrows(AttuneTable, 'StartDate');

%notes(2) is where assign_class_function should be saved, sometimes saved
%as its own variable in class files
if exist('assign_class_function') 
    table_metadata = {assign_class_function; classpath; string(datetime()); 'PE_with_euk_classes:' num2str(PE_euk_celltypes)};
else 
    table_metadata = {notes(2); classpath; string(datetime())};
    if numel(notes) >= 5 %assign_class_v2 class files: class, syn_type and flag code descriptions
        table_metadata = [table_metadata; {notes(1); notes(4); notes(5); ...
            'Crypto_dim/Crypto_bright = classes 3/4 (together = Euk_w_PE); Euk_with_Syn = class 5; Syn_on_particle = class 8 (not in Syn_count); Euk_like_low_chl = class 10 (low-chl Euk-like cells, identity uncertain; not in the Euk columns); Syn type counts by class; flags are class 0 events'; ...
            'Pro_count/biovolume/carbon NaN in files run without Pro settings; Pro_lower_bound_count = class 7 (Pro above the chl trigger) plus Pro-like events (flag 4) there, NaN in Pro-setting files; Pro_near_trigger_fraction = share of Pro within 0.1 log of the SSC or chl trigger edge (Pro-setting files)'; ...
            'Review_exclude = 1: file marked exclude in the class review (review_decisions.csv); rows kept here, left out of the EDI tables'}];
    end
end

save([classpath '..\AttuneTable'],'AttuneTable', 'table_metadata')
disp(['Result file saved:'])
disp([classpath '..\AttuneTable'])

end

