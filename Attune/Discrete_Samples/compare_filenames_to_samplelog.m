
%first load samplelog from here

samplelog = readtable('\\sosiknas1\Lab_data\LTER\LTER_sample_log.xlsx');

samplelog.Properties.VariableNames{1} = 'Cruise';
samplelog.Properties.VariableNames{7} = 'Cast';
samplelog.Properties.VariableNames{8} = 'Niskin';
samplelog.Properties.VariableNames{36} = 'FCMcol'; 

samplelog = samplelog(:, [1 7 8 36]); 
%%

%add spiropa sample log
spiropa_samplelog = readtable('\\sosiknas1\Lab_data\SPIROPA\SPIROPA_Sosik_sample_log.xlsx'); 

spiropa_samplelog = spiropa_samplelog(:, [1 2 7 14]); 
spiropa_samplelog.Properties.VariableNames{4} = 'FCMcol'; 

samplelog = [samplelog; spiropa_samplelog]; 

%%

%make variable for whether FCM sample should exist
FCM_logical = logical(strcmp(samplelog.FCMcol, 'D') | strcmp(samplelog.FCMcol, 'F') | strcmp(samplelog.FCMcol, 'FR')  | strcmp(samplelog.FCMcol, 'Y') | strcmp(samplelog.FCMcol, 'y')) ; 

%cut down sample log only to niskins that have FCM samples
samplelog_small = table(string(samplelog.Cruise), samplelog.Cast, samplelog.Niskin); 

samplelog_small.Properties.VariableNames{1} = 'Cruise';
samplelog_small.Properties.VariableNames{2} = 'Cast';
samplelog_small.Properties.VariableNames{3} = 'Niskin';

samplelog_small = samplelog_small(FCM_logical, :);


%%
%Now get table of processed fcm files 
load('\\sosiknas1\Lab_data\Attune\cruise_data\Attune_Discrete_Table.mat')

ADT = Attune_Discrete_Table; 

%% get total number of cruises to consider

cruiselist = unique([ADT.cruise; samplelog_small.Cruise]); 


%% Finally we want to generate the output table 

%Cruise %Cast %Niskin %FCM sample taken? %File processed on attune?
%%filename

Comparison = table(); 

for c = 1:length(cruiselist)
    %cut down each table to this cruise of interest
    temp_sl = samplelog_small(samplelog_small.Cruise == cruiselist(c), :); 
    temp_att = ADT(ADT.cruise == cruiselist(c), :); 

    %get groups of unique cast and niskin combos
    [G, cast, niskin] = findgroups([temp_sl.Cast; temp_att.cast], [temp_sl.Niskin; temp_att.niskin]);

    %make appropriately sized output columns
    cruisename = repmat(cruiselist(c), max(G), 1); 
    
    sample_taken = zeros(size(cruisename)); 
    processed_on_attune = zeros(size(cruisename)); 
    attune_filenames = strings(length(sample_taken), 3); 
    for g = 1:max(G)
        if cast(g) == 0 %underways aren't in this sample log
            sample_taken(g) = 2; %unevaluated
        elseif sum(temp_sl.Cast == cast(g) & temp_sl.Niskin == niskin(g)) == 1 %did we take a sample from this niskin according to log?
            sample_taken(g) = 1; 
        end
        if sum(temp_att.cast == cast(g) & temp_att.niskin == niskin(g)) == 1 %do we have attune file matching to this niskin? 
            processed_on_attune(g) = 1; 
            ind = find(temp_att.cast == cast(g) & temp_att.niskin == niskin(g)); %go get filenmaes in case that's useful
            attune_filenames(g, 1) = [temp_att.syn_filename{ind}];
            attune_filenames(g, 2) = [temp_att.redeuk_filename{ind}];
            attune_filenames(g, 3) = [temp_att.hetprok_filename{ind}]; 
        end

    end

    cruise_table = table(cruisename, cast, niskin, sample_taken, processed_on_attune, attune_filenames); 

    Comparison = [Comparison; cruise_table]; %combine this cruise table with others

end


Comparison.flag = ones(size(Comparison.sample_taken)).*2; 
Comparison.flag(Comparison.sample_taken == 1 & Comparison.processed_on_attune == 1) = 1; 
Comparison.flag(Comparison.sample_taken == 0 & Comparison.processed_on_attune == 1) = 3; 



writetable(Comparison, '\\sosiknas1\Lab_data\Attune\cruise_data\Compare_attune_discrete_files_to_samplelog.csv')













