
%% 

dirlist = {'\\sosiknas1\Lab_data\Attune\cruise_data\20180131_EN608\preserved'; 
    '\\sosiknas1\Lab_data\Attune\cruise_data\20190820_EN644\preserved';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20190921_AR38\preserved';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20191005_AR39B\preserved';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20200201_EN649\preserved';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20200725_EN655\preserved2';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20201013_EN657\preserved';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20210203_EN661\preserved';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20210716_EN668\preserved';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20211108_AR61B\preserved';
    '\\sosiknas1\Lab_data\Attune\cruise_data\20220216_AT46\preserved'; 
    '\\sosiknas1\Lab_data\Attune\cruise_data\20190705_TN368\preserved'}

Attune_Discrete_Table = table(); 
edipath = ['\\sosiknas1\Lab_data\Attune\EDI_data_packages\Attune_transect_FCMdiscrete\attune-transect-discrete-samples.csv'];


for c = 1:length(dirlist)

    basepath = dirlist{c};

fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'outputs' filesep];
classpath = [outpath 'class' filesep];
awspath = [basepath filesep 'aws\'];

%Load data for this cruise

G = load([outpath '\Gated_Table.mat']);
gated_table = G.gated_table; 
C = load([outpath '\SummaryTable.mat']); 
CNTable = C.CNTable; 


EDI_table = table(CNTable.cruise, CNTable.Cast, CNTable.Niskin, CNTable.latitude, CNTable.longitude, CNTable.depth_m, CNTable.salinity, CNTable.potemp090c); 
EDI_table.Properties.VariableNames = {'cruise'; 'cast'; 'niskin'; 'latitude'; 'longitude'; 'depth_m'; 'salinity'; 'potential_temperature_c';}; 

EDI_table.cruise = string(EDI_table.cruise); %helpful for merging tables when cruisenames are different lengtths 

%reformat dates so they don't suck 
if iscell(CNTable.date_sampled)
dates1 = cell2mat(CNTable.date_sampled); 
EDI_table.date_sampled = datetime(dates1(:, 1:19), 'Format', 'yyyy-MM-dd HH:mm:ss'); 
else
EDI_table.date_sampled = datetime(CNTable.date_sampled, 'Format', 'yyyy-MM-dd HH:mm:ss'); 
end

dates2 = cell2mat(CNTable.date_processed); 
EDI_table.date_processed = datetime(dates2, 'Format', 'yyyy-MM-dd'); 



%% Go back to class files using Gated_table

for i = 1:height(EDI_table)-1
    load([outpath filesep 'EDI_table.mat'], 'EDI_table')

end

for i = height(EDI_table); 
    %first check syn file
    if ~isempty(CNTable.Synfile{i})
    filename = CNTable.Synfile{i}; 
    cfilename = regexprep(filename, '.fcs', '.mat'); 

        if ~exist([[classpath filesep cfilename]])
        EDI_table.syn_cells_per_ml(i) = NaN;
        EDI_table.syn_biovolume_concentration(i) = NaN;
        EDI_table.syn_carbon_concentration(i) = NaN; %divide by 1000 to get micrograms per Liter
        EDI_table.syn_volume_analyzed_ml(i) = NaN; 
        EDI_table.syn_filename{i} = 'NaN';
        continue
        end


    C = load([classpath filesep cfilename]); 

    volume = real(C.volume); 
    carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
    carbon = real(carbon); %having issues with formatting, keeps having valus with + 0i. 


    %match to gated table row 
    gind = find(strcmp(gated_table.fcslist, filename)); 


    EDI_table.syn_cells_per_ml(i) = sum(C.class==2)./gated_table.Vol_analyzed_ml(gind); 
    EDI_table.syn_biovolume_concentration(i) = nansum(volume(C.class==2))./gated_table.Vol_analyzed_ml(gind); 
    EDI_table.syn_carbon_concentration(i) = nansum(carbon(C.class==2))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter
    EDI_table.syn_volume_analyzed_ml(i) = gated_table.Vol_analyzed_ml(gind); 
    EDI_table.syn_filename(i) = CNTable.Synfile(i); 

    else

        EDI_table.syn_cells_per_ml(i) = NaN;
        EDI_table.syn_biovolume_concentration(i) = NaN;
        EDI_table.syn_carbon_concentration(i) = NaN; %divide by 1000 to get micrograms per Liter
        EDI_table.syn_volume_analyzed_ml(i) = NaN; 
        EDI_table.syn_filename{i} = 'NaN'; 
    end

    %% done with Syn, move on to Euks
    if ~isempty(CNTable.Eukfile{i})

    filename = CNTable.Eukfile{i}; 
    cfilename = regexprep(filename, '.fcs', '.mat'); 


     if ~exist([[classpath filesep cfilename]])
        EDI_table.redeuk_leq_2um_cells_per_ml(i) = NaN;
        EDI_table.redeuk_leq_2um_biovolume_concentration(i) =  NaN;
        EDI_table.redeuk_leq_2um_carbon_concentration(i) =  NaN;
    % <= 3
        EDI_table.redeuk_leq_3um_cells_per_ml(i) =  NaN;
        EDI_table.redeuk_leq_3um_biovolume_concentration(i) =  NaN;
        EDI_table.redeuk_leq_3um_carbon_concentration(i) =  NaN;
    % <= 5
        EDI_table.redeuk_leq_5um_cells_per_ml(i) = NaN;
        EDI_table.redeuk_leq_5um_biovolume_concentration(i) =  NaN;
        EDI_table.redeuk_leq_5um_carbon_concentration(i) =  NaN;
        % <= 10
        EDI_table.redeuk_leq_10um_cells_per_ml(i) =  NaN;
        EDI_table.redeuk_leq_10um_biovolume_concentration(i) =  NaN;; 
        EDI_table.redeuk_leq_10um_carbon_concentration(i) =  NaN;
        % <= 20
        EDI_table.redeuk_leq_20um_cells_per_ml(i) =  NaN;
        EDI_table.redeuk_leq_20um_biovolume_concentration(i) =  NaN;
        EDI_table.redeuk_leq_20um_carbon_concentration(i) =  NaN;
    
        EDI_table.redeuk_volume_analyzed_ml(i) =  NaN;
        EDI_table.redeuk_filename{i} = 'NaN';
        continue
     end

    C = load([classpath filesep cfilename]); 

    volume = real(C.volume); 
    carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
    carbon = real(carbon); %having issues with formatting, keeps having valus with + 0i. 

    gind = find(strcmp(gated_table.fcslist, filename)); 


    EukSizes = [0 2 3 5 10 20];
    %size fractions by diameter
    diam = (volume*3/4/pi).^(1/3)*2; %equivalent spherical diam, micrometers
    
    %first < 2 
    bin_particle_ind = find(diam'<=2 & C.class==1)';

    EDI_table.redeuk_leq_2um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
    EDI_table.redeuk_leq_2um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind); 
    EDI_table.redeuk_leq_2um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter

    % <= 3
    bin_particle_ind = find(diam'<=3 & C.class==1)';
    EDI_table.redeuk_leq_3um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
    EDI_table.redeuk_leq_3um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind); 
    EDI_table.redeuk_leq_3um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter

    % <= 5
    bin_particle_ind = find(diam'<=5 & C.class==1)';
    EDI_table.redeuk_leq_5um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
    EDI_table.redeuk_leq_5um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind); 
    EDI_table.redeuk_leq_5um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter

    % <= 10
    bin_particle_ind = find(diam'<=10 & C.class==1)';
    EDI_table.redeuk_leq_10um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
    EDI_table.redeuk_leq_10um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind); 
    EDI_table.redeuk_leq_10um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter

    % <= 20
    bin_particle_ind = find(diam'<=20 & C.class==1)';
    EDI_table.redeuk_leq_20um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
    EDI_table.redeuk_leq_20um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind); 
    EDI_table.redeuk_leq_20um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter


    EDI_table.redeuk_volume_analyzed_ml(i) = gated_table.Vol_analyzed_ml(gind); 
    EDI_table.redeuk_filename(i) = CNTable.Eukfile(i); 

    else
        EDI_table.redeuk_leq_2um_cells_per_ml(i) = NaN;
        EDI_table.redeuk_leq_2um_biovolume_concentration(i) =  NaN;
        EDI_table.redeuk_leq_2um_carbon_concentration(i) =  NaN;
    % <= 3
        EDI_table.redeuk_leq_3um_cells_per_ml(i) =  NaN;
        EDI_table.redeuk_leq_3um_biovolume_concentration(i) =  NaN;
        EDI_table.redeuk_leq_3um_carbon_concentration(i) =  NaN;
    % <= 5
        EDI_table.redeuk_leq_5um_cells_per_ml(i) = NaN;
        EDI_table.redeuk_leq_5um_biovolume_concentration(i) =  NaN;
        EDI_table.redeuk_leq_5um_carbon_concentration(i) =  NaN;
        % <= 10
        EDI_table.redeuk_leq_10um_cells_per_ml(i) =  NaN;
        EDI_table.redeuk_leq_10um_biovolume_concentration(i) =  NaN;; 
        EDI_table.redeuk_leq_10um_carbon_concentration(i) =  NaN;
        % <= 20
        EDI_table.redeuk_leq_20um_cells_per_ml(i) =  NaN;
        EDI_table.redeuk_leq_20um_biovolume_concentration(i) =  NaN;
        EDI_table.redeuk_leq_20um_carbon_concentration(i) =  NaN;
    
        EDI_table.redeuk_volume_analyzed_ml(i) =  NaN;
        EDI_table.redeuk_filename{i} = 'NaN';


    end

    %% finally heterotrophic bacteria
    if ~isempty(CNTable.BacteriaFile{i})

    filename = CNTable.BacteriaFile{i}; 
    cfilename = regexprep(filename, '.fcs', '.mat'); 

    if ~exist([[classpath filesep cfilename]])
        EDI_table.hetprok_cells_per_ml(i) = NaN;
        EDI_table.hetprok_carbon_concentration(i) = NaN;
        EDI_table.hetprok_volume_analyzed_ml(i) = NaN;
        EDI_table.hetprok_filename{i} = 'NaN';
        continue
    end


    C = load([classpath filesep cfilename]); 

    volume = real(C.volume); 
    carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
    carbon = real(carbon); %having issues with formatting, keeps having valus with + 0i. 

    gind = find(strcmp(gated_table.fcslist, filename)); 

   
    EDI_table.hetprok_cells_per_ml(i) = CNTable.bac_per_ml(i); % Do NOT go back to class file, since we had to account for time gate 
    %EDI_table.hetprok_biovolume_concentration(i) = nansum(volume(C.class==3))./gated_table.Vol_analyzed_ml(gind); 
    %EDI_table.hetprok_carbon_concentration(i) = nansum(carbon(C.class==3))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter
    
    
    %doing carbon concentration per cell conversion based on Lee & Furhman 1987 
     EDI_table.hetprok_carbon_concentration(i) = EDI_table.hetprok_cells_per_ml(i).* 20  * 1e-6; %to convert to micrograms per liter

    EDI_table.hetprok_volume_analyzed_ml(i) = gated_table.Vol_analyzed_ml(gind); 
    EDI_table.hetprok_filename(i) = CNTable.BacteriaFile(i); 

    else
        EDI_table.hetprok_cells_per_ml(i) = NaN;
        EDI_table.hetprok_carbon_concentration(i) = NaN; 
        EDI_table.hetprok_volume_analyzed_ml(i) = NaN;
        EDI_table.hetprok_filename{i} = 'NaN';
    end

end
save([outpath filesep 'EDI_table.mat'], 'EDI_table')
disp([outpath filesep 'EDI_table.mat'])

%add cruise to merged table
Attune_Discrete_Table = [Attune_Discrete_Table; EDI_table];

clear EDI_table
end

save('\\sosiknas1\Lab_data\Attune\cruise_data\Attune_Discrete_Table.mat', 'Attune_Discrete_Table')

Attune_Discrete_Table.latitude = round(Attune_Discrete_Table.latitude, 4); 
Attune_Discrete_Table.longitude = round(Attune_Discrete_Table.longitude, 4); 

for j = 6:width(Attune_Discrete_Table)
    if isnumeric(Attune_Discrete_Table{:,j})
        Attune_Discrete_Table{:,j} = round(Attune_Discrete_Table{:,j}, 2); 
    end
end

writetable(Attune_Discrete_Table, edipath)







