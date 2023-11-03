%Working from Attune Table, cut down to EDI ready 


function [] = generate_attune_table_EDI_temp_EP(Attunetablepath)

Attune = load([Attunetablepath filesep 'AttuneVolTable.mat']);
AttuneT = load([Attunetablepath filesep 'AttuneTable.mat']);

cruiseStr = Attunetablepath(50:end-17);
 if startsWith(cruiseStr, 'EN') 
     depth = 5; %5 meters
 elseif startsWith(cruiseStr, 'AR')
     depth = 2.1336; %7 ft
 else 
     keyboard
 end

%all we want is to remove low quality runs, 
% and columns QC_flag, Num_particles, scatter_hv


AttuneTable = Attune.AttuneVolTable; 
AttuneTable(AttuneTable.QC_flag == 0, :) = []; 

Attune2EDI = table;
Attune2EDI.cruise = repmat(cruiseStr,length(AttuneTable{:,1}),1);
Attune2EDI.datetime = datestr(AttuneTable.StartDate, 'yyyy-mm-ddTHH:MM:SS+00:00');
Attune2EDI.latitude = AttuneTable.lat;
Attune2EDI.longitude = AttuneTable.lon;
Attune2EDI.salinity = AttuneTable.salinity;
Attune2EDI.temperature = AttuneTable.temperature;
Attune2EDI.depth_m = depth*ones(size(AttuneTable(:,1)));
%Attune2EDI.Orgpicopro_cells_per_ml = AttuneTable.("Syn_count")./AttuneTable.VolAnalyzed_ml;

Attune2EDI.syn_cells_per_ml = AttuneTable.("Syn_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_2um_cells_per_ml = AttuneTable.("Euk_without_PE_leq2um_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_3um_cells_per_ml = AttuneTable.("Euk_without_PE_leq3um_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_5um_cells_per_ml = AttuneTable.("Euk_without_PE_leq5um_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_10um_cells_per_ml = AttuneTable.("Euk_without_PE_leq10um_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_20um_cells_per_ml = AttuneTable.("Euk_without_PE_leq20um_count")./AttuneTable.VolAnalyzed_ml;

Attune2EDI.syn_biovolume_concentration = AttuneTable.("Syn_biovolume")./AttuneTable.VolAnalyzed_ml; 
Attune2EDI.redeuk_leq_2um_biovolume_concentration = AttuneTable.("Euk_without_PE_leq2um_biovolume")./AttuneTable.VolAnalyzed_ml; 
Attune2EDI.redeuk_leq_3um_biovolume_concentration = AttuneTable.("Euk_without_PE_leq3um_biovolume")./AttuneTable.VolAnalyzed_ml;  
Attune2EDI.redeuk_leq_5um_biovolume_concentration = AttuneTable.("Euk_without_PE_leq5um_biovolume")./AttuneTable.VolAnalyzed_ml;  
Attune2EDI.redeuk_leq_10um_biovolume_concentration = AttuneTable.("Euk_without_PE_leq10um_biovolume")./AttuneTable.VolAnalyzed_ml; 
Attune2EDI.redeuk_leq_20um_biovolume_concentration = AttuneTable.("Euk_without_PE_leq20um_biovolume")./AttuneTable.VolAnalyzed_ml; 


Attune2EDI.syn_carbon_concentration = AttuneTable.("Syn_carbon")./AttuneTable.VolAnalyzed_ml./1000; %divide by 1000 to get micrograms per Liter
Attune2EDI.redeuk_leq_2um_carbon_concentration = AttuneTable.("Euk_without_PE_leq2um_carbon")./AttuneTable.VolAnalyzed_ml./1000; %divide by 1000 to get micrograms per Liter
Attune2EDI.redeuk_leq_3um_carbon_concentration = AttuneTable.("Euk_without_PE_leq3um_carbon")./AttuneTable.VolAnalyzed_ml./1000; 
Attune2EDI.redeuk_leq_5um_carbon_concentration = AttuneTable.("Euk_without_PE_leq5um_carbon")./AttuneTable.VolAnalyzed_ml./1000;   
Attune2EDI.redeuk_leq_10um_carbon_concentration = AttuneTable.("Euk_without_PE_leq10um_carbon")./AttuneTable.VolAnalyzed_ml./1000; 
Attune2EDI.redeuk_leq_20um_carbon_concentration = AttuneTable.("Euk_without_PE_leq20um_carbon")./AttuneTable.VolAnalyzed_ml./1000;  


Attune2EDI.volume_analyzed_ml = AttuneTable.VolAnalyzed_ml;
Attune2EDI.filename = AttuneTable.Filename;



if exist('assign_class_function') 
    table_metadata = {assign_class_function; classpath; string(datetime())};
else 
    table_metadata = AttuneT.table_metadata;
end
AttuneTableEDI = Attune2EDI; 
keyboard

save([Attunetablepath filesep 'AttuneTable_EDI.mat'],'AttuneTableEDI', 'table_metadata')
disp(['Result file saved:'])
disp([Attunetablepath '\AttuneTable_EDI'])

writetable(AttuneTableEDI, [Attunetablepath filesep cruiseStr '_EDI.csv'])

end

