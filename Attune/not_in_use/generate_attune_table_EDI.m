function [AttuneTableEDI] = generate_attune_table_EDI(Attunetablepath)

a = dir([Attunetablepath]); 
a = a(1).folder; 
%Attune = load([a filesep 'AttuneVolTable.mat']);
Attune = load([a filesep 'AttuneTable.mat']);

test = split(a, '\');
test = test{7}; 
test = split(test, '_');
cruiseStr = test{2};
 if startsWith(cruiseStr, 'EN') 
     depth = 5; %5 meters
 elseif startsWith(cruiseStr, 'AR')
     depth = 2.1336; %7 ft
 elseif startsWith(cruiseStr, 'AT')
      depth = 5;
 else
     depth = NaN; 
 end

%all we want is to remove low quality runs, 
% and columns QC_flag, Num_particles, scatter_hv


AttuneTable = Attune.AttuneTable; 
if contains(Attunetablepath, '20180414_AR29')
    AttuneTable.QC_flag(AttuneTable.QC_flag==0 & AttuneTable.QC_flowrates(:,2) < 5) = 2;
else
    AttuneTable(AttuneTable.QC_flag == 0, :) = []; 
end 
Attune2EDI = table;
Attune2EDI.cruise = string(repmat(cruiseStr,length(AttuneTable{:,1}),1));
Attune2EDI.date_time_utc = datestr(AttuneTable.StartDate, 'yyyy-mm-dd HH:MM:SS');
if ismember('lat', AttuneTable.Properties.VariableNames)
    Attune2EDI.latitude = AttuneTable.lat;
    Attune2EDI.longitude = AttuneTable.lon;
    if AttuneTable.lon > 180
        Attune2EDI.longitude = AttuneTable.lon - 360; 
    end
end


Attune2EDI.depth_m = depth*ones(size(AttuneTable(:,1)));

Attune2EDI.syn_cells_per_ml = AttuneTable.("Syn_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_2um_cells_per_ml = AttuneTable.("Euk_without_PE_leq2um_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_3um_cells_per_ml = AttuneTable.("Euk_without_PE_leq3um_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_5um_cells_per_ml = AttuneTable.("Euk_without_PE_leq5um_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_10um_cells_per_ml = AttuneTable.("Euk_without_PE_leq10um_count")./AttuneTable.VolAnalyzed_ml;
Attune2EDI.redeuk_leq_20um_cells_per_ml = AttuneTable.("Euk_without_PE_leq20um_count")./AttuneTable.VolAnalyzed_ml;

Attune2EDI.syn_biovolume = AttuneTable.("Syn_biovolume")./AttuneTable.VolAnalyzed_ml; 
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
if contains(Attunetablepath, '20180414_AR29')
    Attune2EDI.QC_flag = AttuneTable.QC_flag;
end

% if exist('assign_class_function') 
%     table_metadata = {assign_class_function; classpath; string(datetime())};
% else 
    table_metadata = Attune.table_metadata;
%end
AttuneTableEDI = Attune2EDI; 

save([a filesep 'AttuneTable_EDI.mat'],'AttuneTableEDI', 'table_metadata')

disp(['Result file saved:'])
disp([a '\AttuneTable_EDI'])

writetable(AttuneTableEDI, [a filesep cruiseStr '_EDI.csv'])
disp(['Result file saved:'])
disp([a filesep cruiseStr '_EDI.csv'])

end
