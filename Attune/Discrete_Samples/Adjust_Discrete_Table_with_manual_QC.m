%Go get table of processed fcm files 
load('\\sosiknas1\Lab_data\Attune\cruise_data\Attune_Discrete_Table.mat')

ADT = Attune_Discrete_Table; 

j = find(isnan(ADT.syn_cells_per_ml) & isnan(ADT.redeuk_leq_20um_cells_per_ml) & isnan(ADT.hetprok_cells_per_ml)); 
ADT(j, :) = []; 

j = find(ADT.cruise == 'TN368' & isnan(ADT.depth_m))
ADT(j, :) = []; 



%load autoQC results
autoQC = readtable('\\sosiknas1\Lab_data\Attune\cruise_data\Compare_attune_discrete_files_to_samplelog.csv'); 


%load label lookup table
manualQC = readtable('\\sosiknas1\Lab_data\Attune\cruise_data\label_lookup_table_Attune_cruise.csv'); 


%
ADT.samplelog_QC = 2.*ones(height(ADT), 1); 

for i = 1:height(ADT)
    j = find(autoQC.cruisename == ADT.cruise(i) & autoQC.cast == ADT.cast(i) & autoQC.niskin == ADT.niskin(i));
    ADT.samplelog_QC(i) = autoQC.flag(j); %put in quality flag from sample log comparison results


    if autoQC.flag(j) == 3 % if it was identified as low quality, check if there is a recommended label change 
        k = find(manualQC.cruisename == ADT.cruise(i) & manualQC.cast == ADT.cast(i) & manualQC.niskin == ADT.niskin(i));
        
        if ~isempty(k) %there is a recommended change
        
        %find information to change to 
            l = find(ADT.cruise == manualQC.cruisename(k) & ADT.cast == manualQC.new_cast(k) & ADT.niskin == manualQC.new_niskin(k)); 
            
            if ~isempty(l) %there was already an entry for this niskin. Things get confusing
                disp('Existing rows need to be merged! be careful here. Line 37')
                keybaord
            else %just change cast and niskin values in existing row
                ADT.cast(i) = manualQC.new_cast(k); 
                ADT.niskin(i) = manualQC.new_niskin(k); 


                %also have to change environmental metadata 

                ADT.latitude(i) = NaN; 
                ADT.longitude(i) = NaN; 
                ADT.depth_m(i) = NaN;
                ADT.salinity(i) = NaN; 
                ADT.potential_temperature_c(i) = NaN; 
                ADT.date_sampled(i) = NaT; 

            end
                
    
        else %no recommended change
            ADT.samplelog_QC(i) = 4; %bad. Throw it out. 

        end
    end
end


%go back to get list of cruises we need envdata for 


shortstack = ADT(isnan(ADT.depth_m), :); 
cruiselist = unique(shortstack.cruise); 

for c = 1:length(cruiselist)

    if cruiselist(c) == 'EN644'
        restpath = 'https://nes-lter-data.whoi.edu/api/ctd/en644/'; 
    elseif cruiselist(c) == 'AR61B'
        restpath = 'https://nes-lter-data.whoi.edu/api/ctd/ar61b/'; 
    elseif cruiselist(c) == 'EN668'
        restpath = 'https://nes-lter-data.whoi.edu/api/ctd/en668/'; 
    elseif cruiselist(c) == 'TN368'
        restpath = '\\sosiknas1\Lab_data\Attune\cruise_data\20190705_TN368\preserved\tn368_bottle_data_Apr_2020_table.mat'; 
    else
        disp('need new case for cruise restpath')
        keyboard
    end


    if startsWith(restpath, 'https')
        bottledata = webread([restpath 'bottles.csv']); 
        metadata = webread([restpath 'metadata.csv']); 

    
        cruiseind = find(shortstack.cruise == cruiselist(c)); %row indeces that correspond to this cruise
        for i = cruiseind'
                %add cast data
                shortstack.latitude(i) = metadata.latitude(metadata.cast == shortstack.cast(i)); 
                shortstack.longitude(i) = metadata.longitude(metadata.cast == shortstack.cast(i)); 
                shortstack.date_sampled(i) = datetime(metadata.date{metadata.cast == shortstack.cast(i)}, 'InputFormat', 'yyyy-MM-dd HH:mm:SS+00:00');

                %now add info for this specific niskin 
                if sum(bottledata.cast == shortstack.cast(i) & bottledata.niskin == shortstack.niskin(i)) > 0
                    shortstack.depth_m(i) = bottledata.depsm(bottledata.cast == shortstack.cast(i) & bottledata.niskin == shortstack.niskin(i));
                    shortstack.salinity(i) = bottledata.sal00(bottledata.cast == shortstack.cast(i) & bottledata.niskin == shortstack.niskin(i)); 
                    shortstack.potential_temperature_c(i) = bottledata.potemp090c(bottledata.cast == shortstack.cast(i) & bottledata.niskin == shortstack.niskin(i));
                end
        end

    else 
        load(restpath)
        temp = importdata('\\sosiknas1\Lab_data\SPIROPA\20190705_TN368\fromOlga\tn368_niskin_pressure_depth.txt');    
        bottle_depth = temp.data;
        clear temp
        temp = table;
        temp.Cast = bottle_depth(:,1);
        temp.Niskin = bottle_depth(:,2);
        temp.depsm = bottle_depth(:,5);

        for i = 1:height(temp)
            temptemp = BTL(find(BTL.Cast == temp.Cast(i)), :);
            [~,ind] = min(abs(temptemp.TargetDepth_m - bottle_depth(i,3)));
            temp.sal00(i) = temptemp.Sal00(ind);
            temp.potemp090c(i) = temptemp.Potemp090C(ind);
            temp.Latitude_decimalDeg(i) = temptemp.Latitude_decimalDeg(ind);
            temp.Longitude_decimalDeg(i) = temptemp.Longitude_decimalDeg(ind);
            temp.datetime(i) = temptemp.datetime(ind);

        end

        BTL = temp;


        cruiseind = find(shortstack.cruise == cruiselist(c)); %row indeces that correspond to this cruise
        for i = cruiseind'
                %add cast data
                tempmeta = BTL(BTL.Cast == shortstack.cast(i), :); 
                shortstack.latitude(i) = tempmeta.Latitude_decimalDeg(1); 
                shortstack.longitude(i) = tempmeta.Longitude_decimalDeg(1); 
                shortstack.date_sampled(i) = tempmeta.datetime(1); 

                %now add info for this specific niskin 
                if sum(BTL.Cast == shortstack.cast(i) & BTL.Niskin == shortstack.niskin(i)) > 0
                    shortstack.depth_m(i) = BTL.depsm(BTL.Cast == shortstack.cast(i) & BTL.Niskin == shortstack.niskin(i));
                    shortstack.salinity(i) = BTL.sal00(BTL.Cast == shortstack.cast(i) & BTL.Niskin == shortstack.niskin(i)); 
                    shortstack.potential_temperature_c(i) = BTL.potemp090c(BTL.Cast == shortstack.cast(i) & BTL.Niskin == shortstack.niskin(i));
                end
        end

    end 

end


%now replace bad ADT rows with shortstack rows 
ADT(isnan(ADT.depth_m), :) = shortstack; 

ADT.salinity(ADT.salinity == 0) = NaN; 
ADT.depth_m(ADT.depth_m == 0) = NaN; 
ADT.potential_temperature_c(ADT.potential_temperature_c == 0 & isnan(ADT.salinity)) = NaN; 



ADT.hetprok_carbon_concentration(isnan(ADT.hetprok_cells_per_ml)) = NaN; %sometimes this carbon is still 0 


ADT = sortrows(ADT, [1 2 3]);

writetable(ADT, '\\sosiknas1\Lab_data\Attune\cruise_data\Attune_Discrete_Table_postmanualQC.csv');



writetable(ADT, edipath)








