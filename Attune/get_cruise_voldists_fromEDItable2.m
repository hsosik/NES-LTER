%% Get Euk and Syn Volume Data into table 
%% add it to uw table and also standardize columns

%let's go through class files and store Ndist values for syn and euks
% not every volume, but the histogram. Have to make sure we keep track 
% of what volume bins we use. 

%we also need to choose the underway variables that we are going to use for
%our analysis. We will try to update this according to Stace's
%recommendations. Differ between cruises. 

%basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20190820_EN644';

%adding mean fluorescence

function get_cruise_voldists_fromEDItable2(basepath)

euk_volbins = [0 2.^[-5:1/5:8]]; % for euks 
syn_volbins = [0 2.^[-5:0.125:2]]; 

fname = dir([basepath '\bead_calibrated\Attune*uw_match.mat']);

T = load([fname.folder filesep fname.name]);
fname = fname.name;

% First find standard column names and append to underway table. 

    if isfield(T, 'Attune_uw_match')
         T = T.Attune_uw_match; 
    else
        T = T.AttuneTable; 
    end
    
  % Get cruise name from first filename
    cruisename =  split(T.Filename{1}, '_'); 
    cruisename = cruisename{2}; 

    %First case where attune tables are old, set quality flag to 0 if bad, 1 if good 
    Exist_Column = strcmp('QC_flowrate_median',T.Properties.VariableNames); 
    if Exist_Column(Exist_Column==1) 
       T.QC_flag = logical((T.QC_flowrate_median<2 & T.QC_flowrate_std<1.5)); 
    else
        T.QC_flag = logical(T.QC_flag); 
    end
    %now remove low quality files!
    T = T(T.QC_flag == 1, :);
    
    %Find columns with class counts
    Euk_count = T.Euk_without_PE_leq20um_count; 
    T.Euk_count = Euk_count; 
    


    %% OK Now is big env. data parsing 

    %Latitude & Longitude 
    if startsWith(cruisename, 'EN') 
        Exist_Column = strcmp('gps_furuno_latitude', T.Properties.VariableNames); %use gps-furuno if quality flag is not 0
        Exist_Column2 = strcmp('gps_garmin741_latitude', T.Properties.VariableNames); %use gps-furuno if quality flag is not 0

        if Exist_Column(Exist_Column ==1)
            lat = T.gps_furuno_latitude;
            lon = T.gps_furuno_longitude;
            lat(T.gps_furuno_quality == 0) = NaN;
            lon(T.gps_furuno_quality == 0) = NaN;
        elseif Exist_Column2(Exist_Column2 ==1) %super early cruises dont have gps-furuno
            lat = T.gps_garmin741_latitude;
            lon = T.gps_garmin741_longitude;
            lat(T.gps_garmin741_quality == 0) = NaN;
            lon(T.gps_garmin741_quality == 0) = NaN;
        end
       
    elseif startsWith(cruisename, 'AR')
        Exist_Column4 = strcmp('dec_lat',T.Properties.VariableNames); 
        Exist_Column5 = strcmp('Dec_LAT', T.Properties.VariableNames); 
        if Exist_Column4(Exist_Column4==1) 
            lat = T.dec_lat; 
            lon = T.dec_lon; 
        elseif Exist_Column5(Exist_Column5==1) 
            lat = T.Dec_LAT; 
            lon = T.Dec_LON;
        end
    
     elseif startsWith(cruisename, 'HRS')
            lat = T.latitude_deg; 
            lon = T.longitude_deg; 

    else %other cruises
        Exist_Column1 = strcmp('gps_furuno_latitude', T.Properties.VariableNames); %use gps-furuno if quality flag is not 0
        Exist_Column2 = strcmp('dec_lat',T.Properties.VariableNames); 
        Exist_Column3 = strcmp('Dec_LAT', T.Properties.VariableNames); 
        Exist_Column4 = strcmp('lat_flr', T.Properties.VariableNames); 
        Exist_Column5 = strcmp('lat_SAMOS', T.Properties.VariableNames);
        Exist_Column6 = strcmp('lat_tsg',T.Properties.VariableNames); 
        Exist_Column7 = strcmp('lat',T.Properties.VariableNames); 
        if Exist_Column1(Exist_Column1==1) 
           lat = T.gps_furuno_latitude;
           lon = T.gps_furuno_longitude;
        elseif Exist_Column2(Exist_Column2==1) 
            lat = T.dec_lat; 
            lon = T.dec_lon; 
        elseif Exist_Column3(Exist_Column3==1) 
            lat = T.Dec_LAT; 
            lon = T.Dec_LON;
        elseif Exist_Column4(Exist_Column4==1) 
            lat = T.lat_flr; 
            lon = T.lon_flr; 
        elseif Exist_Column5(Exist_Column5==1) 
            lat = T.lat_SAMOS; 
            lon = T.lon_SAMOS;
        elseif Exist_Column6(Exist_Column6==1) 
            lat = T.lat_tsg; 
            lon = T.lon_tsg; 
        elseif Exist_Column7(Exist_Column7==1) 
            lat = T.lat; 
            lon = T.lon; 
        end
    end

    T.lat = lat; 
    T.lon = lon; 
   

    %Temperature & Salinity 
    
     if startsWith(cruisename, 'EN')

         if strcmp(cruisename, 'EN627') %SUPER ODDBALL We use tsg2_temperature
              temperature = T.tsg2_temperature; 
         else
            temperature = T.tsg1_sst;
         end

         if sum(strcmp(cruisename, {'EN627'; 'EN644'; 'EN668'})) 
            salinity = T.tsg2_salinity; 
         else
             salinity = T.tsg1_salinity; 
         end

     elseif startsWith(cruisename, 'AR')
          Exist_Column5 = strcmp('sbe45s', T.Properties.VariableNames);
          Exist_Column6 = strcmp('SBE45S', T.Properties.VariableNames);

          if Exist_Column5(Exist_Column5==1)
                salinity = T.sbe45s; 
                temperature = T.sbe48t; 
          elseif Exist_Column6(Exist_Column6==1)
                salinity = T.SBE45S; 
                temperature = T.SBE48T; 
          end
      elseif startsWith(cruisename, 'HRS')
            salinity = T.salinity_psu; 
            temperature = T.water_temperature_degree_c; 
     else %other cruise vessels 
        Exist_Column = strcmp('sbe45S',T.Properties.VariableNames);  
        Exist_Column2 = strcmp('tsg1_salinity', T.Properties.VariableNames);  %sometimes we want tsg2, depends on cruise
        Exist_Column3 = strcmp('salinity', T.Properties.VariableNames);
        Exist_Column4 = strcmp('s', T.Properties.VariableNames);
        Exist_Column5 = strcmp('sbe45s', T.Properties.VariableNames);
        Exist_Column6 = strcmp('SBE45S', T.Properties.VariableNames);
        if Exist_Column(Exist_Column==1) 
            salinity = T.sbe45S; 
            temperature = T.sbe45T; 
        elseif Exist_Column2(Exist_Column2==1)
            salinity = T.tsg1_salinity; 
            temperature = T.tsg1_sst; %this is good. We don't want tsg1_temperature. Even if we use tsg2 for salinity, we want tsg1_sst for temp. 
        elseif Exist_Column3(Exist_Column3==1)
            salinity = T.salinity; 
            temperature = T.temperature; 
        elseif Exist_Column4(Exist_Column4==1)
            salinity = T.s; 
            temperature = T.t1; 
        elseif Exist_Column5(Exist_Column5==1)
            salinity = T.sbe45s; 
            temperature = T.sbe48t; 
        elseif Exist_Column6(Exist_Column6==1)
            salinity = T.SBE45S; 
            temperature = T.SBE48T; 
        end
     end
        T.salinity = salinity; 
        T.temperature = temperature; 
        

        %Sunlight data

        rad_sw = nan(height(T),1) ;   
        Exist_Column = strcmp('rad1_sw', T.Properties.VariableNames); %we also have a rad2_sw

        if Exist_Column(Exist_Column==1)
            figure
            scatter(T.StartDate, T.rad1_sw, '.k')
            hold on 
            scatter(T.StartDate, T.rad2_sw, '.r')
            legend({'rad1_sw', 'rad2_sw'})
            datetick
            ylabel('Radiation')

            pause
            whichtouse = input('Which to use? 1 or 2'); 

            if whichtouse == 1
                rad_sw = T.rad1_sw; 
            elseif whichtouse ==2
                rad_sw = T.rad2_sw; 
            end

        else %no 1 and 2 to choose from 
        Exist_Column2 = strcmp('rad_sw',T.Properties.VariableNames); 
        Exist_Column3 = strcmp('RAD_SW', T.Properties.VariableNames);  %look for Samos.coaps website to see if armstrong data has a warning

        if Exist_Column2(Exist_Column2==1)
            rad_sw = T.rad_sw; 
        elseif Exist_Column3(Exist_Column3==1)
           rad_sw = T.RAD_SW; 
        end

        if startsWith(cruisename, 'HRS')
            rad_sw = T.qsr_s_n_10367; %pretty sure this sunlight data, not sure what units
        end

        T.rad_sw = rad_sw; 


%% Now go through all class files and get volume data 

classpath = [basepath '\bead_calibrated\class\'];

tempeuk = nan(height(T), length(euk_volbins)-1);
tempsyn = nan(height(T), length(syn_volbins)-1); 
meaneuk = nan(height(T), 1);
meansyn = nan(height(T), 1);
syn_mean_PE = nan(height(T), 1); 

parfor i = 1:height(T)

   S = load([classpath regexprep(T.Filename{i}, 'fcs', 'mat')]);

   if ~isfield(S, 'volume')
       continue
   end

   Svol = real(S.volume); %having some issues with imaginary numbers
   %Svol(S.negA_ind) = NaN; % don't include converted volumes, they aren't great. 
   eukdist = histcounts(Svol(S.class == 1), euk_volbins); 
   syndist = histcounts(Svol(S.class == 2), syn_volbins); 
   
   tempeuk(i, :) = eukdist; 
   tempsyn(i,:) = syndist; 
   
   meaneuk(i) = nanmean(Svol(S.class == 1));
   meansyn(i) = nanmean(Svol(S.class == 2));

end

T.EukDist = tempeuk; 
T.SynDist = tempsyn; 
T.meanEukVol = meaneuk;
T.meanSynVol = meansyn;


AttuneVolTable = T;


notes =  {'EukDist = histcounts(volume(class == 1), euk_volbins)'; 'EDGES(k) <= X(i) < EDGES(k+1)'};


save([basepath '\bead_calibrated\AttuneVolTable.mat'], 'AttuneVolTable', 'euk_volbins', 'syn_volbins', 'notes')


end