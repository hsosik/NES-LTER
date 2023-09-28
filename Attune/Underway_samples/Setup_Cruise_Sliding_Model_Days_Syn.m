%% Starting from Vol Table for a single cruise
%let's make a folder with input .mat files
% for each day we want to apply our division rate model 


%% this script requires input as it runs!
% to see if light data looks ok and dawn is correctly identified. 

%Outputs are 
    % daytable - very helpful table for later analysis of all the files and
    %environmental variables for each sample included in a given day
    % daystarttime - actual time of first sample included in day

    % lastdawn - actual time rounded down to the hour of the last sunrise before the model input
    % starts. So if daystarttime is May 15, 12;32, lastdawn will be
    % something like May 15, 10am. 

    % the rest are model inputs: 
    % ts - inputs for model fit indicating the region of the simulation over which division should be limited.
    %       stands for time start, but that is misleading, since now the
    %       simulation might start any time, but ts is the hour after the
    %       simulation start at which division should start and stop. 

    % Edata - another input for the model fitting, sunlight data 
    % N_dist - distribution of cell sizes used for model fitting, first bin
    % is excluded, so that is different from Vol Table.
    %Vhists - scaled distribution,  Vhists = N_dist./sum(N_dist)

    %cellsperml, Vhists, N_dists, and ts are dfferent between euks and syn 


function Setup_Cruise_Sliding_Model_Days_Syn(AttuneVolTable)

%save them all here 
synoutpath = '\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\SlidingWindow\';
eukoutpath =  '\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\eukaryotes\SlidingWindow\';

% Get cruise name from first filename
cruisename =  split(AttuneVolTable.Filename{1}, '_'); 
cruisename = cruisename{2}; 


%First, we just want to add month and season and some stuff to the table 
AttuneVolTable = AttuneVolTable(AttuneVolTable.QC_flag == 1, :); %double check to make sure low quality files are out
T = sortrows(AttuneVolTable, 'StartDate');

    StartDate = T.StartDate;
    VolAnalyzed_ml = T.VolAnalyzed_ml; 
    Euk_count = T.Euk_without_PE_leq20um_count;
    Syn_count = T.Syn_count;
    lat = T.lat;
    lon = T.lon;
    temperature = T.temperature;
    salinity = T.salinity;
    rad_sw = T.rad_sw; 
    SynDist = T.SynDist;
    EukDist = T.EukDist;
    meanSynVol = T.meanSynVol;
    meanEukVol = T.meanEukVol;
    

Filename = T.Filename; 
cruisenames = repmat(cruisename, [length(Filename), 1]); 

Month = month(T.StartDate); 

seasons = nan(1, length(Filename)) ; 
seasons(Month >= 1 & Month <= 3) = 1; %winter
seasons(Month ==4 | Month == 5 | (Month == 6 & day(T.StartDate) < 15)) = 2; %spring, April 1 to June 15
seasons(Month == 6 & (day(T.StartDate) >= 15) | Month == 7 | Month == 8 | (Month == 9 & day(T.StartDate) < 15)) =3; %summer, June 15 to Sept 15
seasons((Month == 9 & day(T.StartDate) >= 15) | Month == 10 | Month == 11 | Month == 12) = 4; %fall, sept 15 through Dec
syn_seasons = seasons';

Tnew = table(cruisenames, Filename, StartDate, VolAnalyzed_ml, Euk_count, Syn_count, lat, lon, temperature, salinity, rad_sw, SynDist, EukDist, meanSynVol, meanEukVol, syn_seasons); 

Tnew.month = Month; 

%now we actually make the inputs from the table


    mkdir([synoutpath cruisename])
        
    %bin data into hours
    hour_group = hour(Tnew.StartDate); 
    date_group = floor(datenum(Tnew.StartDate)); 
  
    [Groupnum, date_id, hour_id] = findgroups(date_group, hour_group); 
        
    numgroups = max(Groupnum); 
    
    Tnew.GroupNums = Groupnum; 
    
    for dnum = 1:max(Groupnum) %now deal with each overlapping day individually
        disp(dnum)
        
        daystarttime = Tnew.StartDate(find(Groupnum == dnum, 1));
        dayendtime = datenum(daystarttime)+1.0417; %25 hours later 
        
        if dayendtime > datenum(max(Tnew.StartDate))
            continue 
        end
        
        daytable = Tnew((Tnew.StartDate >= daystarttime & datenum(Tnew.StartDate) < dayendtime), :) ;
        numhours = length(unique(hour(daytable.StartDate)));
        if numhours < 22 %missing more than 2 hours of a day 
            continue
        end
        
        t_relstart = datenum(daytable.StartDate) - datenum(daystarttime); 
        hournum = floor(t_relstart.*24)+1;
        
        %format light data
        Edata = [t_relstart.*24 daytable.rad_sw];
        Edata(Edata(:,2) < 0, 2) = 0; %remove negative light values 
        
        %need to keep track of dawn also though, since we don't let syn
        %divide in the first 6 hours
        
        nighttime = find(daytable.rad_sw <= 3);
        %if  %some cruises have nighttime light noise
        %    nighttime = find(daytable.rad_sw <= 100); 
        %end
        dawnind = nighttime((nighttime(2:end) - nighttime(1:end-1)) > 9); % look for periods of darkness where next period of darkness isn't for a while 
        if isempty(dawnind)
            daytime = find(daytable.rad_sw > 10);
            dawnind = daytime(find((daytime(2:end) - daytime(1:end-1))>12)+1); %or look for periods of light where it hasn't been light for while 
        end
        if isempty(dawnind) %if day just happend to start at sunset we wont have more than 1 transition
            dawnind = nighttime(end);
        end

        if length(dawnind) > 1; 
            dawnind = dawnind(1); 
        end
        
        dawn_hr = hournum(dawnind);
        ts = [dawn_hr dawn_hr+6]; %limit division up to 6 hours after dawn
        

        lastdawn = datetime((datenum(daystarttime) + (dawn_hr - 24)./24), 'ConvertFrom', 'datenum'); 
        lastdawn = dateshift(lastdawn, 'start', 'hour'); 

        %check dawn choice
        if rem(dnum, 14) == 1 
        figure(2)
        clf
        plot(daytable.StartDate, daytable.rad_sw)
        hold on 
        scatter(daytable.StartDate(dawnind), dawnind.*0+50, 50, 'filled')
    
        % check = input('look ok? y/n', 's');
        % if strcmp(check, 'n')
        %     keyboard
        % end
        end
        
        synvolbins = 2.^[-5:0.125:2];
        eukvolbins = [0 2.^[-5:1/5:8]]; 

        %format syn data inputs            
        func2 = @(x) nansum(x, 1); 
 
        %much ado about adding empty rows 
        missinghr = find(sum(hournum == [1:25], 1) == 0);
        daytable.Hournum = hournum; 
        daytable.Hournum((height(daytable)+1):(height(daytable)+length(missinghr))) = missinghr; 
        daytable.Syn_count(end-length(missinghr)+1:end) = NaN; %0s are problematic
        daytable.Euk_count(end-length(missinghr)+1:end) = NaN; 
        daytable.VolAnalyzed_ml(end-length(missinghr)+1:end) = NaN; 
        daytable.lat(end-length(missinghr)+1:end) = NaN; 
        daytable.lon(end-length(missinghr)+1:end) = NaN; 
        daytable.temperature(end-length(missinghr)+1:end) = NaN; 
        daytable.salinity(end-length(missinghr)+1:end) = NaN; 
        daytable.rad_sw(end-length(missinghr)+1:end) = NaN; 
        daytable.SynDist(end-length(missinghr)+1:end, :) = NaN; 
        daytable.EukDist(end-length(missinghr)+1:end, :) = NaN; 
        daytable.meanSynVol(end-length(missinghr)+1:end) = NaN; 
        daytable.meanEukVol(end-length(missinghr)+1:end) = NaN; 
        daytable.month(end-length(missinghr)+1:end) = NaN; 
        daytable.syn_seasons(end-length(missinghr)+1:end) = NaN; 
        daytable.GroupNums(end-length(missinghr)+1:end) = NaN; 

        
        %save Synechococcus input files first
        cellsperml = splitapply(func2, daytable.Syn_count, daytable.Hournum)./splitapply(func2, daytable.VolAnalyzed_ml, daytable.Hournum); 
                
      
        N_dist = splitapply(func2, daytable.SynDist, daytable.Hournum)';
        N_dist(1,:) = 0; %smallest size bin is extra large and should not be included in model fitting..... 
        %sometimes 
        N_dist(:,(sum(N_dist)==0)) = NaN;
        Vhists = N_dist./sum(N_dist); 
        
        save([synoutpath cruisename filesep cruisename 'day' num2str(dnum, '%03.f') 'input.mat'], 'daystarttime', 'lastdawn', 'daytable', 'synvolbins', 'N_dist', 'cellsperml', 'Vhists', 'Edata', 'ts')
    
    
        %then do it for Eukaryotes

        cellsperml = splitapply(func2, daytable.Euk_count, daytable.Hournum)./splitapply(func2, daytable.VolAnalyzed_ml, daytable.Hournum); 
                
        N_dist = splitapply(func2, daytable.EukDist, daytable.Hournum)';
        N_dist(1,:) = 0; %smallest size bin is extra large and should not be included in model fitting..... 
        %sometimes 
        N_dist(:,(sum(N_dist)==0)) = NaN;
        Vhists = N_dist./sum(N_dist); 
        
        save([eukoutpath cruisename filesep cruisename 'day' num2str(dnum, '%03.f') 'Euk_input.mat'], 'daystarttime', 'daytable', 'eukvolbins', 'N_dist', 'cellsperml', 'Vhists', 'Edata', 'ts')
  
    
    end

   
    
    end

