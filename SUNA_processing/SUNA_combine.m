% function to read seperate raw SUNA files in and then combine data from
% all files into one single file

% Bofu Zheng
% zhengbofuzju@gmail.com
% USE EN715 as an example


%% generate a combined SUNA file with raw data
dir_root = '/Volumes/Lab_data/SUNA/data/20230111_EN695/CTD_SUNA_NTR1227/rawdata/';  % path to raw SUNA data
dd0 = dir([dir_root '*.csv']);  % find all .csv files
n           = [];  % raw nitrate data
time        = [];
lamp_time   = [];  % suna lamp time
lamp_temp   = [];  % suna lamp temperature
absop_spec  = [];  % all wavelengths - absorption spectrum
dark_value  = [];  % dark value
volt = [];

% *** recommend opening a sample file to make sure data structure is the same ******************

for i = 1:length(dd0)
    data = readtable([dir_root,dd0(i).name]);
    %nitrate uM
    n_temp = data{:,4};
    n=[n; n_temp];
    
    volt_temp = data{:,273};
    volt = [volt; volt_temp];
    
    lamp_time_temp = data{:,271};
    lamp_time=[lamp_time; lamp_time_temp];
    
    lamp_temp_temp = data{:,270};
    lamp_temp=[lamp_temp; lamp_temp_temp];
    
    % absorption spectrum for all data
    as_temp = data{:,12:267};  % need to check the location of this array
    absop_spec = [absop_spec;as_temp];
    
    dark_temp = data{:,10};  % need to check the location of this array
    dark_value = [dark_value; dark_temp];
    
    % time
    time_temp = nan(size(data,1),1);
    for j = 1:size(data,1)  % for each sample   %%% need more atention!!!!!
        
        yearday = num2str(table2array(data(j,2)));
        year = str2num(yearday(1:4));
        day = str2num(yearday(5:7));
        hour = table2array(data(j,3));
        danum = datenum(year,1,1)-1+day+hour/24;
        
        time_temp(j) = danum;
        if j>2 & danum == time_temp(j-1)  % if there are two points having the same time stamp...
            time_temp(j) = time_temp(j-1) + 0.5/86400;  % make a time difference
        end
    end
    time=[time; time_temp];
    
    i
end

disp('finished SUNA')
clear SUNA
SUNA.time        = time;
SUNA.n           = n;
SUNA.absop_spec  = absop_spec;
SUNA.dark_value  = dark_value;
SUNA.lamp_time   = lamp_time;
SUNA.lamp_temp   = lamp_temp;
SUNA.volt = volt;

%% remove SUNA dark count
ind = find(SUNA.n == 0);
SUNA.time(ind)         = [];
SUNA.n(ind)            = [];
SUNA.absop_spec(ind,:) = [];
SUNA.dark_value(ind)   = [];
SUNA.lamp_time(ind)    = [];
SUNA.lamp_temp(ind)    = [];
SUNA.volt(ind) = [];

%% save
save(['/Users/warrbob/Desktop/WHOI/research/sunaQC/',cruise_name,'/SUNA/',cruise_name,'_SUNA.mat'],'SUNA')

