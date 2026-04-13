% function to read seperate CTD files in and then combine data from
% all files into one single file

% Bofu Zheng
% zhengbofuzju@gmail.com
% USE EN715 as an example


%% combine CTD files 
clear all
cruise_name = 'EN715';
cruise_name_low = 'en715';

%%
filepath  =  '/Volumes/Lab_data/LTER/20240503_EN715/ctd/proc/time_bin/';  % filepath to all CTD cast data
dd0 = dir([filepath '*.asc']);
mlat = nan(length(dd0),1);
mlon = nan(length(dd0),1);
castn = nan(length(dd0),1);
for i = 1:length(dd0)
    filename  =  dd0(i).name;   % set file name
    cast      = str2double(filename(7:8));
    filein    =  [filepath filename];
    datain = readtable(filein,'FileType','text');% file to be loaded
    dan = datenum(2024,1,1)-1;  % base time
    %%
    CTD = [];
    for c = 1:length(datain.Properties.VariableNames)
        switch datain.Properties.VariableNames{c}
            case 'TimeQ'
                timeq = datain{:,c}/86400 + dan;
            case 'TimeJ'
                timej = datain{:,c};
            case 'TimeS'
                times = datain{:,c};
            case 'DepSM'
                CTD.depth = datain{:,c};
            case 'PrDM'
                CTD.pres = datain{:,c};
            case 'CStarAt0'
                CTD.beamattenuation = datain{:,c};
            case 'FlECO_AFL'
                CTD.chla = datain{:,c};
            case 'AltM'
                CTD.altimeter = datain{:,c};
            case 'Par'
                CTD.PAR = datain{:,c};
            case 'Spar'
                CTD.surfacePAR = datain{:,c};
            case 'OxsolMm_Kg'
                oxy1 = datain{:,c};
            case 'OxsatMm_Kg'
                oxy2 = datain{:,c};
            case 'T090C'
                temp1 = datain{:,c};
            case 'T190C'
                temp2 = datain{:,c};
            case 'C0S_m'
                CTD.cond = datain{:,c} * 10;
            case 'Latitude'
                CTD.lat = datain{:,c};
            case 'Longitude'
                CTD.lon = datain{:,c};
            otherwise
                disp('others');
        end
    end
    
    
    CTD.time  = dan + timej(1) + times/86400;
%     CTD.pres = data{4};
    CTD.temp  = nanmean([temp1 temp2],2);
%     CTD.cond  = data{9}*10;
    CTD.sali  = gsw_SP_from_C(CTD.cond,CTD.temp,CTD.pres); % need the gsw package: https://www.teos-10.org/software.htm#1
    % CTD.sali  = nanmean([data{32} data{33}],2);
    CTD.dens  = gsw_rho(CTD.sali,CTD.temp,0);
    CTD.oxy  = nanmean([oxy1 oxy2],2);
%     CTD.lat   = data{29};
%     CTD.lon   = data{30};
    CTD.cast  = cast*ones(size(CTD.dens));
    CTD.cruisename = repmat(filename(1:5),size(CTD.dens,1),size(CTD.dens,2));
%     CTD.cruisename = filename(1:5);
    
    %% save
    fname = ['cast',num2str(cast)];
    save(['/Users/warrbob/Desktop/WHOI/research/sunaQC/',cruise_name,'/CTD/mat/',fname,'.mat'],'CTD')  % output folder
    
end



%% combine CTD casts all together and make a CTD structure

dir_root = ['/Users/warrbob/Desktop/WHOI/research/sunaQC/',cruise_name,'/CTD/mat/'];  % filepath to all CTD mat files
dd0 = dir([dir_root '*.mat']);  
CTDr = cell(size(dd0));
CTDrall.time = [];
CTDrall.t    = [];
CTDrall.s    = [];
CTDrall.rho  = [];
CTDrall.d    = [];
CTDrall.lat  = [];
CTDrall.lon  = [];
CTDrall.oxy  = [];
CTDrall.beamatt = [];
CTDrall.chla = [];

for i = 1:length(dd0)
    load([dir_root dd0(i).name]);
    CTDr{i} = CTD;
    CTDrall.time = [CTDrall.time; CTD.time];
    CTDrall.t    = [CTDrall.t; CTD.T_corr];
    CTDrall.s    = [CTDrall.s; CTD.S_corr];
    CTDrall.rho  = [CTDrall.rho; CTD.dens_corr];
    CTDrall.d    = [CTDrall.d; CTD.depth];
    CTDrall.lat  = [CTDrall.lat; CTD.lat];
    CTDrall.lon  = [CTDrall.lon; CTD.lon];
    CTDrall.oxy  = [CTDrall.oxy; CTD.oxy];
    CTDrall.beamatt = [CTDrall.beamatt; CTD.beamattenuation];
    CTDrall.chla = [CTDrall.chla; CTD.chla];
end

%
save(['/Users/warrbob/Desktop/WHOI/research/sunaQC/',cruise_name,'/CTD/matall/',cruise_name,'_CTD.mat'],'CTDrall')  % output folder
