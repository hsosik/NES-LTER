% example file for processing SUNA data using TSP correction and bottle
% data match-ups 

% Bofu Zheng (zhengbofuzju@gmail.com)
% some of the functions were kindly shared by Josh Plant at MBARI.

% this processing script is an example for post-possessing SUNA nitrate data that
% are collected from a SUNA-CTD rosette intrgration. The paper documenting
% these processing steps is:
% Zheng, Bofu, E. Taylor Crockford, Weifeng (Gordon) Zhang, Rubao Ji, and Heidi M. Sosik. Bias-corrected
% high­ resolution nitrate profiles from a CTD rosette­-mounted submersible ultraviolet nitrate
% analyzer. Limnology and Oceanography: Methods, 2024


%% load all necessary files
clear all
% CTD data
% for CTD data, users are required to combine all CTD raw data from
% individual data files into one data structure. 
load(['/Users/warrbob/Desktop/WHOI/code/SUNA_QC_publish/data/EN657_CTD.mat']);
% SUNA data
% for SUNA data, users are required to combine all SUNA raw data from
% individual data files into one data structure. 
load(['/Users/warrbob/Desktop/WHOI/code/SUNA_QC_publish/data/EN657_SUNA.mat']);
% discrete bottle data
load(['/Users/warrbob/Desktop/WHOI/code/SUNA_QC_publish/data/EN657_bottle.mat']);


%% remove outliers in the bottle data
ind_throw = [];  % thrown-away index
for i = 1:length(CTDn.time)
    fake_time = CTDn.time;
    fake_time(i) = 0;
    [~, ind]=min(abs(fake_time(:)-CTDn.time(i)));  % find the closest time to this data point
    % time difference less than 1 min &&  depth difference less than 3 m &&
    % nitrate difference larger than 0.5 -> good replicates; otherwise,
    % bottle data points will be thrown away.
    if abs(CTDn.time(ind)-CTDn.time(i))<1/24/60 & abs(CTDn.d(ind)-CTDn.d(i))<3 & abs(CTDn.n(ind)-CTDn.n(i))>=0.5
        ind_throw = [ind_throw i];
    end
end

% CTDn_new is for comparing results using the TSP correction algorithm
CTDn_new.time = CTDn.time;
CTDn_new.n    = CTDn.n;
CTDn_new.time(ind_throw) = [];
CTDn_new.n(ind_throw)    = [];


%% apply TSP correction onto SUNA data

spec_fp = '/Users/warrbob/Desktop/WHOI/code/SUNA_QC_publish/data/A0000606-SUNA1227.csv';  % one example raw SUNA file
spec    = parse_Bofu_SUNA(spec_fp);  % get parameters from raw SUNA file

cal_path  = ['/Users/warrbob/Desktop/WHOI/code/SUNA_QC_publish/data//SNA1227V_BZ.CAL'];
ncal = parse_SOSIK_NO3cal(cal_path); % parameters from the calibration file

% update spec with SUNA data
spec.DC = SUNA.dark_value*ones(1,size(SUNA.absop_spec,2));
spec.UV_INTEN = SUNA.absop_spec;
[~,I]     = unique(CTDrall.time);
sinpt = interp1(CTDrall.time(I),CTDrall.s(I),SUNA.time);
tinpt = interp1(CTDrall.time(I),CTDrall.t(I),SUNA.time);
pinpt = interp1(CTDrall.time(I),CTDrall.d(I),SUNA.time);
STP = [sinpt tinpt pinpt];
spec.STP = STP;
spec.NO3 = SUNA.n;
spec.SDN = SUNA.time;
spec.spectra_WL_range(1) = ncal.WL(1);
spec.spectra_WL_range(2) = ncal.WL(end);  
Pcor_flag = 1;  % 1 means to apply pressure correction

out = calc_Bofu_NO3(spec, ncal, 1); % TSP correction!
SUNA.n_corr = out.data(:,7); % TSP corrected nitrate data

%% APPLY TEMP RESIDUE CORRECTION

niraw = nan(size(CTDn_new.time));      % raw SUNA measurements at each bottle data point
niterp = nan(size(CTDn_new.time));     % TSP-corrected SUNA measurements at each bottle data point 
niraw_std = nan(size(CTDn_new.time));  % standard deviation of raw SUNA data
niterp_std = nan(size(CTDn_new.time)); % standard deviation of TSP corrected SUNA data
for j = 1:length(CTDn_new.time)
    indtime = find(SUNA.time>=CTDn_new.time(j)-15/86400 & SUNA.time<=CTDn_new.time(j)+15/86400);
    niraw(j) = median(SUNA.n(indtime),'omitnan');
    niterp(j) = median(SUNA.n_corr(indtime),'omitnan');
    niraw_std(j) = std(SUNA.n(indtime),'omitnan');
    niterp_std(j) = std(SUNA.n_corr(indtime),'omitnan');
end

resi0 = niraw - CTDn_new.n;  % the redisue between bottle and raw SUNA data

resi1 = niterp - CTDn_new.n;  % the redisue between bottle and TSP corrected SUNA data

[~,I] = unique(CTDrall.time);
SUNA.t = interp1(CTDrall.time(I),CTDrall.t(I),SUNA.time);
SUNA.p = interp1(CTDrall.time(I),CTDrall.d(I),SUNA.time);

ind = find(SUNA.n_corr<=2.5); % as TSP corrected suna and bottle data are relatively close, so use SUNA data to locate data of less than 2.5 uM
ptemp = [0.0569 -0.7749];   % updated values for temperature-dependent correction
temp_corr = polyval(ptemp,SUNA.t(ind)); % get temperature correction values
SUNA.n_corr2 = SUNA.n_corr;
SUNA.n_corr2(ind) = SUNA.n_corr2(ind) - temp_corr ;  % correct for temperatue effect


niterp2 = nan(size(CTDn_new.time));  % if apply temp correction
for j = 1:length(CTDn_new.time)
    indtime = find(SUNA.time>=CTDn_new.time(j)-15/86400 & SUNA.time<=CTDn_new.time(j)+15/86400);
    niterp2(j) = mean(SUNA.n_corr2(indtime),'omitnan');
end
resi2 = niterp2 - CTDn_new.n;

ind = find(~isnan(niterp2));  % inverse after TSP correction
p = polyfit(CTDn_new.n(ind),niterp2(ind),1);  % bias correction use discrete bottle data
SUNA.n_corr3 = (SUNA.n_corr2-p(2))/p(1);

niterp3 = nan(size(CTDn_new.time));  % TSP + temperature dependent correction + discrete bottle correction
for j = 1:length(CTDn_new.time)
    indtime = find(SUNA.time>=CTDn_new.time(j)-15/86400 & SUNA.time<=CTDn_new.time(j)+15/86400);
    niterp3(j) = mean(SUNA.n_corr3(indtime),'omitnan');
end

resi3 = niterp3 - CTDn_new.n;  % the redisue between bottle and TSP correction /w inverse


%% make a plot to compare results at different stages

MAE(1) = rmse(niraw, CTDn_new.n,'omitnan');   % RMSE of raw data
MAE(2) = rmse(niterp, CTDn_new.n,'omitnan');  % RMSE of TSP corrected data
MAE(3) = rmse(niterp2, CTDn_new.n,'omitnan'); % RMSE of TSP + temp corrected data
MAE(4) = rmse(niterp3, CTDn_new.n,'omitnan'); % RMSE of final inverse data

figure
set(gcf,'position',[25 25 900 300])
subplot(141)
plot(CTDn_new.n,niraw,'k*')
hold on
plot([-5 30],[-5 30],'b')
axis([-5 30 -5 30])
xlabel('bottled nitrate (\muM)')
ylabel('SUNA nitrate (\muM)')
title('Raw SUNA')
text(0,25,['RMSE = ',num2str(MAE(1),'%.2f'),'+/-', num2str(mean(niraw_std,'omitnan'),'%.2f')],'fontsize',10,'color','k')


subplot(142)
plot(CTDn_new.n,niterp,'k*')
hold on
plot([-5 30],[-5 30],'b')
axis([-5 30 -5 30])
xlabel('bottled nitrate (\muM)')
ylabel('SUNA nitrate (\muM)')
title('T/S/P corrected')
text(0,25,['RMSE = ',num2str(MAE(2),'%.2f'),'+/-', num2str(mean(niterp_std,'omitnan'),'%.2f')],'fontsize',10,'color','k')


subplot(143)
plot(CTDn_new.n,niterp2,'k*')
hold on
plot([-5 30],[-5 30],'b')
axis([-5 30 -5 30])
xlabel('bottled nitrate (\muM)')
ylabel('SUNA nitrate (\muM)')
title('T/S/P corrected + temp')
text(0,25,['RMSE = ',num2str(MAE(3),'%.2f')],'fontsize',10,'color','k')


subplot(144)
plot(CTDn_new.n,niterp3,'k*')
hold on
plot([-5 30],[-5 30],'b')
axis([-5 30 -5 30])
xlabel('bottled nitrate (\muM)')
ylabel('SUNA nitrate (\muM)')
title('T/S/P corrected + temp + inverse')
text(0,25,['RMSE = ',num2str(MAE(4),'%.2f')],'fontsize',10,'color','k')

sgtitle([datestr(CTDrall.time(1),'     mmm. yyyy')])


