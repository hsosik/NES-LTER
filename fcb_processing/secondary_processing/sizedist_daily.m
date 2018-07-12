function [ ] = sizedist_daily( year2do )
%sizedist_daily updated from size1micron_daily, now as function and generic for
%all years; April 2017, Heidi

base_path1 = '\\sosiknas1\Lab_data\MVCO\FCB\';
temp = dir([base_path1 'MVCO_*' num2str(year2do)]);
base_path = [base_path1 temp.name filesep];
path_out = [base_path 'data\processed\size\'];

filelist = dir([path_out '*.mat']);

N_dist_crp_all = [];
N_dist_syn_all = [];
N_dist_euk_all = [];
biovol_dist_crp_all = [];
biovol_dist_euk_all = [];
biovol_dist_syn_all = [];
cellresults_all = [];
for count = 1:length(filelist),
    filename = filelist(count).name;
    disp(['loading...' filename])
    eval(['load ' path_out filename])
    
    N_dist_crp_all = [N_dist_crp_all N_dist_crp];
    N_dist_syn_all = [N_dist_syn_all N_dist_syn];
    N_dist_euk_all = [N_dist_euk_all N_dist_euk];
    biovol_dist_crp_all = [biovol_dist_crp_all biovol_dist_crp];
    biovol_dist_euk_all = [biovol_dist_euk_all biovol_dist_euk];
    biovol_dist_syn_all = [biovol_dist_syn_all biovol_dist_syn];
    cellresults_all = [cellresults_all cellresults'];
end;
clear N_dist_crp N_dist_euk N_dist_syn biovol_dist_crp biovol_dist_euk biovol_dist_syn cellresults
clear filename filelist count

day_all = floor(cellresults_all(1,:)');
mdate_day = floor(min(cellresults_all(1,:))):floor(max(cellresults_all(1,:)));

Ndist_crp = NaN(length(diambins),length(mdate_day),1);
Ndist_syn = Ndist_crp;
Ndist_euk = Ndist_crp;
biovoldist_crp = Ndist_crp;
biovoldist_syn = Ndist_crp;
biovoldist_euk = Ndist_crp;
ml_analyzed = NaN(1,length(mdate_day));
for count = 1:length(mdate_day),
    ind = find(day_all == mdate_day(count));
    if ~isempty(ind),
        Ndist_crp(:,count) = nansum(N_dist_crp_all(:,ind),2);
        Ndist_syn(:,count) = nansum(N_dist_syn_all(:,ind),2);
        Ndist_euk(:,count) = nansum(N_dist_euk_all(:,ind),2);
        biovoldist_crp(:,count) = nansum(biovol_dist_crp_all(:,ind),2);
        biovoldist_syn(:,count) = nansum(biovol_dist_syn_all(:,ind),2);
        biovoldist_euk(:,count) = nansum(biovol_dist_euk_all(:,ind),2);
        ml_analyzed(count) = nansum(cellresults_all(3,ind));
    end;
end;
clear biovol_dist* N_dist* ind count day_all cellresults_all

save([base_path1 'summary\' 'sizedist_daily_' num2str(year2do)], 'Ndist*', 'biovoldist*', 'diambins', 'ml_analyzed', 'mdate_day')