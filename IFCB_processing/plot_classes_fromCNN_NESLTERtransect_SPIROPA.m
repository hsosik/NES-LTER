cruises = {'AR22' 'AR24A' 'AR24B' 'AR24C' 'EN608' 'AR28A' 'AR28B' 'EN617'...
    'AR31A' 'AR31B' 'AR31C' 'AR32' 'AR44' 'EN627' 'AR34A' 'AR34B' 'AR38' 'AR39a'...
    'AR39B' 'EN644' 'EN649' 'EN655' 'EN657' 'EN661'};  %'AR16' 'AR48A' 'AR48B'  
cruises2 = {'AR29' 'RB1904' 'TN368' };

%cruises = {'EN608' 'AR28B' 'EN617' 'AR31B'};  %'AR16' 'AR34B' 'AR38' 'AR39'
%cruises = {'EN608' 'EN617' 'EN627' 'EN644'};
%cruises = {'EN617' 'EN644'};

ubase = '\\sosiknas1\IFCB_data\NESLTER_transect\match_up\';
ubase2 = '\\sosiknas1\IFCB_data\SPIROPA\match_up\';
%ubase = 'c:\work\IFCB_products\NESLTER_transect\match_up\';

match_uw = table;
nut = table;
match_cast = table;
slist = {'pid' 'cruise' 'lat' 'lon' 'temperature' 'salinity' 'mdate'};
slist_cast = {'pid' 'cruise' 'lat' 'lon' 't090c' 'sal00' 'depth' 'mdate'};

for count1 = 1:length(cruises)
     disp(cruises{count1})
     if strmatch('EN627', cruises{count1})
         orig_var = {'tsg2_temperature' 'tsg2_salinity'};
     elseif strmatch('EN', cruises{count1}(1:2))
         orig_var = {'tsg1_temperature' 'tsg1_salinity'};
     else
         orig_var = {'sbe48t' 'sbe45s'};
     end
     u = load([ubase 'NESLTER_transect_' cruises{count1} '_uw_match.mat']);
     u.IFCB_match_uw_results.Properties.VariableNames(orig_var) = {'temperature' 'salinity'};
     [~,ia] = ismember(slist, u.IFCB_match_uw_results.Properties.VariableNames);
     match_uw = [match_uw; u.IFCB_match_uw_results(:,ia)];
     if exist([ubase 'NESLTER_transect_' cruises{count1} '_cast_match.mat'], 'file')
        u = load([ubase 'NESLTER_transect_' cruises{count1} '_cast_match.mat']);
        u.IFCB_match_btl_results.mdate = datenum(u.IFCB_match_btl_results.datetime, 'yyyy-mm-dd hh:MM:SS+00:00');
        [~,ia] = ismember(slist_cast, u.IFCB_match_btl_results.Properties.VariableNames);
        match_cast = [match_cast; u.IFCB_match_btl_results(:,ia)];
     end
%  %   n = webread(['https://nes-lter-data.whoi.edu/api/nut/' cruises{count1} '.csv']);
%  %   n.alternate_sample_id = []; %move this column since the type doesn't match between all cruises
%  %   nut = [nut; n];
 end
% %nut.mdate = datenum(nut.date, 'yyyy-mm-dd hh:MM:ss+00:00');
ind = find(strcmp(match_uw.cruise, 'EN627') & match_uw.salinity > 34);
match_uw.temperature(ind) = NaN;
match_uw.salinity(ind) = NaN;
%%
for count1 = 1:length(cruises2)
     disp(cruises2{count1})
     if strmatch('AR29', cruises2{count1})
         orig_var = {'SBE48T' 'SBE45S'};
     else
         orig_var = {'T' 'SSPS'};
     end
     u = load([ubase2 'SPIROPA_' cruises2{count1} '__newdashboard_USEMEuw_match.mat']);
     u.IFCB_match_uw_results.Properties.VariableNames(orig_var) = {'temperature' 'salinity'};
     [~,ia] = ismember(slist, u.IFCB_match_uw_results.Properties.VariableNames);
     match_uw = [match_uw; u.IFCB_match_uw_results(:,ia)];
     if exist([ubase2 'SPIROPA_' cruises2{count1} '__newdashboard_USEMEcast_match.mat'], 'file')
        u = load([ubase2 'SPIROPA_' cruises2{count1} '__newdashboard_USEMEcast_match.mat']);  
        u.IFCB_match_btl_results.Properties.VariableNames = lower(u.IFCB_match_btl_results.Properties.VariableNames);
        if strmatch('AR29', cruises2{count1})
            uind = find(u.IFCB_match_btl_results.cast == 106);
            if isempty(u.IFCB_match_btl_results.datetime{uind(1)})
                u.IFCB_match_btl_results.Depth_m(uind) = 2;
                u.IFCB_match_btl_results.datetime(uind) = {'2018-04-25 04:13:51+00:00'};
            end
        end    
        u.IFCB_match_btl_results.mdate = datenum(u.IFCB_match_btl_results.datetime, 'yyyy-mm-dd hh:MM:SS+00:00');
        [~,ia] = ismember(slist_cast, u.IFCB_match_btl_results.Properties.VariableNames);
        match_cast = [match_cast; u.IFCB_match_btl_results(:,ia)];
     end
 end
%%

s2017 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017.mat');
s2018 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018.mat');
s2019 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019.mat');
s2020 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2020.mat');
%s2021 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2021.mat');
s2018b = load('\\sosiknas1\IFCB_products\SPIROPA\summary\summary_biovol_allHDF_min20_2018.mat');
s2019b = load('\\sosiknas1\IFCB_products\SPIROPA\summary\summary_biovol_allHDF_min20_2019.mat');
tag5 = repmat(cellstr(''),size(s2018b.meta_data,1),1);
s2018b.meta_data = addvars(s2018b.meta_data, tag5, 'After', 'tag4', 'NewVariableNames', 'tag5' );
s2018b.meta_data.cast = cellstr(num2str(s2018b.meta_data.cast));
s2018b.meta_data.tag4 = cellstr(num2str(s2018b.meta_data.tag4));
tag5 = repmat(cellstr(''),size(s2019b.meta_data,1),1);
s2019b.meta_data = addvars(s2019b.meta_data, tag5, 'After', 'tag4', 'NewVariableNames', 'tag5' );
s2019b.meta_data.cast = cellstr(num2str(s2019b.meta_data.cast));
s2019b.meta_data.tag4 = cellstr(num2str(s2019b.meta_data.tag4));
%s2017 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017.mat');
%s2018 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018.mat');
%s2019 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019.mat');

IFCBsum = table;
slist = {'filelist' 'classcount' 'meta_data' 'classbiovol' 'mdate'};
for count = 1:length(slist)
    s = slist{count}; IFCBsum.(s) = [s2017.(s); s2018.(s); s2019.(s) ;s2020.(s) ; s2018b.(s); s2019b.(s)];
end

class2use = s2017.class2use;

match_uw = sortrows(match_uw, 'mdate');
[a,b] = ismember(match_uw.pid, IFCBsum.filelist);
%
match_cast = sortrows(match_cast, 'mdate');
[a_cast,b_cast] = ismember(match_cast.pid, IFCBsum.filelist);


temp = find(a==0);
if ~isempty(temp)
    disp('Skipping unmatched files: ')
    disp(match_uw.pid(temp))
    disp('Skipped above unmatched files')
    match_uw(temp,:) = [];
    b(a==0) = [];
end

temp = find(a_cast==0);
if ~isempty(temp)
    disp('Skipping unmatched files: ')
    disp(match_cast.pid(temp))
    disp('Skipped above unmatched files')
    match_cast(temp,:) = [];
    b_cast(a_cast==0) = [];
end

group_table = readtable('\\sosiknas1\training_sets\IFCB\config\IFCB_classlist_type.csv');
group_table.CNN_classlist(strmatch('Pseudo-nitzschia', group_table.CNN_classlist)) = {'Pseudo_nitzschia'};
[~,ia,ib] = intersect(group_table.CNN_classlist, class2use);
diatom_ind = ib(find(group_table.Diatom(ia)));

dv = datevec(match_uw.mdate);
yd_vec = match_uw.mdate-datenum(dv(:,1),1,0);

Z2 = (sum(IFCBsum.classbiovol(b,diatom_ind),2)./IFCBsum.meta_data.ml_analyzed(b));

Z2all = IFCBsum.classbiovol(b,diatom_ind)./IFCBsum.meta_data.ml_analyzed(b);

[ mdate_mat, y_mat, yearlist, yd ] = timeseries2ydmat( match_uw.mdate, Z2 );

for ii = 1:12, y_month(ii) = nanmean(Z2(dv(:,2)==ii)); end

ilat = 39.7:.1:41.5;
for ii = 1:length(ilat)-1, iii = find(match_uw.lat >ilat(ii) & match_uw.lat<=ilat(ii+1)); Z2latmn(ii) = nanmean(Z2(iii));end

figure
boxplot(Z2,dv(:,2), 'Whisker', 5, 'notch', 'on')
text(1:12,4.5e6*ones(12,1),num2str(histcounts(dv(:,2),0:12)'), 'fontsize', 10)

figure
plot(ilat(1:end-1)+.05, Z2latmn,'-*', 'linewidth', 3)

figure
plot(1:366, nanmean(y_mat,2), 'linewidth', 3)
