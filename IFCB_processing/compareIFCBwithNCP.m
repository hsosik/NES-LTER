cruises = {'AR22' 'AR24A' 'AR24B' 'AR24C' 'EN608' 'AR28A' 'AR28B' 'EN617'...
    'AR31A' 'AR31B' 'AR31C' 'AR32' 'EN627' 'AR34A' 'AR34B' 'EN644'};  %'AR16' 'AR34B' 'AR38' 'AR39A' 'AR39A'

cruises = {'EN617' 'EN644'};

ubase = '\\sosiknas1\IFCB_data\NESLTER_transect\match_up\';
%ubase = 'c:\work\IFCB_products\NESLTER_transect\match_up\';

match_uw = table;
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
 end
ind = find(strcmp(match_uw.cruise, 'EN627') & match_uw.salinity > 34);
match_uw.temperature(ind) = NaN;
match_uw.salinity(ind) = NaN;

s2017 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017.mat');
s2018 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018.mat');
s2019 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019.mat');
%s2017 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017.mat');
%s2018 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018.mat');
%s2019 = load('c:\work\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019.mat');
f2017 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2017lists.mat', 'classFeaList', 'classFeaList_variables');
f2018 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2018lists.mat', 'classFeaList', 'classFeaList_variables');
f2019 = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\summary_biovol_allHDF_min20_2019lists.mat', 'classFeaList', 'classFeaList_variables');

IFCBsum = table;
slist = {'filelist' 'classcount' 'classbiovol' 'mdate' 'meta_data'};
for count = 1:length(slist)
    s = slist{count}; IFCBsum.(s) = [s2017.(s); s2018.(s); s2019.(s)];
end
IFCBsum.classFeaList = [f2017.classFeaList; f2018.classFeaList; f2019.classFeaList];
class2use = s2017.class2use;
% FUDGE for now, already fixed in ifcb-db
IFCBsum.meta_data.skip(strmatch('D20170518T174547_IFCB009', IFCBsum.filelist)) = 1;

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
[~,ia,ib] = intersect(group_table.CNN_classlist, class2use);
diatom_ind = ib(find(group_table.Diatom(ia)));
notalive_ind = [ib(find(group_table.OtherNotAlive(ia))); ib(find(group_table.IFCBArtifact(ia)))];
alive_ind = 1:length(class2use); alive_ind(notalive_ind) = [];
alive_ind(strmatch( 'Phaeocystis', class2use(alive_ind))) = []; %seems to be detritus for en644
alive_ind(strmatch( 'unclassified', class2use(alive_ind))) = [];

fi1 = strmatch('ESD', f2019.classFeaList_variables);
fi2 = strmatch('maxFeretDiameter', f2019.classFeaList_variables);
numfiles = length(IFCBsum.filelist);
alivecount_gt10 = NaN(numfiles,1);
alivebv_gt10 = alivecount_gt10;
diatomcount_gt10 = alivecount_gt10;
diatombv_gt10 = alivecount_gt10;
for ii = 1:numfiles
    if ~rem(ii,100), disp(IFCBsum.filelist(ii)), end
    if ~IFCBsum.meta_data.skip(ii)
        temp = cat(1,IFCBsum.classFeaList{ii,alive_ind}); 
        bv = 4/3*pi*(temp(:,fi1)/2).^3;
        alivecount_gt10(ii,2) = size(temp,1);
        gt10i = find(temp(:,fi1)>10 | temp(:,fi2)>10);
        alivecount_gt10(ii,1) = numel(gt10i); 
        alivebv_gt10(ii,2) = sum(bv);
        alivebv_gt10(ii,1) = sum(bv(gt10i));
        temp = cat(1,IFCBsum.classFeaList{ii,diatom_ind}); 
        bv = 4/3*pi*(temp(:,fi1)/2).^3;
        diatomcount_gt10(ii,2) = size(temp,1);
        gt10i = find(temp(:,fi1)>10 | temp(:,fi2)>10);
        diatomcount_gt10(ii,1) = numel(gt10i); 
        diatombv_gt10(ii,2) = sum(bv);
        diatombv_gt10(ii,1) = sum(bv(gt10i));
    end
end
cc = strmatch('Hemiaulus', class2use);

Hemiaulus_count_gt10 = NaN(numfiles,2);
Hemiaulus_bv_gt10 = Hemiaulus_count_gt10;
for ii = 1:numfiles     
    if ~rem(ii,100), disp(IFCBsum.filelist(ii)), end
    ind = find(IFCBsum.classFeaList{ii,cc}(:,fi1)>10 | IFCBsum.classFeaList{ii,cc}(:,fi2)>10);
    Hemiaulus_count_gt10(ii,1) = numel(ind);
    bv = 4/3*pi*(IFCBsum.classFeaList{ii,cc}(:,fi1)/2).^3;
    Hemiaulus_count_gt10(ii,2) = numel(bv);
    Hemiaulus_bv_gt10(ii,1) = sum(bv(ind));
    Hemiaulus_bv_gt10(ii,2) = sum(bv);
end

ncp617 = load('C:\work\LTER\Stanley_data\ncplterEn617.mat');
ncp644 = load('C:\work\LTER\Stanley_data\ncplterEn644.mat');

en644ind = find(strcmp(IFCBsum.meta_data.cruise, 'EN644') & strcmp(IFCBsum.meta_data.sample_type, 'underway') & ~IFCBsum.meta_data.skip);
en617ind = find(strcmp(IFCBsum.meta_data.cruise, 'EN617') & strcmp(IFCBsum.meta_data.sample_type, 'underway') & ~IFCBsum.meta_data.skip);
%%
ind = en644ind; ncp = ncp644; dr = [datenum(2019,8,20,18,0,0) datenum(2019,8,25,1,0,0)]; yl = [-.3e6 3e6];
%ind = en617ind; ncp = ncp617; dr = [datenum(2018,7,20,12,0,0) datenum(2018,7,25,12,0,0)]; yl = [-.1e6 .5e6]; 

figure('position', [120 340 1120 420])
yyaxis right, plot(ncp.ncplter(:,1), ncp.ncplter(:,9), '-', 'linewidth', 2), ylim([-10 70])
ylabel('Net community production (mmol O_2 m^2 d^{-1})')
yyaxis left
ml = IFCBsum.meta_data.ml_analyzed(ind);
plot(IFCBsum.mdate(ind), alivebv_gt10(ind,1)./ml , 'b.-.', 'linewidth', 1), hold on
%plot(IFCBsum.mdate(ind), diatombv_gt10(ind,1)./ml , 'c.-', 'linewidth', 1)
%legend('Biovolume 10-200 \mum', 'Biovolume, diatoms', 'NCP')
plot(IFCBsum.mdate(ind), Hemiaulus_bv_gt10(ind,1)./ml , 'c.-', 'linewidth', 1)
legend('Biovolume 10-200 \mum', 'Biovolume, \it{Hemiaulus}', 'NCP')
ylabel('Biovolume (\mum^3 ml^{-1})')
yyaxis left, ylim(yl)
set(gca, 'xlim', datenum(dr), 'linewidth', 2)
datetick('x', 1, 'keeplimits')

figure('position', [120 340 1120 200])
yyaxis left, plot(IFCBsum.mdate(ind), IFCBsum.meta_data.latitude(ind), '-', 'linewidth', 2), ylabel('Latitude')
yyaxis right, plot(IFCBsum.mdate(ind), IFCBsum.meta_data.longitude(ind), '-', 'linewidth', 2), ylabel('Longitude')
datetick('x', 1, 'keeplimits')
set(gca, 'xlim', datenum(dr), 'linewidth', 2)

%%
ncpmatch617 = interp1(ncp617.ncplter(:,1),ncp617.ncplter(:,9), IFCBsum.mdate(en617ind));
ncpmatch644 = interp1(ncp644.ncplter(:,1),ncp644.ncplter(:,9), IFCBsum.mdate(en644ind));
figure
plot(alivebv_gt10(en617ind,1)./IFCBsum.meta_data.ml_analyzed(en617ind), ncpmatch617, '.', 'markersize', 8)
axis square
hold on
plot(alivebv_gt10(en644ind,1)./IFCBsum.meta_data.ml_analyzed(en644ind), ncpmatch644, '*')
legend('EN617, summer 2018', 'EN644, summer 2019')
ylabel('Net community production (mmol O_2 m^2 d^{-1})')
xlabel('Biovolume (\mum^3 ml^{-1})')
ylabel('NCP (mmol O_2 m^2 d^{-1})')

figure
plot(diatombv_gt10(en617ind,1)./IFCBsum.meta_data.ml_analyzed(en617ind), ncpmatch617, '.', 'markersize', 8)
axis square
hold on
plot(diatombv_gt10(en644ind,1)./IFCBsum.meta_data.ml_analyzed(en644ind), ncpmatch644, '*')
legend('EN617, summer 2018', 'EN644, summer 2019')
ylabel('Net community production (mmol O_2 m^2 d^{-1})')
xlabel('Biovolume, diatoms (\mum^3 ml^{-1})')
ylabel('NCP (mmol O_2 m^2 d^{-1})')


%% review ESD and maxFeret histograms
for ii = 1:100:numfiles
    if ~rem(ii,100), disp(IFCBsum.filelist(ii)), end
    if ~IFCBsum.meta_data.skip(ii)
        temp = cat(1,IFCBsum.classFeaList{ii,alive_ind}); 
        figure(99), clf
        subplot(1,2,1)
        histogram(temp(:,fi1),0:10:300), set(gca, 'yscale', 'log', 'ylim', [.1 inf])
        subplot(1,2,2)
        histogram(temp(:,fi2),0:10:500), set(gca, 'yscale', 'log', 'ylim', [.1 inf])       
        pause
    end
end