% March 2026, Heidi's reanalysis of March/April 2019 calibration expt
% see "\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\BLF_tying_to_figure_out_where_table_came_from.m" for history 
% Here addressing a few issues:
% 1. histogram fits to Attune SSC were not always resolved appropriately
% with default histcounts binning (e.g., for small signals such as Syn)
% 2. Previous ssc to volume fit was based on ssc with OD2 only (which is not good for the smallest cultures)
% 21 March 2026, change to size_calibration_fit_march2019_v2 to handle
% newly gated FCS files, Emily cleaned up cell and bead gates (e.g., refining to remove long triggers, etc.)

%attune_par = 'FSC-A'; attune_min = -1e4;
%attune_par = 'FSC-H'; attune_min = -1e4;
attune_par = 'SSC-A'; attune_min = -200;
%attune_par = 'SSC-H'; attune_min = -200;

%, SSC-H, FSC-H, FSC-A

load('\\sosiknas1\Lab_data\Attune\Size_calibration_March2019\CCsizecal_table.mat')
%Heidi, Mar 2026, this is diameter in col1 and freq in col2, collated from 'size_calibration_data_March2019.mat'
%CCsize_threshold values are diameters and start with Nano (e.g., Syn and Micromonas not included)
CCsize_threshold = [0 0 CCsize_threshold];
printNow = 0;
saveNow = 0;
%%
figure('WindowState','maximized')
tlcc = tiledlayout(3,4);
cc_mode = table;
cc_mode.files = CCsizecal.Properties.VariableNames';
for i = 1:width(CCsizecal)
%    vname = CCsizecal.Properties.VariableNames{i};
    cc = CCsizecal{:,i};
    cc = cc(cc(:,1) > CCsize_threshold(i),:); % remove CC junk below size threshold
    %cc = 4/3*pi*(cc/2).^3;
    cc(:,1) = 4/3*pi*(cc(:,1)/2).^3; %Heidi: the line above is BF's but it should be this one (very small diff)
    %cc_smooth = smooth(cc(:,2), 10);
    if ismember(cc_mode.files{i}, {'Syn' 'Micromonas'})
        cc_smooth = smooth(medfilt1(cc(:,2),3),4);
        %cc_smooth = smooth(cc(:,2), 10);
  %  elseif ismember(cc_mode.files{i}, {Prymn_pav'})
   %     cc_smooth = smooth(medfilt1(cc(:,2),8),20);
        %cc_smooth = smooth(cc(:,2), 40);
    elseif ismember(cc_mode.files{i}, {'Pptp02'})
        cc_smooth = smooth(medfilt1(cc(:,2),10),30);
        %cc_smooth = smooth(cc(:,2), 40);
    else
        cc_smooth = smooth(medfilt1(cc(:,2),7),3);
        %cc_smooth = smooth(cc(:,2), 20);
    end
    [cc_max, ccind] = max(cc_smooth);
    cc_mode.cc_mode_vol(i) = cc(ccind,1);
    nexttile
    plot(cc(:,1), cc(:,2), cc(:,1), cc_smooth); hold on;
    %set(gca, 'XScale', 'log');
    line([cc_mode.cc_mode_vol(i) cc_mode.cc_mode_vol(i)], [0, cc_max], 'color', 'k'); hold off;
    title(cc_mode.files{i}, 'Interpreter','none')
    xlim([0 cc_mode.cc_mode_vol(i)*10]);
    grid
end
xlabel(tlcc, 'Coulter volume (\mum^3)')
if printNow
    print('\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Analysis_March2026\cc_volume_hist.png')
end
%%
% now od2
attune_bin_num = 120;
attune_smooth_num = 15;
attune_medfilt_num = 6;
str = '\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Gated_FCS_March2026\';
files = dir([str '\\OD2*\\*.fcs']);
fcsOD2T = array2table([{files.name}' {files.folder}'], 'VariableNames',{'file' 'folder'});
runDate = split(fcsOD2T.folder,'_'); 
fcsOD2T.run_date = runDate(:,end); clear runDate
for ii = 1:length(files)
    [fcsdat, fcshdr] = fca_readfcs(fullfile(files(ii).folder, files(ii).name));
    fcsOD2T.dat{ii} = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
end
%%

figure(99)
set(gcf, 'WindowState','maximized')
tlod2 = tiledlayout(4,5);
figure(98)
tlcyto1 = tiledlayout(4,5);
OD2_AttuneSizeCal = table;
for i = 1:height(fcsOD2T)
    dat = fcsOD2T.dat{i}.(attune_par);
    [N, edges] = histcounts(dat,300);
    att_val = edges(1:end-1) + diff(edges)/2;
    varname = regexprep(fcsOD2T.file{i}, '.fcs', '');
    if startsWith(varname, 'beads')
        varname = strcat(varname,'_', fcsOD2T.run_date{i}); 
    end
    figure(99)
    nexttile
    plot(att_val, N./diff(edges)), title(varname, 'Interpreter','none')
    hold on
    %this is how BF did the histograms
    [N, edges] = histcounts(dat);
    %now let's force more resolution around the peak (impt for Syn)
    [att_max, attind] = max(N);
    std_temp = std(dat(dat<edges(attind+1)*3));
    dat_upper = edges(attind+1)+std_temp*3;
    dat_lower = edges(attind+1)-std_temp*3;
    xlim([-inf min(edges(attind+1)+std_temp*3, max(dat))])
    xl = xlim; xl(1) = min([xl(1) -1e4]); xlim(xl), clear xl
    [N, edges] = histcounts(dat(dat<dat_upper & dat>dat_lower & dat>attune_min),attune_bin_num); %dat>-1e4
    if strcmp(varname, 'beads_mcl') & strcmp(attune_par, 'FSC-H')
        N = medfilt1(N,3);
    else
        N = smooth(medfilt1(N,attune_medfilt_num),attune_smooth_num)';
    end
    att_val = edges(1:end-1) + diff(edges)/2;
    plot(att_val,N./diff(edges), 'r', 'LineWidth',2)
    OD2_AttuneSizeCal.(varname) = {[att_val' N']}; %save this final one
    grid
    figure(98)
    nexttile
    loglog(dat, fcsOD2T.dat{i}.('BL3-A'), '.'), title(varname, 'Interpreter','none')
%    varname = regexprep(files(i).name, '.fcs', '');
%    assignin('base', varname, [att_val' N']);
end
xlabel(tlcyto1, attune_par)
ylabel(tlcyto1, 'CHL-A')
figure(99)
xlabel(tlod2,['Attune ' attune_par ', OD2 on SSC'])
ylabel(tlod2,'Number per X-unit')
OD2_AttuneSizeCal = renamevars(OD2_AttuneSizeCal, 'HTriq', 'Htriq'); %match cc
if printNow
    print(['\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Analysis_March2026\attune_' attune_par '_OD2_hist.png'], '-dpng')
end
%%
[bmax, bind] = max(OD2_AttuneSizeCal.beads_1um_20Mar2019{1}(:,2));
beads_1um_mode_OD2a = OD2_AttuneSizeCal.beads_1um_20Mar2019{1}(bind,1); 
[bmax, bind] = max(OD2_AttuneSizeCal.beads_1um_9Apr2019{1}(:,2));
beads_1um_mode_OD2b = OD2_AttuneSizeCal.beads_1um_9Apr2019{1}(bind,1); 
beads_1um_mode_OD2 = mean([beads_1um_mode_OD2a,beads_1um_mode_OD2b]);

OD2_att_mode = table;
OD2_att_mode.files = OD2_AttuneSizeCal.Properties.VariableNames';
OD2_att_mode.file_orig = OD2_att_mode.files;
temp = startsWith(OD2_att_mode.files, 'Syn');
OD2_att_mode.files(temp) = {'20Mar_Syn'}; 
temp = startsWith(OD2_att_mode.files, '1820i');
OD2_att_mode.files(temp) = {'Syn'}; %force this one to match with coulter counter

for i = 1:height(OD2_att_mode)
    att = OD2_AttuneSizeCal.(OD2_att_mode.file_orig{i}){1};
    [att_max, attind] = max(att(:,2));
    OD2_att_mode.OD2_mode(i) = att(attind,1);
end
OD2_att_mode.OD2_mode_bd_norm = OD2_att_mode.OD2_mode./beads_1um_mode_OD2;

ModeTable = join(cc_mode, OD2_att_mode, 'Keys', 'files');

%%
% p = fit(log10(ModeTable.OD2_mode_bd_norm), log10(ModeTable.cc_mode_vol), 'poly1'); % model fit
% 
% figure
% plot(p, log10(ModeTable.OD2_mode_bd_norm), log10(ModeTable.cc_mode_vol), 'predfunc');
% ylabel('log10 [ Coulter Counter Volume (\mum^3)]')
% xlabel(['log10 [Attune ' attune_par ' (OD2), 1 \mum bead normalized]'])
% eqstr = ['y =  ' num2str(p.p1) '*x +' num2str(p.p2)];
% text(-.5,3, eqstr, 'color', 'r')

%%
%No OD case
str = '\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Gated_FCS_March2026\';
files = dir([str '\\NoOD*\\*.fcs']);
fcsNoODT = array2table([{files.name}' {files.folder}'], 'VariableNames',{'file' 'folder'});
runDate = split(fcsNoODT.folder,'_'); 
fcsNoODT.run_date = runDate(:,end); clear runDate
for ii = 1:length(files)
    [fcsdat, fcshdr] = fca_readfcs(fullfile(files(ii).folder, files(ii).name));
    fcsNoODT.dat{ii} = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
end

%%
figure(97)
set(gcf, 'WindowState','maximized')
tlNoOD = tiledlayout(4,5);
figure(96)
tlcyto2 = tiledlayout(4,5);
NoOD_AttuneSizeCal = table;
for i = 1:length(files)
    dat = fcsNoODT.dat{i}.(attune_par);
    [N, edges] = histcounts(dat,300);
    att_val = edges(1:end-1) + diff(edges)/2;
    varname = regexprep(fcsNoODT.file{i}, '.fcs', '');
    if startsWith(varname, 'beads')
        varname = strcat(varname,'_', fcsNoODT.run_date{i}); 
    end
    figure(97)
    nexttile
    plot(att_val, N./diff(edges)), title(varname, 'Interpreter','none')
    hold on
    %this is how BF did the histograms
    [N, edges] = histcounts(dat);
    %att_val = edges(1:end-1) + diff(edges)/2;
    %plot(att_val,N./diff(edges), 'g', 'LineWidth',2)
    %now let's force more resolution around the peak (impt for Syn)
    [att_max, attind] = max(N);
    std_temp = std(dat(dat<edges(attind+1)*3));
    dat_upper = edges(attind+1)+std_temp*3;
    dat_lower = edges(attind+1)-std_temp*3;
    xlim([-inf min(edges(attind+1)+std_temp*3, max(dat))])
    xl = xlim; xl(1) = min([xl(1) -1e4]); xlim(xl), clear xl
    [N, edges] = histcounts(dat(dat<dat_upper & dat>dat_lower & dat>attune_min),attune_bin_num); %dat>-1e4
 %   N = smooth(medfilt1(N,attune_medfilt_num),attune_smooth_num)';
   if startsWith(varname,'SKIPbeads') & strcmp(attune_par, 'SSC-A') %only for SSC-A and No OD
        %N = smooth(medfilt1(N,2),3)';
        N = medfilt1(N,3);
   elseif strcmp(varname, 'beads_mcl') & strcmp(attune_par, 'FSC-H')
        N = medfilt1(N,3);
   else
        N = smooth(medfilt1(N,attune_medfilt_num),attune_smooth_num)';
   end
   att_val = edges(1:end-1) + diff(edges)/2;
   plot(att_val,N./diff(edges), 'r', 'LineWidth',2)
   NoOD_AttuneSizeCal.(varname) = {[att_val' N']}; %save this final one
   grid
   figure(96)
   nexttile
   loglog(dat, fcsNoODT.dat{i}.('BL3-A'), '.'), title(varname, 'Interpreter','none')
end

xlabel(tlcyto2, attune_par)
ylabel(tlcyto2, 'CHL-A')
NoOD_AttuneSizeCal = renamevars(NoOD_AttuneSizeCal, 'HTriq', 'Htriq'); %match cc
xlabel(tlNoOD,['Attune ' attune_par ', No OD on SSC'])
ylabel(tlNoOD,'Number per X-unit')
figure(97)
if printNow
    print(['\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Analysis_March2026\attune_' attune_par '_noOD_hist.png'], '-dpng')
end

[bmax, bind] = max(NoOD_AttuneSizeCal.beads_1um_9Apr2019{1}(:,2));

beads_1um_mode_NoOD = NoOD_AttuneSizeCal.beads_1um_9Apr2019{1}(bind,1); 

NoOD_att_mode = table;
NoOD_att_mode.files = NoOD_AttuneSizeCal.Properties.VariableNames';
NoOD_att_mode.file_orig = NoOD_att_mode.files;
temp = startsWith(NoOD_att_mode.files, 'Syn');
NoOD_att_mode.files(temp) = {'20Mar_Syn'}; 
temp = startsWith(NoOD_att_mode.files, '1820i');
NoOD_att_mode.files(temp) = {'Syn'}; %force this one to match with coulter counter

for i = 1:height(NoOD_att_mode)
    att = NoOD_AttuneSizeCal.(NoOD_att_mode.file_orig{i}){1};
    [att_max, attind] = max(att(:,2));
    NoOD_att_mode.noOD_mode(i) = att(attind,1);
end
NoOD_att_mode.noOD_mode_bd_norm = NoOD_att_mode.noOD_mode./beads_1um_mode_NoOD;

ModeTable = join(ModeTable, NoOD_att_mode, 'Keys', 'files');
%%
figure
loglog(ModeTable.OD2_mode_bd_norm, ModeTable.noOD_mode_bd_norm, '.', 'MarkerSize',10)
xlabel(['log10 [Attune ' attune_par ' (OD2), 1 \mum bead norm]'])
ylabel(['log10 [Attune ' attune_par ' (No OD), 1 \mum bead norm]'])
%axis([.05 120 .05 120])
line(xlim, xlim)
hold on
text(ModeTable.OD2_mode_bd_norm, ModeTable.noOD_mode_bd_norm, ModeTable.files, 'fontsize',8,'interpreter', 'none', 'VerticalAlignment','top')
grid
if printNow
    print(['\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Analysis_March2026\attune_' attune_par '_noODvOD2.png'], '-dpng')
end
%%
ModeTable.merge_mode_bd_norm = ModeTable.OD2_mode_bd_norm;
ind = ModeTable.OD2_mode_bd_norm<1;
ModeTable.merge_mode_bd_norm(ind) = ModeTable.noOD_mode_bd_norm(ind);
%%
%load \\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\FitResults_saved_SSC-A_22Mar2026
[p, stats] = fit(log10(ModeTable.merge_mode_bd_norm), log10(ModeTable.cc_mode_vol), 'poly1'); % model fit
figure
plot(p, log10(ModeTable.merge_mode_bd_norm), log10(ModeTable.cc_mode_vol), 'predfunc');
ylabel('log10 [ Coulter Counter Volume (\mum^3)]')
xlabel(['log10 [Attune ' attune_par ' (merge), 1 \mum bead normalized]'])
eqstr = ['y =  ' num2str(p.p1) '*x +' num2str(p.p2)];
text(-1,3, eqstr, 'color', 'r')
text(-1,2.7, ['r^2 = ' num2str(stats.rsquare)])
text(log10(ModeTable.merge_mode_bd_norm), log10(ModeTable.cc_mode_vol), ModeTable.files, 'fontsize',8,'interpreter', 'none', 'VerticalAlignment','top')
axis([-1.5 2.5 -.5 4])
legend('Location', 'northwest')
grid
%%
if printNow
    print(['\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Analysis_March2026\vol_cal_' attune_par '_merge_' datestr(now, 'ddmmmyyyy')], '-dpng')
end
if saveNow
    notes = {'Created with size_calibration_fit_march2019_v2.m'}
    fout = ['\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\FitResults_saved_' attune_par '_' datestr(now, 'ddmmmyyyy')]; 
    save(fout, 'OD2*', 'NoOD*', 'CC*', 'cc_mode', 'p', 'stats', 'ModeTable', 'beads_1um_mode*', 'notes', 'attune_par')
    disp('Results saved:')
    disp(fout)
end

return

%%
% here plot vol size dists cc and attune
figure
%pstr = 'Pavlova'; s = 20;
%pstr = 'Prymn_pav'; s = 40;
%pstr = 'Dun'; s = 6;
%pstr = 'Nano'; s = 60;
pstr = 'Micromonas'; s = 4;
p.p1
plot(4/3*pi*(CCsizecal.(pstr)(:,1)/2).^3, CCsizecal.(pstr)(:,2), '.-')
hold on
plot(10.^(log10(NoOD_AttuneSizeCal.(pstr){1}(:,1)./beads_1um_mode_NoOD)*p.p1+p.p2), NoOD_AttuneSizeCal.(pstr){1}(:,2)*s, '.-')
plot(10.^(log10(OD2_AttuneSizeCal.(pstr){1}(:,1)./beads_1um_mode_OD2)*p.p1+p.p2), OD2_AttuneSizeCal.(pstr){1}(:,2)*s, '.-')
%plot(10.^(log10(Neutral_AttuneSizeCal.Dun{1}(:,1)./beads_1um_mode_NoOD)*1.1944+1.2049), Neutral_AttuneSizeCal.Dun{1}(:,2)*6, '.-')
%plot(10.^(log10(OD2_AttuneSizeCal.Dun{1}(:,1)./beads_1um_mode_OD2)*1.1944+1.2049), OD2_AttuneSizeCal.Dun{1}(:,2)*6, '.-')
legend('Coulter', 'Attune No OD', 'Attune OD2')
title(pstr)
%%
% Here is a plot to compare new SSC-A cal with old BF cal
new = load('\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\FitResults_saved_SSC-A_22Mar2026');
BF = load('\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\BF_results_saved');
%%
figure
plot(new.p, 'r', log10(new.ModeTable.merge_mode_bd_norm), log10(new.ModeTable.cc_mode_vol),'r', 'predfunc');
eqstr = ['y =  ' num2str(new.p.p1) '*x +' num2str(new.p.p2)];
text(-1,2.5, eqstr, 'color', 'r')
%text(-1,2.7, ['r^2 = ' num2str(new.stats.rsquare)])
text(log10(new.ModeTable.merge_mode_bd_norm), log10(new.ModeTable.cc_mode_vol), new.ModeTable.files, 'fontsize',8,'interpreter', 'none', 'VerticalAlignment','top')
axis([-1.5 2.5 -.5 4])
legend('Location', 'northwest')
grid
hold on
plot(BF.p,'b', log10(BF.OD2_att_mode), log10(BF.cc_mode), '+b')
eqstr = ['y =  ' num2str(BF.p.p1) '*x +' num2str(BF.p.p2)];
text(.5,.5, eqstr, 'color', 'b')
ylabel('log10 [ Coulter Counter Volume (\mum^3)]')
xlabel(['log10 [Attune ' new.attune_par ' (merge), 1 \mum bead normalized]'])
lh = legend;
lh.String([1:2,4:5]) = {'New modes' 'new fit' 'Old modes' 'old fit (BF)'};

%print('\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Analysis_March2026\vol_cal_SSC-A_merge_withBF', '-dpng')

return
%%
%%%%%
%Below here is how to get exactly the calibration line that Bethany used

load('\\sosiknas1\Lab_data\Attune\Size_calibration_March2019\CCsizecal_table.mat')
%Heidi, Mar 2026, this is diameter in col1 and freq in col2, collated from 'size_calibration_data_March2019.mat'
%CCsize_threshold values are diameters and start with Nano (e.g., Syn and Micromonas not included)
CCsize_threshold = [0 0 CCsize_threshold];

%%
figure
for i = 1:width(CCsizecal)
    cc = CCsizecal{:,i};
    cc = cc(cc(:,1) > CCsize_threshold(i),:); % remove CC junk below size threshold
    cc = 4/3*pi*(cc/2).^3;
    %cc(:,1) = 4/3*pi*(cc(:,1)/2).^3; %Heidi: the line above is BF's but it should be this one (very small diff)
    cc_smooth = smooth(cc(:,2), 10);
    [cc_max, ccind] = max(cc_smooth);
    cc_mode(i) = cc(ccind,1);
    subplot(3,4,i)
    plot(cc(:,1), cc(:,2), cc(:,1), cc_smooth); hold on;
    %set(gca, 'XScale', 'log');
    line([cc_mode(i) cc_mode(i)], [0, cc_max], 'color', 'k'); hold off;
end
%%
% now od2

str = '\\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\Gated_FCS\OD2_Filter\';
files = dir('\\sosiknas1\Lab_data\Attune\Size_calibration_March2019\Attune\Gated_FCS\OD2_Filter\*.fcs');

for i = 1:length(files)
    [fcsdat, fcshdr] = fca_readfcs(strcat(str, files(i).name));
    dat = fcsdat(:,3);
    [N, edges] = histcounts(dat);
    diam = edges(1:end-1) + diff(edges)/2;
    varname = regexprep(files(i).name, '.fcs', '');
    assignin('base', varname, [diam' N']);
end

OD2_AttuneSizeCal = {Syn, Micromonas, Nano, Pavlova, I_galbana, Prymn_pav, Dun, Crypto_spmc, HTriq, H_akashiwo, Ambopo, Pptp02};

[bmax, bind] = max(beads_1um(:,2)); 
beads_1um_mode = beads_1um(bind, 1); 

for i = 1:length(files)
    [fcsdat, fcshdr] = fca_readfcs(strcat(str, files(i).name));
    dat = fcsdat(:,3);
    scal_ssca = dat./beads_1um_mode; 
    [N, edges] = histcounts(scal_ssca);
    diam = edges(1:end-1) + diff(edges)/2;
    varname = regexprep(files(i).name, '.fcs', '');
    assignin('base', varname, [diam' N']);
end

OD2_AttuneSizeCal_scaled = {Syn, Micromonas, Nano, Pavlova, I_galbana, Prymn_pav, Dun, Crypto_spmc, HTriq, H_akashiwo, Ambopo, Pptp02};

for i = 1:length(OD2_AttuneSizeCal_scaled)
    att = OD2_AttuneSizeCal_scaled{i};
    [att_max, attind] = max(att(:,2));
    OD2_att_mode(i) = att(attind,1);
end

p = fit(log10(OD2_att_mode)', log10(cc_mode)', 'poly1'); % model fit

figure
plot(p, log10(OD2_att_mode), log10(cc_mode), 'predfunc');
ylabel('log Coulter Counter Volume (\mum^3)')
xlabel('log Scaled Attune SSC-A OD2')
eqstr = ['y =  ' num2str(p.p1) '*x +' num2str(p.p2)];
text(-.5,3, eqstr, 'color', 'r')

%save('\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\BF_results_saved', 'OD2*', 'CC*', 'cc_mode', 'p')

%BF = load('\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\BF_results_saved')