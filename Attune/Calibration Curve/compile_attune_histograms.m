addpath('C:\GitHubRepositories\NES-LTER\Attune')
printNow = 0;
saveNow = 1;
attune_par = {'SSC-A';'GL1-A'; 'FSC-A'; 'SSC-H';'GL1-H';'FSC-H'}; attune_min = -200;
%attune_par = 'GL1-A'

attune_bin_num = 120;
attune_smooth_num = 15;
attune_medfilt_num = 6;
%str = '\\sosiknas1\lab_data\Attune\Size_calibration_March2019\Attune\Gated_FCS_March2026\';
str = '\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Attune\Gated_FCS\';
%str = '\\sosiknas1\lab_data\Attune\size_calibration_July2026\Attune\Gated_FCS\settings_diatom_cocco\';
%files = dir([str '\\OD2*\\*.fcs']);
files = dir([str '*.fcs']);
fcsT = array2table([{files.name}' {files.folder}'], 'VariableNames',{'file' 'folder'});
fcsT.runDate = extractBefore(fcsT.file,'_');
for ii = 1:length(files)
    [fcsdat, fcshdr] = fca_readfcs(fullfile(files(ii).folder, files(ii).name));
    fcsT.dat{ii} = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
end
%%
ibd = contains(fcsT.file, 'bead');
icell = find(~ibd); 
if 0
figure
for ii = 1:length(icell)
    loglog(fcsT.dat{icell(ii)}.("SSC-A"),fcsT.dat{icell(ii)}.("GL1-A"), '.', MarkerSize=2)
    %loglog(fcsT.dat{icell(ii)}.("SSC-H"),fcsT.dat{icell(ii)}.("GL1-H"), '.', MarkerSize=2)
    %loglog(fcsT.dat{icell(ii)}.("FSC-H"),fcsT.dat{icell(ii)}.("GL1-H"), '.', MarkerSize=2)
    hold on, axis([1e3 1e6 1e1 1e6]), grid on,
    legend(fcsT.file(icell), 'location', 'northwest')
    pause
end
end
%%
ibd = find(contains(fcsT.file, 'bead'));
if 0
figure
for ii = 1:length(ibd)
    loglog(fcsT.dat{ibd(ii)}.("SSC-H"),fcsT.dat{ibd(ii)}.("SSC-A"), '.', MarkerSize=2)
    %loglog(fcsT.dat{icell(ii)}.("SSC-H"),fcsT.dat{icell(ii)}.("GL1-H"), '.', MarkerSize=2)
    %loglog(fcsT.dat{icell(ii)}.("FSC-H"),fcsT.dat{icell(ii)}.("GL1-H"), '.', MarkerSize=2)
    hold on, axis([1e3 1e6 1e1 1e6]), grid on,
    legend(fcsT.file(ibd), 'location', 'eastoutside', 'Interpreter','none')
    pause
end
end

%%
%let's pool a few culture runs
warning off
n = height(fcsT);
tt = '28-Jul-2026_T_weissflogii';
ind = find(contains(fcsT.file,tt));
fcsT.file{n+1} = [tt '_pooled']; fcsT.runDate{n+1} = fcsT.runDate{ind(1)};
fcsT.dat{n+1} = cat(1,fcsT.dat{ind}); 
tt = '28-Jul-2026_Nano';
ind = find(contains(fcsT.file,tt));
fcsT.file{n+2} = [tt '_pooled']; fcsT.runDate{n+2} = fcsT.runDate{ind(1)};
fcsT.dat{n+2} = cat(1,fcsT.dat{ind}); 

%%
ibd = contains(fcsT.file, 'bead');
tempT = table;
AttuneSizeCal_smoothhists = table;
for ip = 1:length(attune_par)
    f1 = ip*2-1;
    f2 = ip*2;
    figure(f1)
    set(gcf, 'WindowState','maximized')
    tlod(f1) = tiledlayout(6,6);
    figure(f2)
    set(gcf, 'WindowState','maximized')
    tlod(f2) = tiledlayout(4,6);
    for i = 1:height(fcsT)
        dat = fcsT.dat{i}.(attune_par{ip});
        [N, edges] = histcounts(dat(dat>attune_min & dat<prctile(dat(:),95)*2),300);
        att_val = edges(1:end-1) + diff(edges)/2;
        vname = regexprep(fcsT.file{i}, '.fcs', '');
        if startsWith(vname, 'beads')
            vname = strcat(vname,'_', fcsT.run_date{i});
        end
        if ibd(i)
            figure(f1)
        else
            figure(f2)
        end
        nexttile
        plot(att_val, N./diff(edges)), title(vname, 'Interpreter','none', 'FontSize',6)
        hold on
        %this is how BF did the histograms
        [N, edges] = histcounts(dat);
        %now let's force more resolution around the peak (impt for Syn)
        [att_max, attind] = max(N);
        std_temp = std(dat(dat<edges(attind+1)*3));
        dat_upper = edges(attind+1)+std_temp*3;
        dat_lower = edges(attind)-std_temp*3;
        xlim([-inf min(edges(attind+1)+std_temp*3, max(dat))])
        xl = xlim; xl(1) = min([xl(1) -1e4]); xlim(xl), clear xl
        [N, edges] = histcounts(dat(dat<dat_upper & dat>dat_lower & dat>attune_min),attune_bin_num); %dat>-1e4
        if contains(vname, 'syn') | contains(vname, 'Nano')% & strcmp(attune_par, 'FSC-H')
            N = smooth(medfilt1(N,attune_medfilt_num/2),attune_smooth_num/2)';
        elseif contains(vname, 'McLane') | contains(vname, '1_um')
            N = smooth(medfilt1(N,3),4)';
        elseif contains(vname, "T_weiss") & contains(attune_par{ip}, 'GL1-A')
            N = smooth(medfilt1(N,attune_medfilt_num*4),attune_smooth_num*4)';
        else
            N = smooth(medfilt1(N,attune_medfilt_num*2),attune_smooth_num)';
        end
        att_val = edges(1:end-1) + diff(edges)/2;
        plot(att_val,N./diff(edges), 'r', 'LineWidth',2)
        tempT.(vname) = [att_val' N']; %save this final one
        grid
        %figure(98)
        %nexttile
        %loglog(dat, fcsT.dat{i}.('BL3-A'), '.'), title(varname, 'Interpreter','none')
    end
    AttuneSizeCal_smoothhists.(attune_par{ip}) = tempT;
    figure(f1)
    xlabel(tlod(f1),['Attune ' attune_par{ip}])
    ylabel(tlod(f1),'Number per X-unit')
    figure(f2)
    xlabel(tlod(f2),['Attune ' attune_par{ip}])
    ylabel(tlod(f2),'Number per X-unit')
    if printNow
        figure(f1)
        print(['\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Attune\plots\attune_' attune_par{ip} '_bead_hist.png'], '-dpng')
        figure(f2)
        print(['\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Attune\plots\attune_' attune_par{ip} '_culture_hist.png'], '-dpng')
    end
end

%%
par = AttuneSizeCal_smoothhists.Properties.VariableNames;
files = AttuneSizeCal_smoothhists{:,1}.Properties.VariableNames;
AttuneMode = table;
for ip = 1:length(par)
    for ii = 1:length(files)
        att = AttuneSizeCal_smoothhists.(par{ip}).(files{ii});
        [att_max, attind] = max(att(:,2));
        AttuneMode.(par{ip})(ii) = att(attind,1);
    end
end
AttuneMode.Row = files';
AttuneMode.run_date = datetime(extractBefore(AttuneMode.Row,'_'));

%%
t = AttuneMode.Row';
labelstr(contains(t, '1_um')) = "1_um";
%labelstr(contains(t, '1um')) = "1_um";
labelstr(contains(t,'F8819')) = "F8819";
labelstr(contains(t, '0.5_um')) = "0.5_um";
labelstr(contains(t, 'McLane')) = "6_um";
labelstr(contains(t, 'Dun')) = "Dun";
labelstr(contains(t, 'Micromonas')) = "Micromonas";
labelstr(contains(t, 'Micromonas(2)')) = "C2Micromonas"; % first day, CC large aperature, probably skip later
labelstr(contains(t, 'syn_7335p')) = "Syn7335";
labelstr(contains(t, 'syn_8113p')) = "Syn8113";
labelstr(contains(t, 'CRYPTO')) = "Crypto";
labelstr(contains(t, 'Heterocapsa')) = "Heterocapsa";
labelstr(contains(t, 'Isochrysis')) = "Isochrysis";
labelstr(contains(t, 'syn_WH8018')) = "Syn_WH8018";
labelstr(contains(t, 'I_galbana')) = "I_galbana";
labelstr(contains(t, 'Nano')) = "Nano";
labelstr(contains(t, 'Pavlova')) = "Pavlova";
labelstr(contains(t, 'Rhodomonas')) = "Rhodomonas";
labelstr(contains(t, 'SE62')) = "SE62";
labelstr(contains(t, 'T_weissflogii')) = "T_weissflogii";
AttuneMode.label = labelstr';

%%
if saveNow
    notes = {'Created with compile_attune_histograms.m'};
    fout = ['\\sosiknas1\lab_data\Attune\Size_calibration_July2026\AttuneSummary_saved_' datestr(now, 'ddmmmyyyy')]; 
    save(fout, 'AttuneMode', 'AttuneSizeCal_smoothhists', 'notes', 'fcsT')
    disp('Results saved:')
    disp(fout)
end

%%
%make some bead mode overview plots

ibd = contains(AttuneMode.Row, 'bead');
iculture = (~ibd);
i1um = contains(AttuneMode.label, '1_um');
iF8819 = contains(AttuneMode.label, 'F8819');
ihum = contains(AttuneMode.label, '0.5_um');
i6um = contains(AttuneMode.label, '6_um');

figure, tiledlayout(3,1, 'TileSpacing','none');
lstr = {{'SSC-H' 'SSC-A'} {'GL1-H' 'GL1-A'} {'FSC-H' 'FSC-A'}};
yl = {[4e4 6e4] [800 1200] [1.3e4 2e4]};
for ii = 1:length(lstr)
    nexttile
    bar([AttuneMode.Row(i1um); "_"; AttuneMode.Row(iF8819)], [AttuneMode{i1um,lstr{ii}}; [NaN NaN]; AttuneMode{iF8819,lstr{ii}}])
    grid on
    lh = legend(lstr{ii},'location', 'southeast');    
    if ii<length(lstr)
        set(gca, 'XTickLabel', 'none')
    end
    ylim(yl{ii})
end
set(gca, 'TickLabelInterpreter', 'none');
if printNow
    print(['\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Attune\plots\1_um_F8819.png'], '-dpng')
end
%%
figure, tiledlayout(3,1, 'TileSpacing','none');
lstr = {{'SSC-H' 'SSC-A'} {'GL1-H' 'GL1-A'} {'FSC-H' 'FSC-A'}};
%ix = i1um; yl = {[4e4 6e4] [800 1200] [1.3e4 2e4]}; istr = '1_um';
ix = ihum; yl = {[2.7e3 1e4] [150 225] [1e3 2e3]}; istr = '0.5_um';
%ix = i6um; yl = {[3.5e5 5.5e5] [1.5e4 3e4] [2.5e5 3.5e5]}; istr = '6_um';

for ii = 1:length(lstr)
    nexttile
    bar(AttuneMode.Row(ix), AttuneMode{ix,lstr{ii}})
    grid on
    lh = legend(lstr{ii},'location', 'southeast');    
    if ii<length(lstr)
        set(gca, 'XTickLabel', 'none')
    end
    ylim(yl{ii})
end
set(gca, 'TickLabelInterpreter', 'none');
if printNow
    print(['\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Attune\plots\' istr '.png'], '-dpng')
end
