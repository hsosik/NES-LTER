load \\sosiknas1\lab_data\Attune\Size_calibration_July2026\CCSummary_saved_10Sep2026 %from compile_cc_histograms.m
load \\sosiknas1\lab_data\Attune\Size_calibration_July2026\AttuneSummary_saved_11Sep2026  %from compile_attune_histograms.m
printNow = 1;
saveNow = 1;

%%
ibd = contains(AttuneMode.Row, 'bead');
iculture = (~ibd);
i1um = contains(AttuneMode.label, '1_um');
iF8819 = contains(AttuneMode.label, 'F8819');
ihum = contains(AttuneMode.label, '0.5_um');
i6um = contains(AttuneMode.label, '6_um');

ModeTable= AttuneMode(iculture,:);
beadT = AttuneMode(i1um,:); %1 micron beads
beadT = removevars(beadT, { 'label'});
beadT_day = groupsummary(beadT, 'run_date', {'mean' 'std'});
%%
[~,ib] = ismember(ModeTable.run_date, beadT_day.run_date);
bead_match = beadT_day(ib,:);

%%
figure, tiledlayout(1,2,"TileSpacing","compact");
nexttile
loglog(ModeTable.("GL1-H")./bead_match.("mean_GL1-H"), ModeTable.("SSC-H")./bead_match.("mean_SSC-H"), '*')
ylabel('SSC-H, bd unit'), xlabel('GL1-H, bd unit')
line(xlim, xlim), grid on
axis square
nexttile
loglog(ModeTable.("GL1-A")./bead_match.("mean_GL1-A"), ModeTable.("SSC-A")./bead_match.("mean_SSC-A"), '*')
ylabel('SSC-A, bd unit'), xlabel('GL1-A, bd unit')
line(xlim, xlim), grid on
axis square
if printNow
    print('\\sosiknas1\lab_data\Attune\Size_calibration_July2026\Attune\SSC_vs_GL1_beadunit.png', '-dpng')
end

%%
figure, tiledlayout(1,2,"TileSpacing","compact");
nexttile
loglog(ModeTable.("GL1-H"), ModeTable.("SSC-H"), '*')
hold on
loglog(AttuneMode.("GL1-A")(ibd), AttuneMode.("SSC-A")(ibd), '^')
ylabel('SSC-H'), xlabel('GL1-H')
line(xlim, xlim*50), grid on
axis square
legend('culture', 'bead', '50:1', 'Location','northwest')
nexttile
loglog(ModeTable.("GL1-A"), ModeTable.("SSC-A"), '*')
hold on
loglog(AttuneMode.("GL1-A")(ibd), AttuneMode.("SSC-A")(ibd), '^')
ylabel('SSC-A'), xlabel('GL1-A')
line(xlim, xlim*50), grid on
axis square
legend('culture', 'bead', '50:1', 'Location','northwest')

if printNow
    print('\\sosiknas1\lab_data\Attune\Size_calibration_July2026\Attune\SSC_vs_GL1_raw.png', '-dpng')
end
%%
[ii,ind] = ismember(ModeTable.label, CCMode.Row);
ModeTable.cc_mode_vol(ii) = 4/3*pi*(CCMode.diameter_mode(ind(ind>0))./2).^3;
ModeTable.cc_mode_vol(~ii) = NaN;
ModeTable.cc_mode_diameter(ii) = CCMode.diameter_mode(ind(ind>0));
ModeTable.cc_mode_diameter(~ii) = NaN;

%%

figure, tiledlayout(1,3, 'TileSpacing','compact');
nexttile
loglog(ModeTable.("SSC-H")./bead_match.("mean_SSC-H"), ModeTable.cc_mode_vol,'*')
grid on, axis square

hold on
loglog(ModeTable.("GL1-H")./bead_match.("mean_GL1-H"), ModeTable.cc_mode_vol, '*')
xlabel('SSC-H or GL1-H, bd unit'), ylabel('Coulter volume, \mum^3')
text(ModeTable.("GL1-H")./bead_match.("mean_GL1-H"), ModeTable.cc_mode_vol, ModeTable.label, 'FontSize',6, 'HorizontalAlignment','right', 'VerticalAlignment','bottom')

nexttile
loglog(ModeTable.("SSC-A")./bead_match.("mean_SSC-A"), ModeTable.cc_mode_vol, '*')
grid on, axis square
hold on
loglog(ModeTable.("GL1-A")./bead_match.("mean_GL1-A"), ModeTable.cc_mode_vol, '*')
xlabel('SSC-A or GL1-A, bd unit'), ylabel('Coulter volume, \mum^3')

nexttile
loglog(ModeTable.("FSC-H")./bead_match.("mean_FSC-H"), ModeTable.cc_mode_vol, '*')
grid on, axis square
xlabel('FSC-H, bd unit'), ylabel('Coulter volume, \mum^3')

%%
%attune_par = 'SSC-H'; attune_par2 = 'GL1-H';
attune_par = 'SSC-A'; attune_par2 = 'GL1-A';
%attune_par = 'FSC-H'; attune_par2 = '';
n = sum(iculture)+1; 
ModeTable(n:end,:) = []; %remove means from below if already present
ModeTable.merge_mode_bd_norm = ModeTable.(attune_par)./bead_match.(['mean_' attune_par]);
if ~isempty(attune_par2)
    ind = find(ModeTable.merge_mode_bd_norm>1);
    ModeTable.merge_mode_bd_norm(ind) = ModeTable.(attune_par2)(ind)./bead_match.(['mean_' attune_par2])(ind);
else %FSC case
    ind = find(ModeTable.merge_mode_bd_norm>5);
    ModeTable.merge_mode_bd_norm(ind) = NaN; %only fit to the low values
end

% best not to have two similar Isochrysis galbana points, so average them before fitting
ind = contains(ModeTable.label, {'Iso' 'I_gal'});
ModeTable{n,:} = missing;
ModeTable.Row(n) = {'I_galbana_mean'}; ModeTable.label(end) = "Isochrysis";
ModeTable(n,{'cc_mode_diameter' 'cc_mode_vol' 'merge_mode_bd_norm' }) = mean(ModeTable(ind,{'cc_mode_diameter' 'cc_mode_vol' 'merge_mode_bd_norm' }))

%ind2use = ~isnan(ModeTable.merge_mode_bd_norm) & ~isnan(ModeTable.cc_mode_vol) & ~contains(ModeTable.label, 'C2Micro') & ~contains(ModeTable.Row, 'T_weissflogii_CCMP') & ~(contains(ModeTable.Row, 'Nano') & ~contains(ModeTable.Row, 'pooled'));
ind2use = ~isnan(ModeTable.merge_mode_bd_norm) & ~isnan(ModeTable.cc_mode_vol) & ~contains(ModeTable.Row, {'Micromonas(2)' 'T_weissflogii_CCMP' 'Iso' '28-Jul-2026_I_gal'}) & ~(contains(ModeTable.Row, 'Nano') & ~contains(ModeTable.Row, 'pooled'));

[p, stats] = fit(log10(ModeTable.merge_mode_bd_norm(ind2use)), log10(ModeTable.cc_mode_vol(ind2use)), 'poly1'); % model fit

figure
plot(p, log10(ModeTable.merge_mode_bd_norm(ind2use)), log10(ModeTable.cc_mode_vol(ind2use)), 'predfunc');
ylabel('log10 [ Coulter Counter Volume (\mum^3)]')
if ~isempty(attune_par2)
    xlabel(['log10 [Attune ' attune_par ' (merge), 1 \mum bead normalized]'])
else
    xlabel(['log10 [Attune ' attune_par ', 1 \mum bead normalized]'])
end
eqstr = ['y =  ' num2str(p.p1) '*x +' num2str(p.p2)];
text(-1,3, eqstr, 'color', 'r')
text(-1,2.7, ['r^2 = ' num2str(stats.rsquare)])
text(log10(ModeTable.merge_mode_bd_norm(ind2use)), log10(ModeTable.cc_mode_vol(ind2use)), ModeTable.label(ind2use), 'fontsize',8,'interpreter', 'none', 'VerticalAlignment','top')
axis([-1.5 2.5 -.5 4])
legend('Location', 'northwest')
grid
axis square

if printNow
    print(['\\sosiknas1\lab_data\Attune\Size_calibration_July2026\Attune\vol_cal_' attune_par '_merge_' datestr(now, 'ddmmmyyyy')], '-dpng')
end
if saveNow
    notes = {'Created with size_calibration_fit_jul2026.m'}
    fout = ['\\sosiknas1\lab_data\Attune\Size_calibration_July2026\Attune\FitResults_saved_' attune_par '_' datestr(now, 'ddmmmyyyy')]; 
    save(fout, 'p', 'stats', 'ModeTable', 'notes', 'attune_par', 'ind2use', 'bead_match')
    disp('Results saved:')
    disp(fout)
end
