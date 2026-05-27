%load('\\sosiknas1\Lab_data\Attune\cruise_data\beads\FCB_bead_mix_experiment_settings\between_cruises\outputs\beadstat.mat')
cruise = 'EN706'; %'EN706'; %'EN727';
basepath_temp =  '\\sosiknas1\Lab_data\Attune\cruise_data\';
temp = dir([basepath_temp '*' cruise]);
p.basepath = [basepath_temp temp.name filesep ];
p.fpath = [p.basepath filesep 'FCS' filesep];
p.outpath = [p.basepath filesep 'outputs' filesep];
p.classpath = [p.outpath 'class' filesep];
p.awspath = [p.basepath filesep 'aws\'];
%p.SSCDIM = 'A'; %needed for Step 4 & 5, SSCDIM = 'A' or 'H'

fcslist = dir([p.fpath, '*.fcs']); 
mergeT = table;
mergeT.filename = {fcslist.name}';
%%

for ii = 1:1:10%length(fcslist) %3
   disp(ii)
   figure(100),clf
   set(gcf, 'Position',[100 100 1000 300])
   tl = tiledlayout(1,3);
    filename = [p.fpath, mergeT.filename{ii}];
    [fcsdat,fcshdr] = fca_readfcs(filename);
    fcsdat = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
    file_hv = array2table([NaN [fcshdr.par.hv]], 'VariableNames', {fcshdr.par.name});
    mergeT.SSCA_hv(ii) = file_hv.("SSC-A");
    mergeT.GL1_hv(ii) = file_hv.("GL1-A");
    mergeT.filetime(ii) = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    %temp = fcsdat.("GL1-A")>500 & fcsdat.("GL1-A")<1500;
    temp = fcsdat.("GL1-A")>600 & fcsdat.("GL1-A")<1600;
    X = fcsdat(temp,["GL1-A" "SSC-A" ]);
    %m = fitlm(X.("GL1-A"),X.("SSC-A"),'Intercept',false);
    [fitm, fitstat] = fit(X.("GL1-A"),X.("SSC-A"),fittype('a*x'), 'start', [50]);
    nexttile
    plot(fcsdat.("GL1-A"), fcsdat.("SSC-A"), '.')
    grid on
    ylabel('SSC-A, No OD'), xlabel('GL1-A, OD2')
    axis([-1000 12e5 -1000  12e5])
    nexttile
    plot(fcsdat.("GL1-A"), fcsdat.("SSC-A"), '.')
    grid on, hold on
    axis([-500 10000 -1e4 3e5])
    plot(fitm, xlim,ylim)
    legend('Location', 'southeast', 'FontSize',6)
    ylabel('SSC-A, No OD'), xlabel('GL1-A, OD2')
    nexttile
    plot(fcsdat.("GL1-A"), fcsdat.("SSC-A"), '.')
    grid on, hold on
    axis([-500 3000 -1e4 1.5e5])
    plot(fitm, xlim,ylim)
    legend('Location', 'southeast', 'FontSize',6)
    ylabel('SSC-A, No OD'), xlabel('GL1-A, OD2')
    plot(X{:,1}, X{:,2}, '.r')
    title(['fit slope = ' num2str(coeffvalues(fitm),2)] )
    title(tl, fcshdr.filename, 'Interpreter','none', 'FontSize',12)

    mergeT.slope(ii) = coeffvalues(fitm);
    temp = confint(fitm);
    mergeT.confintLower(ii) = temp(1);
    mergeT.confintUpper(ii) = temp(2);
    mergeT.rsquare(ii) = fitstat.rsquare; 
end
save([basepath_temp cruise '_mergeTb'], 'mergeT', 'cruise')