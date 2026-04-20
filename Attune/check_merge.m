load('\\sosiknas1\Lab_data\Attune\cruise_data\beads\FCB_bead_mix_experiment_settings\between_cruises\outputs\beadstat.mat')
cruise = 'EN727';
basepath_temp =  '\\sosiknas1\Lab_data\Attune\cruise_data\';
temp = dir([basepath_temp '*' cruise]);
    p.basepath = [basepath_temp temp.name filesep 'preserved'];
p.fpath = [p.basepath filesep 'FCS' filesep];
p.outpath = [p.basepath filesep 'outputs' filesep];
p.classpath = [p.outpath 'class' filesep];
p.awspath = [p.basepath filesep 'aws\'];
%p.SSCDIM = 'A'; %needed for Step 4 & 5, SSCDIM = 'A' or 'H'


classlist = dir([p.classpath, '*.mat']); 
%%
figure
mergeT = table;
mergeT.filename = regexprep({classlist.name}, '.mat', '.fcs')';
for ii = 1:length(classlist)
    disp(ii)
    filename = [p.fpath, mergeT.filename{ii}];
    [fcsdat,fcshdr] = fca_readfcs(filename);
    fcsdat = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
    file_hv = array2table([NaN [fcshdr.par.hv]], 'VariableNames', {fcshdr.par.name});
    mergeT.SSCA_hv(ii) = file_hv.("SSC-A");
    mergeT.GL1_hv(ii) = file_hv.("GL1-A");
    mergeT.filetime(ii) = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    [alert,ind1] = min(abs(datenum(B.beadstat.time)-datenum(mergeT.filetime(ii))));
    bead_value1 = [B.beadstat.NoOD2_hv(ind1) B.beadstat.NoOD2centers(ind1,2)]; %bead value on SSC
    bead_value2 = [B.beadstat.OD2_hv(ind1) B.beadstat.OD2centers(ind1,2)]; %bead value on GL1
    bead_valueNoOD = 10.^(0.016588.*(-bead_value1(1)+file_hv.("SSC-A")) + log10(bead_value1(2)));
    bead_valueOD2 = 10.^(0.016659.*(bead_value2(1)-file_hv.("GL1-A")) + log10(bead_value2(2)));
    y = (fcsdat.("GL1-A")./bead_valueOD2-fcsdat.("SSC-A")./bead_valueNoOD)./(fcsdat.("SSC-A")./bead_valueNoOD);
    %plot(fcsdat.("SSC-A")./bead_valueNoOD, y*100, '.')
    %title(fcshdr.filename, 'Interpreter','none', 'FontSize',8)
    ylim([-100 100])
    line(xlim, [0 0], 'linewidth', 2)
    temp = fcsdat.("SSC-A")./bead_valueNoOD>.5 & fcsdat.("SSC-A")./bead_valueNoOD<1.5;
    mergeT.meandiff(ii) = mean(y(temp));
    mergeT.stddiff(ii) = std(y(temp));
    mergeT.bead_valueSSCA(ii,:) = bead_value1;
    mergeT.bead_valueGL1A(ii,:) = bead_value2;
    mergeT.bead_valueNoOD(ii) = bead_valueNoOD;
    mergeT.bead_valueOD2(ii,:) = bead_valueOD2;   
   % pause
end