function [ ] = scattering_merge(p, saverpath, stepsize4plot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%fcslist = dir([p.fpath, '*.fcs']);
fcslist = dir([p.classpath '*.mat']); %just get the ones in the class folder (not the beads, etc.)
fcslist = regexprep({fcslist.name}', '.mat', '.fcs');
mergeT = table;
mergeT.filename = fcslist;
%%
sscstr = strcat("SSC-", p.SSCDIM);
gl1str = strcat("GL1-", p.SSCDIM);
warning off

gl1min = 600;
gl1max = 1600;
sscmax = 1e5;
for ii = 1:length(fcslist) %3
    filename = [p.fpath, mergeT.filename{ii}];
    [fcsdat,fcshdr] = fca_readfcs(filename);
    fcsdat = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
    file_hv = array2table([NaN [fcshdr.par.hv]], 'VariableNames', {fcshdr.par.name});
    mergeT.SSC_hv(ii) = file_hv.(sscstr);
    mergeT.GL1_hv(ii) = file_hv.(gl1str);
    mergeT.filetime(ii) = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    %temp = fcsdat.("GL1-A")>500 & fcsdat.("GL1-A")<1500;
    temp = fcsdat.(gl1str)>gl1min & fcsdat.(gl1str)<gl1max & fcsdat.(sscstr)<sscmax;

    X = fcsdat(temp,[gl1str sscstr ]);
    mergeT.number_pts(ii) = height(X);
    if mergeT.number_pts(ii) < 10 %this might break if the first file doesn't have enough points
        mergeT{ii,setdiff(mergeT.Properties.VariableNames,{'filename' 'filetime' 'GL1_hv' 'SSC_hv'})} = NaN;
    else
        %m = fitlm(X.("GL1-A"),X.("SSC-A"),'Intercept',false);
        [fitm, fitstat] = fit(X.(gl1str),X.(sscstr),fittype('a*x'), 'start', [50]);
        mergeT.slope_median(ii) = median(X.(sscstr)./X.(gl1str));
        mergeT.rmse_median(ii) = rmse(X.(gl1str).*mergeT.slope_median(ii), X.(sscstr));
        temp = corrcoef(X.(gl1str).*mergeT.slope_median(ii), X.(sscstr));
        mergeT.r2_median(ii) = temp(2)^2; clear temp
        if ~rem(ii,stepsize4plot)
            disp(['merging ' num2str(ii) ' of ' num2str(length(fcslist))])

            figure(100),clf
            set(gcf, 'Position',[100 100 1000 300])
            tl = tiledlayout(1,3);
            nexttile
            plot(fcsdat.(gl1str), fcsdat.(sscstr), '.')
            grid on
            ylabel(strcat(sscstr,", No OD")), xlabel(strcat(gl1str," , OD2"))
            axis([-1000 12e5 -1000  12e5])
            nexttile
            plot(fcsdat.(gl1str), fcsdat.(sscstr), '.')
            grid on, hold on
            plot(X{:,1}, X{:,2}, '.r')
            set(gca, 'yscale', 'log', 'xscale', 'log')
            %axis([-500 10000 -1e4 3e5])
            axis([10 20000 1e3 10e5])
            %plot(fitm, xlim,ylim)
            plot(xlim,xlim.*mergeT.slope_median(ii), 'LineWidth',2, 'color', 'm')
            %legend('Location', 'southeast', 'FontSize',6)
            ylabel(strcat(sscstr,", No OD")), xlabel(strcat(gl1str," , OD2"))
            nexttile
            plot(fcsdat.(gl1str), fcsdat.(sscstr), '.')
            hold on
            plot(X{:,1}, X{:,2}, '.r')
            grid on, hold on
            axis([-500 3000 -1e4 1.5e5])
            %plot(fitm, xlim,ylim)
            plot(xlim,xlim.*mergeT.slope_median(ii), 'LineWidth',2, 'Color','m')
            %legend('Location', 'southeast', 'FontSize',6)
            ylabel(strcat(sscstr,", No OD")), xlabel(strcat(gl1str," , OD2"))
            %        title(['fit slope = ' num2str(coeffvalues(fitm),2)] )
            title(['fit slope (median) = ' num2str(mergeT.slope_median(ii),2)] )
            title(tl, fcshdr.filename, 'Interpreter','none', 'FontSize',12)
            print([saverpath regexprep(fcshdr.filename, 'fcs', 'png')],'-dpng')
        end
        mergeT.slope(ii) = coeffvalues(fitm);
        temp = confint(fitm);
        mergeT.confintLower(ii) = temp(1);
        mergeT.confintUpper(ii) = temp(2);
        mergeT.rsquare(ii) = fitstat.rsquare;
    end
end

mergeT = sortrows(mergeT, 'filetime');
save([p.outpath 'mergeTable'], 'mergeT', 'gl1min', 'gl1max', 'sscmax')
figure
yyaxis left
plot(mergeT.filetime, mergeT.slope_median, '.-')
ylabel('SSC merge slope')
grid on
yyaxis right
plot(mergeT.filetime, mergeT.r2_median, '.-')
ylabel('SSC merge r^2')

print([saverpath 'merge_stats'],'-dpng')

end
