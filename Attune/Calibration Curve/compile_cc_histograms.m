p = '\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Coulter\';

f = dir([p '*.0*']);
f(contains({f.name}, 'C47335.001')) = []; %skip a bad sample

fname = {f.name}';
run_date = dateshift(datetime([f.datenum], 'ConvertFrom', 'datenum'), 'start', 'day')
saveNow = 1;
printNow = 1;
%fname = split(fname,'.'); %remove ext
%fname = fname(:,1);

ccdat = table; 
for ii = 1:length(f)
    [dn,n] = ccreadraw_fun(p,f(ii).name);
    tt = table(dn',n, 'VariableNames',["diameter_micron", "count"]);
    ccdat.(fname{ii}) = tt;
end

%%
close all
for ii = 1:length(f)
    if rem(ii,30) == 1
        if printNow && ii > 1
            print(['\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Coulter\plots\cc_hist' num2str(gcf().Number) ' .png'], '-dpng')
        end
        figure
        tiledlayout(5,6,'TileSpacing','compact')
    end
    nexttile
    fn = fname{ii};
    x = ccdat.(fn).diameter_micron;
    plot(x, ccdat.(fn).count, '.')
    hold on
    if contains(fn, {'BD' 'FCB' 'NANO' 'C2MICRO'})
        y = smooth(medfilt1(ccdat.(fn).count,3),4);
        %y(x<1.25) = NaN; %skip noise on micro
        %s = 8;
    else
        %y = smooth(medfilt1(ccdat.(fn).count,7),3);
        y = smooth(medfilt1(ccdat.(fn).count,7),10);
        %s = 20;
    end
    y(1:10) = NaN; %skip any noise at low channels
    %y = smooth(ccdat.(fn).count,s); 
    [m,mi] = max(y);
    plot(ccdat.(fn).diameter_micron, y, '-')
    line([1 1].*x(mi), ylim, 'linestyle', ':', 'color', 'r')
    title([fname{ii} ' mode ' num2str(x(mi))])
end
if printNow
    print(['\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Coulter\plots\cc_hist' num2str(gcf().Number) ' .png'], '-dpng')
end
%% 
tstr = {'C4MICRO' 'NANO' 'GALBAN' 'ISO' 'PAVLOV' 'CRYPTO' 'RHODO' 'DUN' 'TW'   'SE62' 'HETERO' '7335' '8113' 'C2MICRO'};
vname = {'Micromonas' 'Nano' 'I_galbana' 'Isochrysis' 'Pavlova' 'Crypto' 'Rhodomonas' 'Dun' 'T_weissflogii' 'SE62' 'Heterocapsa' 'Syn7335' 'Syn8113' 'C2Micromonas'};
%tstr = {'DUN'};
tempT = table;
CCSizeCal_smoothhists = table;
CCSizeCal_smoothhists.filestr = tstr';
CCSizeCal_smoothhists.Row = vname';
CCMode = CCSizeCal_smoothhists;
figure
tiledlayout(4,4,"TileSpacing","compact")
for ii = 1:length(tstr)
    nexttile
    ind = find(contains(fname, tstr(ii)));
    if all(run_date(ind) == run_date(ind(1)))
        %CCSizeCal_smoothhists.run_date(ii) = run_date(ind(1));
        CCMode.run_date(ii) = run_date(ind(1));
        CCSizeCal_smoothhists.run_date(ii) = run_date(ind(1));
    else
        disp('check: files from different dates!')
        keyboard
    end
%    x = ccdat.(fname{ind(1)}).diameter_micron;
    y = nan(256,length(ind));
    x = y;
    for iii = 1:length(ind)
        y(:,iii) = ccdat.(fname{ind(iii)}).count;
        x(:,iii) = ccdat.(fname{ind(iii)}).diameter_micron;
    end
    check = (all(x == x(:,1)));
    if all(check) %all same diameter vector
        x = x(:,1);
    elseif all(all(x(:,2:end) == x(:,2))) 
        y(:,1) = interp1(x(:,1),y(:,1),x(:,2));
        x = x(:,2);
    end
    plot(x,y,'.')
    hold on
    plot(x,sum(y,2),'.k')
     if contains(tstr(ii), {'BD' 'FCB' 'NANO' 'MICRO' '8113' '7335'})
        ysm = smooth(medfilt1(sum(y,2),3),8); %4);
       %ysm(x<1.25) = NaN; %skip noise on micro
    else
        ysm = smooth(medfilt1(sum(y,2),7),20); %10);
     end
    ysm(1:10) = NaN; %skip any noise at low channels
    %CCSizeCal_smoothhists.histogram(ii) = [x ysm];
    plot(x,ysm, '-k')
    [m,mi] = max(ysm);
    line([1 1].*x(mi), ylim, 'linestyle', ':', 'color', 'r')
    yl = ylim; yl(2) = min([yl(2) m*1.5]); ylim(yl) 
    title([tstr{ii} ' mode ' num2str(x(mi))])
    CCSizeCal_smoothhists.histogram(ii) = {[x ysm]};
    CCMode.diameter_mode(ii) = x(mi);
end

if printNow
     print(['\\sosiknas1\Lab_data\Attune\size_calibration_July2026\Coulter\plots\cc_hist_culture_pooled.png'], '-dpng')
end

if saveNow
    notes = {'Created with compile_cc_histograms.m'};
    fout = ['\\sosiknas1\lab_data\Attune\Size_calibration_July2026\CCSummary_saved_' datestr(now, 'ddmmmyyyy')]; 
    save(fout, 'CCMode', 'CCSizeCal_smoothhists', 'notes', 'ccdat')
    disp('Results saved:')
    disp(fout)
end
