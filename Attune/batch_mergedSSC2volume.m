function batch_mergedSSC2volume(p)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%get the SSC merge stats
load([p.outpath 'mergeTable'])

if exist([p.outpath, 'beadstat_2026.mat'])
    load([p.outpath, 'beadstat_2026.mat'])
else
    beadtype = 'PT only'; 
end

%get the 1 micron bead means for SSC-A and SSC-H for cases starting on cruise day and later
ind = beadstat_2026.time>dateshift(min(mergeT.filetime),'start', 'day');
beadT = array2table([beadstat_2026.NoOD2centersA(ind,2) beadstat_2026.NoOD2centersH(ind,2) beadstat_2026.NoOD2_hv(ind)], 'VariableName',{'SSCA_1micron' 'SSCH_1micron' 'SSC_hv'});
beadSSCmean = groupsummary(beadT,'SSC_hv','mean'); %compute mean grouped by HV
beadT = array2table([beadstat_2026.OD2centersA(ind,2) beadstat_2026.OD2centersH(ind,2) beadstat_2026.OD2_hv(ind)], 'VariableName',{'GL1A_1micron' 'GL1H_1micron' 'GL1_hv'});
beadGL1mean = groupsummary(beadT,'GL1_hv','mean'); %compute mean grouped by HV

sscstr = strcat("SSC-", p.SSCDIM);
gl1str = strcat("GL1-", p.SSCDIM);
SSCAmin = 500; %use SSC-H below this value

for filecount = 1:height(mergeT)
   if ~rem(filecount,50)
        disp(['volume calculating ' num2str(filecount) ' of ' num2str(height(mergeT))])
   end
    if ~isnan(mergeT.slope_median(filecount))
        filename = [p.fpath, mergeT.filename{filecount}];
        [fcsdat,fcshdr] = fca_readfcs(filename);
        fcsdat = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
        bdrow = find(mergeT.SSC_hv(filecount)==beadSSCmean.SSC_hv);  %FIX THIS LATER IF KEEPING
        ssc_bdnorm = fcsdat.(sscstr)./beadSSCmean.mean_SSCA_1micron(bdrow);
        gl1_bdnorm = fcsdat.(gl1str)./beadGL1mean.mean_GL1A_1micron(bdrow);
        sscmerge_bdnorm = ssc_bdnorm;
        itemp = ssc_bdnorm>2;
        sscmerge_bdnorm(itemp) = gl1_bdnorm(itemp);
        itemp = ssc_bdnorm>0.5 & ssc_bdnorm<2;
        sscmerge_bdnorm(itemp) = mean([ssc_bdnorm(itemp) gl1_bdnorm(itemp)],2);
        volume_cubic_micron = NaN(size(sscmerge_bdnorm)); 
        volume_cubic_micronH = volume_cubic_micron;
        ii = sscmerge_bdnorm>0;
        [volume_cubic_micron(ii) vol_func_string] = Attune_SC2vol(sscmerge_bdnorm(ii),strcat('SSC-',p.SSCDIM));
        ii = fcsdat.(sscstr)<SSCAmin; %for small signals use SSC-H instead (A goes negative)
        volume_cubic_micron(ii) = Attune_SC2vol(fcsdat.("SSC-H")(ii)./beadSSCmean.mean_SSCH_1micron(bdrow),strcat('SSC-H'));
        % sscmerge = fcsdat.(sscstr);
        % itemp = fcsdat.(gl1str)>gl1max | fcsdat.(sscstr)>sscmax; %use merged GL1 for large signals
        % sscmerge(itemp) = fcsdat.(gl1str)(itemp)*mergeT.slope_median(filecount); 
        % %now use merged GL1 for every other point the overlap
        % itemp = find(fcsdat.(gl1str)>gl1min & fcsdat.(gl1str)<gl1max & fcsdat.(sscstr)<sscmax);
        % sscmerge(itemp(1:2:end)) = fcsdat.(gl1str)(itemp(1:2:end)) *mergeT.slope_median(filecount);
        % bdrow = find(mergeT.SSC_hv(filecount)==beadSSCmean.SSC_hv);
        % sscmerge_bdnorm = sscmerge/beadSSCmean.mean_SSCA_1micron(bdrow);
        % volume_cubic_micron = NaN(size(sscmerge_bdnorm)); 
        % volume_cubic_micronH = volume_cubic_micron;
        % ii = sscmerge>0;
        % [volume_cubic_micron(ii) vol_func_string] = Attune_SC2vol(sscmerge_bdnorm(ii),strcat('SSC-',p.SSCDIM));
        % ii = sscmerge<SSCAmin; %for small signals use SSC-H instead (A goes negative)
        % volume_cubic_micron(ii) = Attune_SC2vol(fcsdat.("SSC-H")(ii)./beadSSCmean.mean_SSCH_1micron(bdrow),strcat('SSC-H'));
        %%ctemp = load([p.classpath regexprep(mergeT.filename{filecount},'.fcs', '.mat')]);
        %%c.class = ctemp.class;  clear ctemp %TEMP get rid of old stuff
        c = load([p.classpath regexprep(mergeT.filename{filecount},'.fcs', '.mat')]);
        c.vol_notes = {strcat('calibrated: ', string(datetime())); strcat(' using SSC-', p.SSCDIM); vol_func_string};
        c.ssc_merge_bdnorm = sscmerge_bdnorm; c.beadSSCmean = beadSSCmean; c.beadGL1mean = beadGL1mean; c.file_hv.SSC = mergeT.SSC_hv(filecount); c.file_hv.GL1 = mergeT.GL1_hv(filecount);
        c.volume_cubic_microns = volume_cubic_micron;
        c.merge_info = array2table([mergeT.slope_median(filecount) gl1min gl1max sscmax SSCAmin], 'variablenames', {'slope' 'GL1min' 'GL1max' 'SSCmax' 'SSCAmin'}); 
        save([p.classpath regexprep(mergeT.filename{filecount},'.fcs', '.mat')],'-struct', "c")
    end
end