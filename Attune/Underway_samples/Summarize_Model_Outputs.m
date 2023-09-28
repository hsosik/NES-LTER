

% Quality Control, and create summary table for days that start at dawn


% inputspath = '\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\SlidingWindow\HRS2303';
% outpath = '\\sosiknas1\Lab_data\Attune\cruise_data\20230429_HRS2303\bead_calibrated'; 



function Summarize_Model_Outputs(inputspath, outpath)

inputlist = dir([inputspath filesep '*input.mat']); 

if isempty(inputlist)
    return
end

ModelOutputs = table; 

cruisename = split(inputspath, '\');
if isempty(cruisename{end}) %in case of filesep at end of modelpath 
    cruisename = cruisename{end-1}; 
else
    cruisename = cruisename{end}; 
end



for i = 1:length(inputlist)

    %load data
    load([inputspath filesep inputlist(i).name])
    outputname = inputlist(i).name; 
    outputname = regexprep(outputname, 'input', 'output'); 
    if exist([inputspath filesep outputname])
        load([inputspath filesep outputname])
    else
        continue
    end

    ModelOutputs.cruisename(i) = string(cruisename); 
    ModelOutputs.outputfile(i) = string(outputname); 
    ModelOutputs.daystarttime(i) = daystarttime; 
    ModelOutputs.sampledates{i} = daytable.StartDate;
 

    %Big Pain to go get this variable if we haven't already saved it
     if ~exist('lastdawn')
        t_relstart = datenum(daytable.StartDate) - datenum(daystarttime); 
        hournum = floor(t_relstart.*24)+1;

        nighttime = find(daytable.rad_sw <= 3);
        dawnind = nighttime((nighttime(2:end) - nighttime(1:end-1)) > 9); % look for periods of darkness where next period of darkness isn't for a while 
        if isempty(dawnind)
            daytime = find(daytable.rad_sw > 10);
            dawnind = daytime(find((daytime(2:end) - daytime(1:end-1))>12)+1); %or look for periods of light where it hasn't been light for while 
        end
        if isempty(dawnind) %if day just happend to start at sunset we wont have more than 1 transition
            dawnind = nighttime(end);
        end

        if length(dawnind) > 1; 
            dawnind = dawnind(1); 
        end
        
        dawn_hr = hournum(dawnind);
        

        lastdawn = datetime((datenum(daystarttime) + (dawn_hr - 24)./24), 'ConvertFrom', 'datenum'); 
        lastdawn = dateshift(lastdawn, 'start', 'hour'); 


     end

    ModelOutputs.lastdawn(i) = lastdawn; 
   
    ModelOutputs.Einterps{i} = Einterp; 
    ModelOutputs.modelresults{i} = modelresults; 
    ModelOutputs.divrate(i) = modelresults(17); 
    ModelOutputs.lat(i) = {daytable.lat}; 
    ModelOutputs.lon(i) = {daytable.lon}; 
    ModelOutputs.temperature(i) = {daytable.temperature}; 
    ModelOutputs.salinity(i) = {daytable.salinity}; 
    ModelOutputs.rad_sw(i) = {daytable.rad_sw}; 
    ModelOutputs.Euk_conc(i) = {daytable.Euk_count./daytable.VolAnalyzed_ml};
    ModelOutputs.Syn_conc(i) = {daytable.Syn_count./daytable.VolAnalyzed_ml};


    if hour(daystarttime) == hour(lastdawn); %this is a day that starts at dawn
        ModelOutputs.dawnstart(i) = 1;
    else
        ModelOutputs.dawnstart(i) = 0;
    end
        figure(2)
        subplot(1,2,1)
        pcolor(Vhists)
        title('Data')
        subplot(1,2,2)
        pcolor(simPROPS)
        title('Simulation')

        drawnow 

        ModelOutputs.QC_check(i) = 1; %input('Look good? (1= yes. 0 = no)');

    clear lastdawn
end

if exist('synvolbins')
    SynModelOutputs = ModelOutputs;
    if ~exist([outpath filesep 'ModelOutputs.mat'])
        save([outpath filesep 'ModelOutputs.mat'], 'SynModelOutputs')
    else
        save([outpath filesep 'ModelOutputs.mat'], 'SynModelOutputs', '-append')
    end
else
    EukModelOutputs = ModelOutputs;
    if ~exist([outpath filesep 'ModelOutputs.mat'])
        save([outpath filesep 'ModelOutputs.mat'], 'EukModelOutputs')
    else
        save([outpath filesep 'ModelOutputs.mat'], 'EukModelOutputs', '-append')
    end
end




%now let's merge the model outputs and find sliding window variability 

load([outpath filesep 'ModelOutputs.mat']) %only if both summaries are done already 
if exist('EukModelOutputs') 
    
if ~exist('SynModelOutputs'); 
    SynModelOutputs = table; 
end

MergedOutputs = outerjoin(EukModelOutputs, SynModelOutputs, 'Keys', {'daystarttime', 'daystarttime'}, 'MergeKeys', true); 
MergedOutputs(isnat(MergedOutputs.daystarttime), :) = []; 

ModelResults = table; 
ModelResults.cruisename = MergedOutputs.cruisename_EukModelOutputs; 
ModelResults.cruisename(ismissing(ModelResults.cruisename)) =  MergedOutputs.cruisename_EukModelOutputs(ismissing(ModelResults.cruisename)) ; 

ModelResults.daystarttime = MergedOutputs.daystarttime; 


eukind = find(MergedOutputs.QC_check_EukModelOutputs == 1); 
synind = find(MergedOutputs.QC_check_SynModelOutputs == 1); 

ModelResults.dawnstart(eukind) = MergedOutputs.dawnstart_EukModelOutputs(eukind); 
ModelResults.dawnstart(synind) = MergedOutputs.dawnstart_SynModelOutputs(synind); 

ModelResults.latitude(eukind) = MergedOutputs.lat_EukModelOutputs(eukind); 
ModelResults.latitude(synind) = MergedOutputs.lat_SynModelOutputs(synind); 


ModelResults.longitude(eukind) = MergedOutputs.lon_EukModelOutputs(eukind); 
ModelResults.longitude(synind) = MergedOutputs.lon_SynModelOutputs(synind);


ModelResults.temperature(eukind) = MergedOutputs.temperature_EukModelOutputs(eukind); 
ModelResults.temperature(synind) = MergedOutputs.temperature_SynModelOutputs(synind);

ModelResults.salinity(eukind) = MergedOutputs.salinity_EukModelOutputs(eukind); 
ModelResults.salinity(synind) = MergedOutputs.salinity_SynModelOutputs(synind);

ModelResults.rad_sw(eukind) = MergedOutputs.rad_sw_EukModelOutputs(eukind); 
ModelResults.rad_sw(synind) = MergedOutputs.rad_sw_SynModelOutputs(synind);


ModelResults.euk_per_ml(eukind) = MergedOutputs.Euk_conc_EukModelOutputs(eukind); 
ModelResults.euk_per_ml(synind) = MergedOutputs.Euk_conc_SynModelOutputs(synind);


ModelResults.syn_per_ml(eukind) = MergedOutputs.Syn_conc_EukModelOutputs(eukind); 
ModelResults.syn_per_ml(synind) = MergedOutputs.Syn_conc_SynModelOutputs(synind);


ModelResults.euk_divrate(eukind) = MergedOutputs.divrate_EukModelOutputs(eukind); 

ModelResults.syn_divrate(synind) = MergedOutputs.divrate_SynModelOutputs(synind); 


save([outpath filesep 'ModelOutputs.mat'], 'ModelResults', '-append')

end


end
