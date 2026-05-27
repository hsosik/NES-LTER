function Process_Preserved_AllSteps_2026(cruise, steps2do)
% function Process_Preserved_AllSteps_2026(cruise, steps2do)
% e.g.,
%   Process_Preserved_AllSteps_2026('EN661', [1 3]) %for just steps 1 and 3
% or
%   Process_Preserved_AllSteps_2026('EN661', [1:7]) %for all steps at once
%
% February 2026 (Heidi M. Sosik, WHOI)
% update from Process_Preserved_AllSteps (from Bethany);
% streamline as function
% handle defaults for all steps

%Trying to pull together all steps into one:
% Process_preserved_pt1
% Proces_Preserved_samples_with_AWS
% Get_metadata_for_Gated_Table
% Classify_from_Gated_Tables
% create class files with volumes
% Convert_Gated_Table_to_CNTable

%% INPUTS
fclose('all');
addpath(fileparts(pwd)) %add path one folder up for shared functions from underway code
cruise = lower(cruise);
%make step vector
step = zeros(1,7);
step(steps2do) = 1;

basepath_temp =  '\\sosiknas1\Lab_data\Attune\cruise_data\';
temp = dir([basepath_temp '*' cruise]);
if ~isempty(temp)
    p.basepath = [basepath_temp temp.name filesep 'preserved'];
else
    disp('No directory for cruise:')
    disp([basepath_temp '*' cruise])
    keyboard
end
clear temp

%Set all steps to 1 if starting from beginning
% Step1 = 0; %make FCSList
% Step2 = 0; %go look at AWS files to find gate assignments
% Step3 = 0; % add metadata to gated table
% Step4 = 0; % classify using gate_table
% Step5 = 0; %size calibrate and create class files
% Step6 = 1; %convert gated table to Summary table
% Step7 = 0; %Reformat Summary Table to have EDI headers


%% Set up

% some file structure setup
p.fpath = [p.basepath filesep 'FCS' filesep];
p.outpath = [p.basepath filesep 'outputs' filesep];
p.classpath = [p.outpath 'class' filesep];
p.awspath = [p.basepath filesep 'aws\'];

if ~exist(p.outpath, 'dir')
    mkdir(p.outpath)
end
if ~exist(p.classpath, 'dir')
    mkdir(p.classpath)
end

save([p.outpath '\Processing_variables.mat'])

p.uw_fullname = '';
p.elogpath = '';
switch cruise
    case 'AR29'
        p.uw_fullname = '\\sosiknas1\Lab_data\SPIROPA\20180414_AR29\underway\procAR29_underway_all.mat';
    case 'SR1812'
        p.elogpath = '\\sosiknas1\Lab_data\EXPORTS\SR1812\share\R2R_ELOG_SR1812_FINAL_EVENTLOG_20180913_022931_edited_Sosik_postcruise.csv';
        p.bottlefile = '\\sosiknas1\Lab_data\EXPORTS\SallyRideSIOBottleFiles_v6.csv';
    case 'RB1904'
        p.uw_fullname = '\\sosiknas1\Lab_data\SPIROPA\20190503_RB1904\compiled_underway\rb1904_uw_compiled.mat';
    case 'TN368'
        p.uw_fullname = '\\sosiknas1\Lab_data\SPIROPA\20190705_TN368\compiled_underway\tn368_uw_compiled.mat';
        p.bottlefile =  '\\sosiknas1\Lab_data\Attune\cruise_data\20190705_TN368\preserved\tn368_bottle_data_Apr_2020_table.mat';
    case 'HB1907'
        p.bottlefile = "\\sosiknas1\lab_data\OTZ\20190725_HB1907\HB1907_ctd_Deepseefolder\HB1907_bottles.csv";
        %this file does not have temperature and salinity but does have nutrients
        %p.bottlefile = '\\sosiknas1\Lab_data\Attune\cruise_data\20190725_HB1907\preserved\bottle_environmental_data_partial.csv';
    case 'AR43'
        p.uw_fullname = '\\sosiknas1\Lab_data\OTZ\20200311_AR43\underway\proc\ar43_underway.csv';
        p.bottlefile = '\\sosiknas1\Lab_data\OTZ\20200311_AR43\ctd\ar43_ctd_bottles.csv';
    case 'DY131'
        p.elogpath = '\\sosiknas1\Lab_data\EXPORTS\DY131\r2r\R2R_ELOG_dy131_FINAL_EVENTLOG_20210601_094434_edited_Sosik_postcruise.csv';
        p.bottlefile =  '\\sosiknas1\Lab_data\EXPORTS\DY131\CTDbottle\EXPORTS2021_DY131_BottleFile_R0d_20220615.csv';
    case 'SG2105'
        p.bottlefile = '\\sosiknas1\Lab_data\Attune\cruise_data\20210512_SG2105\preserved\EXPORTS2021_SDG2105_BottleFile_R0_20210720T124833.csv';
        % case 'EN657'
        %     p.uw_fullname = ['https://nes-lter-data.whoi.edu/api/underway/' lower(cruise)];
        %     p.bottlefile = ['https://nes-lter-api.whoi.edu/api/ctd/bottles/' lower(cruise)];
        %     p.elogpath = '\\sosiknas1\Lab_data\LTER\20201013_EN657\eLog\R2R_ELOG_EN657_FINAL_EVENTLOG_20201018_134037.csv';
    otherwise %DEFAULT, NES LTER api
        p.uw_fullname = ['https://nes-lter-api.whoi.edu/api/underway/' cruise];
        p.bottlefile = ['https://nes-lter-api.whoi.edu/api/ctd/bottles/' cruise];
        p.elogpath = ['https://nes-lter-api.whoi.edu//api/events/' cruise];
end

% Step 1 - make FCSlist defaults
p.filetype2exclude = {'Rinses', 'skip', 'xxx', 'qwater';}; % 'Sample(',
%p.beadfiles2include = {'FCB_bead'};  %%%%WHERE ARE THE BEAD RUNS????
% Step 5 - %size calibrate and create class files
p.SSCDIM = 'A'; %needed for Step 4 & 5, SSCDIM = 'A' or 'H'

%% Step 1 - make FCSlist
if step(1)
    [FCSfileinfo] = FCS_DateTimeList(p.fpath);
    FCSfileinfo(contains(FCSfileinfo.fcslist,p.filetype2exclude),:) = []; %remove entries for excluded cases
    %extract cast and niskin from filename
    temp = regexp(FCSfileinfo.fcslist,'_C\d*N\d*', 'match');
    temp_ind = ~cellfun(@isempty, temp);
    FCSfileinfo.cast(temp_ind) = str2double(extractBetween(string(temp(temp_ind)), 'C', 'N'));
    FCSfileinfo.niskin(temp_ind) = str2double(extractAfter(string(temp(temp_ind)), 'N'));
    FCSfileinfo = movevars(FCSfileinfo, {'cast' 'niskin'}, 'After', 'fcslist');
    save([p.outpath 'FCSfileinfo.mat'], 'FCSfileinfo');
else
    load([p.outpath filesep 'FCSfileinfo.mat'])
end
clearvars -except cruise step p

%         Vol_analyzed_ml(i) = .8.*fcshdr.VOL ./ 1e6; %changed to 80% based on uniform time gate
%         Vol_ml(i) = fcshdr.VOL ;

%% Step 2 - go look at AWS files to find gate assignments

if step(2)
    %T = load([p.outpath '\FCSfileinfo.mat']);
    %T = T.FCSfileinfo;
    load([p.outpath 'FCSfileinfo']);
    gated_table = FCSfileinfo;
    gated_table = renamevars(gated_table, 'vol_analyzed', 'vol_run'); %save vol_analyzed for later to account for time gate

    % go through files, for each file, find matching aws file
    %
    % Make one figure for each niskin Bottle, reporting depth and total counts
    %
    %Make a list of fcs files that didn%t have an aws files

    %first get runtype directory names
    runtypes = dir(p.awspath); runtypes = string({runtypes.name})';
    runtypes = runtypes(~startsWith(runtypes, '.'));
    fcslist_CN = regexp(FCSfileinfo.fcslist,'_C\d*N\d*.*fcs', 'match');
    fcslist_CN = extractBetween(string(fcslist_CN), '_','.fcs');

    for ii = 1:length(runtypes)
        tind = find(contains(gated_table.fcslist, runtypes(ii)));
        awslist = dir(strcat(p.awspath, runtypes(ii),filesep, '*.aws'));
        awslist2 = regexprep({awslist.name}, '.aws', '')';
        [~,ia,ib] = intersect(fcslist_CN(tind), awslist2);
        gated_table.awsfilename(tind(ia)) = strcat(runtypes(ii),filesep, {awslist(ib).name});
        figpath = strcat(p.outpath , 'figs\', runtypes(ii)) ;
        if ~exist(figpath, 'dir')
            mkdir(figpath)
        end
        %%%NEED SPECIAL CASES FOR SR1812 and DY131 and uw samples
    end
    gated_table = movevars(gated_table,'awsfilename', 'After', 'fcslist');
    no_aws_files = gated_table.fcslist(ismissing(gated_table.awsfilename));

    %                 if strcmp(cruise,'SR1812')
    %                     temp = split(filename, '_');
    %                     ind = find(awslist == [temp{end-1} '.aws']);
    %                 elseif strcmp(cruise, 'DY131') && T.Cast(i) == 0
    %                     temp = split(filename, 'Ex');
    %                     ind = find(awslist == ['Ex' temp{end}(1:end-4) '.aws']);
    %                 elseif T.Cast(i) == 0 && T.niskin(i) == 0 %need case for UW data
    %                     uwname = split(filename, '_');
    %                     uwname = regexprep(uwname{end}, '.fcs', '.aws');
    %                     ind = find(awslist == uwname);
    %                 else
    %                     ind = find(awslist == strcat("C", num2str(T.Cast(i), '%02.f'), 'N', num2str(T.niskin(i), '%02.f'), '.aws'));
    %                 end


    for i = 1:height(gated_table)
        if ~ismissing(gated_table.awsfilename(i))
            %if we have an aws file proceed with gating
            [fcsdat, fcshdr]  = fca_readfcs([p.fpath gated_table.fcslist{i}]);
            [gate_assignments, polygon_names, polygon_vars, polygon_vals, gate_names, gate_logic_legible, parent_logic, time_gate_fraction] = ApplyAWSgates_hierarchical(strcat(p.awspath, gated_table.awsfilename{i}), fcsdat, fcshdr);

            gated_table.gate_names{i} = gate_names;
            gated_table.gate_assignments{i} = gate_assignments;
            gated_table.gate_logic{i} = gate_logic_legible;
            gated_table.polygon_names{i} = polygon_names;
            gated_table.polygon_vars{i} = polygon_vars;
            gated_table.polygon_vals{i} = polygon_vals;
            make_figure_aws(fcsdat, fcshdr, gate_assignments, polygon_vars, polygon_vals, gate_names, figpath);
        end
    end
    gated_table.vol_analyzed_ml = gated_table.vol_run*time_gate_fraction; %account for time gate

    save([p.outpath 'Gated_Table.mat'], 'gated_table', 'no_aws_files')

    clearvars -except cruise step p
end

%% Step 3 - add metadata to gated table

if step(3)

    load([p.outpath '\Gated_Table.mat']);

    switch cruise
        case {'sr1812' 'dy131' 'sg2105'}
            bottledata = readtable(p.bottlefile);
            bottledata = renamevars(bottledata, {'Cast' 'BottleNo' 'Lat' 'Lon', 'depth' 'potTemp1' 'sal1' 'sdate'},...
                {'cast' 'niskin' 'latitude' 'longitude' 'depth_m' 'potemp090c' 'sal00' 'date_sampled'});
            %       bottledata.date_sampled = datetime(bottledata.sdate, 'InputFormat', 'dd/MM/yyyy HH:mm');
        case 'tn368'
            load(p.bottlefile)         %use mat file for SPIROPA Cruises not the same format >:(
            temp = importdata('\\sosiknas1\Lab_data\SPIROPA\20190705_TN368\fromOlga\tn368_niskin_pressure_depth.txt');
            bottle_depth = array2table(temp.data, 'VariableNames', temp.colheaders); clear temp
            bottledata = outerjoin(bottle_depth,BTL, 'RightKeys', {'Cast' 'TargetDepth_m'}, 'LeftKeys', {'%Cast' 'TargetDepth(m)'}, 'RightVariables', { 'datetime' 'Latitude_decimalDeg' 'Longitude_decimalDeg' 'Potemp090C' 'Sal00'}, 'Type', 'left');
            bottledata = unique(bottledata,'rows'); %some rows are repeated in BTL because of Productivity entries
            bottledata = renamevars(bottledata, {'%Cast' 'Bottle', 'Depth(m)' 'Potemp090C' 'Sal00' 'Latitude_decimalDeg' 'Longitude_decimalDeg' 'datetime'},...
                {'cast' 'niskin' 'depth_m' 'potemp090c' 'salinity' 'latitude' 'longitude' 'date_sampled'});
        case 'hb1907'
            bottledata = readtable(p.bottlefile);
            %%%%%BE CAREFUL THIS ONE ISN'TREALLY POTENTIAL TEMPERATURE (renamed for consistency with downstream code)
            bottledata = renamevars(bottledata, {'depth' 'date' 't090c'},{'depth_m' 'date_sampled' 'potemp090c'});
        case 'ar43'
            p.bottlefile = '\\sosiknas1\Lab_data\OTZ\20200311_AR43\ctd\ar43_ctd_bottles.csv';
            bottledata = renamevars(bottledata, {'depth' 'date'},{'depth_m' 'date_sampled'});
        otherwise %DEFAULT, NES LTER api
            bottledata = readtable(p.bottlefile, 'delimiter', ',');
            cast_meta = readtable(regexprep(p.bottlefile, 'bottles', 'metadata'), 'delimiter', ',');
            bottledata = outerjoin(bottledata, cast_meta, 'keys', 'cast', 'RightVariables', 'nearest_station', 'Type','left');
            bottledata = renamevars(bottledata, {'date' 'depsm' 'sal00'}, {'date_sampled' 'depth_m' 'salinity'});
    end

    switch cruise
        case {'sr1812' 'dy131' 'sg2105'}
            gated_table = outerjoin(gated_table, bottledata, 'keys', {'cast' 'niskin'}, 'RightVariables', {'latitude' 'longitude' 'date_sampled' 'depth_m' 'potemp090c' 'salinity' 'r2r_event'}, 'Type','left');
        otherwise
            if ~ismember('nearest_station',bottledata.Properties.VariableNames)
                bottledata.nearest_station(:) = '';
            end
            gated_table = outerjoin(gated_table, bottledata, 'keys', {'cast' 'niskin'}, 'RightVariables', {'latitude' 'longitude' 'date_sampled' 'nearest_station' 'depth_m' 'potemp090c' 'salinity'}, 'Type','left');
    end

    save([p.outpath '\Gated_Table.mat'], 'gated_table', 'no_aws_files');


    % now underways
    %gated_table.fcslist(gated_table.Cast==0)
    %EN657 - one underway with FCM in elog comment field
    %SR1812 & DY131 - uw marked with Cast = -9999 (after step3)
    %TN368 - 7 Cast=0 but fcs marked with 'maybe' throwing off matchup, also 'SPIROPA_TN368_Jul2019_preserved(5)_phyto_PE_SSC_15N24(1).fcs'
    %AR38 - B01, also 'NESLTER_AR38_Sep2019_preserved_phyto_CHL_SSC_CxxN07.fcs', 'bucket' in elog comment field

    if sum(gated_table.Cast==0) %|| sum(gated_table.Cast==-9999) %cases for underway or other samples
        elog = readtable(p.elogpath, 'delimiter', ',');
        disp('No cast matches:')
        disp(gated_table.fcslist(gated_table.cast==0))
        disp(gated_table.fcslist(gated_table.cast==-9999))
        switch cruise
            case {'sr1812' 'dy131'}  %no underway temp and salinity yet for sr1812
                ind = gated_table.cast== 0; %-9999;
                if startsWith(cruise, 'sr')
                    temp = (regexprep(regexprep(string(regexp(gated_table.fcslist(ind), '_\d\d\d\d\d_UW', 'match')),'_UW', ''), '_', 'EX'));
                    [tt, ib] = ismember(lower(temp), lower(elog.RoeslerEvent));
                    %matching date-time to format from bottle files for EXPORTS
                    gated_table.date_sampled(ind) = datetime(elog.dateTime8601(ib), 'InputFormat','uuuu-MM-dd''T''HH:mm:ss+0000', 'Format','dd/MM/yyyy HH:mm');
                elseif startsWith(cruise, 'dy')
                    temp = string(regexp(gated_table.fcslist(ind),'Ex\d\d\d\d\d', 'match'));
                    [tt, ib] = ismember(lower(temp), lower(elog.RoeslerEvent));
                    %matching date-time to format from bottle files for EXPORTS
                    gated_table.date_sampled(ind) = datetime(elog.dateTime8601(ib), 'InputFormat','uuuu-MM-dd''T''HH:mm:ss+00:00', 'Format','dd/MM/yyyy HH:mm');
                end
                gated_table(ind,{'latitude' 'longitude' 'r2r_event'}) = elog(ib,{'Latitude' 'Longitude' 'R2R_Event'});
                gated_table.depth_m(ind) = 6;
                gated_table.cast(ind) = -9999;
                gated_table.niskin(ind) = -9999;
            case {'en657' 'ar38'} %one case wth FCM or bucket in elog comment
                ind = gated_table.cast==0;
                ii = contains(elog.Comment, 'FCM') | contains(elog.comment, 'bucket');
                gated_table(ind,{'date_sampled' 'latitude' 'longitude'}) = repmat(elog(ii,{'dateTime8601' 'Latitude' 'Longitude'}),1,sum(ind));
                uw = readtable(p.uw_fullname, 'delimiter', ',');
                uw.datetime = datetime(uw.date, 'InputFormat', 'uuuu-MM-dd HH:mm:ss+00:00');
                [dmin mind] = min(abs(uw.datetime-datetime(elog.dateTime8601(ii),'InputFormat', 'uuuu-MM-dd HH:mm:ss+00:00')));
                if startsWith(cruise, 'en')
                    gated_table(ind,{'potemp090c' 'salinity'}) = repmat(uw(mind,{'tsg1_temperature' 'tsg1_salinity'}),1,sum(ind));
                else
                    gated_table(ind,{'potemp090c' 'salinity'}) = repmat(uw(mind,{'sbe48t' 'sbe45s'}),1,sum(ind));
                end
                gated_table.depth_m(ind) = 5;
        end
    end
    %overwrite saved gated table with metadata
    save([p.outpath '\Gated_Table.mat'], 'gated_table', 'no_aws_files');

    clearvars -except cruise step p

end

%% Step 4 - classify using gate_table

if step(4)

    load([p.outpath 'Gated_Table.mat']);

    cut_off_pro_pop = [];

    % Here is where we decide which gates we are interested in and how we will find them

    % gates_of_interest = {'syn'; 'Syn'; 'Euks'; 'euks'; 'pro'; 'Pro'; 'bacteria'; 'Bacteria'};

    classnames = {'Euks == 1'; 'Syn == 2'; 'Bacteria == 3'; 'Pro == 4'; 'LowPE_Euks = 5'; 'HighPE_Euks = 6'};

    class = cell(1,height(gated_table));
    syn_conc = nan(height(gated_table),1);
    euk_conc = syn_conc;
    bact_incl_pro_conc = syn_conc;
    pro_conc = syn_conc;
    lp_euk_conc = syn_conc;
    hp_euk_conc = syn_conc;

    for i = 1:height(gated_table)
        filename = gated_table.fcslist{i};
        [fcsdat, fcshdr]  = fca_readfcs([p.fpath filename]);

        a = size(gated_table.gate_assignments{i});

        if max(a) == 0
            class{i} = [];
            continue
        end

        %first go through and generate a single class assingment for each
        %particle
        gate_assign_i =gated_table.gate_assignments{i};
        class_i = zeros(1, max(size(gate_assign_i)));
        gate_names = gated_table.gate_names{i};


        % look for euks
        if sum(strcmp(gate_names, 'euk'))
            gate_num_1 = strcmp(gate_names, 'euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_1, :)) = 1;
        end
        if sum(strcmp(gate_names, 'Euk'))
            gate_num_1 = strcmp(gate_names, 'Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_1, :)) = 1;
        end

        if sum(strcmp(gate_names, 'LowP_Euk'))
            gate_num_5 = strcmp(gate_names, 'LowP_Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_5, :)) = 5;
        end

        if sum(strcmp(gate_names, 'LowPE_Euk'))
            gate_num_5 = strcmp(gate_names, 'LowPE_Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_5, :)) = 5;
        end
        if sum(strcmp(gate_names, 'HighP_Euk'))
            gate_num_6 = strcmp(gate_names, 'HighP_Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_6, :)) = 6;
        end
        if sum(strcmp(gate_names, 'HighPE_Euk'))
            gate_num_6 = strcmp(gate_names, 'HighPE_Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_6, :)) = 6;
        end

        %next syn
        if sum(strcmp(gate_names, 'Syn'))
            gate_num_2 = strcmp(gate_names, 'Syn');
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2;
        end
        if sum(strcmp(gate_names, 'syn'))
            gate_num_2 = strcmp(gate_names, 'syn');
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2;
        end
        if sum(strcmp(gate_names, 'SYN'))
            gate_num_2 = strcmp(gate_names, 'SYN');
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2;
        end


        %next bacteria, may include prochlorococcus, but those will be
        %reassigned below
        if sum(strcmp(gate_names, 'Bacteria'))
            gate_num_b = strcmp(gate_names, 'Bacteria');
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3;
        end
        if sum(strcmp(gate_names, 'bacteria'))
            gate_num_b = strcmp(gate_names, 'bacteria');
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3;
        end
        if sum(strcmp(gate_names, 'BacPro'))
            gate_num_b = strcmp(gate_names, 'BacPro');
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3;
        end



        %next assign pro, allowed to overwrite bacteria label
        if sum(strcmp(gate_names, 'pro'))
            gate_num_4 = strcmp(gate_names, 'pro');
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4;
        end
        if sum(strcmp(gate_names, 'Pro'))
            gate_num_4 = strcmp(gate_names, 'Pro');
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4;
        end
        if sum(strcmp(gate_names, 'Euk_sm'))
            gate_num_4 = strcmp(gate_names, 'Euk_sm');
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4;
        end
        if sum(strcmp(gate_names, 'Proc'))
            gate_num_4 = strcmp(gate_names, 'Proc');
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4;
        end


        %ok now save class assignments

        class{i} = class_i;


        %onto calculating concentrations
        concent_i = nan(1, 6);

        sizefrac = gated_table.gate_logic{i};

        if exist('gate_num_1')
            concent_i(1) = sum(class_i == 1)./gated_table.vol_analyzed_ml(i);
        end

        if exist('gate_num_2')
            concent_i(2) = sum(class_i == 2)./gated_table.vol_analyzed_ml(i);
        end

        if exist('gate_num_b')
            concent_i(3) = sum(class_i == 3 | class_i ==4)./gated_table.vol_analyzed_ml(i); %here is where we combine Bacteria and Syn class assignments!
        end

        if exist('gate_num_4')
            concent_i(4) = sum(class_i == 4)./gated_table.vol_analyzed_ml(i);

            %Here is where we check if pro gate is cut off
            if ~contains(filename, 'pro', 'IgnoreCase', true) %if it's not a pro run specifically but there is a pro gate

                x_ch = strmatch('SSC-H', {fcshdr.par.name});
                y_ch = strmatch('BL3-H', {fcshdr.par.name});

                figure(2)
                [N, X] = hist(fcsdat(class_i==4, y_ch)); %check to see if cutoff on BL3
                hist(fcsdat(class_i==4, y_ch))
                if N(1) < N(2) & N(1) < N(3) & N(1) < N(4) & N(4) < N(3)
                    title('good')
                elseif N(1) < N(2) & N(2) < N(3) & N(3) < N(4)
                    title('good')
                else
                    %l = input('population cut off? y/n', 's');
                    %if l == 'y'
                    concent_i(4) = NaN;
                    title('cut off')
                    cut_off_pro_pop = [cut_off_pro_pop; string(filename)];
                end

            end

        end

        if exist('gate_num_5')
            concent_i(5) = sum(class_i == 5)./gated_table.vol_analyzed_ml(i);
        end

        if exist('gate_num_6')
            concent_i(6) = sum(class_i == 6)./gated_table.vol_analyzed_ml(i);
        end

        disp(gated_table(i,:))
        %disp(concent_i)

        %now save concentrations
        euk_conc(i) = concent_i(1);
        syn_conc(i) = concent_i(2);
        bact_incl_pro_conc(i) = concent_i(3) ;
        pro_conc(i) = concent_i(4);
        lp_euk_conc(i) = concent_i(5);
        hp_euk_conc(i) = concent_i(6);

        clear gate_num_1 gate_num_2 gate_num_4 gate_num_b gate_assign_i class_i gate_num_5 gate_num_6

        notes = "'Euks == 1'; 'Syn == 2'; 'Bacteria & Pro == 3'; 'Pro == 4'; 'LowPE_Euks = 5'; 'HighPE_Euks - 6'";
        save([p.classpath regexprep(gated_table.fcslist{i}, '.fcs', '')], 'class', 'notes') %I think class files wont have volume yet if we don't have bead statistics to calibrate 

    end

    gated_table.class = class';
    gated_table.Euk_conc = euk_conc;
    gated_table.LowP_Euk_conc = lp_euk_conc;
    gated_table.HighP_Euk_conc = hp_euk_conc;

    gated_table.Syn_conc = syn_conc;
    gated_table.Bact_incl_pro_conc = bact_incl_pro_conc;
    gated_table.Pro_conc = pro_conc;

    save([p.outpath '\Gated_Table.mat'], 'gated_table', 'no_aws_files', 'cut_off_pro_pop');

    clearvars -except cruise step p

end

%% Step 5 - size calibrate and create class files

if step(5)

   %this part Bethany calls calibration is SSC and GL1 merging
   %if strcmp(p.OD2setting, 'GL1')  %%%%%%%%DOUBLE check that OD2 setting is always GL1 for discretes
        get_calibration_stats_linear_2021(p.outpath, p.classpath, 1, p.SSCDIM) %A means ssch_ch_num is ssc-a
   %end

    %%%CHECK WHAT NEEDS TO BE DONE TO UPDATE THIS?? if anything
    B = load('\\sosiknas1\Lab_data\Attune\cruise_data\beads\FCB_bead_mix_experiment_settings\between_cruises\outputs\beadstat.mat')

    load([p.outpath '\Gated_Table.mat']);
    
    if ~exist('cut_off_pro_pop', 'var')
        cut_off_pro_pop = NaN;
    end

    median_volumes = nan(height(gated_table), 6);

%    DIM = 'A';
%    saverpath = [p.classpath 'calibration'];
%    if ~exist([p.classpath 'calibration'], 'dir')
%        mkdir([p.classpath 'calibration'])
%    end
%
% 
% %    joint_table = [];
% 
%     for counti = 1:height(gated_table)
%         qc_warning = 0;
% 
%         if isempty(gated_table.class{counti})
%             continue
%         end
% 
%         %load corresponding fcs file
%         filename = gated_table.fcslist{counti};
%         [fcsdat,fcshdr] = fca_readfcs([p.fpath filename]);

   %load back joint_table now saved by get_calibration_stats_linear_2021
        joint_table = load([p.classpath '/calibration/table.mat']);
        joint_table = joint_table.joint_table;

        for counti = 1:height(joint_table)
        classfilename = [p.classpath regexprep(joint_table.filename{counti}, '.fcs', '.mat')];
        class = gated_table.class{counti};

        [fcsdat,fcshdr] = fca_readfcs([p.fpath joint_table.filename{counti}]);

         ssc_ch_num = strmatch(['SSC-' DIM], {fcshdr.par.name});
         gl1_ch_num = strmatch(['GL1-' DIM], {fcshdr.par.name});
        % %bl3_ch_num = strmatch(['BL3-' DIM], {fcshdr.par.name});
        % 
        % file_hv = fcshdr.par(ssc_ch_num).hv; %heidi

        %if file_hv == 220 %useful for working with bacteria size calibration
            %keyboard
            %  continue
        %end

         gl1_vals = fcsdat(:, gl1_ch_num);
         ssc_value = fcsdat(:,ssc_ch_num);
        % 
        % ssch_ch_num = strmatch(['SSC-H'], {fcshdr.par.name});
        % ssc_value(ssc_value<=0) = fcsdat(ssc_value<=0, ssch_ch_num); %replace negative A values with H value as proxy
        % 
        % l_bound = 1e3;
        % r_bound = 9e3; %fixed
        % 
        % if r_bound < l_bound;
        %    qc_warning =1;
        %end
        % 
        % ind_to_fit = find(gl1_vals>l_bound & gl1_vals<r_bound); %not a super robust way to choose the range to fit
        % 
        % LM = fitlm(gl1_vals(ind_to_fit), ssc_value(ind_to_fit));
        % intercept = LM.Coefficients.Estimate(1);
        % slope = LM.Coefficients.Estimate(2);
        % 
        % if isnan(intercept) %if data is bad, linear model doesn't work
        %     qc_warning = 1;
        %     scatter_value = fcsdat(:,ssc_ch_num);
        % elseif LM.Rsquared.Adjusted < .8
        %     qc_warning = 1;
        % end
        % 
        % new_ssc_vals = ssc_value;
        % new_ssc_vals(gl1_vals>r_bound) = 10.^[intercept + slope.*(gl1_vals(gl1_vals>r_bound))];
        % 
        % %dither according to heidi's example
        % t = find((gl1_vals <= r_bound & ssc_value >= 10.^(intercept + slope.*(r_bound))) | (gl1_vals >= r_bound & ssc_value <= 10.^(intercept + slope.*(r_bound))));
        % new_ssc_vals(t(1:2:end)) = 10.^[intercept + slope.*(gl1_vals(t(1:2:end)))];
        % 

        % calibrate_info = table;
        % calibrate_info.filename = string(filename);
        % calibrate_info.ssc_ch_num = ssc_ch_num;
        % calibrate_info.qc = qc_warning;
        % calibrate_info.intercept = intercept;
        % calibrate_info.slope = slope;
        % calibrate_info.rsquared = LM.Rsquared.Adjusted;
        % calibrate_info.numpoints = length(ind_to_fit);
        % calibrate_info.rightbound = r_bound;
        % 
        % joint_table = [joint_table; calibrate_info];
        % 
        % figure(98)
        % plot(gl1_vals, ssc_value, '.')
        % xlabel('GL1 - low sensitivity')
        % ylabel('SSC - high sensitivity')
        % hold on
        % plot([1 max(gl1_vals)], [intercept (intercept + slope*(max(gl1_vals)))])
        % plot(gl1_vals(t), ssc_value(t), '.b')
        % plot(gl1_vals, new_ssc_vals, 'g.')
        % plot(gl1_vals(ind_to_fit), ssc_value(ind_to_fit), 'r.')
        % title({filename; [num2str(qc_warning)]}, 'Interpreter', 'none')
        % axis([0 1e5 0 1e6])
        % 
        % if ~contains(filename, 'hbac') & ~contains(filename, 'pro')
        %     print(figure(98), [saverpath filesep regexprep(filename, '.fcs', '.png')], '-dpng')
        % end
        % clf(98)

     
        
        %now use calibration to get volume data
        negA_ind = ssc_value<=0; %keep track of particles for which area is neagtive.
        ssch_ch_num = strmatch(['SSC-H'], {fcshdr.par.name});
        ssc_value(negA_ind) = fcsdat(negA_ind, ssch_ch_num); %replace negative A values with H value as proxy

        file_hv = fcshdr.par(ssc_ch_num).hv;

        %get median 1micron bead center for a bead run that is nearby in time
        %in the lab.
        %use OD2 measurements to project to NoOD2 values

        filetime = datetime([fcshdr.date, ' ', fcshdr.starttime]);
        beadstat = beadstat(beadstat.QC_flag ==0,:);
        [alert,ind1] = min(abs(datenum(beadstat.time)-datenum(filetime)));
        if alert > 30
            disp('more than a month between bead run and file run. Update beadstat.')
            keyboard
        end

        bead_file = beadstat.filename(ind1);
        if beadstat.QC_flag(ind1) == 1
            keyboard
        else
            bead_value = [beadstat.NoOD2_hv(ind1) beadstat.NoOD2centers(ind1,2)]; %bead value on SSC
            bead_value_to_convert = [beadstat.OD2_hv(ind1) beadstat.OD2centers(ind1,2)]; %bead value on GL1

        end

        %use results of bead hv experiment to convert to appropriate value for file_hv
        bead_value = 10.^(0.016588.*(-bead_value(1)+file_hv) + log10(bead_value(2)));
        bead_value_to_convert = 10.^(0.016659.*(bead_value_to_convert(1)-file_hv) + log10(bead_value_to_convert(2)));


        %if using SSC to GL1 calibration
        l_bound = 5e2;

        %replace high values with estimates from low sensitivity channel
        new_ssc_vals = ssc_value;
        intercept = joint_table.intercept(counti); slope = joint_table.slope(counti);
        r_bound = joint_table.rightbound(counti);
        if ~contains(filename, 'hbac') & ~contains(filename, 'pro')
            new_ssc_vals(gl1_vals>r_bound) = [intercept + slope.*(gl1_vals(gl1_vals>r_bound))];

            %dither according to heidi's example
            t = find((gl1_vals <= r_bound & ssc_value >= intercept + slope.*(r_bound)) | (gl1_vals >= r_bound & ssc_value <= (intercept + slope.*(r_bound))));
            new_ssc_vals(t(1:2:end)) = [intercept + slope.*(gl1_vals(t(1:2:end)))];
        end

        scatter_value = real(new_ssc_vals);

        %also convert bead value if it is large
        if bead_value_to_convert>r_bound
            %bead_value = [intercept + slope.*(bead_value_to_convert)];
        end %if not, bead_value is SSC no filter measurement


        %now convert modified scatter_values to volumes & save results
        volume = 10.^(1.24*log10(scatter_value./bead_value) + 1.064); %based on linear fit to scaled ssch on OD2 filter March 2019
        volumestring = '10.^(1.24*log10(scatter_value./bead_value) + 1.064';

        %treat H values differently from A values
        volume(negA_ind) = 10.^(1.4225*log10(scatter_value(negA_ind)./bead_value) + 1.1432);
        volume = real(volume);

        vol_notes = {strcat('calibrated: ', string(datetime())); strcat('using SSC-', DIM, ' and GL1 Linear Scale Fit: right bound ', num2str(r_bound), 'intercept ', num2str(intercept), 'slope ', num2str(slope));
            volumestring};

        save([classfilename], 'class', 'notes', 'volume', 'ssc_value', 'vol_notes', 'bead_file', 'bead_value', 'negA_ind', 'file_hv')


        for c = 1:6
            median_volumes(counti, c) = nanmedian(volume(class == c));
        end

    end

    save([saverpath '/table.mat'], 'joint_table')


    gated_table.median_volumes_euk = median_volumes(:,1);
    gated_table.median_volumes_syn = median_volumes(:,2);
    gated_table.median_volumes_bact = median_volumes(:,3);
    gated_table.median_volumes_pro = median_volumes(:,4);
    gated_table.median_volumes_low_pe_euk = median_volumes(:,5);
    gated_table.median_volumes_high_pe_euk = median_volumes(:,6);


    save([p.outpath '\Gated_Table.mat'], 'gated_table', 'no_aws_files', 'cut_off_pro_pop');

    clearvars -except cruise step p

end

%% Step 6 - convert gated table to Summary table

if step(6)

    load([p.outpath '\Gated_Table.mat']);
    
    %gated_table(contains(gated_table.fcslist, 'lower_thresh'), :) = [];
    %gated_table(gated_table.cast == 0, :) = [];

    %first some counting of underway samples
    gind = (gated_table.cast==0 | gated_table.cast==-9999);
    temp = gated_table(gind, :);
    uw_list = unique(string(temp.date_sampled));
    %uw_list = string(temp.date_sampled);
    uw_num = zeros(height(gated_table), 1);
    for i = 1:length(uw_list)
        uw_num((gated_table.cast==0 | gated_table.cast==-9999) & strcmp(string(gated_table.date_sampled), uw_list(i))) = i;
    end

    if strcmp(cruise, 'dy131') && ~contains(basepath, 'thawed')
        uw_num = zeros(height(gated_table), 1);
        t = split(temp.fcslist, 'Ex');
        exnum = extractBefore(t(:,2),'_');
        sizefrac = cell(size(exnum));
        sizefrac(contains(t(:,2), 'total')) = {'total'};
        sizefrac(contains(t(:,2), 'less5')) = {'less5'};
        sizefrac(contains(t(:,2), 'less20')) = {'less20'};
        [uw_num(gind), g1,g2] = findgroups(exnum, sizefrac);
    end

    [G, C, N, uw] = findgroups(gated_table.cast, gated_table.niskin, uw_num);
    ia = 1:max(G);

    CNTable = table;
    cruisevec = repmat(cruise, length(C), 1) ;
    CNTable.cruise = cruisevec;
    CNTable.cast = C;
    CNTable.niskin = N;

    %preallocate columns so we can get nans rather than zeros
    baccol = nan(height(CNTable), 1);
    syncol = baccol;
    procol = syncol;
    eukcol = syncol;
    lp_eukcol = syncol;
    hp_eukcol = syncol;
    vols = nan(height(CNTable), 6);


    for g = 1:max(G)


        temp = gated_table(G == g, :);

        CNTable.latitude(g) = temp.Latitude(1);
        CNTable.longitude(g) = temp.Longitude(1);
        CNTable.nearest_station(g) = temp.nearest_station(1);
        CNTable.salinity(g) = temp.salinity(1);
        CNTable.potemp090c(g) = temp.potemp090c(1);
        CNTable.depth_m(g) = temp.depth_m(1);
        CNTable.date_sampled(g) = temp.date_sampled(1);
        CNTable.date_processed(g) = temp.Date_processed(1);
        if ismember('r2r_event', temp.Properties.VariableNames)
            CNTable.r2r_event(g) = temp.r2r_event(1);
        end

        %pick file to count Syn
        ind = find(contains(temp.fcslist, 'phyto_PE') & ~cellfun(@isempty, temp.awsfilename));
        use_euk_for_syn = 0;
        if length(ind) == 1 %if only one fcs file of this type, use that.
            CNTable.Synfile(g) = temp.fcslist(ind);
            syncol(g) = temp.Syn_conc(ind);
            vols(g, 2) = temp.median_volumes_syn(ind);
        elseif ~isempty(ind) %if more than 1, get most recent
            timesince = [];
            for f = 1:length(ind)
                time = dir(fullfile(p.fpath, temp.fcslist{ind(f)}));
                time = datetime(time.date);
                timesince = [timesince datetime()-time];
            end
            [~,truind] = min(timesince);
            CNTable.Synfile(g) = temp.fcslist(ind(truind));
            syncol(g) = temp.Syn_conc(ind(truind));
            vols(g, 2) = temp.median_volumes_syn(ind(truind));

        elseif isempty(ind) %no special PE runs
            use_euk_for_syn = 1;
        end


        %now pick file to count Euks
        ind = find(contains(temp.fcslist, 'phyto_CHL') & ~contains(temp.fcslist, 'pro', 'IgnoreCase', true) & ~cellfun(@isempty, temp.awsfilename));
        if length(ind) == 1 %if only one fcs file of this type, use that.
            if use_euk_for_syn
                CNTable.Synfile(g) = temp.fcslist(ind);
                syncol(g) = temp.Syn_conc(ind);
                vols(g, 2) = temp.median_volumes_syn(ind);

            end
            CNTable.Eukfile(g) = temp.fcslist(ind);
            eukcol(g) = temp.Euk_conc(ind);
            lp_eukcol(g) = temp.LowP_Euk_conc(ind);
            hp_eukcol(g) = temp.HighP_Euk_conc(ind);
            vols(g, 1) = temp.median_volumes_euk(ind);
            vols(g, 5) = temp.median_volumes_low_pe_euk(ind);
            vols(g, 6) = temp.median_volumes_high_pe_euk(ind);
        elseif ~isempty(ind) %if more than 1, get most recent
            timesince = [];
            for f = 1:length(ind)
                time = dir(fullfile(p.fpath, temp.fcslist{ind(f)}));
                time = datetime(time.date);
                timesince = [timesince datetime()-time];
            end
            [~,truind] = min(timesince);
            if use_euk_for_syn
                CNTable.Synfile(g) = temp.fcslist(ind(truind)) ;
                syncol(g) = temp.Syn_conc(ind(truind));
                %vols(g, 2) = temp.median_volumes_syn(ind(truind));
            end
            CNTable.Eukfile(g) = temp.fcslist(ind(truind));
            eukcol(g) = temp.Euk_conc(ind(truind));
            lp_eukcol(g) = temp.LowP_Euk_conc(ind(truind));
            hp_eukcol(g) = temp.HighP_Euk_conc(ind(truind));
            vols(g, 1) = temp.median_volumes_euk(ind(truind));
            vols(g, 5) = temp.median_volumes_low_pe_euk(ind(truind));
            vols(g, 6) = temp.median_volumes_high_pe_euk(ind(truind));

        end

        %pick file to count Prochlorococcus
        ind = find(contains(temp.fcslist, 'pro', 'IgnoreCase', true) & ~cellfun(@isempty, temp.awsfilename));
        noprofile = 1;
        if isempty(ind)
            if strcmp('Eukfile',CNTable.Properties.VariableNames)
                CNTable.ProFile{g} = CNTable.Eukfile{g}; %if no pro run, check to see if there is a pro gate in euk run
                if ~isnan(temp.Pro_conc(ind))
                    procol(g) = temp.Pro_conc(ind);
                end
            end
            noprofile = 1;
        elseif length(ind) == 1 %if only one fcs file of this type, use that.
            CNTable.ProFile(g) =  temp.fcslist(ind);
            procol(g) = temp.Pro_conc(ind);
            vols(g, 4) = temp.median_volumes_pro(ind);
            if isnan(procol(g)) %if there was no pro gate, make conc zero
                procol(g) = 0;
            end
        else %if more than 1, get most recent
            timesince = [];
            for f = 1:length(ind)
                time = dir(fullfile(p.fpath, temp.fcslist{ind(f)}));
                time = datetime(time.date);
                timesince = [timesince datetime()-time];
            end
            [~,truind] = min(timesince);
            CNTable.ProFile(g) =  temp.fcslist(ind(truind));
            procol(g) = temp.Pro_conc(ind(truind));
            vols(g, 4) = temp.median_volumes_pro(ind(truind));
            if isnan(procol(g)) %if there was no pro gate, make conc zero
                procol(g) = 0;
            end
        end




        %pick file to count bacteria, and subtract prochlorococs
        %if a pro file exists and there is no pro gate in bacteria file, then
        %subtract concentration
        ind = find(contains(temp.fcslist, 'hbac') & ~cellfun(@isempty, temp.awsfilename));
        if length(ind) == 1 %if only one fcs file of this type, use that.
            CNTable.BacteriaFile(g) =  temp.fcslist(ind);

            %we want to add prochloro
            baccol(g) = temp.Bact_incl_pro_conc(ind);

            vols(g, 3) = temp.median_volumes_bact(ind);

        elseif ~isempty(ind) %if more than 1, get most recent
            timesince = [];
            for f = 1:length(ind)
                time = dir(fullfile(p.fpath, temp.fcslist{ind(f)}));
                time = datetime(time.date);
                timesince = [timesince datetime()-time];
            end
            [~,truind] = min(timesince);
            CNTable.BacteriaFile(g) =  temp.fcslist(ind(truind));

            baccol(g) = temp.Bact_incl_pro_conc(ind(truind));

            vols(g, 3) = temp.median_volumes_bact(ind(truind));

        end


    end


    CNTable.euk_per_ml = eukcol;
    CNTable.syn_per_ml = syncol;
    CNTable.pro_per_ml = procol;
    CNTable.bac_per_ml = baccol;

    CNTable.low_pe_euk_per_ml = lp_eukcol;
    CNTable.high_pe_euk_per_ml = hp_eukcol;


    CNTable.median_volumes_euk = vols(:,1);
    CNTable.median_volumes_syn = vols(:,2);
    CNTable.median_volumes_bact = vols(:,3);
    CNTable.median_volumes_pro = vols(:,4);
    CNTable.median_volumes_low_pe_euk = vols(:,5);
    CNTable.median_volumes_high_pe_euk = vols(:,6);


    save([p.outpath 'SummaryTable.mat'], 'CNTable')

    clearvars -except cruise step

end


%% Step 7 - Now we want to reformat the Summary table to be like EDI table with carbon counts etc

if step(7)

    %Load data for this cruise

    G = load([p.outpath 'Gated_Table.mat']);
    gated_table = G.gated_table;


    C = load([p.outpath 'SummaryTable.mat']);
    CNTable = C.CNTable;

    if strcmp(cruise, 'tn368')
        %if cruise == 'TN368'
        CNTable(CNTable.cast == 0, :) = []; %two bad entries, not matched to casts
    end

    EDI_table = table(CNTable.cruise, CNTable.cast, CNTable.niskin, CNTable.latitude, CNTable.longitude, CNTable.depth_m, CNTable.salinity, CNTable.potemp090c); %, CNTable.r2r_event);
    EDI_table.Properties.VariableNames = {'cruise'; 'cast'; 'niskin'; 'latitude'; 'longitude'; 'depth_m'; 'salinity'; 'potential_temperature_c'; 'r2r_event'};

    EDI_table.cruise = string(EDI_table.cruise); %helpful for merging tables when cruises are different lengtths

    %reformat dates so they aren't terrible
    if iscell(CNTable.date_sampled)
        dates1 = cell2mat(CNTable.date_sampled);
        EDI_table.date_sampled = datetime(dates1(:, 1:19), 'Format', 'yyyy-MM-dd HH:mm:ss');

    else
        EDI_table.date_sampled = CNTable.date_sampled;
    end


    dates2 = cell2mat(CNTable.date_processed);
    EDI_table.date_processed = datetime(dates2, 'Format', 'yyyy-MM-dd');


    % Go back to class files using Gated_table


    for i = 1:height(EDI_table);
        %first check syn file
        if ~isempty(CNTable.Synfile{i})
            filename = CNTable.Synfile{i};
            cfilename = regexprep(filename, '.fcs', '.mat');

            if ~exist([[p.classpath filesep cfilename]])
                EDI_table.syn_cells_per_ml(i) = NaN;
                EDI_table.syn_biovolume_concentration(i) = NaN;
                EDI_table.syn_carbon_concentration(i) = NaN; %divide by 1000 to get micrograms per Liter
                EDI_table.syn_volume_analyzed_ml(i) = NaN;
                EDI_table.syn_filename{i} = 'NaN';
                continue
            end

            C = load([p.classpath filesep cfilename]);

            volume = real(C.volume);
            carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
            carbon = real(carbon); %having issues with formatting, keeps having valus with + 0i.


            %match to gated table row
            gind = find(strcmp(gated_table.fcslist, filename));


            EDI_table.syn_cells_per_ml(i) = sum(C.class==2)./gated_table.Vol_analyzed_ml(gind);
            EDI_table.syn_biovolume_concentration(i) = nansum(volume(C.class==2))./gated_table.Vol_analyzed_ml(gind);
            EDI_table.syn_carbon_concentration(i) = nansum(carbon(C.class==2))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter
            EDI_table.syn_volume_analyzed_ml(i) = gated_table.Vol_analyzed_ml(gind);
            EDI_table.syn_filename(i) = CNTable.Synfile(i);

        else

            EDI_table.syn_cells_per_ml(i) = NaN;
            EDI_table.syn_biovolume_concentration(i) = NaN;
            EDI_table.syn_carbon_concentration(i) = NaN; %divide by 1000 to get micrograms per Liter
            EDI_table.syn_volume_analyzed_ml(i) = NaN;
            EDI_table.syn_filename{i} = 'NaN';
        end

        % done with Syn, move on to Euks
        if ~isempty(CNTable.Eukfile{i})

            filename = CNTable.Eukfile{i};
            cfilename = regexprep(filename, '.fcs', '.mat');

            if ~exist([[p.classpath filesep cfilename]])
                EDI_table.redeuk_leq_2um_cells_per_ml(i) = NaN;
                EDI_table.redeuk_leq_2um_biovolume_concentration(i) =  NaN;
                EDI_table.redeuk_leq_2um_carbon_concentration(i) =  NaN;
                % <= 3
                EDI_table.redeuk_leq_3um_cells_per_ml(i) =  NaN;
                EDI_table.redeuk_leq_3um_biovolume_concentration(i) =  NaN;
                EDI_table.redeuk_leq_3um_carbon_concentration(i) =  NaN;
                % <= 5
                EDI_table.redeuk_leq_5um_cells_per_ml(i) = NaN;
                EDI_table.redeuk_leq_5um_biovolume_concentration(i) =  NaN;
                EDI_table.redeuk_leq_5um_carbon_concentration(i) =  NaN;
                % <= 10
                EDI_table.redeuk_leq_10um_cells_per_ml(i) =  NaN;
                EDI_table.redeuk_leq_10um_biovolume_concentration(i) =  NaN;;
                EDI_table.redeuk_leq_10um_carbon_concentration(i) =  NaN;
                % <= 20
                EDI_table.redeuk_leq_20um_cells_per_ml(i) =  NaN;
                EDI_table.redeuk_leq_20um_biovolume_concentration(i) =  NaN;
                EDI_table.redeuk_leq_20um_carbon_concentration(i) =  NaN;

                EDI_table.redeuk_volume_analyzed_ml(i) =  NaN;
                EDI_table.redeuk_filename{i} = 'NaN';
                continue
            end

            C = load([p.classpath filesep cfilename]);

            volume = real(C.volume);
            carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
            carbon = real(carbon); %having issues with formatting, keeps having valus with + 0i.

            gind = find(strcmp(gated_table.fcslist, filename));


            EukSizes = [0 2 3 5 10 20];
            %size fractions by diameter
            diam = (volume*3/4/pi).^(1/3)*2; %equivalent spherical diam, micrometers

            %first < 2
            bin_particle_ind = find(diam'<=2 & C.class==1)';

            EDI_table.redeuk_leq_2um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
            EDI_table.redeuk_leq_2um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind);
            EDI_table.redeuk_leq_2um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter

            % <= 3
            bin_particle_ind = find(diam'<=3 & C.class==1)';
            EDI_table.redeuk_leq_3um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
            EDI_table.redeuk_leq_3um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind);
            EDI_table.redeuk_leq_3um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter

            % <= 5
            bin_particle_ind = find(diam'<=5 & C.class==1)';
            EDI_table.redeuk_leq_5um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
            EDI_table.redeuk_leq_5um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind);
            EDI_table.redeuk_leq_5um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter

            % <= 10
            bin_particle_ind = find(diam'<=10 & C.class==1)';
            EDI_table.redeuk_leq_10um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
            EDI_table.redeuk_leq_10um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind);
            EDI_table.redeuk_leq_10um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter

            % <= 20
            bin_particle_ind = find(diam'<=20 & C.class==1)';
            EDI_table.redeuk_leq_20um_cells_per_ml(i) = length(bin_particle_ind)./gated_table.Vol_analyzed_ml(gind); %counts over volume
            EDI_table.redeuk_leq_20um_biovolume_concentration(i) = nansum(volume(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind);
            EDI_table.redeuk_leq_20um_carbon_concentration(i) = nansum(carbon(bin_particle_ind))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter


            EDI_table.redeuk_volume_analyzed_ml(i) = gated_table.Vol_analyzed_ml(gind);
            EDI_table.redeuk_filename(i) = CNTable.Eukfile(i);

        else
            EDI_table.redeuk_leq_2um_cells_per_ml(i) = NaN;
            EDI_table.redeuk_leq_2um_biovolume_concentration(i) =  NaN;
            EDI_table.redeuk_leq_2um_carbon_concentration(i) =  NaN;
            % <= 3
            EDI_table.redeuk_leq_3um_cells_per_ml(i) =  NaN;
            EDI_table.redeuk_leq_3um_biovolume_concentration(i) =  NaN;
            EDI_table.redeuk_leq_3um_carbon_concentration(i) =  NaN;
            % <= 5
            EDI_table.redeuk_leq_5um_cells_per_ml(i) = NaN;
            EDI_table.redeuk_leq_5um_biovolume_concentration(i) =  NaN;
            EDI_table.redeuk_leq_5um_carbon_concentration(i) =  NaN;
            % <= 10
            EDI_table.redeuk_leq_10um_cells_per_ml(i) =  NaN;
            EDI_table.redeuk_leq_10um_biovolume_concentration(i) =  NaN;;
            EDI_table.redeuk_leq_10um_carbon_concentration(i) =  NaN;
            % <= 20
            EDI_table.redeuk_leq_20um_cells_per_ml(i) =  NaN;
            EDI_table.redeuk_leq_20um_biovolume_concentration(i) =  NaN;
            EDI_table.redeuk_leq_20um_carbon_concentration(i) =  NaN;

            EDI_table.redeuk_volume_analyzed_ml(i) =  NaN;
            EDI_table.redeuk_filename{i} = 'NaN';


        end

        % finally heterotrophic bacteria
        if ~isempty(CNTable.BacteriaFile{i})

            filename = CNTable.BacteriaFile{i};
            cfilename = regexprep(filename, '.fcs', '.mat');

            if ~exist([[p.classpath filesep cfilename]])
                EDI_table.hetprok_cells_per_ml(i) = NaN;
                EDI_table.hetprok_carbon_concentration(i)
                EDI_table.hetprok_volume_analyzed_ml(i) = NaN;
                EDI_table.hetprok_filename{i} = 'NaN';
                continue
            end

            C = load([p.classpath filesep cfilename]);

            volume = real(C.volume);
            carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
            carbon = real(carbon); %having issues with formatting, keeps having valus with + 0i.

            gind = find(strcmp(gated_table.fcslist, filename));


            EDI_table.hetprok_cells_per_ml(i) = CNTable.bac_per_ml(i); % Do NOT go back to class file, since we had to account for time gate
            %EDI_table.hetprok_biovolume_concentration(i) = nansum(volume(C.class==3))./gated_table.Vol_analyzed_ml(gind);
            %EDI_table.hetprok_carbon_concentration(i) = nansum(carbon(C.class==3))./gated_table.Vol_analyzed_ml(gind)./1000; %divide by 1000 to get micrograms per Liter


            %doing carbon concentration per cell conversion based on Lee & Furhman 1987
            EDI_table.hetprok_carbon_concentration(i) = EDI_table.hetprok_cells_per_ml(i).* 20  * 1e-6; %to convert to micrograms per liter

            EDI_table.hetprok_volume_analyzed_ml(i) = gated_table.Vol_analyzed_ml(gind);
            EDI_table.hetprok_filename(i) = CNTable.BacteriaFile(i);

        else
            EDI_table.hetprok_cells_per_ml(i) = NaN;
            EDI_table.hetprok_carbon_concentration(i)
            EDI_table.hetprok_volume_analyzed_ml(i) = NaN;
            EDI_table.hetprok_filename{i} = 'NaN';
        end

    end
    save([p.outpath 'EDI_table.mat'], 'EDI_table')
    disp([p.outpath 'EDI_table.mat'])


end


%% function needed for Step 2

    function make_figure_aws(fcsdat, fcshdr, gate_assignments, polygon_vars, polygon_vals, gate_names, figpath);

        timegateind = find(contains(gate_names, 'time_to_include'));  %remove this gate
        gate_names(timegateind) = [];
        gate_assignments(timegateind, :) = [];


        %and remove time polygon
        timepolyind = find(contains(polygon_vars(1,:), 'Time'));
        polygon_vars(:, timepolyind) = [];
        polygon_vals(:, timepolyind) = [];

        vars_plotted = [];
        clf
        for i = 1:length(polygon_vars) %make an axis for each polygon
            x_ch = strmatch(polygon_vars{1, i}, {fcshdr.par.name});

            if isempty(polygon_vars{2,i})
                y_ch = 12; %just need a var to plot for time gates
            else
                y_ch = strmatch(polygon_vars{2,i}, {fcshdr.par.name});
            end

            if isempty(vars_plotted)
                vars_plotted = [vars_plotted [x_ch y_ch]'];
            elseif sum(ismember(vars_plotted, [x_ch y_ch]')) ~=2
                vars_plotted = [vars_plotted [x_ch y_ch]'];
            end

        end

        numpanels = size(vars_plotted, 2);

        for i = 1:height(gate_assignments); %plot variabel combos for each gate

            for j = 1:numpanels
                x_ch = vars_plotted(1, j);
                y_ch = vars_plotted(2, j);

                subplot(height(gate_assignments), numpanels, ((i-1).*numpanels+j))
                scatter(fcsdat(:, x_ch), fcsdat(:, y_ch), '.')
                hold on
                scatter(fcsdat(logical(gate_assignments(i,:)), x_ch), fcsdat(logical(gate_assignments(i,:)), y_ch), '.')

                xlabel(fcshdr.par(x_ch).name)
                ylabel(fcshdr.par(y_ch).name)

                set(gca, 'YScale', 'log')
                set(gca, 'XScale', 'log')

            end

            title(gate_names{i}, 'interpreter', 'none')

        end

        %add a title for the figure
        subplot(height(gate_assignments), numpanels, 1)
        title(fcshdr.filename, 'interpreter', 'none')
        figname = regexprep(fcshdr.filename, '.fcs', '.jpg');
        print(strcat(figpath, '\', figname), '-djpeg')

    end

end
