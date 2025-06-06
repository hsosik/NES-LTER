 
%Trying to pull together all steps into one:
    % Process_preserved_pt1
    % Proces_Preserved_samples_with_AWS
    % Get_metadata_for_Gated_Table
    % Classify_from_Gated_Tables
    % create class files with volumes 
    % Convert_Gated_Table_to_CNTable

%% INPUTS 
clear all 
fclose('all');

% % Manually choose cruise to process

%basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20210501_DY131\';
basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\MVCO\preserved\';
%cruisename = 'EN688';
cruisename = 'MVCO';
hierarchical_gates = 'False';  %set to 'True' or 'False'; 

%%
%restpath = '\\sosiknas1\Lab_data\Attune\cruise_data\20190725_HB1907\preserved\bottle_environmental_data_partial.csv'; 
%restpath = 'https://nes-lter-data.whoi.edu/api/ctd/ar61b/'
%only relevat elog  path if discrete underway or bucket samples were taken,
%their position is from the elog
restpath = '\\sosiknas1\Lab_data\MVCO\mvcodata\matlab_files\FCM_sample_list.mat';
elogpath = '';%'\\sosiknas1\Lab_data\LTER\20201013_EN657\eLog\R2R_ELOG_EN657_FINAL_EVENTLOG_20201018_134037.csv'; %set to '' if there are no discrete underway samples 

% uw_fullname also only needed if discret UW or bucket samples were taken
uw_fullname = ''; %'https://nes-lter-data.whoi.edu/api/underway/en657.csv';


%Set all steps to 1 if starting from beginning 
Step1 = 0; %make FCSList
Step2 = 0; %go look at AWS files to find gate assignments 
Step3 = 0; % add metadata to gated table
Step4 = 0; % classify using gate_table
Step5 = 1; %size calibrate and create class files 
Step6 = 1; %convert gated table to Summary table 
Step7 = 1; %Reformat Summary Table to have EDI headers


%% Set up 

% some file structure setup  
fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'outputs' filesep];
classpath = [outpath 'class' filesep];
awspath = [basepath filesep 'aws\'];

if ~exist(outpath, 'dir')
    mkdir(outpath)
end
if ~exist(classpath, 'dir')
    mkdir(classpath)
end


save([outpath '\Processing_variables.mat'])

%% Step 1 - make FCSlist

if Step1

filetype2exclude = {'Rinses', 'skip', 'xxx', 'qwater', 'dilution', 'autoclaved', '2024';}; % 'Sample(',
beadfiles2include = {'FCB_bead'};

%look at files in FCS folder
fcslist= dir(fullfile(fpath, '*fcs'));
fcslist = {fcslist.name}';
samplelist = fcslist;
%remove filetype2exclude from filelist to generate samplelist
for iii = 1:length(filetype2exclude)
    t = contains(samplelist, filetype2exclude{iii});
    if ~isempty(t)
        samplelist(t, :) = [];
    end
end
clear t iii fcslist


%headers = {'fcslist', 'Event', 'Niskin'};
headers = {'fcslist', 'Event', 'Niskin'};

if ~exist([outpath '/FCSList.mat']) %if no table has been started
    %initialize table
    FCSList = cell2table(cell(0,3), 'VariableNames', headers);

else
    load([outpath '/FCSList.mat'])

    %find list of samples not yet included in table
    samplelist = setdiff(samplelist, FCSList.fcslist);
end
Table_Add = table(samplelist, 'VariableNames', {'fcslist'});

%go through files and parse filenames
%to extract Event and niskin numbers

for ii = 1:length(samplelist)
    filename = samplelist{ii};
    Table_Add.fcslist(ii) = {filename};
    %s = strfind(filename, '_C');
    s = strfind(filename, '_MV');
    if isempty(s)
        continue
    end
    s = s(end); %get index of last _ in filename
    %filenameshort = filename(s+1:end-3);
    filenameshort = filename(s+1:end-4); %for MVCO processing
    n = strfind(filenameshort, 'N');
    if ~isempty(str2num(filenameshort(3:n-1)))
        Table_Add.Event(ii) = str2num(filenameshort(3:n-1));
    end
    if ~isempty(str2num(filenameshort(n+1:n+2)))
        Table_Add.Niskin(ii) = str2num(filenameshort(n+1:n+2));
    elseif strcmp(filenameshort(n+2), '(')
        Table_Add.Niskin(ii) = str2num(filenameshort(n+1));
    end
    clear n s filenameshort filename
end

FCSList = [FCSList; Table_Add];

save([outpath 'FCSList.mat'], 'FCSList')


% add trigger settings and info to FCSlist

Vol_analyzed_ml = nan(height(FCSList), 1);
Vol_ml = nan(height(FCSList), 1); 
Date_processed = cell(height(FCSList), 1);
trigger_1 = Date_processed;
trigger_2 = Date_processed;
trigger_hv1 = Vol_analyzed_ml;
trigger_hv2 = Vol_analyzed_ml;

for i = 1:height(FCSList)
    filename = FCSList{i, 1};
    filename = [fpath filename{:}];
    [fcsdat, fcshdr] = fca_readfcs(filename);

    Date_processed{i} = fcshdr.date;
    Vol_analyzed_ml(i) = .8.*fcshdr.VOL ./ 1e6; %changed to 80% based on uniform time gate
    Vol_ml(i) = fcshdr.VOL ; 

    trigger_1{i} = fcshdr.trigger1;
    trigger_2{i} = fcshdr.trigger2;

    %parse to find trigger ch numbers
    s = strfind(fcshdr.trigger1, '_');
    ch1 = fcshdr.trigger1(s+1:s+3);
    chnum = strmatch([ch1 '-H'], {fcshdr.par.name});
    trigger_hv1(i) = fcshdr.par(chnum).hv;

    s = strfind(fcshdr.trigger2, '_');
    ch2 = fcshdr.trigger2(s+1:s+3);
    chnum = strmatch([ch2 '-H'], {fcshdr.par.name});
    if ~isempty(chnum)
        trigger_hv2(i) = fcshdr.par(chnum).hv;
    end

end

FCSList.Date_processed = Date_processed;
FCSList.Vol_analyzed_ml = Vol_analyzed_ml;
FCSList.trigger_1 = trigger_1;
FCSList.trigger_2 = trigger_2;
FCSList.trigger_hv1 = trigger_hv1;
FCSList.trigger_hv2 = trigger_hv2;

save([outpath 'FCSList.mat'], 'FCSList')

clearvars -except basepath restpath fpath outpath classpath awspath cruisename elogpath uw_fullname cruisename hierarchical_gates Step1 Step2 Step3 Step4 Step5 Step6 Step7

end

%% Step 2 - go look at AWS files to find gate assignments 

if Step2
T = load([outpath '\FCSlist.mat']); 
T = T.FCSList; 
gated_table = T; 
no_aws_files = []; 

% go through files, for each file, find matching aws file
% 
% Make one figure for each Niskin Bottle, reporting depth and total counts
%
%Make a list of fcs files that didn%t have an aws files

%first get runtype directory names 
runtypes = dir(awspath); runtypes = struct2table(runtypes); runtypes = string(runtypes.name); 
runtypes = runtypes(~startsWith(runtypes, '.')); 

for i = 1:height(T) 
    filename = T.fcslist{i} 
    %step 1 find, aws file

    awsfile = [];
    awsfilename = [];
    filetypevec = ""; 

    for j = 1:length(runtypes)
        if contains(filename, runtypes(j))
            if length(runtypes{j}) > length(filetypevec) %check for things like CHLS_SSC_pro_hi and pick the longest match
                filetypevec = runtypes{j}; 
                k = j; 

            awslist = dir(strcat(awspath, '', runtypes(j), '\*.aws'));
            awslist = struct2table(awslist); awslist = string(awslist.name);
            %right now this only works if there are only two digit Events
            %and niskins

            if strcmp(cruisename,'SR1812')
                temp = split(filename, '_');
                ind = find(awslist == [temp{end-1} '.aws']);
            elseif strcmp(cruisename, 'DY131') && T.Event(i) == 0
                temp = split(filename, 'Ex');
                ind = find(awslist == ['Ex' temp{end}(1:end-4) '.aws']);
            elseif T.Event(i) == 0 && T.Niskin(i) == 0 %need case for UW data
                uwname = split(filename, '_');
                uwname = regexprep(uwname{end}, '.fcs', '.aws');
                ind = find(awslist == uwname);
            else
                ind = find(awslist == strcat("C", num2str(T.Event(i), '%02.f'), 'N', num2str(T.Niskin(i), '%02.f'), '.aws'));
            end

            awsfilename = awslist(ind) ;
            awsfile = strcat(awspath, '', runtypes(j), '\', awslist(ind));

            figpath = strcat(outpath , 'figs\', runtypes(j)) ;
            if ~exist(figpath, 'dir')
                mkdir(figpath)
            end

            end
        end
    end


        %track whether no aws file chosen 
        if isempty(awsfilename)
            
            temp = split(filename, '_'); 
            temp = temp{end}; 

            awsfilename = regexprep(temp, 'fcs', 'aws'); 
            if isfile(strcat(awspath, '', runtypes(k), '\', awsfilename)); 
                  awsfile = strcat(awspath, '', runtypes(k), '\', awsfilename); 
            else
                awsfilename = []; 
                %if you want manual inputs
                %awsfilename = input(['Please choose aws file for this fcs: ' filename 'If none or does not matter hit enter'], 's');
                %awsfile = strcat(awspath, '', runtypes(k), '\', awsfilename); 
            end 
        end
        if isempty(awsfilename) %still empty 
            no_aws_files = [no_aws_files; string(filename)]; 
        else
            %if we have an aws file proceed with gating
            [fcsdat, fcshdr]  = fca_readfcs([fpath filename]);
            %if isempty(awsfile)%I truly don't understand how this can still be empty
             %               awsfile = strcat(awspath, '', runtypes(k), '\', awslist(ind));
            %end
            if strcmp(hierarchical_gates, 'True')
                 [gate_assignments, polygon_names, polygon_vars, polygon_vals, gate_names, gate_logic_legible, parent_logic] = ApplyAWSgates_hierarchical(awsfile, fcsdat, fcshdr);
            elseif strcmp(hierarchical_gates, 'False')
                [gate_assignments, polygon_names, polygon_vars, polygon_vals, gate_names, gate_logic_legible] = ApplyAWSgates2(awsfile, fcsdat, fcshdr); 
            end

            gated_table.awsfilename{i} = strcat(filetypevec, '/', awsfilename); 
            gated_table.gate_names{i} = gate_names; 
            gated_table.gate_assignments{i} = gate_assignments; 
            gated_table.gate_logic{i} = gate_logic_legible; 
            gated_table.polygon_names{i} = polygon_names; 
            gated_table.polygon_vars{i} = polygon_vars; 
            gated_table.polygon_vals{i} = polygon_vals; 
%%
            make_figure_aws(fcsdat, fcshdr, gate_assignments, polygon_vars, polygon_vals, gate_names, figpath);


        end
    
    end

    save([outpath 'Gated_Table.mat'], 'gated_table', 'no_aws_files', 'hierarchical_gates')

clearvars -except basepath restpath fpath outpath classpath awspath cruisename elogpath uw_fullname cruisename hierarchical_gates Step1 Step2 Step3 Step4 Step5 Step6 Step7

end

%% Step 3 - add metadata to gated table

if Step3

T = load([outpath '\FCSlist.mat']);
T = T.FCSList;

G = load([outpath '\Gated_Table.mat']);
gated_table = G.gated_table;
no_aws_files = G.no_aws_files; 

load(restpath)
temp = table;
ev = char(data(:,2));
temp.Event = str2num(ev(:,6:8));
for i = 1:length(temp.Event)
if strcmp(ev(i,10), '0')
temp.Niskin(i) = str2num(ev(i, 10:11));
else
    temp.Niskin(i) = -9999;
end %clear ev;
end;
temp.depsm = str2num(char(string(data(:,7))));
temp.Latitude_decimalDeg = str2num(char(string(data(:,5))));
temp.Longitude_decimalDeg = str2num(char(string(data(:,6))));
temp.datetime= datetime(datevec(matdate), 'Format', 'uuuu-MM-dd HH:mm:ss');




BTL = temp; 

Eventlist = unique(gated_table.Event);

for f = 1:height(gated_table)

    if gated_table.Event(f) == 0 || gated_table.Event(f) == -9999
        continue
    end


    temp = BTL(BTL.Event == gated_table.Event(f), :); 

    %first get Event metdata
    gated_table.Latitude(f) = temp.Latitude_decimalDeg(1);
    gated_table.Longitude(f) = temp.Longitude_decimalDeg(1);
    gated_table.date_sampled(f) = temp.datetime(1);

   
    %no station names for SPIROPA cruises
    gated_table.nearest_station{f} = '';


    %now add depths for each niskin number
    if sum(temp.Event == gated_table.Event(f) & temp.Niskin == gated_table.Niskin(f)) > 0
        %gated_table.salinity(f) = temp.sal00(temp.Event == gated_table.Event(f) & temp.Niskin == gated_table.Niskin(f));
        %gated_table.potemp090c(f) = temp.potemp090c(temp.Event == gated_table.Event(f) & temp.Niskin == gated_table.Niskin(f));
        %EP reprocessing HB1907 for quick check of problem with PE runs May
        %2024
        %gated_table.depth_m(f) = temp.depth(f);
        %temp2 = temp.depsm(temp.Event == gated_table.Event(f) & temp.Niskin
        %== gated_table.Niskin(f)); 
        %gated_table.depth_m(f) = temp2(1) 
        %these two lines are ueful only if you are
        %working with a metadata table that is actually a nutrient table
        %and has more than one sample per Event and niskin. 
        gated_table.depth_m(f) = unique(temp.depsm(temp.Event == gated_table.Event(f) & temp.Niskin == gated_table.Niskin(f)));
    end


end

%     gated_table.Longitude(gated_table.Latitude == 0) = NaN;
%     gated_table.salinity(gated_table.Latitude == 0) = NaN;
%     gated_table.potemp090c(gated_table.Latitude == 0) = NaN;
%     gated_table.depth_m(gated_table.Latitude == 0) = NaN;
%     gated_table.Latitude(gated_table.Latitude == 0) = NaN;




save([outpath '\Gated_Table.mat'], 'gated_table', 'no_aws_files', 'hierarchical_gates'); 


no_aws_files = G.no_aws_files; 


clearvars -except basepath restpath fpath outpath classpath awspath cruisename elogpath uw_fullname  hierarchical_gates Step1 Step2 Step3 Step4 Step5 Step6 Step7
end


%% Step 4 - classify using gate_table

if Step4

G = load([outpath '\Gated_Table.mat']);
gated_table = G.gated_table; 
no_aws_files = G.no_aws_files; 

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
    [fcsdat, fcshdr]  = fca_readfcs([fpath filename]);

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
            gate_num_1 = find(strcmp(gate_names, 'euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_1, :)) = 1; 
    end
    if sum(strcmp(gate_names, 'Euk'))
            gate_num_1 = find(strcmp(gate_names, 'Euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_1, :)) = 1; 
    end

    if sum(strcmp(gate_names, 'LowP_Euk'))
            gate_num_5 = find(strcmp(gate_names, 'LowP_Euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_5, :)) = 5; 
    end

    if sum(strcmp(gate_names, 'LowPE_Euk'))
            gate_num_5 = find(strcmp(gate_names, 'LowPE_Euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_5, :)) = 5; 
    end
    if sum(strcmp(gate_names, 'HighP_Euk'))
            gate_num_6 = find(strcmp(gate_names, 'HighP_Euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_6, :)) = 6; 
    end
    if sum(strcmp(gate_names, 'HighPE_Euk'))
            gate_num_6 = find(strcmp(gate_names, 'HighPE_Euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_6, :)) = 6; 
    end

    %next syn 
    if sum(strcmp(gate_names, 'Syn'))
            gate_num_2 = find(strcmp(gate_names, 'Syn')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2; 
    end
    if sum(strcmp(gate_names, 'syn'))
            gate_num_2 = find(strcmp(gate_names, 'syn')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2; 
    end
   if sum(strcmp(gate_names, 'SYN'))
            gate_num_2 = find(strcmp(gate_names, 'SYN')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2; 
    end


    %next bacteria, may include prochlorococcus, but those will be
    %reassigned below
    if sum(strcmp(gate_names, 'Bacteria'))
            gate_num_b = find(strcmp(gate_names, 'Bacteria')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3; 
    end
    if sum(strcmp(gate_names, 'bacteria'))
            gate_num_b = find(strcmp(gate_names, 'bacteria')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3; 
    end
    if sum(strcmp(gate_names, 'BacPro'))
            gate_num_b = find(strcmp(gate_names, 'BacPro')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3; 
    end



    %next assign pro, allowed to overwrite bacteria label 
    if sum(strcmp(gate_names, 'pro'))
            gate_num_4 = find(strcmp(gate_names, 'pro')); 
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4; 
    end
    if sum(strcmp(gate_names, 'Pro'))
            gate_num_4 = find(strcmp(gate_names, 'Pro')); 
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4; 
    end
    if sum(strcmp(gate_names, 'Euk_sm'))
            gate_num_4 = find(strcmp(gate_names, 'Euk_sm')); 
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4; 
    end
    if sum(strcmp(gate_names, 'Proc'))
            gate_num_4 = find(strcmp(gate_names, 'Proc')); 
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4; 
    end


    %ok now save class assignments 

    class{i} = class_i;

    
    %onto calculating concentrations
    concent_i = nan(1, 6);

    sizefrac = gated_table.gate_logic{i}; 

    % this is deleted since now all gates have post hoc 80% time window
    % added. 
    % if sum(gates_that_use_time) ~= 0 %go get volume included in time_to_include
    %     p = gated_table.polygon_vars{i}; 
    %     emptyCells = cellfun(@isempty,p) ;
    %     [a,b] = find(emptyCells); 
    %     p{a,b} = '';
    %     [a,b] = find(cellfun(@(x) contains(x, 'Time'), p));
    %     pv = gated_table.polygon_vals{i} ;
    %     vals = pv{b};
    %     t_val = vals(:,1); %this is the actual value range used for the Time variable


        % matdate1 = [fcshdr.date ' ' fcshdr.starttime];
        % matdate1 = datenum(matdate1, 'dd-mmm-yyyy HH:MM:SS');
        % matdate2 = [fcshdr.date ' ' fcshdr.stoptime];
        % matdate2 = datenum(matdate2, 'dd-mmm-yyyy HH:MM:SS');
        % 
        % fulltime = 1000*etime(datevec(matdate2), datevec(matdate1)); %this should be 60 seconds * 1000 but just to make sure
        % 
        % fraction_of_time = (t_val(2) - t_val(1))./fulltime;

    %volume_considered = fraction_of_time.*gated_table.Vol_analyzed_ml(i);
    %Vol_analyzed_ml column has already been adjusted to 80% total sample volume


    if exist('gate_num_1')
            concent_i(1) = sum(class_i == 1)./gated_table.Vol_analyzed_ml(i);
    end

    if exist('gate_num_2')
            concent_i(2) = sum(class_i == 2)./gated_table.Vol_analyzed_ml(i);
    end

    if exist('gate_num_b')
            concent_i(3) = sum(class_i == 3 | class_i ==4)./gated_table.Vol_analyzed_ml(i); %here is where we combine Bacteria and Syn class assignments!
    end

    if exist('gate_num_4')
            concent_i(4) = sum(class_i == 4)./gated_table.Vol_analyzed_ml(i);

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
            concent_i(5) = sum(class_i == 5)./gated_table.Vol_analyzed_ml(i);
    end

    if exist('gate_num_6')        
            concent_i(6) = sum(class_i == 6)./gated_table.Vol_analyzed_ml(i);
    end

    
    disp(gated_table(i,:))
    disp(concent_i)

    %now save concentrations
    euk_conc(i) = concent_i(1);
    syn_conc(i) = concent_i(2);
    bact_incl_pro_conc(i) = concent_i(3) ;
    pro_conc(i) = concent_i(4);
    lp_euk_conc(i) = concent_i(5); 
    hp_euk_conc(i) = concent_i(6); 
   
        clear gate_num_1 gate_num_2 gate_num_4 gate_num_b gate_assign_i class_i gate_num_5 gate_num_6
end

gated_table.class = class'; 
gated_table.Euk_conc = euk_conc; 
gated_table.LowP_Euk_conc = lp_euk_conc; 
gated_table.HighP_Euk_conc = hp_euk_conc; 

gated_table.Syn_conc = syn_conc; 
gated_table.Bact_incl_pro_conc = bact_incl_pro_conc; 
gated_table.Pro_conc = pro_conc; 

save([outpath '\Gated_Table.mat'], 'gated_table', 'no_aws_files', 'cut_off_pro_pop', 'hierarchical_gates'); 


clearvars -except basepath restpath fpath outpath classpath awspath  cruisename hierarchical_gates Step1 Step2 Step3 Step4 Step5 Step6 Step7


end

%% Step 5 - size calibrate and create class files 

if Step5

load('\\sosiknas1\Lab_data\Attune\cruise_data\beads\FCB_bead_mix_experiment_settings\between_cruises\outputs\beadstat.mat')


G = load([outpath '\Gated_Table.mat']);
gated_table = G.gated_table; 
no_aws_files = G.no_aws_files; 

if isfield(G, 'cut_off_pro_pop')
    cut_off_pro_pop = G.cut_off_pro_pop; 
else 
    cut_off_pro_pop = NaN; 
end

median_volumes = nan(height(gated_table), 6); 

DIM = 'A';
saverpath = [classpath 'calibration']; 
if ~exist([classpath 'calibration'], 'dir')
    mkdir([classpath 'calibration'])
end

 
 joint_table = []; 
    
     for counti = 1:height(gated_table)
         qc_warning = 0; 

         if isempty(gated_table.class{counti})
             continue
         end

        %load corresponding fcs file
        filename = gated_table.fcslist{counti}; 
        [fcsdat,fcshdr] = fca_readfcs([fpath filename]);
         
        classfilename = [classpath regexprep(filename, '.fcs', '.mat')];
        class = gated_table.class{counti};
         notes = "'Euks == 1'; 'Syn == 2'; 'Bacteria & Pro == 3'; 'Pro == 4'; 'LowPE_Euks = 5'; 'HighPE_Euks - 6'";

        ssc_ch_num = strmatch(['SSC-' DIM], {fcshdr.par.name});
        gl1_ch_num = strmatch(['GL1-' DIM], {fcshdr.par.name});
        %bl3_ch_num = strmatch(['BL3-' DIM], {fcshdr.par.name});

        file_hv = fcshdr.par(ssc_ch_num).hv; %heidi    
    
        if file_hv == 220 %useful for working with bacteria size calibration
            %keyboard
          %  continue
        end
 
        gl1_vals = fcsdat(:, gl1_ch_num);
        ssc_value = fcsdat(:,ssc_ch_num); 
   
        ssch_ch_num = strmatch(['SSC-H'], {fcshdr.par.name});
        ssc_value(ssc_value<=0) = fcsdat(ssc_value<=0, ssch_ch_num); %replace negative A values with H value as proxy
                   
        l_bound = 1e3;
        r_bound = 9e3; %fixed
    
        if r_bound < l_bound; 
            qc_warning =1; 
        end
    
        ind_to_fit = find(gl1_vals>l_bound & gl1_vals<r_bound); %not a super robust way to choose the range to fit
   
        LM = fitlm(gl1_vals(ind_to_fit), ssc_value(ind_to_fit));     
        intercept = LM.Coefficients.Estimate(1); 
        slope = LM.Coefficients.Estimate(2); 
    
        if isnan(intercept) %if data is bad, linear model doesn't work 
            qc_warning = 1; 
            scatter_value = fcsdat(:,ssc_ch_num); 
        elseif LM.Rsquared.Adjusted < .8
            qc_warning = 1; 
        end
    
    new_ssc_vals = ssc_value;
    new_ssc_vals(gl1_vals>r_bound) = 10.^[intercept + slope.*(gl1_vals(gl1_vals>r_bound))]; 
    
    %dither according to heidi's example
    t = find((gl1_vals <= r_bound & ssc_value >= 10.^(intercept + slope.*(r_bound))) | (gl1_vals >= r_bound & ssc_value <= 10.^(intercept + slope.*(r_bound)))); 
    new_ssc_vals(t(1:2:end)) = 10.^[intercept + slope.*(gl1_vals(t(1:2:end)))];
    

    calibrate_info = table; 
    calibrate_info.filename = string(filename);
    calibrate_info.ssc_ch_num = ssc_ch_num; 
    calibrate_info.qc = qc_warning; 
    calibrate_info.intercept = intercept; 
    calibrate_info.slope = slope; 
    calibrate_info.rsquared = LM.Rsquared.Adjusted; 
    calibrate_info.numpoints = length(ind_to_fit); 
    calibrate_info.rightbound = r_bound; 
    
    joint_table = [joint_table; calibrate_info];

        figure(98)
        plot(gl1_vals, ssc_value, '.')
        xlabel('GL1 - low sensitivity')
        ylabel('SSC - high sensitivity')
        hold on 
        plot([1 max(gl1_vals)], [intercept (intercept + slope*(max(gl1_vals)))])
        plot(gl1_vals(t), ssc_value(t), '.b')
        plot(gl1_vals, new_ssc_vals, 'g.')
        plot(gl1_vals(ind_to_fit), ssc_value(ind_to_fit), 'r.')
        title({filename; [num2str(qc_warning)]}, 'Interpreter', 'none')
        axis([0 1e5 0 1e6])

        if ~contains(filename, 'hbac') & ~contains(filename, 'pro')
            print(figure(98), [saverpath filesep regexprep(filename, '.fcs', '.png')], '-dpng')
        end
        clf(98)
                

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


    save([outpath '\Gated_Table.mat'], 'gated_table', 'no_aws_files', 'cut_off_pro_pop', 'hierarchical_gates');


clearvars -except basepath restpath fpath outpath classpath awspath cruisename hierarchical_gates Step1 Step2 Step3 Step4 Step5 Step6 Step7

end

%% Step 6 - convert gated table to Summary table 

if Step6
 
G = load([outpath '\Gated_Table.mat']);
gated_table = G.gated_table; 


%gated_table(contains(gated_table.fcslist, 'lower_thresh'), :) = []; 
%gated_table(gated_table.Event == 0, :) = []; 


%first some counting of underway samples
    gind = (gated_table.Event==0 | gated_table.Event==-9999);
    %gind = (gated_table.Event==-9999);%comment out after, Emily using for EN688 live discretes to ignore these samples 
    temp = gated_table(gind, :);
    uw_list = unique(string(temp.date_sampled)); 
    %uw_list = string(temp.date_sampled); 
    uw_num = zeros(height(gated_table), 1); 
    for i = 1:length(uw_list)
        uw_num((gated_table.Event==0 | gated_table.Event==-9999) & strcmp(string(gated_table.date_sampled), uw_list(i))) = i; 
    end

    if strcmp(cruisename, 'DY131') && ~contains(basepath, 'thawed')
        uw_num = zeros(height(gated_table), 1); 
        t = split(temp.fcslist, 'Ex');
        exnum = extractBefore(t(:,2),'_');
        sizefrac = cell(size(exnum));
        sizefrac(contains(t(:,2), 'total')) = {'total'};
        sizefrac(contains(t(:,2), 'less5')) = {'less5'};
        sizefrac(contains(t(:,2), 'less20')) = {'less20'};
        [uw_num(gind), g1,g2] = findgroups(exnum, sizefrac);
    end

[G, C, N, uw] = findgroups(gated_table.Event, gated_table.Niskin, uw_num);
ia = 1:max(G); 

CNTable = table; 
cruisevec = repmat(cruisename, length(C), 1) ;
CNTable.cruise = cruisevec;
CNTable.Event = C; 
CNTable.Niskin = N; 

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
%    CNTable.salinity(g) = temp.salinity(1);
 %   CNTable.potemp090c(g) = temp.potemp090c(1);
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
            time = dir(fullfile(fpath, temp.fcslist{ind(f)}));
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
            time = dir(fullfile(fpath, temp.fcslist{ind(f)}));
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
            time = dir(fullfile(fpath, temp.fcslist{ind(f)}));
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
            time = dir(fullfile(fpath, temp.fcslist{ind(f)}));
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


save([outpath 'SummaryTable.mat'], 'CNTable')

clearvars -except basepath restpath fpath outpath classpath awspath cruisename hierarchical_gates Step1 Step2 Step3 Step4 Step5 Step6 Step7

end


%% Step 7 - Now we want to reformat the Summary table to be like EDI table with carbon counts etc 

if Step7

%Load data for this cruise

G = load([outpath 'Gated_Table.mat']);
gated_table = G.gated_table; 


C = load([outpath 'SummaryTable.mat']); 
CNTable = C.CNTable; 

if strcmp(cruisename, 'TN368')
%if cruisename == 'TN368'
    CNTable(CNTable.Event == 0, :) = []; %two bad entries, not matched to Events
end

EDI_table = table(CNTable.cruise, CNTable.Event, CNTable.Niskin, CNTable.latitude, CNTable.longitude, CNTable.depth_m); %, CNTable.salinity, CNTable.potemp090c, CNTable.r2r_event); 
EDI_table.Properties.VariableNames = {'cruise'; 'Event'; 'niskin'; 'latitude'; 'longitude';  'depth_m';}; %'salinity'; 'potential_temperature_c'; 'r2r_event'}; 

EDI_table.cruise = string(EDI_table.cruise); %helpful for merging tables when cruisenames are different lengths 

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

    if ~exist([[classpath filesep cfilename]])
        EDI_table.syn_cells_per_ml(i) = NaN;
        EDI_table.syn_biovolume_concentration(i) = NaN;
        EDI_table.syn_carbon_concentration(i) = NaN; %divide by 1000 to get micrograms per Liter
        EDI_table.syn_volume_analyzed_ml(i) = NaN; 
        EDI_table.syn_filename{i} = 'NaN';
        continue
    end

    C = load([classpath filesep cfilename]); 

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

    if ~exist([[classpath filesep cfilename]])
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

    C = load([classpath filesep cfilename]); 

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

    if ~exist([[classpath filesep cfilename]])
        EDI_table.hetprok_cells_per_ml(i) = NaN;
        EDI_table.hetprok_carbon_concentration(i)
        EDI_table.hetprok_volume_analyzed_ml(i) = NaN;
        EDI_table.hetprok_filename{i} = 'NaN';
        continue
    end

    C = load([classpath filesep cfilename]); 

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
save([outpath 'EDI_table.mat'], 'EDI_table')
disp([outpath 'EDI_table.mat'])


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
