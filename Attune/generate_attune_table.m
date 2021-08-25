%Generates the attune table from Class files and FCSfileinfo.mat only
%Created 5/3/21 as parth of steps to make attune processing more modular. 

%Saves a table of stats about the fcs files as well as a variable 
%table_metadata which saves date generated and assign class function used

%classpath should end in \

function [] = generate_attune_table(classpath, FCSfileinfopath)

Attune = load([FCSfileinfopath]);

% identify bead run files and determine mean bead SSC-H for cruise
% may need to adjust epsilon and/or minpts depending on cruise

% read in and process Attune data files
%classlist = dir([classpath, '*.mat']); 
filelist = regexprep(Attune.FCSfileinfo.fcslist,'fcs', 'mat'); 
AttuneTable = table(Attune.FCSfileinfo.fcslist, datetime(Attune.FCSfileinfo.matdate_start, 'ConvertFrom', 'datenum'), datetime(Attune.FCSfileinfo.matdate_stop, 'ConvertFrom', 'datenum'), Attune.FCSfileinfo.vol_analyzed/1e6, 'VariableNames', {'Filename' 'StartDate' 'StopDate' 'VolAnalyzed_ml'});

% Creating the variables
numClass = 6;
diamEdges = [0 2 5 10 20 50 inf];
numBins = length(diamEdges)-1;
Count = NaN(length(filelist),numClass);
Biovol = Count;
Carbon = Count;
CountBin = NaN(length(filelist),numBins);
BiovolBin = CountBin;
CarbonBin = CountBin;
rem_ind = []; 

for count = 1:length(filelist) %go through each of the files in the FCSfileinfo
    if ~rem(count,10)
        disp([num2str(count) ' of ' num2str(length(filelist))])
    end
    filename = [classpath filelist{count}]; %load class file 
    if exist(filename)
    
    load(filename)
    if ~exist('volume', 'var') %% if we haven't done size calibration yet, function should still run 
        %eval('classvec = class;'); %again issues with class as variable name
        eval('volume = NaN.*class;'); 
    end
    
    carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
    
    eval(['rename_class = class;']) %issues with class as a variable name since its a function in matlab 
    clear class 
    for ii = 1:numClass
        Count(count,ii) = sum(rename_class==ii);
        Biovol(count,ii) = nansum(volume(rename_class==ii));
        Carbon(count,ii) = nansum(carbon(rename_class==ii));
    end
    
    diam = (volume*3/4/pi).^(1/3)*2; %equivalent spherical diam, micrometers
    for ii = 1:length(diamEdges)-1
        ind = find(diam>=diamEdges(ii) & diam<diamEdges(ii+1) & rename_class~=0);
        CountBin(count,ii) = size(ind,1);
        BiovolBin(count,ii) = sum(volume(ind));
        CarbonBin(count,ii) = nansum(carbon(ind));
    end
    else
        rem_ind = [rem_ind count]; %indeces to remove from final table
        disp(['skipped ', filelist{count}])
    end
    
    clear volume 
end

for ii = 1:numBins
    binlabel{ii} = ['X' num2str(diamEdges(ii)) 'to' num2str(diamEdges(ii+1))];
end


%back classnames out of notes string, sometimes a string array... 
if isstring(notes)
    classnames = split(notes(1), '='); 
else
    classnames = split(notes, '='); 
end
    classnames = classnames(2:numClass+1); 
    classnames = regexprep(classnames, 'Class \d', '');
    classnames = regexprep(classnames, '_euk_coincident', 'Euk');

AttuneTable = [AttuneTable array2table(Count, 'VariableNames', regexprep(classnames, ', ', '_count'))];
AttuneTable = [AttuneTable array2table(Biovol, 'VariableNames', regexprep(classnames, ', ', '_biovolume'))];
AttuneTable = [AttuneTable array2table(Carbon, 'VariableNames', regexprep(classnames, ', ', '_carbon'))];
AttuneTable = [AttuneTable array2table(CountBin, 'VariableNames', regexprep(binlabel, 'X', 'count_'))];
AttuneTable = [AttuneTable array2table(BiovolBin, 'VariableNames', regexprep(binlabel, 'X', 'biovolume'))];
AttuneTable = [AttuneTable array2table(CarbonBin, 'VariableNames', regexprep(binlabel, 'X', 'carbon'))];

AttuneTable.QC_flag = Attune.FCSfileinfo.QC_flag; 

AttuneTable(rem_ind,:) = []; 
AttuneTable = sortrows(AttuneTable, 'StartDate');

%notes(2) is where assign_class_function should be saved, sometimes saved
%as its own variable in class files
if exist('assign_class_function') 
    table_metadata = {assign_class_function; classpath; string(datetime())};
else 
    table_metadata = {notes(2); classpath; string(datetime())};
end

save([classpath '..\AttuneTable'],'AttuneTable', 'table_metadata')
disp(['Result file saved:'])
disp([classpath '..\AttuneTable'])

end

