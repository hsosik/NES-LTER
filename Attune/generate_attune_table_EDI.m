%Generates the attune table from Class files and FCSfileinfo.mat only
%Created 5/3/21 as parth of steps to make attune processing more modular. 

%Saves a table of stats about the fcs files as well as a variable 
%table_metadata which saves date generated and assign class function used

%classpath should end in \

function [] = generate_attune_table_EDI(classpath, FCSfileinfopath)

Attune = load([FCSfileinfopath]);

PE_euk_celltypes = input('which class numbers are PE containing Euks for this cruise?');

% identify bead run files and determine mean bead SSC-H for cruise
% may need to adjust epsilon and/or minpts depending on cruise

% read in and process Attune data files
%classlist = dir([classpath, '*.mat']); 
filelist = regexprep(Attune.FCSfileinfo.fcslist,'fcs', 'mat'); 
AttuneTable = table(Attune.FCSfileinfo.fcslist, datetime(Attune.FCSfileinfo.matdate_start, 'ConvertFrom', 'datenum'), datetime(Attune.FCSfileinfo.matdate_stop, 'ConvertFrom', 'datenum'), Attune.FCSfileinfo.vol_analyzed/1e6, 'VariableNames', {'Filename' 'StartDate' 'StopDate' 'VolAnalyzed_ml'});

% Creating the variables
%we want Syn, Euk<=2, Euk<=3, Euk<=5, Euk<=10, Euk<=20, PEeuk<=2, PEeuk<=3,
%PEeuk<5, PEeuk<10, PEeuk<20


EukSizes = [0 2 3 5 10 20];
numBins = 1+2*(length(EukSizes)-1);
Count = NaN(length(filelist),numBins);
Biovol = Count;
Carbon = Count;
rem_ind = []; 
Scatter_hv = NaN(length(filelist), 1); 

for count = 1:length(filelist) %go through each of the files in the FCSfileinfo
    if ~rem(count,10)
        disp([num2str(count) ' of ' num2str(length(filelist))])
    end
    filename = [classpath filelist{count}]; %load class file 
    if exist(filename)
    
    load(filename)
    if ~exist('volume', 'var') %% if we haven't done size calibration yet, function should still run 
        eval('volume = NaN.*class;'); 
    end
    
    volume = real(volume); 
    carbon = biovol2carbon(volume, 0); % carbon, picograms per cell
    carbon = real(carbon); %having issues with formatting, keeps having valus with + 0i. 

    if exist('file_hv', 'var')
        Scatter_hv(count) = file_hv;
    end
    
    eval(['rename_class = class;']) %issues with class as a variable name since its a function in matlab 
    clear class 
    
    
    %SynFirst
        Count(count,1) = sum(rename_class==2);
        Biovol(count,1) = nansum(volume(rename_class==2));
        Carbon(count,1) = nansum(carbon(rename_class==2));

    %now Euks
    celltype = 1; 
    diam = (volume*3/4/pi).^(1/3)*2; %equivalent spherical diam, micrometers
    for ii = 1:length(EukSizes)-1
        ind = find(diam<=EukSizes(ii+1) & rename_class==celltype);
        Count(count,ii+1) = size(ind,1);
        Biovol(count,ii+1) = nansum(volume(ind));
        Carbon(count,ii+1) = nansum(carbon(ind));
    end

    %now PE containing Euks
    celltype = PE_euk_celltypes;
    for ii = 3:length(EukSizes)-1
        ind = find(diam<=EukSizes(ii+1) & ismember(rename_class, celltype));
        Count(count,ii+length(EukSizes)) = size(ind,1);
        Biovol(count,ii+length(EukSizes)) = nansum(volume(ind));
        Carbon(count,ii+length(EukSizes)) = nansum(carbon(ind));
    end

    
    else
        rem_ind = [rem_ind count]; %indeces to remove from final table
        disp(['skipped ', filelist{count}])
    end
    
    clear volume 
end



%define column names
classnames = {'SynX'};
for ii = 2:length(EukSizes)
    classnames = [classnames; ['Euk_without_PE_leq' num2str(EukSizes(ii)) 'umX']]; 
end
for ii = 4:length(EukSizes)
    classnames = [classnames; ['Euk_w_PE_leq' num2str(EukSizes(ii)) 'umX']] ;
end

AttuneTable = [AttuneTable array2table(Count, 'VariableNames', regexprep(classnames, 'X', '_count'))];
AttuneTable = [AttuneTable array2table(Biovol, 'VariableNames', regexprep(classnames, 'X', '_biovolume'))];
AttuneTable = [AttuneTable array2table(Carbon, 'VariableNames', regexprep(classnames, 'X', '_carbon'))];


AttuneTable.QC_flag = Attune.FCSfileinfo.QC_flag; 
AttuneTable.Scatter_hv = Scatter_hv; 

AttuneTable(rem_ind,:) = []; 
AttuneTable = sortrows(AttuneTable, 'StartDate');

%notes(2) is where assign_class_function should be saved, sometimes saved
%as its own variable in class files
if exist('assign_class_function') 
    table_metadata = {assign_class_function; classpath; string(datetime())};
else 
    table_metadata = {notes(2); classpath; string(datetime())};
end

save([classpath '..\AttuneTable_EDI'],'AttuneTable', 'table_metadata')
disp(['Result file saved:'])
disp([classpath '..\AttuneTable_EDI'])

end

