%Analyze preserved discrete samples from cruise CTD Casts

% We will generate 2 tables:

% A) one with cast, niskin, settings, and counts for
% each file, automatically generated from names

% B) a table with files to look at for counts of each class for each Cast
% and Niskin. This one can be edited if there are duplicates etc. that need
% to be manually dealt with. 
% Both A & B will only ever add new files or casts to the tables so as not
% to overwrite any manual work. 

% We'll make another script to classify and generate a final table with 
% Cast Niskin metadata and counts reported from the
% files indicated by B. 


%% Manually choose cruise to process
basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20210716_EN668\preserved'; 
restpath = 'https://nes-lter-data.whoi.edu/api/ctd/en668/';

%Choose which steps to do 
maketable1 = 0
maketable2 = 1

% Vars for parsing filenames
filetype2exclude = {'Rinses', 'skip', 'xxx';}; 
beadfiles2include = {'FCB_bead'};


%% some file structure setup  
fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'outputs' filesep];
classpath = [outpath 'class' filesep];

if ~exist(outpath, 'dir')
    mkdir(outpath)
end
if ~exist(classpath, 'dir')
    mkdir(classpath)
end

%% Table 1 

if maketable1 

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


headers = {'fcslist', 'Cast', 'Niskin'};

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
     %to extract cast and niskin numbers
     
     for ii = 1:length(samplelist) 
            filename = samplelist{ii};
            Table_Add.fcslist(ii) = {filename}; 
            s = strfind(filename, '_C');
            s = s(end); %get index of last _ in filename
            if ~isempty(str2num(filename(s+6:s+7)))
                Table_Add.Cast(ii) = str2num(filename(s+2:s+4));
            end
            if ~isempty(str2num(filename(s+6:s+7)))
                Table_Add.Niskin(ii) = str2num(filename(s+6:s+7));
            end
     end
     
    FCSList = [FCSList; Table_Add];

 save([outpath 'FCSList.mat'], 'FCSList')

end

%% Make Table 2 

if maketable2

if ~exist('FCSList')
    P = load([outpath '/FCSList.mat']);
    FCSList = P.FCSList; 
end

if ~exist([outpath '/FilesToUse.csv'])
    FilesToUse = table(); 
    [G, C, N] = findgroups(FCSList.Cast, FCSList.Niskin);
    ia = 1:length(G); 
else
    FilesToUse = readtable([outpath 'FilesToUse.csv']);
    if ~iscell(FilesToUse.ProFile) %have to make it a cell, if it was read as NaNs
        FilesToUse.ProFile = cell(height(FilesToUse), 1);
    end
    FilesToUse.Flags = num2cell(FilesToUse.Flags);
    [G, C, N] = findgroups(FCSList.Cast, FCSList.Niskin);
    [G2, C2, N2] = findgroups(FilesToUse.Cast, FilesToUse.Niskin);
    [~, ia] = setdiff([C N], [C2 N2], 'rows');
end

FilesToADD = table(C(ia), N(ia), 'VariableNames', {'Cast'; 'Niskin'});

for g = 1:length(ia)% go through unique casts, which have not already been included
    
    flag = ''; 

    if C(ia(g)) == 0 
        continue %parsing error, no cast number
    end

    Subset = FCSList(G==ia(g), :); %limit list to just this Niskin

    %Get best guess filename for each cell type 
    %Syn First
    ind = find(contains(Subset.fcslist, 'phyto_PE'));
    if length(ind) == 1 %if only one fcs file of this type, use that. 
        FilesToADD.SynFile(g) =  Subset.fcslist(ind);
    elseif ~isempty(ind) %if more than 1, get most recent
        flag = strcat(flag, '1');
        timesince = []; 
        for f = 1:length(ind)
            time = dir(fullfile(fpath, Subset.fcslist{ind(f)}));
            time = datetime(time.date); 
            timesince = [timesince datetime()-time];
        end
        [~,truind] = min(timesince);
        FilesToADD.SynFile(g) =  Subset.fcslist(ind(truind));
    end

     %Picoeuks
    ind = find(contains(Subset.fcslist, 'phyto_CHL'));
    if length(ind) == 1 %if only one fcs file of this type, use that. 
        FilesToADD.EukFile(g) =  Subset.fcslist(ind);
    elseif ~isempty(ind) %if more than 1, get most recent
        flag = strcat(flag, '2');
        timesince = []; 
        for f = 1:length(ind)
            time = dir(fullfile(fpath, Subset.fcslist{ind(f)}));
            time = datetime(time.date); 
            timesince = [timesince datetime()-time];
        end
        [~,truind] = min(timesince);
        FilesToADD.EukFile(g) =  Subset.fcslist(ind(truind));
    end

     %Heterotrophic bacteria 
    ind = find(contains(Subset.fcslist, 'hbac'));
    if length(ind) == 1 %if only one fcs file of this type, use that. 
        FilesToADD.BacteriaFile(g) =  Subset.fcslist(ind);
    elseif ~isempty(ind) %if more than 1, get most recent
        flag = strcat(flag, '3');
        timesince = []; 
        for f = 1:length(ind)
            time = dir(fullfile(fpath, Subset.fcslist{ind(f)}));
            time = datetime(time.date); 
            timesince = [timesince datetime()-time];
        end
        [~,truind] = min(timesince);
        FilesToADD.BacteriaFile(g) =  Subset.fcslist(ind(truind));
    end

     %Prochlorococus if present
    ind = find(contains(Subset.fcslist, 'pro', 'IgnoreCase', true));
    if isempty(ind)
        FilesToADD.ProFile{g} = '';
    elseif length(ind) == 1 %if only one fcs file of this type, use that. 
        FilesToADD.ProFile(g) =  Subset.fcslist(ind);
    else %if more than 1, get most recent
        flag = strcat(flag, '4');
        timesince = []; 
        for f = 1:length(ind)
            time = dir(fullfile(fpath, Subset.fcslist{ind(f)}));
            time = datetime(time.date); 
            timesince = [timesince datetime()-time];
        end
        [~,truind] = min(timesince);
        FilesToADD.ProFile(g) =  Subset.fcslist(ind(truind));
    end

    FilesToADD.Flags{g} = flag; 

end


%% add metadata

bottledata = webread([restpath 'bottles.csv']); 
metadata = webread([restpath 'metadata.csv']); 

castlist = unique(C);
FilesToADD_w_meta = table(); 

for c = 1:length(castlist)

    if castlist(c) == 0 
        continue
    end

    tempP = FilesToADD(FilesToADD.Cast == castlist(c), :); 
    %first add depths for each niskin number
    tempB = bottledata(bottledata.cast == castlist(c), :); 
    Niskin = tempB.niskin; 
    Depth = tempB.depsm; 
    Salinity = tempB.sal00;
    poTemp = tempB.potemp090c;
    par = tempB.par; 

    tempB = table(Niskin, Depth, Salinity, poTemp, par); 
    tempP = join(tempP, tempB); 
    
    %now get overall cast info 
    tempM = metadata(metadata.cast == castlist(c), :); 
    tempM = repelem(tempM, height(tempP), 1); 
    
    tempP = [tempP tempM]; 
    
    FilesToADD_w_meta = [FilesToADD_w_meta; tempP]; 
end
    

FilesToUse = [FilesToUse; FilesToADD_w_meta];

writetable(FilesToUse, [outpath 'FilesToUse.csv'])

end

end
