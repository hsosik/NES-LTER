
%Script to process CTD casts into raw .cnv files with Sea-Bird SBE Data
%Processing software:

CTDpath=fullfile('\\maddie\TaylorF\from_Samwise\data\MVCO\'); % path to folders with raw CTD data...
%mastersourcepath=fullfile('/Volumes/TaylorF/from_Samwise/data/MVCO/'); % path to folders with raw CTD data...

%%
%It appears that the CTD casts are stored in a few files:
%Those with the Tioga designation
%Survery Cruises
%todownload

d = dir(CTDpath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders

%%%% Taylor comment Nov 30, 2017 - Tioga cruise #s are now over 1,000 thus
%%%% 4 digit cruise number which will cause the list to be out of sequential order in time
temp=regexp(foldernames,'Tioga_\d{3,4}'); %find all folders with this naming scheme
%%%%

datafolders = foldernames(cellfun('isempty',temp)==0);

%%% Taylor Nov 30, 2017 directory "SurveyCruises" is now obsolete after
%%% directory re-structuring
%Survey Cruises:
%{
d = dir(fullfile([CTDpath '/SurveyCruises/']));
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'Tioga_\d{3}'); %find
ii=find(cellfun('isempty',temp)==0);
datafolders =[datafolders;  strcat(cellstr(repmat('/SurveyCruises/',length(ii),1)),foldernames(ii))];
%}

%obsolete: no more "todownload" subfolder
%{
%todownload folder:
d = dir(fullfile([CTDpath '/todownload/']));
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'(Tioga_\d{3,4})|(ti\d{3,4})'); 
ii=find(cellfun('isempty',temp)==0);
datafolders =[datafolders;  strcat(cellstr(repmat('/todownload/',length(ii),1)),foldernames(ii))];
%}

datafiles={};
%% Gather a complete list of data files by checking each of these folder locations for raw data files!

for j=1:length(datafolders)
    
    %check to see what files/folders are in each one of the data folders:
    %a bit more complicated because Matlab does not appear to have a
    %recursive search function....arghhh...
    
    filepath=fullfile(CTDpath,datafolders{j});
    filelist=dir(filepath);
    disp(['Searching folder: ' filepath])
    
    % keep on searching until have found a .hex or a .dat file...
    flag=1; %datafiles={};
    
    while flag
        
        rawind=find(cellfun('isempty',regexp({filelist(:).name}','(\w*\.hex)|(\w*_\w*\.dat)'))==0);
        
        if isempty(rawind)
            %disp('No data files here...looking in subdirectories...')
            
            subfolders=regexpi({filelist([filelist(:).isdir]).name}','(\w*ctd)|(original_files)','match');
            subind=find(cellfun('isempty',subfolders)==0);
            if length(subind)==1
                %disp(['found a ' char(subfolders{subind}) ' folder'])
                filepath=fullfile(filepath,subfolders{subind});
                filelist=dir(char(filepath));
            elseif length(subind)>1
                disp('Uh-Oh! Two subfolders?')
                keyboard
            else %Taylor added 6Dec17 - some cruise dir have underway data but no rosette casts were done
                disp(['Skipping cruise, no CTD data: ' char(filepath)])
                flag = 0;
            end    
        else
            %disp(['I found ' num2str(length(rawind)) ' file(s)...adding to list'])
            datafiles=[datafiles; fullfile(filepath,{filelist(rawind).name}')];
            flag=0;
        end
    end    
end

%%
%%%%%HAD TO MANUALLY PROCESS 4 CASTS FROM DISCOVERY AND R+R
% RR_21May10, RR_09Jun10, Discovery_27Jul10, Discovery_11Dec10 - used
% handcast_config_file:SBE19plus_4393_factoryPARfluor.con


%This portion of the script calls SBE Data-processing software. If the
%software cannot find the variables that you requested, you will need to
%manually click continue (doesn't take too long to process the ~200 or so
%files)
handcasts = {'Tioga_475_27Aug10','Tioga_487_11Oct10','Tioga_490_24Oct10'};

for q=1:length(datafiles)
    
    disp(['Processing raw data file for: ' datafiles{q}])
    
    infile=datafiles{q};
    
    outputdir=fullfile('\\sosiknas1\Lab_data\MVCO\','processed_CTD_casts_20171206');
    outtemp=strsplit(datafiles{q},'\');
    outfile=char(outtemp{end}(1:end-4));
    
%get associated config file. Some cruises use handcasts with departhment seabird, only 1 config file for those in different dir
if ~isempty(find(strcmp(outtemp{7},handcasts),1))
    confile = '\\MADDIE\TaylorF\from_Samwise\data\MVCO\handcast_config_file\SBE19plus_4393_factoryPARfluor.con';
else
    confile=regexprep(datafiles{q},'(\.hex)|(\.dat)','\.CON'); %use the corresponding .con file - this should be with the data file!
end
    psafile='C:\Users\heidi\AppData\Local\Sea-Bird\SBEDataProcessing-Win32\DatCnv.psa';
    
    %/s = start processing now
    %/m = minimize window
    
    cmdstr = ['"C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32\datcnvw.exe" ' '/s /c' confile ' /i' infile ' /o' outputdir ' /f' outfile ' /p' psafile];
    [s,w] = system(cmdstr);
    
%     if ~regexp(w,'XML document loaded successfully')
%         keyboard
%     end
    
end

%% and save the list of datafiles for futre info processing:

save(fullfile(outputdir,'list_and_location_of_raw_ctd_files'),'datafiles')
