function [ ] = sizedist_save2( year2do )
%sizedist_save2 updated from sizedist_save, now as function and generic for
%all years; April 2017, Heidi

base_path1 = '\\sosiknas1\Lab_data\MVCO\FCB\';
temp = dir([base_path1 'MVCO_*' num2str(year2do)]);
base_path = [base_path1 temp.name filesep];

mergedpath = [base_path 'data\processed\grouped\merged\'];
groupedpath = [base_path 'data\processed\grouped\'];
mergedpath2 = mergedpath;
groupedpath2 = groupedpath;

g = load([groupedpath 'groupsum']);
nind = find(~isnan(g.cellresultsall(:,1)));
whos nind

path_out = [base_path 'data\processed\size\'];
if ~exist(path_out, 'dir')
    mkdir(path_out)
end

if year2do <=2005,
    filelist = dir([groupedpath '*.mat']);
    filelist(strmatch('groupsum.mat',{filelist.name}')) = [];
    filelist(strmatch('cellC_summary.mat',{filelist.name}')) = [];
else
    filelist = dir([groupedpath 'FCB*.mat']);
end
if year2do <=2004
    classcol = 7;
else
    classcol = 8;
end
if year2do <= 2005
    upn = 1;
else 
    upn = 3;
end

%filelist = dir([groupedpath 'FCB*.mat']);
%diambins = [0:1:15];
diambins = [-inf 2.^[-3:.15:4] inf];
%binwidth = diff(diambins); binwidth(1) = diambins(2); 
%diambins = [2.^[-3:.25:6]];
%cellclass = 4;  %cell group to model (1 = syn, 4 = euks)        

for count = 1:length(filelist),
   filename = filelist(count).name;
   disp(['loading...' filename])
   eval(['load ' groupedpath filename])
   [~, ia, ib] = intersect(cellresults(:,1), g.cellresultsall(:,1));
   if length(ib) ~= size(cellresults,1)
       disp('missing times in groupedsum')
       whos ia cellresults
       beadmatchSSCsmooth = interp1(g.cellresultsall(nind,1), g.beadmatchSSCsmooth(nind), cellresults(:,1));
%       keyboard
   else
       beadmatchSSCsmooth = g.beadmatchSSCsmooth(ib);
   end
   upos = findstr('_', filename); 
%            eval(['load ' mergedpath filename(1:end-6) 'merged' filename(end-5:end-4)])
   eval(['load ' mergedpath2 filename(1:upos(upn)-1) 'merged' filename(upos(upn):end-4)])
    for class = 1:3,
        if class == 1,
            cellclass = 4; %euk
        elseif class == 2
            cellclass = 1; %syn
        else
            cellclass = 3; %case for crypto
        end;
        N_dist = NaN(length(diambins),length(mergedwithclass));
        biovol_dist = N_dist;          
        for count2 = 1:length(mergedwithclass),
            temp = mergedwithclass{count2};
            if size(temp,1) > 1,
                if cellclass ~=3
                    ind = find(temp(:,classcol) == cellclass);
                else
                    ind = find(temp(:,classcol) == 2 | temp(:,classcol) == 6);
                end;
                %cellvol = cytosub_SSC2vol(temp(ind,5)./beadmatch(count2,5));  %SSC bu values converted to volume
                cellvol = cytosub_SSC2vol(temp(ind,5)./beadmatchSSCsmooth(count2));  %SSC bu values converted to volume
                celldiam = (cellvol.*3./4./pi).^(1/3)*2; %clear cellvol
                N_dist(:,count2) = histc(celldiam, diambins);  %size distribution
                for bincount = 1:length(diambins)-1,
                    ind = find(celldiam >= diambins(bincount) & celldiam < diambins(bincount+1));
                    if length(ind) ~= N_dist(bincount,count2), keyboard, end;
                    if ~isempty(ind),
                        biovol_dist(bincount,count2) = sum(cellvol(ind));
                    else
                        biovol_dist(bincount,count2) = 0;
                    end;
                end;
                %N_dist2(:,count2) = N_dist(1:end-1,count2)./cellresults(count2,3)./binwidth';  %cells/ml/micron    
            end;
        end;
        if class == 1,
            N_dist_euk = N_dist;
            biovol_dist_euk = biovol_dist;
            %N_dist2_euk = N_dist2;
        elseif class == 2,
            N_dist_syn = N_dist;
            biovol_dist_syn = biovol_dist;
            %N_dist2_syn = N_dist2;    
        else
            N_dist_crp = N_dist;
            biovol_dist_crp = biovol_dist;
            %N_dist2_crp = N_dist2;   
        end;
        clear N_dist biovol_dist cellvol celldiam %N_dist2
    end; %for class = 1:2    
        [~,filename] = fileparts(filename);
        %eval(['save ' filename ' N* biovol* diambins cellresults'])  
        save([path_out filename],'N*', 'biovol*', 'diambins', 'cellresults')
        clear N* cellresults
end;  %for count = 1:length(filelist)

