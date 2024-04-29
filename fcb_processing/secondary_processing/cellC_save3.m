function [ ] = cellC_save3( year2do )
%cellC_save2 updated from cellC_save to include cyrptos and to be more
%efficient, 3/19/2010, Heidi
%cellC_save3, update from cellC_save to be generic for all years

base_path1 = '\\sosiknas1\Lab_data\MVCO\FCB\';
temp = dir([base_path1 'MVCO_*' num2str(year2do)]);
base_path = [base_path1 temp.name filesep];

mergedpath = [base_path 'data\processed\grouped\merged\'];
groupedpath = [base_path 'data\processed\grouped\'];
mergedpath2 = mergedpath;
groupedpath2 = groupedpath;

g = load([groupedpath 'groupsum']); %get the summary file with smoothed bead values

if year2do <=2005
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

cellresults_all = [];
sumvol_all = [];
sumC_all = [];
sumnum_all = [];
for count = 1:length(filelist)
    filename = filelist(count).name;
    disp(['loading...' filename])
    groupedpath = groupedpath2;
    eval(['load ' groupedpath filename])
    mergedpath = mergedpath2;
    if year2do <= 2005
        eval(['load ' mergedpath filename(1:end-6) 'merged' filename(end-5:end-4)])
    else
        eval(['load ' mergedpath filename(1:11) 'merged' filename(12:end-4)])
    end
    cellresults_all = [cellresults_all; cellresults];
    sumvol = NaN(length(mergedwithclass),8);
    sumC = sumvol;
    sumnum = sumvol;
    beadmatchSSCsmooth = NaN(size(beadmatch,1),1);
    nnan = ~isnan(beadmatch(:,1));
    [~,ia] = ismember(beadmatch(nnan,1), g.beadmatchall(:,1));
    beadmatchSSCsmooth(nnan) = g.beadmatchSSCsmooth(ia); %smoothed SSC bead value for this file
    for count2 = 1:length(mergedwithclass)
        temp = mergedwithclass{count2};    
        for class = 1:4
            if size(temp,1) > 1
                if class == 1
                    ind = find(temp(:,classcol) == 1);  %syn
                elseif class == 2
                    ind = find(temp(:,classcol) == 4 | temp(:,classcol) == 2 | temp(:,classcol) == 6);  %euk + sm&lg crypto
                elseif class == 3
                    ind = find(temp(:,classcol) == 6);  %dim crypto
                elseif class == 4
                    ind = find(temp(:,classcol) == 2);  %bright crypto
                end
                %cellvol = cytosub_SSC2vol(temp(ind,5)./beadmatch(count2,5));  %SSC bu values converted to volume
                %use smoothed SSC bead value as saved in groupsum file after running bead_smooth.m
                cellvol = cytosub_SSC2vol(temp(ind,5)./beadmatchSSCsmooth(count2));  %SSC bu values converted to volume
                cellC = 10.^(-.665+.939.*log10(cellvol)); %Menden-Deuer & Lessard 2000, protist excl. diatoms
                celldiam = (cellvol.*3./4./pi).^(1/3)*2;
                sumvol(count2,class) = sum(cellvol);
                sumC(count2,class) = sum(cellC);
                sumnum(count2,class) = length(ind);
                if class == 2
                    ind2 = find(celldiam <=2);
                    sumvol(count2,5) = sum(cellvol(ind2));
                    sumC(count2,5) = sum(cellC(ind2));
                    sumnum(count2,5) = length(ind2);
                    ind2 = find(celldiam <=3);
                    sumvol(count2,6) = sum(cellvol(ind2));
                    sumC(count2,6) = sum(cellC(ind2));
                    sumnum(count2,6) = length(ind2);
                    ind2 = find(celldiam <=5);
                    sumvol(count2,7) = sum(cellvol(ind2));
                    sumC(count2,7) = sum(cellC(ind2));
                    sumnum(count2,7) = length(ind2);
                    ind2 = find(celldiam <= 10);
                    %ind2 = find(celldiam > 2 & celldiam <= 10);
                    sumvol(count2,8) = sum(cellvol(ind2));
                    sumC(count2,8) = sum(cellC(ind2));
                    sumnum(count2,8) = length(ind2);
                end
                %else  %empty
                %    sumvol(count2,class) = NaN;
                %    sumC(count2,class) = NaN;
                %    sumnum(count2,class) = NaN;
            end
        end
        clear cellvol cellC celldiam
    end %for class = 1:2
    sumvol_all = [sumvol_all; sumvol(:,[1:2,5:8,3:4])];
    sumC_all = [sumC_all; sumC(:,[1:2,5:8,3:4])];
    sumnum_all = [sumnum_all; sumnum(:,[1:2,5:8,3:4])];
    clear sumvol sumC sumnum
    clear mergedwithclass beadmatch cellresults
end  %for count = 1:length(filelist)

%sumtitles = {'Syn', 'Euk total', 'Euk < 2 microns', 'Euk 2-10 microns' 'Dim cryptos', 'Bright cryptos'};
sumtitles = {'Syn', 'Euk_total', 'Eukleq2microns', 'Eukleq3microns', 'Eukleq5microns', 'Eukleq10microns', 'Dim_cryptos', 'Bright_cryptos'};
clear beadmatch* cellCHL* cellSSC* cellFLS* cellNUM* cellPE* merged* class* count* grouped* ind* file* temp
save([base_path1 '\summary\cellC_summary_' num2str(year2do)])
