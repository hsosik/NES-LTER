%mergebatch2 from mergebatch, modified for new file format, 3/7/03 Heidi
%created from leobatch3, Heidi 1/10/02
%low and high gain merging 

global special_flag

close all
if filetype == 'cell' 
    allslopes = [];
    allintercepts = [];
    allcount = [];
    alldate = [];
    allr = [];
elseif filetype == 'bead'
    eval(['load ' procpath 'mergefitsum'])%load fit results from cell analysis to use with bead files
end;

if year2do <= 2004,
    dattitles = {'event' 'PE' 'FLS' 'CHL' 'SSC' 'CHLpeak'};
    fittitles = {'Slope (PE,FLS,CHL,SSC)'; 'Intercept (PE,FLS,CHL,SSC)'; 'Count (PE,FLS,CHL,SSC1,SSC2)'};
else
    dattitles = {'event' 'PE' 'FLS' 'CHL' 'SSC1_2' 'CHLpeak1' 'CHLpeak2'};
    fittitles = {'Slope (PE,FLS,CHL,SSC1,SSC2,SSC1_2)'; 'Intercept (PE,FLS,CHL,SSC1,SSC2)' ; 'Count (PE,FLS,CHL,SSC1,SSC2)'; 'r (PE,FLS,CHL,SSC1,SSC2)'};
end;
special_flag = 0;
for typenum = 1:size(filetypelist,1),
    if ismember(filetypelist(typenum,:), {'FCB2_2009_1' 'FCB2_2009_2'})
        special_flag = 1; %handle mergeproc case with bad SSC2 data
    end;
    filelist = dir([datapath filetypelist(typenum,:) '*']);
    if year2do <= 2006,
       date = datenum(cat(1,filelist.date));
       [temp, fileorder] = sort(date);
       filelistmain = filelist(fileorder); %consequecutive order
       clear temp date
    else
        filelistmain = filelist;
    end;
    if year2do == 2008 %special case dealing with day of mixed local and UTC time stamps (22 Oct 2008)
        ii1 = strmatch('FCB1_2008_296_092206', {filelistmain.name});
        ii2 = strmatch('FCB1_2008_296_130826',{filelistmain.name});
        ii3 = strmatch('FCB1_2008_296_114008',{filelistmain.name});
        filelistmain(sort([ii1,ii2,ii3])) = filelistmain([ii1,ii2,ii3])
    end;
    filesections = ceil(length(filelistmain)/setsize);
    for sectcount = 1:filesections,
        if sectcount < filesections,
            filelist = filelistmain((sectcount-1)*setsize+1:sectcount*setsize);
        else
            filelist = filelistmain((sectcount-1)*setsize+1:end);
        end;
        nextevent = 0;
        for filenum = 1:length(filelist),       
            filename = filelist(filenum).name;  
            filedate = filelist(filenum).date;
            disp(filename)
            eval(['[header, dat] = ' readrawstr '(datapath, filename);'])
            %for ii = 1:size(dat,1)/100, z(ii) = length(find(dat(ii*100-99:ii*100,6)<=10)); end
            %figure(1), clf, plot(z/100), ylim([0 1])
            %pause(.1)
            if filetype == 'bead',
                [junk,ind] = min(abs(datenum(filedate)-alldate)); %find closest date
                fittouse = [allslopes(:,ind) allintercepts(:,ind)];
                clear junk ind
            else
                fittouse = NaN;
            end;
            
            eval(['[datmerged, fit] = ' mergeprocstr '(dat, fittouse);'])
            %            clear fittouse
            allslopes = [allslopes fit(:,1)];
            allintercepts = [allintercepts fit(:,2)];
            allcount = [allcount fit(:,3)];
            allr = [allr fit(:,4)];
            alldate = [alldate datenum(filedate)];
            lastevent = nextevent+size(datmerged,1);
            datmerged = [(nextevent+1:lastevent)' datmerged];
            nextevent = lastevent;
            datmerged = uint32(datmerged);  %saves disk and memory
            if year2do <=2004,
                datmerged(:,6) = dat(:,9);
            else
                %these are chl peak, unmerged, currentlyoverwrites SSC1 and
                %SSC2 merge products that come from fcbmergeproc2 , 4-13-05
                datmerged(:,6:7) = dat(:,11:12);  
            end;
            clear dat
            if year2do <= 2005
                dotpos = findstr('.', filename);
                eval(['save ' procpath filetypelist(typenum,:) filename(dotpos+1:end) ' dattitles fittitles datmerged fit'])
            else %no extension on files any more, 6/6/06 heidi
                eval(['save ' procpath filename ' dattitles fittitles datmerged fit'])
            end;
            clear datmerged dotpos fit
        end; %for filenum = 1:length(filelist)
    end; %for sectcount
end; %for typenum = 1:length(filetypelist)
if size(allslopes,2) == 4, 
    allslopes(:,5:6) = NaN;
    allintercepts(:,5:6) = NaN;
    allr(:,5:6) = NaN;
end;
if filetype == 'cell'
    eval(['save ' procpath 'mergefitsum all* fittitles'])%save fit results to use with bead files
    if plotflag,
        figure(3)
        subplot(311)
        plot(allslopes(1:6,:)', '.-')
        legend('PE', 'FLS', 'CHL', 'SSC', 'SSC2', 'SSC1/2')
        title('slopes')
        subplot(312)
        plot(allintercepts(1:6,:)', '.-')
        legend('PE', 'FLS', 'CHL', 'SSC', 'SSC2', 'SSC1/2')
        title('intercepts')
        subplot(313)
        plot(allcount(1:6,:)', '.-')
        legend('PE', 'FLS', 'CHL', 'SSC', 'SSC2', 'SSC1/2')
        title('counts')
   		disp('pause for graph...')
        pause
    end;
end;

clear typenum alldate allslopes allintercepts allcount filename filelist fileorder filedate filenum
clear fittitles dattitles lastevent nextevent header