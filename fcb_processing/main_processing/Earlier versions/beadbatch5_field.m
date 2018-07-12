%beadbatch5_field modified to correct bug in loop to do a double section if
%bead sequence split between them; Feb 2014 (previous version wasn't really
%working before now)
%beadbatch5_field modified to deal with low bead numbers encountered part
%way through 2006 and badly in early 2007; changed method for finding bead
%window to use histmulti5 and low resolution bins (64 channels instead of
%1024); Also edited to correct bad beadtitles last entry; April 2007 Heidi
%beadbatch4_field modified to deal with problem after addition of 1 micron
%beads, need to deal with cases where these become the mode...force mode
%picking to get the large beads, 12/9/03 Heidi
%beadbatch3_field modified to deal with 1 hr sections and
%subsections of partialdatmerged (instead of accumulating all files for
%each typenum (solved memory probelm at ~50 files on beans), 5/16/03 Heidi
%beadbatch3 modified from beadbatch_lab for new file format (current
%version for separate bead filpes as in lab), 3/03 Heidi
%beadbatch for oc382, modified from beadreproc for leo-15 2001
%beadreproc: based on heiditemp5, change to more channels and histc to match cellreproc (at crypto stage)
%       also add smoothing of hist before mode est, and mean and std dev for singlet beads only, 2/18/02
%process for beads

timeinterval = 1/24;  %sec (1/24 = 1 hr), resolution for final values

for typenum = 1:size(filetypelist,1),
    filelist = dir([datapath filetypelist(typenum,:) '*.mat']);
  %  date = datenum(cat(1,filelist.date));
  %  [temp, fileorder] = sort(date);
    filelistmain = filelist;
    if year2do <= 2005
        [~,fileorder] = sort(str2num(char(regexprep(regexprep({filelist.name}, '.mat', ''), filetypelist(typenum,:), ''))));
        filelistmain = filelist(fileorder); %consecutive order
    end;
    if year2do == 2008 %special case dealing with day of mixed local and UTC time stamps (22 Oct 2008)
        ii1 = strmatch('FCB1_2008_296_092206', {filelistmain.name});
        ii2 = strmatch('FCB1_2008_296_130826',{filelistmain.name});
        ii3 = strmatch('FCB1_2008_296_114008',{filelistmain.name});
        filelistmain(sort([ii1,ii2,ii3])) = filelistmain([ii1,ii2,ii3])
    end;
    clear temp date fileorder

    filesections = ceil(length(filelistmain)/setsize);
    for sectcount = 1:filesections,
        if sectcount < filesections,
            filelist = filelistmain((sectcount-1)*setsize+1:sectcount*setsize);
        else
            filelist = filelistmain((sectcount-1)*setsize+1:end);
        end;
        eval(['load ' timepath filetypelist(typenum,:) 'time_' num2str(sectcount)])  %load file with processed time
        eval(['totaltime = ' filetypelist(typenum,:) 'time; clear ' filetypelist(typenum,:) 'time'])

        %    timesectionendbin = [1 size(totaltime,1)];
        timesectionendbin = find(diff(floor(totaltime(:,2)/timeinterval)));  %location of hour changes, last point in hr
        timesectionendbin =  [timesectionendbin; size(totaltime,1)];
        timestartind = 1;
        filenum = 1;
        filename = filelist(filenum).name;
        disp(filename)
        eval(['load ' datapath filename])
        partialdatmerged = double(datmerged);
        filenum = 2;
        beadsection = 1;
        %    for sectionnum = 1:length(timesectionendbin) - 1,
        sectionnum = 1;
        while sectionnum <= length(timesectionendbin),
            %        disp(['sectionnum: ' num2str(sectionnum)])
            timeendind = timesectionendbin(sectionnum); %new 5/16/03, heidi
            if sectionnum > 1, timestartind = timesectionendbin(sectionnum-1) + 1; end;
            if totaltime(timeendind,6) == beadport & sectionnum ~= length(timesectionendbin), %do a double "section" if one ends in the middle of bead analysis
                sectionnum = sectionnum + 1;
                timeendind = timesectionendbin(sectionnum);
            end;
            datendind = totaltime(timeendind,1); %this line moved after do a double loop, Feb 2014 - otherwise the loop doesn't actually do anything with more data
            while (partialdatmerged(end,1) < datendind) & filenum <= length(filelist),  %keep adding on files until get full hour
                %        while totaltime(timesectionendbin(sectionnum),1) > partialdatmerged(end,1),  %get to first file
                filename = filelist(filenum).name;
                disp(filename)
                eval(['load ' datapath filename])
                partialdatmerged = [partialdatmerged; double(datmerged)]; %may 2006 heidi changed to concatenate instead of replace partialdatmerged, how did this work before??
                %                partialdatmerged = double(datmerged);
                filenum = filenum + 1;
            end;
            clear datendind
            %timebeadbins = timesectionendbin(sectionnum):timesectionendbin(sectionnum+1); %all
            %timebeadbins = find(totaltime(timesectionendbin(sectionnum):timesectionendbin(sectionnum+1),6) == 1); %syringe port = 1
            timebeadbins = find(totaltime(timestartind:timeendind,6) == beadport); %syringe port = 1, new timestartind...5/16/03 heidi, make sure to change "do a double..." above if port number changes for beads
            timebeadbins = timebeadbins + timestartind - 1;
            datbeadbins = [];
            for count = 1:length(timebeadbins),
                if count == 1 & timebeadbins(1) == 1,
                    datbeadbins = [datbeadbins 1:totaltime(timebeadbins(1),1)];
                else
                    datbeadbins = [datbeadbins totaltime(timebeadbins(count)-1,1)+1:totaltime(timebeadbins(count),1)];
                end;
            end;
            %datbeadbins = datbeadbins-double(partialdatmerged(1,1)) + 1;  %index into existing partialdatmerged, new 5/16/03, heidi
            [junk, junk, datbinstouse] = intersect(datbeadbins, partialdatmerged(:,1));
            clear junk
            maxvalue = 1e6;  %is this too high?
            %            bins = 10.^(0:log10(maxvalue)/1023:log10(maxvalue));  %make 1024 log spaced bins
            bins = 10.^(0:log10(maxvalue)/63:log10(maxvalue));  %make 1024 log spaced bins
            if ~isempty(datbinstouse),  %skip cases where no beads in hour
                [nmergedhist,x,nbins] = histmulti5(partialdatmerged(datbinstouse,[2,4:5]),[bins' bins' bins']);
                %add ad hoc criteria to force bead mode to be found, may
                %need more later like in beadbatch4_field...or maybe this is better and don't need so many?
                %force CHL and SSC to be > 2e4
                %keyboard  %SSC needs to be lower for 18 nov 2010; needs to be 2e4 for low bd num times in July 2010
                mind = find(bins < 2e4); %default
                if year2do == 2010 & strmatch('FCB2', filetypelist(typenum,1:4))
                    mind = find(bins < 1e4);
                end;
                if year2do > 2004,
                    nmergedhist(:,:,mind) = 0; %SSC
                end;
                mind = find(bins < 2e4);  
                nmergedhist(:,mind,:) = 0; %CHL
                %force PE > 1e2
                mind = find(bins < 1e2);
                nmergedhist(mind,:,:) = 0;
                [y,ind] = max(nmergedhist(:));
                [i,j,k] = ind2sub(size(nmergedhist),ind);
                modepos(:,[1,3:4]) = [i,j,k];
                ulim = 1.5;
                ulim = 1.8;
                llim = .5;
                beadstocount = find((partialdatmerged(datbinstouse,2) > bins(modepos(1))*llim) & (partialdatmerged(datbinstouse,2) < bins(modepos(1))*ulim));
                temp = datbinstouse(beadstocount);
                beadstocount = find((partialdatmerged(temp,4) > bins(modepos(3))*llim) & (partialdatmerged(temp,4) < bins(modepos(3))*ulim));
                temp = temp(beadstocount);
                beadstocount = find((partialdatmerged(temp,5) > bins(modepos(4))*llim) & (partialdatmerged(temp,5) < bins(modepos(4))*ulim));
                beadstocount = temp(beadstocount);
                
                if length(beadstocount) > 10,  %otherwise skip
                    %get better resolved mode positions for beads only points
                    numbins = 1024-1;
                    bins = 10.^(0:log10(maxvalue)/numbins:log10(maxvalue));
                    mergedhist = smooth(histc(partialdatmerged(beadstocount,2:5), bins),4);  %smooth over 2 with 512 ch, and 4 with 1024
                    [junk, modepos] = max(mergedhist);
                    bin2 = 10.^(0:log10(2^14)/numbins:log10(2^14));
                    temp = smooth(histc(partialdatmerged(datbinstouse,6), bin2),4);  %chl peak
                    [junk, modepos(:,5)] = max(temp);

                    beadresults(beadsection,1:2) = totaltime(timebeadbins(1),2:3);  % start & end time, days
                    beadresults(beadsection,3) = sum(totaltime(timebeadbins,4));  %acqtime, seconds
                    beadresults(beadsection,4) = length(beadstocount);  %number of beads
                    beadresults(beadsection,5:8) = bins(modepos(1:4));  %bead parameter values
                    beadresults(beadsection, 9) = bin2(modepos(5));  %chl peak
                    beadresults(beadsection,10:14) = mean(partialdatmerged(beadstocount,2:6));  %means
                    beadresults(beadsection,15:19) = std(partialdatmerged(beadstocount,2:6));  %bead std dev.
                    %beadresults(sectionnum,20) = NaN; %syringeinterval;
                    beadresults(beadsection,20) = sum(totaltime(timebeadbins,5)); %analyzed volume (ml)
                    beadsection = beadsection + 1;

                    %keyboard  %stop here and look at hist(partialdatmerged(datbinstouse,4), bins) to check for pre-peak triggers
                    if plotflag,
                        figure(1)
                        clf
                        bar(bins, mergedhist, 'histc')
                        legend('PE','FLS', 'CHL','SSC')
                        line(repmat([bins(modepos(1)); bins(modepos(1))*ulim; bins(modepos(1))*llim],1,2)', repmat([0 3500],3,1)', 'color', 'b', 'linestyle', '--', 'linewidth', 2)
                        line(repmat([bins(modepos(2)); bins(modepos(2))*ulim; bins(modepos(2))*llim],1,2)', repmat([0 3500],3,1)', 'color', 'c', 'linestyle', '--', 'linewidth', 2)
                        line(repmat([bins(modepos(3)); bins(modepos(3))*ulim; bins(modepos(3))*llim],1,2)', repmat([0 3500],3,1)', 'color', 'y', 'linestyle', '--', 'linewidth', 2)
                        line(repmat([bins(modepos(4)); bins(modepos(4))*ulim; bins(modepos(4))*llim],1,2)', repmat([0 3500],3,1)', 'color', 'r', 'linestyle', '--', 'linewidth', 2)
                        set(gca, 'xscale', 'log')
                        axis([1e3 1e6 0 20])
                        figure(2)
                        clf
                        subplot(121)
                        loglog(partialdatmerged(datbinstouse,3),partialdatmerged(datbinstouse,4),'.', 'markersize', .5)
                        axis('square')
                        xlabel('FLS')
                        ylabel('CHL')
                        line(repmat([bins(modepos(2))*ulim  bins(modepos(2))*llim],2,1), [1 1e6], 'color', 'k', 'linestyle', ':')
                        line([1 1e6], repmat([bins(modepos(3))*ulim  bins(modepos(3))*llim],2,1),  'color', 'k', 'linestyle', ':')
                        hold on
                        loglog(partialdatmerged(beadstocount,3),partialdatmerged(beadstocount,4),'r.', 'markersize', .5)
                        axis([1 1e6 1 1e6])
                        title(datestr(beadresults(beadsection-1,1)))
                        subplot(122)
                        loglog(partialdatmerged(datbinstouse,5),partialdatmerged(datbinstouse,2),'.', 'markersize', .5)
                        axis('square')
                        xlabel('SSC')
                        ylabel('PE')
                        line(repmat([bins(modepos(4))*ulim  bins(modepos(4))*llim],2,1), [1 1e6], 'color', 'k', 'linestyle', ':')
                        line([1 1e6], repmat([bins(modepos(1))*ulim  bins(modepos(1))*llim],2,1),  'color', 'k', 'linestyle', ':')
                        hold on
                        loglog(partialdatmerged(beadstocount,5),partialdatmerged(beadstocount,2),'r.', 'markersize', .5)
                        axis([1 1e6 1 1e6])
                        disp('pause for graph')
                        pause(0)
                    end;  %if 0 to plot or not to plot...
                else
                    disp(['Too few beads?:' num2str(length(beadstocount))])
                end; % if length(beadstocount) > 10),
            end; %if ~isempty(datbinstouse)
            partialdatmerged = double(datmerged); %reset partialdat with file partly completed
            sectionnum = sectionnum + 1;
        end;  %while sectionnum...%sectionnum = 1:length(timesectionendbin) - 1
        if beadsection > 1,
            beadtitles = {'start (day)' 'end (day)' 'acq time (sec)' 'bead number' 'beadmodePE' 'beadmodeFLS' 'beadmodeCHL' 'beadmodeSSC' 'beadmodeCHLpk' 'beadmeanPE' 'beadmeanFLS' 'beadmeanCHL' 'beadmeanSSC' 'beadmeanCHLpk' 'beadstdPE' 'beadstdFLS' 'beadstdCHL' 'beadstdSSC' 'beadstdCHLpk' 'analyzed volume (ml)'};
            eval(['save ' savepath filetypelist(typenum,:) 'beads_' num2str(sectcount) ' beadtitles beadresults'])
        end;
        clear beadresults %added April 2007 to fix big with repeating bead results in later short sections
    end %for sectcount
    clear beadresults  %added 10/18/05 to prevent extra rows from previous files, Heidi
end; %for typenum = 1:length(filelist)

beadsummary