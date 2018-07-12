%timebatch2 created from timebatch_labalt to read new file format,
%calls fcbtimeproc, 3/03 Heidi
%timebatch_labalt created from timebatch to call cytosubproc2_labalt, 3/7/03 Heidi
%timebatch for oc382 modified from leobatch2.m, 10/16/02 heidi

close all

for typenum = 1:size(filetypelist,1),  
    disp('initializing...')
    filelist = dir([datapath filetypelist(typenum,:) '*']);
    date = datenum(cat(1,filelist.date));
    [temp, fileorder] = sort(date);
    clear temp date
    filelistmain = filelist(fileorder); %consecutive order
    if year2do == 2008 %special case dealing with day of mixed local and UTC time stamps (22 Oct 2008)
        ii1 = strmatch('FCB1_2008_296_092206', {filelistmain.name});
        ii2 = strmatch('FCB1_2008_296_130826',{filelistmain.name});
        ii3 = strmatch('FCB1_2008_296_114008',{filelistmain.name});
        filelistmain(sort([ii1,ii2,ii3])) = filelistmain([ii1,ii2,ii3])
    end;
    filesections = ceil(length(filelist)/setsize);
    for sectcount = 1:filesections;
        if sectcount < filesections,
            filelist = filelistmain((sectcount-1)*setsize+1:sectcount*setsize);
        else
            filelist = filelistmain((sectcount-1)*setsize+1:end);
        end;
        totalstartsec = zeros(1e5*length(filelist),1);
        totalendsec = totalstartsec;
        syrpumpinfo = zeros(1e5*length(filelist),8);
        nextevent = 1;
        disp('reading...')
        for filenum = 1:length(filelist),       
            filename = filelist(filenum).name;  %fileorder makes sure that files are viewed in date order
            filedate = filelist(filenum).date;
            disp(filename)
            
            %[header] = fcbreadraw2(datapath, filename);  %just get header for times
            eval(['[header] = ' readrawstr '(datapath, filename);'])  %just get header for times
            if year2do <= 2005, %old file naming system requires time stamp from file info
                t = datevec(filedate);
                year = datenum(['1-0-' num2str(t(1))]);
                realdaysfile = datenum(filedate) - year; %get the real time (in days after Jan 1) each file was stored
                if year2do == 2003 %special case to handle some local time stamps in 2003
                    if realdaysfile < 300 %days before 27-Oct-2003
                        realdaysfile = realdaysfile + 4/24;
                    end
                end
                clear t
                timesec = header(:,4).*3600 + header(:,5) + header(:,6)./1000;  %time (sec) for all events in whole file
                timesec(find(timesec<timesec(1)))=timesec(find(timesec<timesec(1)))+2^32/1000;% fix rollovers - straighten up time
                endreference=timesec(end);
                timeendsec=realdaysfile*3600*24-(endreference-timesec);   %determine real time (days) for each event            
                timesec = header(:,1).*3600 + header(:,2) + header(:,3)./1000;  %time (sec) for all events in whole file
                timesec(find(timesec<timesec(1)))=timesec(find(timesec<timesec(1)))+2^32/1000;% fix rollovers - straighten up time
                timestartsec=realdaysfile*3600*24-(endreference-timesec);   %determine real time (days) for each event
            else %new files naming with time stamp
                % this new section uses filename for start time (rather than file closing date as previously), starting may 2006 at MVCO, heidi 6/7/06
                year = datenum(str2num(filename(6:9)),0,0,0,0,0);
                realdaysfile = str2num(filename(11:13)) + str2num(filename(15:16))/24 + str2num(filename(17:18))/60/24 + str2num(filename(19:20))/60/60/24;  % file start
                if year2do == 2008  %special case to handle some local time stamps in Oct/Nov 2008
                    if realdaysfile > 307,
                        realdaysfile = realdaysfile + 5/24; 
                    elseif realdaysfile > 297 | ismember(filename, {'FCB1_2008_296_114008' 'FCB1_2008_296_141136' 'FCB1_2008_296_190118' 'FCB1_2008_296_233516'}),
                        realdaysfile = realdaysfile + 4/24;
                    end;
                end;
                timesec = header(:,1).*3600 + header(:,2) + header(:,3)./1000;  %time (sec) for all events in whole file
                if numel(timesec) == 0, keyboard, end;
                timesec(find(timesec<timesec(1)))=timesec(find(timesec<timesec(1)))+2^32/1000;% fix rollovers - straighten up time
                startreference = timesec(1);
                timestartsec = realdaysfile*3600*24 + timesec - startreference;
                timesec = header(:,4).*3600 + header(:,5) + header(:,6)./1000;  %time (sec) for all events in whole file
                timesec(find(timesec<timesec(1)))=timesec(find(timesec<timesec(1)))+2^32/1000;% fix rollovers - straighten up time
                timeendsec = realdaysfile*3600*24 + timesec - startreference;
            end;
            lastevent = nextevent+length(timestartsec)-1;
            %next loop needed for a few files on May 11, 2007 when acq PC clock seemed to be slipping back in time by 30-60 seconds
            %between a few files, after that we reset clock and then on May 12 Rob replaced the battery, Heidi 14 May 2007
            %if year2do == 2007 || year2do == 2010,
            %if year == datenum('31-Dec-06') || strncmp(filename, 'FCB2_2010_235_234230', length(filename)) %only consider in 2007 and problem date in 2010
                if nextevent > 1 & totalstartsec(nextevent-1) > timestartsec(1) & totalstartsec(nextevent-2) < totalstartsec(nextevent-1),  %last part fudge for bad file ends
                    %keyboard
                    %estimate clock offset as time for nextevent-1 plus diff between start event 1 and 2
                    offset = (totalstartsec(nextevent-1)-timestartsec(1)) + (timestartsec(2) - timestartsec(1));
                    disp(['time offset fudged: ' num2str(round(totalstartsec(nextevent-1)-timestartsec(1) + (timestartsec(2) - timestartsec(1)))) ' seconds'])
                    %timestartsec = timestartsec + (totalstartsec(nextevent-1)-timestartsec(1)) + (timestartsec(2) - timestartsec(1));
                    timestartsec = timestartsec + offset;
                    timeendsec = timeendsec + offset;  %add same offset to end, June 2014   
                end;
            %end;
            totalstartsec(nextevent:lastevent) = timestartsec;
            totalendsec(nextevent:lastevent) = timeendsec;
            tt = find(timeendsec<timestartsec);
            if ~isempty(tt)
                keyboard
            end
            %change to end-1 for March 16, 2006 addition of extra header value at end, change to 7:14 so stable before and after extra header, Mar 2014
            syrpumpinfo(nextevent:lastevent,:) = header(:,7:14); %1=temp,2=humidity,3=start port,4=end port,5=start syr#,6=end syr#,7=start syr pos,8=end syr pos.
            nextevent = lastevent+1;
         % if strmatch(filename, 'FCB1_2011_241_173414')
         %     keyboard
         % end;
        end; %for filenum = 1:length(filelist)
                
        if length(filelist) > 0,
            disp('processing...')
            totalstartsec = totalstartsec(1:lastevent);
            totalendsec = totalendsec(1:lastevent);
            syrpumpinfo = syrpumpinfo(1:lastevent,:);
            
            %fcbtimeproc100;  %new format 4/2005 with 100 events per record (and new channels for SSC2 PMT)
            eval(timeprocstr)
            %reset start and stop times to be in days (instead of sec), 6.5.03 heidi
            outmatrix(:,2:3) = outmatrix(:,2:3)/3600/24 + year;
            
            if plotflag,    
                subplot(211)
                plot(outmatrix(:,2),outmatrix(:,4), '.')
                ylabel('Acquisition time (sec)')
                xlabel('Record start time (day of year)')
                title([filetypelist(typenum,:) '.*'])
                hold on
                %plot port numbers for data to use, with different colors
                colorstr = 'kkrkkg';
                for count = 1:6,
                    t = find(outmatrix(:,6) == count);
                    if ~isempty(t)
                        eval(['text(outmatrix(t,2), outmatrix(t,4), num2str(outmatrix(t,6)), ''color'', ''' num2str(colorstr(count)) ''')'])    
                    end;
                    t = find(outmatrix(:,6) == 99);  %start of gaps
                    text(outmatrix(t,2), outmatrix(t,4), 'G', 'color', 'k')        
                end;
                datetick('x', 6, 'keeplimits', 'keepticks')
                disp('Check plot...')
                subplot(212)
                for count = 1:6,
                    t = find(outmatrix(:,6) == count);
                    if ~isempty(t)
                        eval(['plot(outmatrix(t,2), outmatrix(t,4), ''' num2str(colorstr(count)) '.'')'])    
                    end;
                    hold on    
                end;
                datetick('x', 6, 'keeplimits', 'keepticks')
                ylabel('Syringe volume analyzed (  \ml)')
                xlabel('Record start time (date)')
                pause
                clf
            end;
            clear colorstr t count
            
            eval([filetypelist(typenum,:) 'time = outmatrix;']);       
            %        timetitles = {'endrec index number' 'rec start (sec)' 'rec end (sec)' 'acq time (sec)' 'median acq interval (sec)' 'flag'}; 
            timetitles = {'endrec index number' 'rec start (days)' 'rec end (days)' 'acq time (sec)' 'volume analyzed (ml)' 'flag'}; 
            eval(['save ' procpath filetypelist(typenum,:) 'time_' num2str(sectcount) ' timetitles syrpumpinfo ' filetypelist(typenum,:) 'time']) %datatransinterval'])
            eval(['clear ' filetypelist(typenum,:) 'time'])
            clear time*sec timetitles totalstartsec totalendsec outmatrix *event endreference header realdaysfile syrpumpinfo
        end;  %if length(filelist) > 0
    end; %for sectcount
end; %for typenum = 1:length(filetypelist)

clear filelist filename filedate fileorder filenum typenum