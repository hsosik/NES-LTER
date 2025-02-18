% A carbon copy of beadbatch6 (as of 11 July 2016) from Heidi, with some
% added lines to capture the plots into structure movies to be played in
% MATLAB with implay:

Q=0;
bead_date=[];
figure, set(gcf,'Position',[3         120        1583         846])
cmap=parula(4); %for line colors to match histograms
cmap(4,:)=[ 0.9       0.5       0.1]; %better orange color instead of yellow

timeinterval = 1/24;  %sec (1/24 = 1 hr), resolution for final values

for typenum = 1:size(filetypelist,1)
    filelist = dir([datapath filetypelist(typenum,:) '*.mat']);
    %  date = datenum(cat(1,filelist.date));
    %  [temp, fileorder] = sort(date);
    filelistmain = filelist;
    if year2do <= 2005
        [~,fileorder] = sort(str2num(char(regexprep(regexprep({filelist.name}, '.mat', ''), filetypelist(typenum,:), ''))));
        filelistmain = filelist(fileorder); %consecutive order
    end
    if year2do == 2008 %special case dealing with day of mixed local and UTC time stamps (22 Oct 2008)
        ii1 = strmatch('FCB1_2008_296_092206', {filelistmain.name});
        ii2 = strmatch('FCB1_2008_296_130826',{filelistmain.name});
        ii3 = strmatch('FCB1_2008_296_114008',{filelistmain.name});
        filelistmain(sort([ii1,ii2,ii3])) = filelistmain([ii1,ii2,ii3])
    end
    clear temp date fileorder
    
    filesections = ceil(length(filelistmain)/setsize);
    for sectcount = 1:filesections
        if sectcount < filesections
            filelist = filelistmain((sectcount-1)*setsize+1:sectcount*setsize);
        else
            filelist = filelistmain((sectcount-1)*setsize+1:end);
        end
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
        beadresults = NaN(1,33);
        %    for sectionnum = 1:length(timesectionendbin) - 1,
        sectionnum = 1;
        while sectionnum <= length(timesectionendbin)
            %        disp(['sectionnum: ' num2str(sectionnum)])
            timeendind = timesectionendbin(sectionnum); %new 5/16/03, heidi
            if sectionnum > 1, timestartind = timesectionendbin(sectionnum-1) + 1; end;
            if totaltime(timeendind,6) == beadport & sectionnum ~= length(timesectionendbin), %do a double "section" if one ends in the middle of bead analysis
                sectionnum = sectionnum + 1;
                timeendind = timesectionendbin(sectionnum);
            end
            datendind = totaltime(timeendind,1); %this line moved after do a double loop, Feb 2014 - otherwise the loop doesn't actually do anything with more data
            while (partialdatmerged(end,1) < datendind) & filenum <= length(filelist),  %keep adding on files until get full hour
                %        while totaltime(timesectionendbin(sectionnum),1) > partialdatmerged(end,1),  %get to first file
                filename = filelist(filenum).name;
                disp(filename)
                eval(['load ' datapath filename])
                partialdatmerged = [partialdatmerged; double(datmerged)]; %may 2006 heidi changed to concatenate instead of replace partialdatmerged, how did this work before??
                %                partialdatmerged = double(datmerged);
                filenum = filenum + 1;
            end
            clear datendind
            %timebeadbins = timesectionendbin(sectionnum):timesectionendbin(sectionnum+1); %all
            %timebeadbins = find(totaltime(timesectionendbin(sectionnum):timesectionendbin(sectionnum+1),6) == 1); %syringe port = 1
            timebeadbins = find(totaltime(timestartind:timeendind,6) == beadport); %syringe port = 1, new timestartind...5/16/03 heidi, make sure to change "do a double..." above if port number changes for beads
            timebeadbins = timebeadbins + timestartind - 1;
            datbeadbins = [];
            for count = 1:length(timebeadbins)
                if count == 1 & timebeadbins(1) == 1
                    datbeadbins = [datbeadbins 1:totaltime(timebeadbins(1),1)];
                else
                    datbeadbins = [datbeadbins totaltime(timebeadbins(count)-1,1)+1:totaltime(timebeadbins(count),1)];
                end
            end
            %datbeadbins = datbeadbins-double(partialdatmerged(1,1)) + 1;  %index into existing partialdatmerged, new 5/16/03, heidi
            [~, ~, datbinstouse] = intersect(datbeadbins, partialdatmerged(:,1));
            %clear junk
            maxvalue = 1e6;  %is this too high?
            %            bins = 10.^(0:log10(maxvalue)/1023:log10(maxvalue));  %make 1024 log spaced bins
            bins = 10.^(0:log10(maxvalue)/63:log10(maxvalue));  %make 1024 log spaced bins
            if ~isempty(datbinstouse)  %skip cases where no beads in hour
                [nmergedhist1,x,nbins] = histmulti5(partialdatmerged(datbinstouse,[2,4:5]),[bins' bins' bins']);
                %add ad hoc criteria to force bead mode to be found, may
                %need more later like in beadbatch4_field...or maybe this is better and don't need so many?
                %force CHL and SSC to be > 2e4
                %keyboard  %SSC needs to be lower for 18 nov 2010; needs to be 2e4 for low bd num times in July 2010
                nmergedhist = nmergedhist1;
                mind = find(bins < 2e4); %default
                if year2do == 2010 & strmatch('FCB2', filetypelist(typenum,1:4))
                    mind = find(bins < 1e4);
                elseif ismember(cellstr(filename(1:13)),strcat('FCB2_2023_0', num2str([35:45]')))
                    %handle case 4-14 Feb 2023 when SSC drops very low
                    mind = find(bins < 200);
                end
                if year2do > 2004
                    nmergedhist(:,:,mind) = 0; %SSC
                end
                mind = find(bins < 2e4);
                nmergedhist(:,mind,:) = 0; %CHL
                %force PE > 1e2
                mind = find(bins < 1e2);
                nmergedhist(mind,:,:) = 0;
                [y,ind] = max(nmergedhist(:));
                [i,j,k] = ind2sub(size(nmergedhist),ind);
                modepos(:,[1,3:4]) = [i,j,k];
                ulim = 1.5;
                llim = .5;
                beadstocount = find((partialdatmerged(datbinstouse,2) > bins(modepos(1))*llim) & (partialdatmerged(datbinstouse,2) < bins(modepos(1))*ulim));
                temp = datbinstouse(beadstocount);
                beadstocount = find((partialdatmerged(temp,4) > bins(modepos(3))*llim) & (partialdatmerged(temp,4) < bins(modepos(3))*ulim));
                temp = temp(beadstocount);
                beadstocount = find((partialdatmerged(temp,5) > bins(modepos(4))*llim) & (partialdatmerged(temp,5) < bins(modepos(4))*ulim));
                beadstocount = temp(beadstocount);
                %now get the smaller beads
                nmergedhist = nmergedhist1; %reset to full 3D histogram
                %mind = find(bins < 1e3); %default
                mind = find(bins < bins(modepos(4))*.12 | bins > bins(modepos(4))*.5); %force larger smaller than large beads, but not too small on SSC
                %if year2do == 2010 & strmatch('FCB2', filetypelist(typenum,1:4))
                %    mind = find(bins < 1e4);
                %end;
                %if year2do > 2004,
                nmergedhist(:,:,mind) = 0; %SSC
                %end;
                mind = find(bins > 1e4);
                nmergedhist(:,mind,:) = 0; %CHL
                %force PE > 1e2
                mind = find(bins < 3e2); %original 1e2, but some high noise cases in 2008 - 300 okay? 12 July 2016 KRHC
                nmergedhist(mind,:,:) = 0;
                mind = find(bins > bins(modepos(1))); % force PE < large beads
                nmergedhist(mind,:,:) = 0;
                [y,ind] = max(nmergedhist(:));
                [i,j,k] = ind2sub(size(nmergedhist),ind);
                modepos(:,[1,3:4]) = [i,j,k];
                ulim = 1.5;
                llim = .7; %adjust below to use 0.7 to minimize clipping into syn
                beadstocount2 = find((partialdatmerged(datbinstouse,2) > bins(modepos(1))*llim) & (partialdatmerged(datbinstouse,2) < bins(modepos(1))*ulim));
                temp = datbinstouse(beadstocount2);
                beadstocount2 = find((partialdatmerged(temp,4) > bins(modepos(3))*llim) & (partialdatmerged(temp,4) < bins(modepos(3))*ulim));
                temp = temp(beadstocount2);
                beadstocount2 = find((partialdatmerged(temp,5) > bins(modepos(4))*llim) & (partialdatmerged(temp,5) < bins(modepos(4))*ulim));
                beadstocount2 = temp(beadstocount2);
                %
                bead1flag = 0;
                if length(beadstocount) > 10  %otherwise skip
                    bead1flag = 1;
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
                    bead2flag = 0;
                    if length(beadstocount2) > 10
                        mergedhist2 = smooth(histc(partialdatmerged(beadstocount2,2:5), bins),4);  %smooth over 2 with 512 ch, and 4 with 1024
                        [junk, modepos2] = max(mergedhist2);
                        
                        % if (bins(modepos2(4))./bins(modepos(4)) > .125) & (bins(modepos2(3))./bins(modepos(3)) < .015) &  %otherwise not bead mode;
                        if bins(modepos2(4))./bins(modepos(4)) > .125 & bins(modepos2(3))./bins(modepos(3)) < .015 %otherwise not bead mode, likely because too few
                            
                            % need to tighten bead window, as can have some Syn carry over:
                            ulim2=1.5;
                            llim2=0.75;
                            
                            tightbeads2 = find((partialdatmerged(beadstocount2,2) > bins(modepos2(1))*llim2) & (partialdatmerged(beadstocount2,2) < bins(modepos2(1))*ulim2)...
                                & (partialdatmerged(beadstocount2,4) > bins(modepos2(3))*llim2) & (partialdatmerged(beadstocount2,4) < bins(modepos2(3))*ulim2)...
                                & (partialdatmerged(beadstocount2,5) > bins(modepos2(4))*llim2) & (partialdatmerged(beadstocount2,5) < bins(modepos2(4))*ulim2));
                            
                            %bin2 = 10.^(0:log10(2^14)/numbins:log10(2^14));
                            %temp = smooth(histc(partialdatmerged(datbinstouse,6), bin2),4);  %chl peak
                            %[junk, modepos2(:,5)] = max(temp);
                            beadstocount2=beadstocount2(tightbeads2);
                            bead2flag = 1;
                            %then after tightening up, does the distribution still have a spread out nature?
                            %If so, likely just noise that got past the other flags...
                            mm=nanmean(partialdatmerged(beadstocount2,2:5),1); %mean of PE FSC CHL SSC
                            ss=nanstd(partialdatmerged(beadstocount2,2:5),0,1); %stddeviation of PE FSC CHL SSC
                            coeffvar=mm./ss; %mean/std deviation -> use a scaled SSC CV based on beadnum for cutoff?
                            
                            if coeffvar(4) > 7 && coeffvar(3) > 6 %further check to make sure there is not too much noise on ssc..%(length(beadstocount2)/coeffvar(4) < 10
                                beadresults(beadsection,21) = length(beadstocount2);  %number of beads
                                beadresults(beadsection,22:25) = bins(modepos2(1:4));  %bead parameter values
                                %beadresults(beadsection, 26) = bin2(modepos2(5));  %chl peak
                                beadresults(beadsection,26:29) = mean(partialdatmerged(beadstocount2,2:5));  %means
                                beadresults(beadsection,30:33) = std(partialdatmerged(beadstocount2,2:5));  %bead std dev.
                            end
                        end
                    end
                    beadsection = beadsection + 1;                 
                else
                    disp(['Too few beads?:' num2str(length(beadstocount))])
                    %keyboard
                end; % if length(beadstocount) > 10),
                
                
                %---------------------PLOTTING AND MOVIE------------------------------------------------------------------------------------------
                
                % keyboard  %stop here and look at hist(partialdatmerged(datbinstouse,4), bins) to check for pre-peak triggers
                if plotflag || beadmovieflag,
                    
                    Q=Q+1;
                    if bead1flag %made counts of beads :)
                        subplot(2,3,1,'replace')
                        bar(bins, mergedhist, 'histc'), colormap(cmap)
                        legend('PE','FLS', 'CHL','SSC')
                        line(repmat([bins(modepos(1)); bins(modepos(1))*ulim; bins(modepos(1))*llim],1,2)', repmat([0 3500],3,1)', 'color', cmap(1,:), 'linestyle', '--', 'linewidth', 2)
                        line(repmat([bins(modepos(2)); bins(modepos(2))*ulim; bins(modepos(2))*llim],1,2)', repmat([0 3500],3,1)', 'color', cmap(2,:), 'linestyle', '--', 'linewidth', 2)
                        line(repmat([bins(modepos(3)); bins(modepos(3))*ulim; bins(modepos(3))*llim],1,2)', repmat([0 3500],3,1)', 'color', cmap(3,:), 'linestyle', '--', 'linewidth', 2)
                        line(repmat([bins(modepos(4)); bins(modepos(4))*ulim; bins(modepos(4))*llim],1,2)', repmat([0 3500],3,1)', 'color', cmap(4,:), 'linestyle', '--', 'linewidth', 2)
                        set(gca, 'xscale', 'log')
                        axis([1e3 1e6 0 20])
                        title('Large bead histogram')
                        
                        if length(beadstocount2) > 10 %then go ahead and make histogram...
                            
                            subplot(2,3,2,'replace')
                            bar(bins, mergedhist2, 'histc'), colormap(cmap) %does not have tightbeads....
                            legend('PE','FLS', 'CHL','SSC')
                            line([.015*bins(modepos(3)) .015*bins(modepos(3))],[0 3500],'color','k','linestyle','--')
                            line([.125*bins(modepos(4)) .125*bins(modepos(4))],[0 3500],'color','k','linestyle','--')
                            %box to exclude:
                            line(repmat([bins(modepos2(1))*ulim2; bins(modepos2(1))*llim2],1,2)', repmat([0 3500],2,1)', 'color', cmap(1,:), 'linestyle', '--', 'linewidth', 2)
                            line(repmat([bins(modepos2(3))*ulim2; bins(modepos2(3))*llim2],1,2)', repmat([0 3500],2,1)', 'color', cmap(3,:), 'linestyle', '--', 'linewidth', 2)
                            line(repmat([bins(modepos2(4))*ulim2; bins(modepos2(4))*llim2],1,2)', repmat([0 3500],2,1)', 'color', cmap(4,:), 'linestyle', '--', 'linewidth', 2)
                            set(gca, 'xscale', 'log')
                            axis([1e3 1e6 0 20])
                            
                            title('Small bead histogram')
                            
                            subplot(2,3,3,'replace')
                            set(gca,'visible','off')
                            text(0.2, 0.7,{['Num of large beads: '  num2str(length(beadstocount))];
                                ['Num of small beads: '  num2str(length(beadstocount2))]})
                            if bins(modepos2(4))./bins(modepos(4)) > 0.125
                                text(0.2,0.5,['SSC Critera: ' num2str(bins(modepos2(4))./bins(modepos(4))) ' | should be: >' num2str(0.125)],'color',[0 0.7 0])
                            else
                                text(0.2,0.5,['SSC Critera: ' num2str(bins(modepos2(4))./bins(modepos(4))) ' | should be: >' num2str(0.125)],'color','r')
                            end
                            
                            if bins(modepos2(3))./bins(modepos(3)) < 0.015
                                text(0.2,0.4,['Chl Criteria: ' num2str(bins(modepos2(3))./bins(modepos(3)))  ' | should be: <' num2str(0.015)],'color',[0 0.7 0])
                            else
                                text(0.2,0.4,['Chl Criteria: ' num2str(bins(modepos2(3))./bins(modepos(3)))  ' | should be: <' num2str(0.015)],'color',[1 0 0])
                            end
                            
                            if bead2flag
                                if coeffvar(4) > 7
                                % {['CV scaled: ' num2str(length(beadstocount2)/coeffvar(4))  ' | should be: <' num2str(10) '  or... ']
                                    text(0.2,0.3,['CV: ' num2str(coeffvar(4))  ' | should be: >' num2str(7)],'color',[0 0.7 0])
                                else
                                %{['CV scaled: ' num2str(length(beadstocount2)/coeffvar(4))  ' | should be: <' num2str(10) '  or... '];
                                    text(0.2,0.3,['CV: ' num2str(coeffvar(4))  ' | should be: >' num2str(7)],'color',[1 0 0])
                                end
                                if coeffvar(3) > 6
                                % {['CV scaled: ' num2str(length(beadstocount2)/coeffvar(4))  ' | should be: <' num2str(10) '  or... ']
                                    text(0.2,0.2,['CV: ' num2str(coeffvar(3))  ' | should be: >' num2str(6)],'color',[0 0.7 0])
                                else
                                %{['CV scaled: ' num2str(length(beadstocount2)/coeffvar(4))  ' | should be: <' num2str(10) '  or... '];
                                    text(0.2,0.2,['CV: ' num2str(coeffvar(3))  ' | should be: >' num2str(6)],'color',[1 0 0])
                                end
                            end
                        else
                            subplot(2,3,2,'replace')
                            set(gca,'visible','off')
                            text(0.3,0.5,['Too few small beads: ' num2str(length(beadstocount2))])
                            
                            subplot(2,3,3,'replace')
                            set(gca,'visible','off')
                        end
                                             
                   else
                        subplot(2,3,1,'replace')
                        set(gca,'visible','off')
                        text(0.3,0.5,['Too few beads:' num2str(length(beadstocount)) '?'])
                        
                        subplot(2,3,2,'replace')
                        set(gca,'visible','off')
                        text(0.3,0.5,['Didn''t process any small beads: ' num2str(length(beadstocount2)) '?'])
                        
                        subplot(2,3,3,'replace')
                        set(gca,'visible','off')
                        
                    end

                    subplot(2,3,4,'replace')
                    loglog(partialdatmerged(datbinstouse,5),partialdatmerged(datbinstouse,4),'.', 'markersize', .5)
                    axis('square')
                    %xlabel('FLS')
                    xlabel('SSC')
                    ylabel('CHL')
                    if bead1flag %add the lines :)
                        line(repmat([bins(modepos(4))*ulim  bins(modepos(4))*llim],2,1), [1 1e6], 'color', 'k', 'linestyle', ':')
                        line([1 1e6], repmat([bins(modepos(3))*ulim  bins(modepos(3))*llim],2,1),  'color', 'k', 'linestyle', ':')
                        line([1 1e6],[.015*bins(modepos(3)) .015*bins(modepos(3))],'color','g','linestyle',':')
                        line([.125*bins(modepos(4)) .125*bins(modepos(4))],[1 1e6],'color','g','linestyle',':')
                    end
                    hold on
                    loglog(partialdatmerged(beadstocount,5),partialdatmerged(beadstocount,4),'r.', 'markersize', .5)
                    %for small bead cutoffs:
                    
                    if bead2flag
                        loglog(partialdatmerged(beadstocount2,5),partialdatmerged(beadstocount2,4),'.g', 'markersize', 1)
                    end
                    
                    axis([1 1e6 1 1e6])
                    
                    if bead1flag
                        title(datestr(beadresults(beadsection-1,1)))
                    else
                        title(datestr(totaltime(timebeadbins(1),2)))% is this correct for time?
                    end
                    
                    subplot(2,3,5,'replace')
                    loglog(partialdatmerged(datbinstouse,5),partialdatmerged(datbinstouse,2),'.', 'markersize', .5)
                    axis('square')
                    xlabel('SSC')
                    ylabel('PE')
                    if bead1flag
                        line(repmat([bins(modepos(4))*ulim  bins(modepos(4))*llim],2,1), [1 1e6], 'color', 'k', 'linestyle', ':')
                        line([1 1e6], repmat([bins(modepos(1))*ulim  bins(modepos(1))*llim],2,1),  'color', 'k', 'linestyle', ':')
                        line([.125*bins(modepos(4)) .125*bins(modepos(4))],[1 1e6],'color','g','linestyle',':')
                    end
                    hold on
                    loglog(partialdatmerged(beadstocount,5),partialdatmerged(beadstocount,2),'r.', 'markersize', .5)
                    %for small bead cutoffs:
                    if bead2flag
                        loglog(partialdatmerged(beadstocount2,5),partialdatmerged(beadstocount2,2),'g.', 'markersize', 1)
                    end
                    axis([1 1e6 1 1e6])
                    
                    subplot(2,3,6,'replace')
                    loglog(partialdatmerged(datbinstouse,4),partialdatmerged(datbinstouse,2),'.', 'markersize', .5)
                    axis('square')
                    xlabel('CHL')
                    ylabel('PE')
                    if bead1flag
                        line(repmat([bins(modepos(3))*ulim  bins(modepos(3))*llim],2,1),[1 1e6], 'color', 'k', 'linestyle', ':')                       
                        line([1 1e6],repmat([bins(modepos(1))*ulim  bins(modepos(1))*llim],2,1), 'color', 'k', 'linestyle', ':')                       
                        line([.125*bins(modepos(4)) .125*bins(modepos(4))],[1 1e6],'color','g','linestyle',':')
                    end
                    hold on
                    loglog(partialdatmerged(beadstocount,4),partialdatmerged(beadstocount,2),'r.', 'markersize', .5)
                    %for small bead cutoffs:
                    if bead2flag
                        loglog(partialdatmerged(beadstocount2,4),partialdatmerged(beadstocount2,2),'g.', 'markersize', 1)
                    end
                    axis([1 1e6 1 1e6])
                    
                    %                         mm=nanmean(partialdatmerged(beadstocount2,2:5),1);
                    %                         ss=nanstd(partialdatmerged(beadstocount2,2:5),0,1);
                    %                         disp(['Num of beads2: ' num2str(length(beadstocount2))])
                    %                         disp(['CV: ' num2str(mm./ss)])
                    %                          disp(['3rd criteria < 10: ' num2str(length(beadstocount2)/(mm(4)./ss(4)))])
                    %
                    if plotflag
                        disp('pause for graph')
                        if length(beadstocount2) > 10
                            disp(['CV: ' num2str(mm./ss)])
                            disp(['STD: ' num2str(ss)])
                            disp(['Scaled STD: ' num2str(1/length(beadstocount2)*ss)])
                        end
                        pause%(0)
                    elseif beadmovieflag %capture frame for movie
                        F1(Q)=getframe(gcf);
                        if bead1flag
                            bead_date=[bead_date; beadresults(beadsection-1,1)];
                        else
                            bead_date=[bead_date; totaltime(timebeadbins(1),2)];
                        end
                    end
                    
                    %  F1(Q) = im2frame(zbuffer_cdata(fig1))
                end  %if 0 to plot or not to plot...
                
            end %if ~isempty(datbinstouse)
            partialdatmerged = double(datmerged); %reset partialdat with file partly completed
            sectionnum = sectionnum + 1
        end  %while sectionnum...%sectionnum = 1:length(timesectionendbin) - 1
        if beadsection > 1
            beadtitles = {'start (day)' 'end (day)' 'acq time (sec)' 'bead number' 'beadmodePE' 'beadmodeFLS' 'beadmodeCHL' 'beadmodeSSC' 'beadmodeCHLpk' 'beadmeanPE' 'beadmeanFLS' 'beadmeanCHL' 'beadmeanSSC' 'beadmeanCHLpk' 'beadstdPE' 'beadstdFLS' 'beadstdCHL' 'beadstdSSC' 'beadstdCHLpk' 'analyzed volume (ml)' 'bead2 number' 'bead2modePE' 'bead2modeFLS' 'bead2modeCHL' 'bead2modeSSC' 'bead2meanPE' 'bead2meanFLS' 'bead2meanCHL' 'bead2meanSSC' 'bead2stdPE' 'bead2stdFLS' 'bead2stdCHL' 'bead2stdSSC'};
            eval(['save ' savepath filetypelist(typenum,:) 'beads_' num2str(sectcount) ' beadtitles beadresults'])
        end
        clear beadresults %added April 2007 to fix big with repeating bead results in later short sections
    end %for sectcount
    clear beadresults  %added 10/18/05 to prevent extra rows from previous files, Heidi
end %for typenum = 1:length(filelist)

if beadmovieflag
    [temp ii]=sort(bead_date);
    eval(['beads_' num2str(year2do) '=F1(ii);'])
    eval(['save ' savepath 'beadmovie' num2str(year2do) ' beads_' num2str(year2do)])
end

beadsummary