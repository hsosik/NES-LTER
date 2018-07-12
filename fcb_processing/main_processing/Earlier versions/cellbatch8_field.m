%cellbatch7_field created from cellbatch6_field, now use bead modes to
%eliminate large and small beads only in 2 hours right after bead run
%cellbatch5_field created from cellbatch4_field for MVCO data, corrected
%logic error in pe separation of Syn, junk, cryto (junkind with OR not
%AND), April 2007 Heidi
%cellbatch4_field created from cellbatch3_field for MVCO May 2003 data (new
%PE separation since signals for large euks are above baseline), 5/11/03 Heidi
%cellbatch3_field created from cellbatch3 to handle single sample syringe port (6)
%cellbatch3 created from cellbatch2_labalt, for new file format, 3/03 Heidi
%cellbatch2 for MVCO, modified from cellbatch.m to add simple euk junk
%discrimination and simple crypto separations
%created from cellproc, separate with and without pe, no clustering (for micro&8018 exp't) heidi 1/10/02
%no bead match yet since no beads in data

%process for cells
%heiditemp4 - separate with and without pe, then cluster

%timeinterval = 3600;  %sec (3600 = 1 hr), resolution for final values

timeinterval = 1/24;

if year2do <= 2007,
    mergedtitles = {'rec number' 'PE' 'FLS' 'CHL' 'SSC' 'CHLpk' 'Class'};
else
    mergedtitles = {'rec number' 'PE' 'FLS' 'CHL' 'SSC' 'CHLpk1' 'CHLpk2' 'Class'};
end;

%beadmatchtitles = {'start time (matlab days)' 'beadPE' 'beadFLS' 'beadCHL' 'beadSSC' 'beadnumber' 'bead acq time (min)' 'bead pump rate (ml/min)'};
beadmatchtitles = {'start time (matlab days)' 'beadPE' 'beadFLS' 'beadCHL' 'beadSSC' 'beadnumber' 'bead acq time (min)' 'analyzed volume (ml)'};  %april 2007
classnotes = 'Class numbers in last column of mergedwithclass: 1 = syn; 2=bright cryptos; 3 = junk w/pe; 4 = euks (i.e., no pe); 5 = euk junk; 6 = dim cryptos';
numcluster = 6;

eval(['load ' beadpath 'beadresults'])  %load bead result file
%beadresults = ones(1,10);

culture = cellport;  %port for samples
clear cell*
cellrestitles = {'mean start time (matlab days)' 'acq time (min)' 'volume analyzed (ml)'};
for typenum = 1:size(filetypelist,1),  
    filelist = dir([procpath filetypelist(typenum,:) '*.mat']);
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
    for sectcount = 1:filesections %3
        if sectcount < filesections,
            filelist = filelistmain((sectcount-1)*setsize+1:sectcount*setsize);
        else
            filelist = filelistmain((sectcount-1)*setsize+1:end);
        end;
        eval(['load ' timepath filetypelist(typenum,:) 'time_' num2str(sectcount)])  %load file with processed time
        eval(['totaltime = ' filetypelist(typenum,:) 'time; clear ' filetypelist(typenum,:) 'time']) 
        %    pumprate = 0.05;  %assume for all of OC382
        
        timesectionendbin = find(diff(floor(totaltime(:,2)/timeinterval)));  %location of hour changes, last point in hr
        timesectionendbin =  [timesectionendbin; size(totaltime,1)];
        %timesectionendbin = find(totaltime(:,6) == 99);  %these are gaps between sets (as in lab)
        %timesectionenedbin = timesectionendbin(2:2:end);  %skip odd gaps for 2 culture alternation -- NOT GENERAL
        %    timesectionendbin =  [timesectionendbin; size(totaltime,1)];
        
        timestartind = 1;
        %    timestartind = timesectionendbin(1);
        filenum = 1;  %15
        filename = filelist(filenum).name;
        disp(filename) 
        eval(['load ' procpath filename])
        fit1 = fit; clear fit %reserved function name in matlab now
        partialdatmerged = datmerged;
        cellPE = NaN(length(timesectionendbin),numcluster);
        cellFLS = NaN;
        cellCHL = NaN;
        cellSSC = NaN;
        cellPEmode = NaN;
        cellFLSmode = NaN;
        cellCHLmode = NaN;
        cellSSCmode = NaN;
        cellNUM = NaN;
        cellresults = NaN(length(timesectionendbin),3);
        beadmatch = NaN(length(timesectionendbin),13);
        for sectionnum = 1:length(timesectionendbin) %7
            % for sectionnum = 1:length(timesectionendbin)-1,
            %disp(['sectionnum ' num2str(sectionnum)])
            timeendind = timesectionendbin(sectionnum);
            %timeendind = timesectionendbin(sectionnum+1);
            %following for case where start at sectionnum > 1, otherwise could reinit timestartind at end of loop with partialdatmerged
            if sectionnum > 1, timestartind = timesectionendbin(sectionnum-1) + 1; end;
            %if sectionnum > 1, timestartind = timesectionendbin(sectionnum) + 1; end;
            datendind = totaltime(timeendind,1);      
            while (partialdatmerged(end,1) < datendind) & filenum < length(filelist),  %keep adding on files until get full hour
                filenum = filenum + 1;    
                filename = filelist(filenum).name;
                disp(filename)
                eval(['load ' procpath filename])
                fit1 = fit; clear fit %reserved function name in matlab now
                partialdatmerged = [partialdatmerged; datmerged];    
            end;
            clear datendind
            %        goodtimebins = find(totaltime(timestartind:timeendind,6) == S | totaltime(timestartind:timeendind,6) == 0);
            goodtimebins = find(totaltime(timestartind:timeendind,6) == culture);
            goodtimebins = goodtimebins + timestartind - 1;
            
            if ~isempty(goodtimebins),  %added 3/11/03 to handle sections with no good data
                datbins = [];
                for count = 1:length(goodtimebins),
                    if count == 1 & goodtimebins(1) == 1,
                        datbins = [datbins 1:totaltime(goodtimebins(1),1)];
                    else
                        datbins = [datbins totaltime(goodtimebins(count)-1,1)+1:totaltime(goodtimebins(count),1)];  
                    end;
                end;
                datbins = datbins-double(partialdatmerged(1,1)) + 1;  %index into existing partialdatmerged
           
                cellresults(sectionnum,1) = mean(totaltime(goodtimebins,2));  %mean start time (days)
                cellresults(sectionnum,2) = sum(totaltime(goodtimebins,4))/60;  %acqtime (min)
                cellresults(sectionnum,3) = sum(totaltime(goodtimebins,5));  %vol analyzed from syringe positions (ml)
                [junk, beadind] = min(abs(cellresults(sectionnum,1) - beadresults(:,1)));  
                
                beadmatch(sectionnum,1) = beadresults(beadind,1);
                beadmatch(sectionnum,2:7) = beadresults(beadind,[5:8,4,3]);
                beadmatch(sectionnum,8) = beadresults(beadind,20); %April 2007, now analyzed ml
                beadmatch(sectionnum,9) = beadresults(beadind,21); %bead2 number
                beadmatch(sectionnum,10:13) = beadresults(beadind,22:25); %bead2 modes
                clear junk beadind
                
                partialdatmerged = double(partialdatmerged);
                partialdatmerged(:,2:5) = partialdatmerged(:,2:5) + 1;
                
                a = find(partialdatmerged(datbins,4) ~= 1 & partialdatmerged(datbins,5) ~= 1);  %zero chl is not allowed...also skip 0 SSC
                datbins2 = datbins(a);
                clear a
                
maxvalue = 1e6;
bins = 10.^(0:log10(maxvalue)/15:log10(maxvalue));
binc = [mean([bins(1:end-1); bins(2:end)]) maxvalue];
maxvaluessc = 1e7;
bins_ssc = 10.^(1:(log10(maxvaluessc)-1)/15:log10(maxvaluessc));
binc_ssc = [mean([bins_ssc(1:end-1); bins_ssc(2:end)]) maxvaluessc];
ii2 = find(partialdatmerged(datbins2,4)./partialdatmerged(datbins2,2)<=4 & partialdatmerged(datbins2,2)>100 & partialdatmerged(datbins2,4)./partialdatmerged(datbins2,2)>.05);
h = histmulti5([partialdatmerged(datbins2(ii2),4),partialdatmerged(datbins2(ii2),2)],[bins; bins]');
%flipud(h(1:end-2,:)')
[n,p] = max(h);
t = find(n>=5);
fitobject = fit(log10(binc(t))', log10(binc(p(t)))', 'poly1');
%ci = confint(fitobject);
c = coeffvalues(fitobject);
figure(99), clf
loglog(partialdatmerged(datbins2,4), partialdatmerged(datbins2,2), '.')
hold on
loglog(partialdatmerged(datbins2(ii2),4), partialdatmerged(datbins2(ii2),2), '.r')
plot(binc(p(t)),binc(t), 'g*')
hold on
yest = binc(t).^(c(1)).*10^(c(2));
plot(yest, binc(t), 'g-')
%yest = (bins(t)).^(c(1)).*10^(ci(2,2));
%plot(log10(yest), log10(bins(t)), 'g-')


%ii = find(partialdatmerged(datbins2,4)>= partialdatmerged(datbins2,2).^c(1).*10^(c(2)+.65));
ii = find(partialdatmerged(datbins2,4)>= partialdatmerged(datbins2,2).^c(1).*10^(c(2)+.5) | partialdatmerged(datbins2,2) < 100);
loglog(partialdatmerged(datbins2(ii),4), partialdatmerged(datbins2(ii),2), '.g')

h = histmulti5([partialdatmerged(datbins2(ii),5),partialdatmerged(datbins2(ii),2)],[bins_ssc; bins]');
[nt,pt] = max(h');
mask = ones(size(h));
for idx = 1:length(bins), mask(idx,1:pt(idx)) = NaN; end;
h2 = h.*mask;
[nt2,pt2] = max(h2');
ind = find(nt2>nt/10 & nt>=20);
if ~isempty(ind)
    for ii = 1:length(ind)
        t1 = diff(h(ind(ii),pt(ind(ii)):end));
        t2 = find(t1>0);
        if ~isempty(t2),
            setnow = h(ind(ii),pt(ind(ii))+t2(1):end);
            [nt2, pt2] = max(setnow);
            if nt2 > nt/10
                pt2 = pt2+pt(ind(ii))+t2(1)-1;
                pt(ind) = pt2;
            end
        end;
        %if min(h(ind(ii),pt(ind(ii)):pt2(ind(ii)))) < nt2(ind(ii))
           
        %end
    end
    mask = ones(size(h));
    for idx = 1:length(bins), mask(idx,1:pt(idx)) = NaN; end;
end;
%now use everything to find lows above the rough no PE peaks
ii = find(partialdatmerged(datbins2,4)>= partialdatmerged(datbins2,2).^c(1).*10^(c(2)+.5) & partialdatmerged(datbins2,2) < 1000);
h2 = histmulti5([partialdatmerged(datbins2(ii),5),partialdatmerged(datbins2(ii),2)],[bins_ssc; bins]');
h2 = h2.*mask;
h3 = zeros(size(h2));
%minv = max([repmat(10,length(bins),1)'; min(h2')*1.2])';
%v = repmat(minv,1, length(bins));
%h3(h2<v) = 1;
%h3(h2>=v) = 0;
h3 = diff(h2')';
h3(h3<=0) = NaN;
h3(h3>0) = 0;
h3(h2(:,1:end-1)<=3)=0;
[~,ph] = min(h3');
s = find(nt>=5);
figure(98), clf
loglog((partialdatmerged(datbins2,5)), (partialdatmerged(datbins2,2)), '.c')
hold on
loglog((partialdatmerged(datbins2(ii),5)), (partialdatmerged(datbins2(ii),2)), '.b')
plot(binc_ssc,binc(pt),'m^', 'markerfacecolor', 'm')
plot(binc_ssc(s), binc(ph(s)), 'sk', 'markerfacecolor', 'k', 'markersize', 16)

%yint = interp1(binc_ssc(s), mean([binc(ph(s)); binc(ph(s)-1)]), binc, 'pchip');
yint = interp1(binc_ssc(s), binc(ph(s)), binc_ssc, 'pchip');
fitobject = fit(log10(binc_ssc)', log10(yint*1)', 'smoothingspline', 'SmoothingParam', 0.99);
%fitobject = fit(log10(bins(s))', log10(bins(ph(s)))', 'smoothingspline', 'SmoothingParam', 0.99);
y = feval(fitobject, log10(binc_ssc));
plot(binc_ssc, 10.^y, 'r-')

a = find(partialdatmerged(datbins2,2) > 10.^(feval(fitobject,log10(partialdatmerged(datbins2,5)))));
datbins2pe = datbins2(a);
loglog((partialdatmerged(datbins2pe,5)), (partialdatmerged(datbins2pe,2)), 'or')

ii = find(partialdatmerged(datbins2,4)< partialdatmerged(datbins2,2).^c(1).*10^(c(2)+.5) & partialdatmerged(datbins2,2)>1);
b = .01:.01:1;
hc = histc(partialdatmerged(datbins2(ii),2)./partialdatmerged(datbins2(ii),5),b);
[~,p] = max(hc(2:end)); p = p + 1;
figure(1), clf, bar(b,hc), hold on, plot(b(p),hc(p), '*r'), disp(b(p))
ii2 = find(partialdatmerged(datbins2,2)./partialdatmerged(datbins2,5) <b(p)/100);
figure(2), clf
loglog((partialdatmerged(datbins2(ii),5)), (partialdatmerged(datbins2(ii),2)), '.b')
hold on
loglog((partialdatmerged(datbins2(ii2),5)), (partialdatmerged(datbins2(ii2),2)), '.r')
figure(3), clf
loglog(partialdatmerged(datbins2,4), partialdatmerged(datbins2,2)./partialdatmerged(datbins2,5), '.')
hold on
loglog(partialdatmerged(datbins2(ii),4), partialdatmerged(datbins2(ii),2)./partialdatmerged(datbins2(ii),5), '.')

pause
                
                
 %%               
%                 iitemp = find(partialdatmerged(datbins2,2)<200 & partialdatmerged(datbins2,5)<1e5);
%                 [pecutoff, f1] = mode(partialdatmerged(datbins2(iitemp),2));
%                 if pecutoff == 1
%                     iitemp2 = find(partialdatmerged(datbins2,2)>1 & partialdatmerged(datbins2,2)<200 & partialdatmerged(datbins2,5)<1e5);
%                     [pecutoff2, f2] = mode(partialdatmerged(datbins2(iitemp2),2));
%                 end
%                 if pecutoff > 1 | f2 > f1/50
%                     iitemp = find(partialdatmerged(datbins2,2)>pecutoff & partialdatmerged(datbins2,2)<800 & partialdatmerged(datbins2,5)<1e4);
%                     bins = logspace(1,log10(800));
%                     bins = bins(bins>pecutoff);
%                     h = histc(partialdatmerged(datbins2(iitemp),2),bins);
%                     [m,p] = min(h(1:end-1));
%                     madd = 0;
%                     b = find(bins<250);
%                     [maxh,maxhp] = max(h(b));
%                     minh = min(h(maxhp:end-1));
%                     %if maxh ~= minh
%                     if maxh/minh > 2 & minh<250, m = minh*1.2; else m = m*1.2; end %>2
%                     if maxh/minh > 1.5 & minh<250, m = minh*1.1; else m = m*1.2; end %>2
%                     %low = find(h(maxhp:end-1) <= m + max([1 madd]))';
%                     low = find(h(maxhp:end-1) <= m )';
%                     if isempty(low),
%                         low = find(h(1:end-1) <= m )';
%                     else
%                         %if maxh == minh
%                         [maxh_test,maxhp_test] = max(h(maxhp-1+low(1):maxhp-1+low(end)));
%                         if maxh_test > maxh/2 %case for secondard baseline peak (e.g., baseline jump within hr)
%                             maxh = maxh_test;
%                             maxhp = maxhp-1+low(1)+maxhp_test-1;
%                         end
%                         %end
%                         bins = bins(maxhp:end);
%                         low = find(h(maxhp:end-1) <= m )';
%                     end
%                     if low(1) == 1
%                         pecutoff = 1;
%                     else
%                         pecutoff = bins(low(1));
%                     end
%                     %end
%                 end
%                 %keyboard
%                 %disp(pecutoff)
%                 iitemp2 = find(partialdatmerged(datbins2,2)>=100 & partialdatmerged(datbins2,5)<1e4);
%                 [ssccutoff, f1] = mode(partialdatmerged(datbins2(iitemp2),5));
%                 iitemp2 = find(partialdatmerged(datbins2,2)>=100 & partialdatmerged(datbins2,5)>ssccutoff/2 & partialdatmerged(datbins2,5)<ssccutoff*2);
%                 bins = logspace(2,log10(800));
%                 h = histc(partialdatmerged(datbins2(iitemp2),2),bins);
%                 [m,p] = min(h(1:end-1));
%                 if bins(p) > pecutoff
%                     pecutoff = bins(p);
%                 end;
%                 %disp(pecutoff)
%                 %keyboard
%%              

%                %consider PE/CHL and PE/SSC criteria to handle really bright euks when they have some signal on PE channel
%                %synSSCmin = 100;
%                %a = find(partialdatmerged(datbins2,4).^.8./partialdatmerged(datbins2,2) < 1); %rough cut
%                a = find(partialdatmerged(datbins2,4).^.8./partialdatmerged(datbins2,2) < 2); %rough cut
%                datbins2pe = datbins2(a);
%                clear a
%                a = find(partialdatmerged(datbins2pe,2) > pecutoff & (partialdatmerged(datbins2pe,5)./partialdatmerged(datbins2pe,2) < SSC2PE_cutoff/4 )); %| partialdatmerged(datbins2pe,5) < 2e3));                    
%                datbins2pe = datbins2pe(a);
%                clear a
                
                datbins2nope = setdiff(datbins2, datbins2pe);
                classpe = (ones(size(datbins2pe))*1)';  %syn = class 1
                classnope = (ones(size(datbins2nope))*4)'; %euks = class 4
                if ~isempty(datbins2pe)
                    temp = partialdatmerged(datbins2pe,[4:5]);  %chl and ssc
                    %tempind = find(temp(:,2) < 5e3);  %crude SSC screening to cut out junk   
                    %tempind = find(temp(:,2) < 5e3 & temp(:,2) > 2e3 & temp(:,1) > 200);  %crude SSC screening to cut out junk   
                    tempind = find(temp(:,2) < 5e3 & temp(:,2) > 500 & temp(:,1) > 200);  %crude SSC screening to cut out junk   
                    maxvalue = 1e6;  
                    bins = 10.^(0:log10(maxvalue)/63:log10(maxvalue));  %make 64 log spaced bins
                    [nmergedhist,x,nbins] = histmulti5(temp(tempind,:),[bins' bins']);
                    [y,ind] = max(nmergedhist(:));
                    [i,j] = ind2sub(size(nmergedhist),ind);
                    tempmode = [i,j];
                    tempmode = bins(tempmode);
                    tempind = find(temp(:,2) <= tempmode(2)*5 & temp(:,2) >= tempmode(2)/5 & temp(:,1) < tempmode(1)*5 & temp(:,1) > tempmode(1)/5);% & temp(:,2) > synSSCmin); %May 2015, add hard SSC lower limit, deal with noise bursts in April 2013
                    if size(tempind,1) < 2, %special crude case for very few points
                        tempind = find(temp(:,2) < 5e3);  %crude SSC screening to cut out junk   
                        tempmode = mean(temp(tempind,:),1); 
                        tempind = find(temp(:,2) <= tempmode(2)*5 & temp(:,2) >= tempmode(2)/5 & temp(:,1) < tempmode(1)*5 & temp(:,1) > tempmode(1)/5);
                    end;
                    if length(tempind) > 1, %~isempty(tempind) %chl vs ssc
                        pedist = mahal(log10(temp),log10(temp(tempind,:)));  %distance of each point from cluster
                        threshhold = pedist_thre; %changed from 8 to 5 for 2006 and 2007, April 2007
                        junkind = find(pedist > threshhold);  %don't take all really small stuff *** -8?
                        classpe(junkind) = 3; %reassign class to PE junk
                        clear temp
                        %now conservatively only keep Syn that are tight on PE for now
                        temp = partialdatmerged(datbins2pe,[2,5]);
                        pedist = mahal(log10(temp),log10(temp(classpe==1,:)));  %distance of each point from cluster
                        %threshhold = 1; %changed from 8 to 5 for 2006 and 2007, April 2007
                        junkind = find(pedist > threshhold);  %don't take all really small stuff
                        classpe(junkind) = 3; %reassign class to PE junk
                    end;
               end;
               %figure(99), clf, loglog(partialdatmerged(datbins2pe,5), partialdatmerged(datbins2pe,2), 'r.'), hold on, loglog(partialdatmerged(datbins2pe(classpe==1),5), partialdatmerged(datbins2pe(classpe==1),2), 'g.')
               %keyboard
                %now have first cut syn, find pe vs. chl slope
%                 temp = partialdatmerged(datbins2pe(classpe==1),[2,4]);
%                 lfit = polyfit(log10(temp(:,2)), log10(temp(:,1)),1);
%                 tempcoeff = 10.^(lfit(2)-.9); temppower = lfit(1); %0.9 with log10 shift for cutoff line for euks with pe
%                 iitemp = find(partialdatmerged(datbins2nope,5)< 5e4);
%                 highpe_baseline_index = length(find(partialdatmerged(datbins2nope(iitemp),2)>pecutoff))./length(iitemp); clear iitemp
%                 %disp(highpe_baseline_index)
%                 if mode(partialdatmerged(datbins2pe(classpe == 1),2)) > pecutoff*50 & (highpe_baseline_index > .3 | pecutoff >1) ;%& pecutoff > 1 %syn mode far above baseline and baseline is noisy
%                      SynPE_lower = prctile(partialdatmerged(datbins2pe(classpe == 1),2),.2); %.15
%                  else
%                      SynPE_lower = pecutoff; %
%                 end
                %repeat steps above with better euk pe cutoff
                %consider PE/CHL and PE/SSC criteria to handle really bright euks when they have some signal on PE channel
                %THIS ONE 3.24.16: a = find((partialdatmerged(datbins2,4).^temppower./partialdatmerged(datbins2,2) < 1/tempcoeff & partialdatmerged(datbins2,2) > pecutoff) | (partialdatmerged(datbins2,4) < tempmode(1)*2 & partialdatmerged(datbins2,2) > mode(temp(:,1))/1000));  
                %clear temp
                %%a = find(partialdatmerged(datbins2,4).^temppower./partialdatmerged(datbins2,2) < 1/tempcoeff & partialdatmerged(datbins2,5) > 1e3); %chl > 200 or 300? May 2015, add chl limit to avoid junk with tiny chl
                %datbins2pe = datbins2(a);
                clear a
                %keyboard
                %disp([pecutoff SynPE_lower/6 mode(temp(:,1))/1000])
                %%CONSIDER WHETHER CAN GO TO SynPE_lower/5 might help early 2009
%%                a = find(partialdatmerged(datbins2,2) > max([pecutoff SynPE_lower/5]) & (partialdatmerged(datbins2,5)./partialdatmerged(datbins2,2) < SSC2PE_cutoff )); %| partialdatmerged(datbins2pe,5) < 2e3));                    
%                a = find(partialdatmerged(datbins2,2) > pecutoff & (partialdatmerged(datbins2,5)./partialdatmerged(datbins2,2) < SSC2PE_cutoff )); %| partialdatmerged(datbins2pe,5) < 2e3));                    
%                datbins2pe = datbins2(a);
%                clear a
%                datbins2nope = setdiff(datbins2, datbins2pe);
%                classpe = (ones(size(datbins2pe))*1)';  %syn = class 1
%                classnope = (ones(size(datbins2nope))*4)'; %euks = class 4
%                temp = partialdatmerged(datbins2pe,[4:5]);  %chl and ssc
%               tempind = find(temp(:,2) < 5e3 & temp(:,2) > 500 & temp(:,1) > 200);  %crude SSC screening to cut out junk   
%                maxvalue = 1e6;  
%                bins = 10.^(0:log10(maxvalue)/63:log10(maxvalue));  %make 64 log spaced bins
%                [nmergedhist,x,nbins] = histmulti5(temp(tempind,:),[bins' bins']);
%                [y,ind] = max(nmergedhist(:));
%                [i,j] = ind2sub(size(nmergedhist),ind);
%                tempmode = [i,j];
%                tempmode = bins(tempmode);
%                tempind = find(temp(:,2) <= tempmode(2)*5 & temp(:,2) >= tempmode(2)/4 & temp(:,1) < tempmode(1)*4 & temp(:,1) > tempmode(1)/5);% & partialdatmerged(datbins2pe,2)>1.5*pecutoff);% & temp(:,2)) > synSSCmin); %May 2015, add hard SSC lower limit; Oct 2015 change to *4 (from *5)

                if length(tempind) < 2, %special crude case for very few points (e.g., first hr in FCB1_2009_099_213424)
                    tempind = find(temp(:,2) < 5e3);  %crude SSC screening to cut out junk   
                    tempmode = mean(temp(tempind,:),1);
                    tempind = find(temp(:,2) <= tempmode(2)*4 & temp(:,2) >= tempmode(2)/5 & temp(:,1) < tempmode(1)*4 & temp(:,1) > tempmode(1)/5);
                end;                
                
                if length(tempind) > 1, %~isempty(tempind), %chl vs ssc
                    pedist = mahal(log10(temp),log10(temp(tempind,:)));  %distance of each point from cluster
                    threshhold = pedist_thre - 6; %changed from 8 to 5 for 2006 and 2007, April 2007
                    junkind = find(pedist > threshhold);  %don't take all really small stuff
                    classpe(junkind) = 3; %reassign class to PE junk
                    clear temp
                    %end repeat

                    %now try to find regular cryptos (i.e., the ones bigger and brighter than syn)
                    tind = find(classpe == 1);  %SYN
                    % new crypto scheme 5-10-03, must be larger than syn mean on SSC and have high PE/CHL (beads are rel. low on PE/CHL, Heidi
                    tempmean = mean(partialdatmerged(datbins2pe(tind),5));  %mean ssc of syn
                    tempmean2 = mean(partialdatmerged(datbins2pe(tind),2));  %mean pe of syn
                    tind = find(partialdatmerged(datbins2pe(junkind),5) > tempmean*1.5 & partialdatmerged(datbins2pe(junkind),2) > tempmean2*.1 & partialdatmerged(datbins2pe(junkind),2)./partialdatmerged(datbins2pe(junkind),4) > 5e-2); 
                    classpe(junkind(tind)) = 6; %reassign class to "dim" cryptos
                    %next two lines added for "lg cryptos", Heidi 6/2/03
                    tind = find(partialdatmerged(datbins2pe(junkind),2) > 5e4); %PE above cutoff
                    classpe(junkind(tind)) = 2; %reassign class to "bright" cryptos
                    clear tind tempmean junkind tempmedian pedist threshhold chlcutoff tempind               

                    %now consider cluster of syn points on PE vs. SSC and add back any within threshold
                    temp = partialdatmerged(datbins2pe,[2,5]);  %pe and ssc
                    %classpe(temp(:,2) < synSSCmin) = 3; %reassign class to PE junk
                    if ~isempty(find(temp(:,2)<1)), keyboard, end
                    tempind = find(classpe == 1);
                    if length(tempind) > 1%~isempty(tempind) %May 2015 fixed indexing mistakes here
                        pedist = mahal(log10(temp), log10(temp(tempind,:)));  
                        %May 2015 increase threshold applied to pe vs ssc, final step, but make sure syn SSC not too much higher than initial mean
                        %%synind = find(pedist < pedist_thre + 6 & temp(:,2) > synSSCmin & temp(:,2) < mean(partialdatmerged(datbins2pe(tempind),2))*4);
                        %synind = find(pedist < pedist_thre + 6 & temp(:,2) < mean(partialdatmerged(datbins2pe(tempind),2))*4);
                        [nmergedhist,x,nbins] = histmulti5(temp(tempind,:),[bins' bins']);
                        [y,ind] = max(nmergedhist(:));
                        [i,j] = ind2sub(size(nmergedhist),ind);
                        tempmode = [i,j];
                        tempmode = bins(tempmode);
                        ssc_cut = prctile(temp(classpe==1,2),95);
%                         synind = find(pedist < pedist_thre + 10 & (temp(:,1) <= tempmode(1) | temp(:,2) <= tempmode(2)));
%                         classpe(synind) = 1;
%                         synind = find(pedist < pedist_thre + 6 & temp(:,1) > tempmode(1) & temp(:,2) > tempmode(2));
%                         classpe(synind) = 1;
%                         junkind = find(pedist > pedist_thre + 6 & temp(:,1) > tempmode(1) & temp(:,2) > tempmode(2));
%                         classpe(junkind) = 3;
                        synind = find(pedist < pedist_thre + 8); %12
                        classpe(synind) = 1;
                        junkind = find(classpe == 1 & pedist > pedist_thre + 8 & temp(:,2) > ssc_cut);
                        classpe(junkind) = 3;
                        ssc_cut = prctile(temp(classpe==1,2),98);
                        %junkind = find(classpe == 1 & pedist > pedist_thre + 10 & temp(:,2) > ssc_cut); %12
                        junkind = find(classpe == 1 & temp(:,2) > ssc_cut); %12
                        classpe(junkind) = 3;
                         %ssc_synlow = prctile(temp(classpe==1,2),1);
                         ssc_synlow = min(temp(classpe==1,2));
                         ssc_synhigh = prctile(temp(classpe==1,2),98);
                         pe_synhigh =  prctile(temp(classpe==1,1),98);
                         pe_synlow = prctile(temp(classpe==1,1),1);
                         synind = find(temp(:,1) > pe_synlow & temp(:,1) < pe_synhigh & temp(:,2) > ssc_synlow/10 & temp(:,2) < ssc_synlow);
                         junkind = find(temp(:,1) > pe_synhigh & temp(:,2) > ssc_synlow/10 & temp(:,2) < ssc_synlow/2);
                         %junkind = find((temp(:,1) > pe_synhigh | temp(:,1) < pe_synlow) & temp(:,2) > ssc_synlow/10 & temp(:,2) < ssc_synlow/2);
                         %[ length(junkind) length(synind)/4]
                         %keyboard
                         if length(junkind) < length(synind)/10 % | length(junkind)<40
                             classpe(synind) = 1; 
                         end; %otherwise too much noise to extend syn into very low ssc (for same pe)
                         synind = find(temp(:,2) > ssc_synlow & temp(:,2) < ssc_synhigh & temp(:,1) < pe_synhigh *10);
                         classpe(synind) = 1;
                    end
                    clear pedist threshhold temp tempmode

                    %now do cryptos one more time, after final syn
                    tind = find(classpe == 1);  %SYN
                    junkind = find(classpe == 3);
                    % new crypto scheme 5-10-03, must be larger than syn mean on SSC and have high PE/CHL (beads are rel. low on PE/CHL, Heidi
                    tempmean = mean(partialdatmerged(datbins2pe(tind),5));  %mean ssc of syn
                    tempmean2 = mean(partialdatmerged(datbins2pe(tind),2));  %mean pe of syn
                    tind = find(partialdatmerged(datbins2pe(junkind),5) > tempmean*1.5 & partialdatmerged(datbins2pe(junkind),2) > tempmean2*.1 & partialdatmerged(datbins2pe(junkind),2)./partialdatmerged(datbins2pe(junkind),4) > 5e-2); 
                    classpe(junkind(tind)) = 6; %reassign class to "dim" cryptos
                    %next two lines added for "lg cryptos", Heidi 6/2/03
                    tind = find(partialdatmerged(datbins2pe(junkind),2) > 5e4); %PE above cutoff
                    classpe(junkind(tind)) = 2; %reassign class to "bright" cryptos
                    clear tind tempmean junkind tempmedian pedist threshhold chlcutoff tempind               

                end;
                
                %synSSCmode = mode(temp(classpe == 1,2)); %Nov 2015, try forcing small euk SSC mode to be greater or near to syn SSC mode
                temp = partialdatmerged(datbins2pe,[4,5]);  %chl and ssc
                %synCHLmode = mode(temp(classpe == 1,1)); %Dec 2015, try forcing small euk CHL to be not far below syn CHL mod
                synSSCpop_min = min(temp(classpe == 1,2)); %prctile(temp(classpe == 1,2),.5); %
                synCHLpop_min = prctile(temp(classpe == 1,1),5); %min(temp(classpe == 1,1)); %
                synSSCpop_upper = prctile(temp(classpe == 1,2),75); 
                clear temp
                %now do junk elimination for euks        
                %one special case set in late 2005 with oddly stretched out SSC signals
                if ismember(filename(1:7), {'oc2105a' 'oc2405a' 'oc2905a' 'no1505a'})
                    chljunk_power = .6; 
                end
                coeff = chljunk_coeff; %change from .04 to .01, 7-8-03, back to .04, 11/6/03, back to .01 for lab VolCal, 4-13-05; 4-18-06 for dock change to 0.1 (from .05)
                    
                power = chljunk_power; %change from 1 to 1.1 from MVCO_May2003, 5-11-03 Heidi; change from 1.1 to .8 for lab VolCal, 4-10-05 heidi          
                %Mar 2007 - add chl threshold to eliminate new noise on chl baseline
                temp = partialdatmerged(datbins2nope,[4:5]);  %chl and ssc
                if ~isempty(temp)
                    % tempind = find(temp(:,1) >  coeff.*temp(:,2).^power & temp(:,1) < 1e4);  %line on chl v. ssc for junk cutoff, plus only consider chl < 1e4 for smallest euk peak
                    %tempind = find(temp(:,1) >  coeff.*temp(:,2).^power & temp(:,1) > 300 & temp(:,1) < 4e3);  %line on chl v. ssc for junk cutoff, plus only consider 100 < chl < 3e3 for smallest euk peak
                    %tempind = find(temp(:,1) >  coeff.*temp(:,2).^power & temp(:,1) > 500 & temp(:,1) < 4e3 & temp(:,2) < 2e4 & temp(:,2) > 1e3);  %line on chl v. ssc for junk cutoff, plus only consider 100 < chl < 3e3 for smallest euk peak; May 2015 also add SSC < 2e4 for smallest euks; Oct 2015 add SSC > 1e3 for smallest euks
                    %tempind = find(temp(:,1) >  coeff.*temp(:,2).^power & temp(:,1) > 1000 & temp(:,1) < 4e3 & temp(:,2) < 2e4 & temp(:,2) > synSSCmode);  %line on chl v. ssc for junk cutoff, plus only consider 100 < chl < 3e3 for smallest euk peak; May 2015 also add SSC < 2e4 for smallest euks; Oct 2015 add SSC > 1e3 for smallest euks
                    %tempind = find(temp(:,1) >  coeff.*temp(:,2).^power & temp(:,1) > max([2000 synCHLmode/2]) & temp(:,1) <  4e3 & temp(:,2) < 2e4 & temp(:,2) > synSSCmode); 
                    tempind = find(temp(:,1) >  coeff.*temp(:,2).^power & temp(:,1) > min([1000 synCHLpop_min]) & temp(:,1) < 4e3 & temp(:,2) < 2e4 & temp(:,2) > synSSCpop_upper); 
                    maxvalue = 1e6;  %is this too high?
                    %bins = 10.^(0:log10(maxvalue)/63:log10(maxvalue));  %make 1024 log spaced bins
                    bins = 10.^(0:log10(maxvalue)/31:log10(maxvalue));  %make 1024 log spaced bins
                    [nmergedhist,x,nbins] = histmulti5(temp(tempind,:),[bins' bins']);
                    [y,ind] = max(nmergedhist(:));
                    [i,j] = ind2sub(size(nmergedhist),ind);
                    tempmode = [i,j];
                    tempmode = bins(tempmode);
                    %tempind = find(temp(:,2) <= tempmode(2)*5 & temp(:,2) >= tempmode(2)/5 & temp(:,1) < tempmode(1)*5 & temp(:,1) > tempmode(1)/5);
                    tempind = find(temp(:,2) <= tempmode(2)*3 & temp(:,2) >= tempmode(2)/3 & temp(:,1) < tempmode(1)*3 & temp(:,1) > tempmode(1)/3); %change to tighten up euk cluster (handling too much debris taken in Nov 2013, etc.)
                    if length(tempind) > 1, %~isempty(tempind),
                        coeff = 10.^(log10(tempmode(1))-power*log10(tempmode(2))-0.9); %disp(coeff)%Oct 2015, -0.9 (from -0.5)                        
                        nopedist = mahal(log10(temp),log10(temp(tempind,:)));  %distance of each point from cluster
                        %              threshhold = 8;     %2/25/03 switched back to 10 (from 5) for syn_lab, 6/2/03 changed from 4 to 5; 7/8/03 changed from 5 to 10
                        threshhold = 12; %change from 6 to 5, Jan 2015 trying to address too much debris in euk cluster at some times in late 2013
                        %junkind = find((nopedist > threshhold  &  ((temp(:,2) < tempmode(2)/2 & temp(:,1) < tempmode(1))  |  (temp(:,1) < tempmode(1)/2) | temp(:,1) <  coeff.*temp(:,2).^power)) | temp(:,2) <  synSSCpop_min ); % | temp(:,1) < synCHLpop_min/2);  %don't take all really small stuff
                        junkind = find((nopedist > threshhold  &  ((temp(:,2) < synSSCpop_min)  |  (temp(:,1) < tempmode(1)/2 & temp(:,2) > tempmode(2)) | (temp(:,1) < tempmode(1) & temp(:,2) < tempmode(2)) | temp(:,1) <  coeff.*temp(:,2).^power))); % | temp(:,1) < synCHLpop_min/2);  %don't take all really small stuff
                    else
                        junkind = find(temp(:,1) <  coeff.*temp(:,2).^power);
                    end;
                    classnope(junkind) = 5; %reassign class to euk junk
                    %tempind = find(temp(2) > tempmode(2) & temp(1) > tempmode(1) & tempind = find(temp(:,1) >  coeff.*temp(:,2).^power,
                end;
       
            bd2cell_time = cellresults(sectionnum,1) - beadmatch(sectionnum,1);
            modeflag = 0;
            if bd2cell_time < 2/24 &  bd2cell_time > 0, %if cell section within 2 hours of a bead run
                modeflag = 1;
                %this is to eliminate large beads
                tdata = partialdatmerged(datbins2nope,2:5); tbd = beadmatch(sectionnum,2:5); %PE FLS CHL SSC
                %temp = find(tdata(:,1)>tbd(1)/2 & tdata(:,1)<tbd(1)*6 & tdata(:,2)>tbd(2)/2 & tdata(:,2)<tbd(2)*6 & tdata(:,3)>tbd(3)/2 & tdata(:,3)<tbd(3)*6 &tdata(:,4)>tbd(4)/2 & tdata(:,4)<tbd(4)*6);
                temp = find(tdata(:,1)>tbd(1)*.7 & tdata(:,1)<tbd(1)*1.3 & tdata(:,3)>tbd(3)*.7 & tdata(:,3)<tbd(3)*2 &tdata(:,4)>tbd(4)*.7 & tdata(:,4)<tbd(4)*1.3);
                classnope(temp) = 5; %disp(length(temp))
                tdata = partialdatmerged(datbins2pe,2:5); tbd = beadmatch(sectionnum,2:5);
                %temp = find(tdata(:,1)>tbd(1)/2 & tdata(:,1)<tbd(1)*6 & tdata(:,2)>tbd(2)/2 & tdata(:,2)<tbd(2)*6 & tdata(:,3)>tbd(3)/2 & tdata(:,3)<tbd(3)*6 &tdata(:,4)>tbd(4)/2 & tdata(:,4)<tbd(4)*6);
                temp = find(tdata(:,1)>tbd(1)*.7 & tdata(:,1)<tbd(1)*1.3 & tdata(:,3)>tbd(3)*.7 & tdata(:,3)<tbd(3)*1.3 &tdata(:,4)>tbd(4)*.7 & tdata(:,4)<tbd(4)*1.3);
                classpe(temp) = 5; %disp(length(temp))
                
                %now use the same steps to eliminate small beads (after analysis in beadbatch5_field that now includes their stats)
                %tdata = partialdatmerged(datbins2nope,2:5); tbd = beadmatch(sectionnum,10:13); %PE FLS CHL SSC 
                %temp = find(tdata(:,1)>tbd(1)/2 & tdata(:,1)<tbd(1)*6 & tdata(:,3)>1e4 &tdata(:,4)>tbd(4)/2 & tdata(:,4)<tbd(4)*6);
                %disp(length(temp))
                %classnope(temp) = 5; %disp(length(temp))
                tdata = partialdatmerged(datbins2pe,2:5); tbd = beadmatch(sectionnum,10:13);
                temp = find(tdata(:,1)>tbd(1)*.7 & tdata(:,1)<tbd(1)*1.3 & tdata(:,3)>tbd(3)*.7 & tdata(:,3)<tbd(3)*1.3 & tdata(:,4)>tbd(4)*.7 & tdata(:,4)<tbd(4)*1.3);
                classpe(temp) = 5; %disp(length(temp))
                %disp(length(temp))
            end;
                clear temp %power coeff
                %if intersect(filename(1:4), ['my15'; 'my16'; 'my17'; 'my18'], 'rows')
                %    keyboard
                %    temp = find(partialdatmerged(datbins2nope,4) < 100);
                %    classnope(temp) = 5;                
                %end;
                if plotflag %& ~mod(sectionnum,10)%, %mod(sectionnum,6) == 1,   %make surf plots
                    maxvalue = 1e6;  %is this too high?
                    bins = 10.^(0:log10(maxvalue)/127:log10(maxvalue));  %make 256 log spaced bins
                    maxvalueSSC = 1e7;  %is this too high?
                    binsSSC = 10.^(0:log10(maxvalueSSC)/127:log10(maxvalueSSC));  %make 256 log spaced bins
                    figure(1)
                    clf
                    subplot(221)
                    [n,x] = histmulti5(partialdatmerged(datbins,[3:4]), [bins' bins']);
                    ind = find(n == 0); n(ind) = NaN;
                    surf(x(:,1), x(:,2), log10(n)')
                    ylabel('CHL'), xlabel('FLS')
                    set(gca, 'xscale', 'log', 'yscale', 'log')
                    shading flat
                    axis([1 1e6 1 1e6])
                    view(2)
                    title(datestr(cellresults(sectionnum,1)))
                    subplot(222)        
                    [n,x] = histmulti5([partialdatmerged(datbins,5), partialdatmerged(datbins,2)], [binsSSC' bins']);
                    ind = find(n == 0); n(ind) = NaN;
                    surf(x(:,1), x(:,2), log10(n)')
                    ylabel('PE'), xlabel('SSC')
                    set(gca, 'xscale', 'log', 'yscale', 'log')
                    shading flat
                    view(2)
                    axis([1 1e7 1 1e6])
                    subplot(223)
                    [n,x] = histmulti5([partialdatmerged(datbins,5), partialdatmerged(datbins,4)], [binsSSC' bins']);
                    ind = find(n == 0); n(ind) = NaN;
                    surf(x(:,1), x(:,2), log10(n)')
                    ylabel('CHL'), xlabel('SSC')
                    set(gca, 'xscale', 'log', 'yscale', 'log')
                    shading flat
                    view(2)
                    axis([1 1e7 1 1e6])
                    subplot(224)
                    %        [n,x] = histmulti5([partialdatmerged(datbins2pe,5), partialdatmerged(datbins2pe,4)], [bins' bins']);
                    [n,x] = histmulti5([partialdatmerged(datbins2pe,4), partialdatmerged(datbins2pe,2)], [bins' bins']);
                    ind = find(n == 0); n(ind) = NaN;
                    surf(x(:,1), x(:,2), log10(n)')
                    hold on
                    loglog([1:10:1e5], 10.^(log10([1:10:1e5])*fit1(1) + fit1(2)), 'r')
                    %ylabel('CHL'), xlabel('SSC')
                    ylabel('PE'), xlabel('CHL')
                    set(gca, 'xscale', 'log', 'yscale', 'log')
                    shading flat
                    view(2)
                    axis([1 1e6 1 1e6])
                    %        title('PE containing cells only')
                    clear maxvalue bins
                end;
                
                colorstr = ['r', 'b', 'k', 'g', 'y', 'c', 'm'];
                
                %        numcluster = numclusterpe + numclusternope;
                %        numcluster = 5;
                
                mergedwithclass = [partialdatmerged NaN*ones(size(partialdatmerged,1),1)];
                mergedwithclass(datbins2pe,end) = classpe;
                %    mergedwithclass(datbins2nope,end) = classnope + numclusterpe;  %euks class nums start after pe cells
                mergedwithclass(datbins2nope,end) = classnope;  
                mergedwithclass = mergedwithclass(datbins,:);
                clear classpe classnope
                if plotflag %& ~mod(sectionnum,2), %~mod(sectionnum+2,4), %mod(sectionnum,6) == 1, %make cluster plots
                    figure(2)
                    clf, 
                    subplot(221)
                    hold on
                    ylabel('CHL'), xlabel('FLS')
                    for c = 1:numcluster,
                        ind = find(mergedwithclass(:,end) == c);
                        eval(['loglog(mergedwithclass(ind,3),mergedwithclass(ind,4), ''' colorstr(c) 'o'', ''markersize'', 1)'])
                    end;
                    set(gca, 'xscale', 'log', 'yscale', 'log')
                    axis([1 1e6 1 1e6])
                    title(datestr(cellresults(sectionnum,1)))
                    %X = 1:1000:1e6; plot(X,X*5+50000, 'k-') 
                    %line([1 5000], [5e4 5e4], 'color', 'k'), line([5000 5000], [5e4 1e6], 'color', 'k')
                    subplot(222)        
                    hold on
                    ylabel('PE'), xlabel('SSC')
                    if modeflag
                        plot(beadmatch(sectionnum,5), beadmatch(sectionnum,2), 'ob', 'markersize', 15)
                        plot(beadmatch(sectionnum,13), beadmatch(sectionnum,10), 'ob', 'markersize', 15)
                    end
                    for c = 1:numcluster,
                        ind = find(mergedwithclass(:,end) == c);
                        eval(['loglog(mergedwithclass(ind,5),mergedwithclass(ind,2), ''' colorstr(c) 'o'', ''markersize'', 1)'])
                    end;
                    set(gca, 'xscale', 'log', 'yscale', 'log')  
                    fplot(['x/' num2str(SSC2PE_cutoff)], [1 1e6], 'linestyle', '--')
                    axis([1 1e7 1 1e6])
                    subplot(223)
                    hold on
                    ylabel('CHL'), xlabel('SSC')
                    for c = 1:numcluster,
                        ind = find(mergedwithclass(:,end) == c);
                        eval(['loglog(mergedwithclass(ind,5),mergedwithclass(ind,4), ''' colorstr(c) 'o'', ''markersize'', 1)'])
                    end;
                    c = 1;%overlay syn again on top
                    ind = find(mergedwithclass(:,end) == c);
                    eval(['loglog(mergedwithclass(ind,5),mergedwithclass(ind,4), ''' colorstr(c) 'o'', ''markersize'', 1)'])
                    set(gca, 'xscale', 'log', 'yscale', 'log')
                    set(gca, 'Ygrid', 'on')
                    axis([1 1e7 1 1e6])
                    plot(tempmode(2), tempmode(1), '^k', 'markerfacecolor', 'k', 'markersize',8)
                    fplot([num2str(coeff) '*x.^' num2str(power)], [10 1e6], 'linestyle', '--')
                    axis([1 1e7 1 1e6])
                    subplot(224)
                    hold on
                    ylabel('PE'), xlabel('CHL')
                    for c = 1:numcluster,
                        ind = find(mergedwithclass(:,end) == c);
                        eval(['loglog(mergedwithclass(ind,4),mergedwithclass(ind,2), ''' colorstr(c) 'o'', ''markersize'', 1)'])
                    end;
                    loglog([1:10:1e5], 10.^(log10([1:10:1e5])*fit1(1) + fit1(2)), 'r')
                    set(gca, 'xscale', 'log', 'yscale', 'log')
 %                   fplot([num2str(tempcoeff) '*x.^' num2str(temppower)], [10 1e6], 'linestyle', '--')
                    axis([1 1e6 1 1e6])
                    disp('pause for graphs...')
                    pause %(0.1)
                    disp('reading next...')
                end;  %if 1, (to plot)
                
                maxvalue = 1e6;  %is this too high?
                bins = 10.^(0:log10(maxvalue)/255:log10(maxvalue));  %make 256 log spaced bins
                maxvalueSSC = 1e7;  %is this too high?
                binsSSC = 10.^(0:log10(maxvalueSSC)/255:log10(maxvalueSSC));  %make 256 log spaced bins                
                for c = 1:numcluster,
                    ind = find(mergedwithclass(:,end) == c);
                    cellNUM(sectionnum,c) = length(ind);
                    if length(ind) > 1,
                        n = hist(mergedwithclass(ind,2:4), bins);
                        [junk, maxind] = max(n);
                        n2 = hist(mergedwithclass(ind,5), binsSSC);
                        [junk, maxind2] = max(n2);
                        cellPEmode(sectionnum,c) = bins(maxind(1));  
                        cellFLSmode(sectionnum,c) = bins(maxind(2));  
                        cellCHLmode(sectionnum,c) = bins(maxind(3));  
                        cellSSCmode(sectionnum,c) = binsSSC(maxind2); 
                        cellPE(sectionnum,c) = mean(mergedwithclass(ind,2));  %mean params
                        cellFLS(sectionnum,c) = mean(mergedwithclass(ind,3));  %mean params
                        cellCHL(sectionnum,c) = mean(mergedwithclass(ind,4));  %mean params
                        cellSSC(sectionnum,c) = mean(mergedwithclass(ind,5));  %mean params
                        clear n junk maxind
                    else
                        cellPE(sectionnum,c) = NaN;  
                        cellFLS(sectionnum,c) = NaN;  
                        cellCHL(sectionnum,c) = NaN;  
                        cellSSC(sectionnum,c) = NaN;
                        cellPEmode(sectionnum,c) = NaN;  
                        cellFLSmode(sectionnum,c) = NaN;  
                        cellCHLmode(sectionnum,c) = NaN;  
                        cellSSCmode(sectionnum,c) = NaN;
                    end;  %if ~isempty(ind)
                end; %for c = 1:numcluster
                clear c ind maxvalue bins 
                
                allmergedwithclass{sectionnum} = mergedwithclass;
            else
                allmergedwithclass{sectionnum} = NaN;
                cellPE(sectionnum,1:numcluster) = NaN;
                cellFLS(sectionnum,1:numcluster) = NaN;
                cellCHL(sectionnum,1:numcluster) = NaN;
                cellSSC(sectionnum,1:numcluster) = NaN;
                cellPEmode(sectionnum,1:numcluster) = NaN;
                cellFLSmode(sectionnum,1:numcluster) = NaN;
                cellCHLmode(sectionnum,1:numcluster) = NaN;
                cellSSCmode(sectionnum,1:numcluster) = NaN;
                cellNUM(sectionnum,1:numcluster) = NaN;
                cellresults(sectionnum,1:3) = NaN;
                beadmatch(sectionnum,1:13) = NaN;
            end; %if ~isempty(goodtimebins)
            %get ready for next loop           
            partialdatmerged = datmerged;  %reset partialdat with file partly completed
        end; %sectionnum = 1:length(timesectionendbin)
        mergedwithclass = allmergedwithclass;

        eval(['save ' groupedpath filetypelist(typenum,:) '_' num2str(sectcount) ' beadmatch* cell* classnotes'])
        eval(['save ' mergedpath filetypelist(typenum,:) 'merged_' num2str(sectcount) ' merged*'])
        clear beadmatch cellresults mergedwithclass link* allmergedwithclass cellNUM cellPE cellFLS cellCHL cellSSC cell*mode
        clear datbins datbins2* goodtimebins fit1 fittitles sectionnum
    end; %for sectcount
end; %for typenum = 1:length(filetypelist)

%end; %for count = 1:2

clear count culture typenum numcluster dattitles datmerged partialdatmerged cell*mode
clear timeinterval time*ind timesectionendbin totaltime timetitles mergedtitles cellrestitles classnotes colorstr
clear beadresults beadtitles beadmatchtitles
clear filelist filenum filename ans