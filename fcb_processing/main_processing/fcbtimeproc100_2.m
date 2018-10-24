%fcbtimeproc = modified from fcbtimeproc to work for 12 channel setup in
%lab, now with 100 events per record (instead of 200 as before) 2/24/05
%3/7/03 fcbtimeproc = modified from cytosubproc2_labalt to handle new file format (with
%header), should work with any combination of port switching during analysis, Heidi
%3/1/03 modified from cytosubproc2 to handle files from SYN expt's in lab
%with two cultures analyzed (switch every half hour); Heidi

%In this script:
% Goal is to calculate time and volume analyzed for each record of each syringe:
% ----Go through and remove syringes that are being refilled, etc. during acquitions
% ----Calculate an accurate syringe pump speed in order to remove volume that hasn't been analyzed
% ----Remove syringes where acquisition time of triggers do not conform to standard
%       expection - indicating a possible clog or flow problem or sheath intrusion

%flag summary:
%99 - fast syringes
%98 - syr refill at start of acqusition
%97 - syr refill during acquistion
%60 - outlier record (based on acq time) right after refill or just before end of syr
%61 - acq times within syr do not resemble a normal distribution
%62 - cell syr with just 1 record

% initial info and setup
flag = zeros(size(totalstartsec));

querytime = 0.0908; %0.0930; %actually this is query+deadtime for 100 tirggers; querytime referrring to time it takes to query syringe position
acqtime = totalendsec - totalstartsec - querytime; %for ts-acq-st  %July 2016, from high noise bead runs
%acqtime = totalendsec - totalstartsec - .405;  %for tq-acq-qt?
%acqtime = totalendsec - totalstartsec - .1326; %for vnts-acq-stnv
%acqtime = totalendsec - totalstartsec; %for qt-acq-tq (mr07)
%%
start = syrpumpinfo(:,7);
stop = syrpumpinfo(:,8);
totalvol = .25;  %ml vol of syringe
maxpos = 48000; %48000 changed to 48200 for 2014 test, Heidi June 2016
analvol = stop*NaN;

t = find(start - stop >= 0);  %start > stop
analvol(t) = (start(t)-stop(t))/maxpos*totalvol;


%syrpumpinfo legend:
%1=temp,2=humidity,3=start port,4=end port,5=start syr#,6=end syr#,7=start syr pos,8=end syr pos.

%% CHECK FOR FLAGS:

%good records for cells or beads, skipping cases where syringe refills or
%port (valve) changes, also skipping any fast syringes (i.e., syr# = 0)
%non-zero port, same port at start and end, same syringe # at start and end (i.e. syringe did not refill during acquisition

flag(:) = 0;  %default all records to Not use
%ind = find((syrpumpinfo(:,3) == syrpumpinfo(:,4)) & syrpumpinfo(:,5) & syrpumpinfo(:,6) & (syrpumpinfo(:,5) == syrpumpinfo(:,6))));
ind = 1 + find((syrpumpinfo(2:end,3) == syrpumpinfo(2:end,4)) & syrpumpinfo(2:end,5) & syrpumpinfo(2:end,6) & (syrpumpinfo(2:end,5) == syrpumpinfo(2:end,6)) & (syrpumpinfo(1:end-1,6) == syrpumpinfo(2:end,5)));
%last part to skip cases where syringe starts to refill during preceeding data transmission (start syringe # ~= previous end syringe #)
%these records often seem to have flow problem or something that makes apparent volume too high
flag(ind) = syrpumpinfo(ind,3);  %set good records to syringe port # (3/03 6=culture 1, 3=culture2, ?=beads)
flag = flag(1:length(acqtime)); %eliminate zeros added past end

%follwing finds syringe refills starting mid-record
%ind = find((syrpumpinfo(:,3) == 3 & syrpumpinfo(:,4) == 3) & syrpumpinfo(:,5) & syrpumpinfo(:,6) & (syrpumpinfo(:,5) ~= syrpumpinfo(:,6)));
ind = find(syrpumpinfo(1:end-1,5) & syrpumpinfo(2:end,5) == 0);  %transitions to fast syringes
flag(ind) = 99;
% ind = find(acqtime < 0);
% flag(ind) = 97; %record acquired too fast near?; indicator for high noise conditions in 2014, unreliable analvol, etc.

t = find(syrpumpinfo(:,5) < syrpumpinfo(:,6)); %syringe refilled during acquisition
analvol(t) = (start(t) + maxpos-stop(t))/maxpos*totalvol;
flag(t) = 97; %add specific flag for these syringes KRHC, 7/13/16

t = find(syrpumpinfo(:,5) == syrpumpinfo(:,6) & start - stop <= 0); %syringe in middle of refilling at start; June 2015 <=0 (not < 0) for cases in 2014 with stopped syr at top or bottom during v. noisy fast acq
analvol(t) = (maxpos-stop(t))/maxpos*totalvol;
flag(t) = 98; %syringes refilling in middle

%last record before end of set (switch to new valve) when syringe is at end, sometimes these have really long acq times (renegade triggers?)
flag(find(syrpumpinfo(2:end,6) < syrpumpinfo(1:end-1,6) & syrpumpinfo(1:end-1,8) == 10)) = 98;
if totalstartsec(end)/3600/24 + year > datenum('8-20-2011') &&  totalstartsec(end)/3600/24 + year < datenum('2-11-2012'), %datenum('1-10-2014 17:50'),
    %if year == datenum('31-Dec-11') && totalstartsec(end)/3600/24 < datenum('4-01-12')-year, %case for 2012 bad records after cleaning with clorox and not enough rinsing??
    addgap_sec = 60*60; %minutes * 60 sec/min
    t = find(flag == 99); % & totalstartsec/3600/24 > datenum('8-19-11')-year);
    for count = 1:length(t),
        a = totalstartsec(t(count));
        b = a + addgap_sec;
        tt = find(totalstartsec > a & totalstartsec < b & flag == 3); %find syr 3 = seawater in add gap interval
        flag(tt) = 0; %set these records to 'not use'
    end;
end;


%% FIND VOLUME ANALYZED...turns out, not exactly straightforward...

%Notes:
%There are syringe movements (and hence volume) not yet accounted for in
%the measurements. The time is recorded, then syringe position queried. %
%Volume is being pushed through while we are querying for pump position and
%hence this volume is not accounted for. We are also not able to measure cells
% directly after a trigger (dead time), such that volume is being pushed through that
% we cannot analyze (and this should be subtracted from the final volume
% measurement). These numbers really matter when there is a crazy amount of
% noise in the system, such as years 2014 and 2017.

% From Alexi's measurements (7/12/16) of full noise, we see
%that the average dead time for 100 triggers is .0275s (27.5ms) (without query time)

% From Alexi's measurements, at a pump speed of 160 steps/sec,
% avg steps taken is 14.89 per 100 records with query time
% At 40 steps/sec, avg steps taken is 3.72.
% In these cases, since the triggers could be considered instantaneous, these steps/time are
% essentially due to deadtime and querytime, which leaves:

%       14.89 steps = (160 steps/s) * (dead time + query time)
%       (dead+query) = 0.0931 sec

%       3.72 steps = (40 steps/s) * (dead time + query time)
%       (dead+query) = 0.0930 sec !!!!!

%So, the lost time in query (total - deadtime = 0.0930 - 0.0275) ends up being 0.0655...
%....but, we only really care about the total time lost (dead +query), so we just
%multiply pump speed by 0.0930 sec (total) to get steps lost, which gives volume lost

%older comments, but perhaps useful?
%analvol = analvol - 1.692e-4;  %volume offset for query, for mr07 set (3/03)
%analvol = analvol - 5.48e-5;  %older, corresponds to P2 speed for two syringe queries? July 2016

% FIND PUMP SPEED

%3 different speeds of pump:
%P3 - 40 steps/sec -> 20min per syr -> syrnum rolls over at 50
%P2 - 80 steps/sec -> 10min per syr -> syrnum rolls over at 100
%P1 - 160 steps/sec -> 5min per syr -> syrnum rolls over at 200

% can find pump speed by a few different methods:
% ---- looking at steps per second
% ---- avg syring time
% ---- or the most robust metric - slope of records/time:

stepdist=syrpumpinfo(:,7)-syrpumpinfo(:,8);
timediff=totalendsec-totalstartsec;
pump_speed=stepdist./timediff; %rough speed

%We will use the slope of records over time. To do this, we need to find
%the individual syringes and when the number rolls over:

%legned for syrchangeinfo=[syr_starttime  syr_endtime  start_index end_index syrnum avgsyrtime syrnum_slope]

% a slightly clunkier way to find syringe changes, but perhaps more straight forward:
if length(unique(syrpumpinfo(:,5))) > 1
    count=1;
    ii=1;
    syrnum=syrpumpinfo(1,5);
    syrchangeinfo=[totalstartsec(1) NaN 1 NaN syrpumpinfo(1,5) syrpumpinfo(1,3)];
    while ii < size(syrpumpinfo,1)
        ii=ii+1;
        if syrnum ~= syrpumpinfo(ii,5) %moved onto new syringe
            syrchangeinfo(count,2)=totalendsec(ii-1); %fill in previous record
            syrchangeinfo(count,4)=ii-1; %record ending index
            count=count+1; %advance syr count
            syrchangeinfo(count,1)=totalstartsec(ii); %start new syringe record
            syrnum=syrpumpinfo(ii,5);
            syrchangeinfo(count,5)=syrnum;
            syrchangeinfo(count,3)=ii; %record beginning index
            syrchangeinfo(count,6)=syrpumpinfo(ii,3);
        elseif syrnum ~= syrpumpinfo(ii,6) %moved onto new syringe      %this logic needs to come 2nd after 1st piece
            syrchangeinfo(count,2)=totalstartsec(ii); %fill in record
            syrchangeinfo(count,4)=ii; %record ending index
            syrnum=syrpumpinfo(ii,6);
            count=count+1;  %advance syr count
            syrchangeinfo(count,1)=totalendsec(ii);
            syrchangeinfo(count,5)=syrnum;
            syrchangeinfo(count,3)=ii; %record beginning index
            syrchangeinfo(count,6)=syrpumpinfo(ii,4);  %type of syringe
        end
    end
    %to end this matrix:
    syrchangeinfo(count,2)=totalendsec(end);
    syrchangeinfo(count,4)=size(syrpumpinfo,1);
    
    %calculate average syringe time:
    syrchangeinfo=[syrchangeinfo (1/60)*(syrchangeinfo(:,2)-syrchangeinfo(:,1))];
    
    %find the slope of syringe numbers over time:
    celli=find(syrchangeinfo(:,6)==3); %use only the cell syringes for more reliable slopes
    ii=find(diff(syrchangeinfo(celli,5))<0);
    ro1=celli(ii); %rollover indexes
    ro2=celli(ii+1); %start of new cell syringe
    
    ro=[1; ro1; ro2; size(syrchangeinfo,1)];
    ro=sort(ro);
    
    % occassionally, FCB is down and then resumes in the middle of a syringe.
    % This causes issues for finding the slope. Below is a check for these 'split syringes'
    % use 3hr gap ~ 10000 sec as cutoff to split the syringe -> case where acq stopped, but resumes syr num?
    ind2add=[];
    for q=1:2:length(ro)-1
        if any(diff(syrchangeinfo(ro(q):ro(q+1),1)) > 10000) %3hr gap ~ 10000 sec - split the syringe -> case where acq stopped, but resumes syr num
            figure, plot(totalstartsec,syrpumpinfo(:,5),'.-')
            zz=find(diff(syrchangeinfo(ro(q):ro(q+1),1)) > 10000);
            line([syrchangeinfo(ro(q)+zz-1) syrchangeinfo(ro(q)+zz-1)],ylim,'color','r')
            %keyboard
            disp(['Discontinuous syringe?' num2str(q)])
            jj=ro(q):ro(q+1); %easier to handle indexes
            kk=find(diff(syrchangeinfo(jj,1)) > 10000);
            ind2add=[ind2add; jj(kk)'; jj(kk+1)']; %add this split into the rollover indexes
        end
    end
    
    ro=unique(sort([ro; ind2add])); %unique prevents duplicated indices at beginning or end of data batch
    
    if length(ro)>2
        syrslope=(syrchangeinfo(ro(2:2:end),2)-syrchangeinfo(ro(1:2:end-1)+1,1))./(syrchangeinfo(ro(2:2:end),5)-syrchangeinfo(ro(1:2:end-1)+1,5)); %two indexes to avoid some 0 syringes on restart
    else %case where no rollover detected
        syrslope=(syrchangeinfo(ro(end),2)-syrchangeinfo(ro(1)+1,1))./(syrchangeinfo(ro(end),5)-syrchangeinfo(ro(1)+1,5)); %two indexes to avoid some 0 syringes on restart
    end
    
    %% fill in these slopes for all syringes:
    for q=1:2:length(ro)-1
        
        if length(ro(q):ro(q+1)) < 4 & length(ro)>2
            syrchangeinfo(ro(q)+1:ro(q+1),8)=0; %will be caught later - > unreliable slope estimates if low syrnum and many 0's and 1' for 'little' restarts
            
        else %Typical case:
            syrchangeinfo(ro(q):ro(q+1),8)=syrslope((q+1)/2);
        end
    end
    
    %use the syrchangeinfo matrix to populate a matrix for each record:
    avgsyrtime=nan(size(syrpumpinfo,1),2); %syringe times and slopes
    for q=1:size(syrchangeinfo,1)
        avgsyrtime(syrchangeinfo(q,3):syrchangeinfo(q,4),1)=(1/60)*(syrchangeinfo(q,2)-syrchangeinfo(q,1));
        avgsyrtime(syrchangeinfo(q,3):syrchangeinfo(q,4),2)=syrchangeinfo(q,8);
    end
    
    %assign actually syringe speed based on slope:
    P3=find(avgsyrtime(:,2) > 1000 & avgsyrtime(:,2) < 1800); %slope = ~1250 syringes/time ~ rollover
    P2=find(avgsyrtime(:,2) > 500 & avgsyrtime(:,2) < 700); %slope = ~615 syringes/time ~ rollover
    P1=find(avgsyrtime(:,2) > 250 & avgsyrtime(:,2) < 400); %slope = ~314 syringes/time ~ rollover
    
    %% any cells or bead runs unaccounted for?
    % If cannot find a speed based on the above intervals, use the closest in time sped:
    ind1=find(flag==3 | flag ==6); %cells and beads
    test0=ismember(ind1,[P1;P2;P3]);
    goodrate=find(test0==1); goodrate=ind1(goodrate); %index still into syrpumpinfo
    test=find(test0==0);
    
    if ~isempty(test)
        jj=ind1(test); %easier to handle indexes
        disp(['bead rate and other rates unaccounted for: ' num2str(length(test)) ' out of ' num2str(length(avgsyrtime)) ' ...using closest speed'])
        
        for i=1:length(jj)
            qq=find(goodrate <= jj(i)); %for looking back in time....
            [~, im]=min(totalstartsec(jj(i))-totalstartsec(goodrate(qq)));  %find closest good bead or cell run in (backwards) time that has a pumprate
            
            if ismember(goodrate(qq(im)),P3) %use this index speed
                P3=[P3; jj(i)]; %add the test index to a speed index group
            elseif ismember(goodrate(qq(im)),P2)
                P2=[P2; jj(i)];
            elseif ismember(goodrate(qq(im)),P1)
                P1=[P1; jj(i)];
            else %then can try a forward time look:
                [~, im]=min(abs(totalstartsec(jj(i))-totalstartsec(goodrate)));
                if ismember(goodrate(im),P3) %use this index speed
                    P3=[P3; jj(i)]; %add the test index to a speed index group
                elseif ismember(goodrate(im),P2)
                    P2=[P2; jj(i)];
                elseif ismember(goodrate(im),P1)
                    P1=[P1; jj(i)];
                else
                    disp('Cannot find a speed for these measurements! Help!')
                    keyboard
                end
            end
        end
    end
    
    pumprate=zeros(size(syrpumpinfo,1),1);
    pumprate(P3)=40;
    pumprate(P2)=80;
    pumprate(P1)=160;
    
else %for rare case of only 1 syringe being processed...
    
    P1=find((pump_speed > 70 & pump_speed < 85));
    P2=find((pump_speed > 150 & pump_speed < 165));
    P3=find((pump_speed > 30 & pump_speed < 45));
    
    pumprate=zeros(size(syrpumpinfo,1),1);
    pumprate(P3)=40;
    pumprate(P2)=80;
    pumprate(P1)=160;
end

%% SCREEN SYRINGES FOR EXPECTED DISTRIBUTION OF ACQ TIMES

% We now would like to identify syrignes that may be partially clogged, have sheath
% intrusions or are not operating as expected for smooth, continuous flow.

%We can do this by looking at the acquisition time over time and see if these values
% follow a normal distribution (they should), and if not, flag for removal

% To check for a normal distribution, for each syringe:
% ---- from Andy's stats notes, plot: y(1) < y(2) < y(3)  (where y's are the acquistion times of all the records within a syringe)
% ---- against F^-1(1/(n+1), F^-1(2/n+1), F^-1(3/n+1), where F^-1 is the norm inv parameterized from teh mean and variance of the syringe data
% ---- sum of squares between obs and expected is a good pre-check before testing for linearity
% ---- if ss is high, perform linear regression between ordered data and the normal inverse of data positions and check how good fit is

for q=1:length(syrchangeinfo) %for each syringe...
    
    %syrchangeinfo has the indexes per syringe at columns 3 and 4:
    inds=syrchangeinfo(q,3):syrchangeinfo(q,4);
    tempacq=acqtime(inds); %dataslices
    tempflag=flag(inds);
    tf=find(tempflag==3); %only evaluate cell records (ask of these, should any records be removed?)
    n=length(tempacq(tf));
    
    if ~isempty(tf) && n~=1 %can't be empty and length can't be 1
        
        finv = (norminv(((1:n)/(n+1)),mean(tempacq(tf)),std(tempacq(tf))))'; %expected value given an underling normal distribution
        [ordered_acq, oind] = sort(tempacq(tf)); %ordered acquisition times
        ss=sum(sqrt((ordered_acq-finv).^2)); %difference between expected and observed
 
        %Hmmm...maybe we should just start with linearity and then see if removing outliers fixes it?
        
        [b,~,~,~,stats]=regress(finv,[ones(length(ordered_acq),1) ordered_acq]);
        b0=b; stats0=stats;
        counter=0;
        while stats(1) < 0.875 && counter < 3
            counter=counter+1; %disp(num2str(counter))
            if counter ~= 3
                %disp('outlier problem?')
                finv0 = (norminv(((1:(n-counter))/((n-counter)+1)),mean(ordered_acq(1:end-counter)),std(ordered_acq(1:end-counter))))'; %expected value given an underling normal distribution
                [b,~,~,~,stats]=regress(finv0,[ones(length(ordered_acq(1:end-counter)),1) ordered_acq(1:end-counter)]);
                disp([num2str(q) ': linear regression count: ' num2str(counter) ' : ' num2str(stats(1)) ' : ' num2str(stats(3))])
            end
        end
        
        switch counter
            case 1 %exclusion of one record was enough!
                disp([num2str(q) ': Excluded a starting or ending outlier record from syringe'])
                flag(inds(tf(oind(end))))=60; %flag it!
            case 2 %two records excluded
                disp([num2str(q) ': Excluded two starting or ending outlier records from syringe'])
                flag(inds(tf(oind(end))))=60; %flag it!
                flag(inds(tf(oind(end-1))))=60; %flag it!
            case 3
                disp([num2str(q) ': trying to exclude outliers did not fix data distribution, flagging syringe'])
                flag(inds(tf))=61;
        end
        
        if counter ~=0
            subplot(1,2,1,'replace')
            plot(finv,ordered_acq,'o'), hold on
            ff=find(flag(inds(tf(oind)))==60);
            plot(finv(ff),ordered_acq(ff),'ro','markerface','r') %identify outliers           
            line([ordered_acq(1) ordered_acq(end)],[ordered_acq(1) ordered_acq(end)]) %expectation
            title(['q: ' num2str(q) 'sum of squares:' num2str(ss)])
            plot(finv, b(1)+b(2)*finv,'r.-')
            
            subplot(1,2,2), hold on
            line([totalstartsec(inds(1)) totalstartsec(inds(1))],ylim,'linestyle','-','color',[0.6 0.6 0.6],'linewidth',2)
            line([totalstartsec(inds(end)) totalstartsec(inds(end))],ylim,'linestyle','-','color',[0.6 0.6 0.6],'linewidth',2)            
            plot(totalstartsec(max(1,inds(1)-500):min(inds(end),inds(end)+500)),acqtime(max(1,inds(1)-500):min(inds(end),inds(end)+500)),'.-','color',[0 0 0]), hold on
            xlim([totalstartsec(inds(1))-500 totalstartsec(inds(end))+500])
            plot(totalstartsec(inds(tf(oind(ff)))),acqtime(inds(tf(oind(ff)))),'ro','markerface','r') %identify outliers 
            
            keyboard
        end
         
    elseif n==1 %if you are a cell syringe of length 1 - ignore
        flag(inds(tf))=62; %length of 1
    end
end

%% if plotflag is on, plot a final summary for sanity check:
if syrplotflag
    
    figure(16)
    subplot(3,1,1,'replace')
    plot(totalstartsec,syrpumpinfo(:,5),'k.-'); %syringe number and time
    hold on
    if ~isempty(test)
        h1=plot(totalstartsec(ind1(test)),syrpumpinfo(ind1(test),5),'r.');
        legend(h1(1),'assigned pump rate from closest syr','location','NorthOutside')
    end
    
    ylabel('Syringe Number')
    xlabel('Time (sec)')
    
    subplot(3,1,2,'replace'), hold on
    h1=plot(totalstartsec,pump_speed,'.-','color',[0.6 0.6 0.6]); %time and speed; rough calc of speed for all syr
    h2=plot(totalstartsec(ind1),pump_speed(ind1),'.','color',[0 0.5 1]); %rough calc of speed for cell/bead syr
    h3=plot(totalstartsec(ind1),pumprate(ind1),'.','color',[0 0 1]); %assgined rate of cell/bead syr
    legend([h1(1); h2(1); h3(1)],'calc rough speed','cell syr rough speed','assigned speed for cell syr','location','NorthOutside','Orientation','Horizontal')
    if ~isempty(test)
        h4(1)=plot(totalstartsec(ind1(test)),pumprate(ind1(test)),'c.');
        legend([h1(1); h2(1); h3(1); h4(1)],'calc rough speed','cell syr rough speed','assigned speed for cell syr','assigned pump rate from closest syr','location','NorthOutside','Orientation','Horizontal')
    end
    ylabel('Calc and assigned pump rate (steps/sec)')
    xlabel('Time (sec)')
    ylim([-10 180])
    xlabel('Time (sec)')
    
    subplot(3,1,3,'replace')
    plot(totalstartsec,acqtime,'.-') %acquistion time vs. time
    hold on
    ii=find(flag==3);
    h1=plot(totalstartsec(ii),acqtime(ii),'.');
    jj=find(flag==60 | flag==61 | flag ==62);
    h2=plot(totalstartsec(jj),acqtime(jj),'ko'); %cell syringes excluded
    legend([h1(1); h2(1)],'cell syr.','excl. syr','location','NorthOutside')
    ylabel('Acquisition time (sec)')
    xlabel('Time (sec)')
    
    pause(0.5)
    
end %syrplotflag


%% and now for the magic:

offset = pumprate*querytime*(totalvol/maxpos); %offset = steps/sec*(q+d time)*(totalvol/maxpos)
analvol = analvol - offset;

%%
outmatrix = [100*(1:length(totalstartsec))' totalstartsec totalendsec acqtime analvol flag];  %2/24/05 heidi, 100 event records for 12 channels
clear acqtime flag t ind analvol maxpos start stop totalvol

