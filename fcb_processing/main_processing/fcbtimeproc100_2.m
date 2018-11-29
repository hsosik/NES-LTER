%fcbtimeproc = modified from fcbtimeproc to work for 12 channel setup in
%lab, now with 100 events per record (instead of 200 as before) 2/24/05
%3/7/03 fcbtimeproc = modified from cytosubproc2_labalt to handle new file format (with
%header), should work with any combination of port switching during analysis, Heidi
%3/1/03 modified from cytosubproc2 to handle files from SYN expt's in lab
%with two cultures analyzed (switch every half hour); Heidi

%In this script:
% Goal is to calculate time and volume analyzed for each record of each syringe:
% ----Go through and remove syringes that are being refilled, etc. during acquitions
% ----Calculate an accurate syringe pump speed in order to remove correct volume that hasn't been analyzed
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

%good records for cells or beads, skipping cases where syringe refills or port (valve) changes,
%also skipping any fast syringes (i.e., syr# = 0)
%non-zero port, same port at start and end, same syringe # at start and end (i.e. syringe did not refill during acquisition)

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
%hence this volume is not measured, and needs to be substracted. We are also not able to measure cells
% directly after a trigger (dead time), such that volume is being pushed through that
% we cannot analyze (and this should be subtracted from the final volume
% measurement). These numbers really matter when there is a large amount of
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
% ---- or the most robust metric: slope of records/time:

% Note that pump speed is changed manually - it is on a given speed, until
% the instrument is turned off and reset....

stepdist=syrpumpinfo(:,7)-syrpumpinfo(:,8);
timediff=totalendsec-totalstartsec;
pump_speed=stepdist./timediff; %rough speed

%% identify when syringe number advances:

if length(unique(syrpumpinfo(:,5))) > 1
    
    j0=find(diff(totalstartsec) > 3600); j0=j0+1; %include any time gap over an hour...
    j1=find(syrpumpinfo(:,5)~=syrpumpinfo(:,6));
    j2 =find(syrpumpinfo(2:end,5)~=syrpumpinfo(1:end-1,6)); j2=j2+1;
    
    syr_ind=unique([1;j0;j1;j2]);
    
    %legend for syrchangeinfo=[syr_starttime  syr_endtime  start_index end_index syrnum syrtype]
    syrchangeinfo=nan(length(syr_ind),6);
    
    for j=1:length(syr_ind)-1
        
        %finding when a syringe starts is easy:
        syrchangeinfo(j,:)=[totalstartsec(syr_ind(j)) NaN syr_ind(j) NaN syrpumpinfo(syr_ind(j),6) syrpumpinfo(syr_ind(j),4)];
        
        %finding when it ends...
        if ismember(syr_ind(j+1),j2)  %must look at next syringe change:
            syrchangeinfo(j,[2 4])=[totalendsec(syr_ind(j+1)-1) syr_ind(j+1)-1];
        else
            syrchangeinfo(j,[2 4])=[totalendsec(syr_ind(j+1)) syr_ind(j+1)];
        end
    end
    
    % add end entry:
    syrchangeinfo(end,:)=[totalstartsec(syr_ind(end)) totalendsec(end) syr_ind(end) size(syrpumpinfo,1) syrpumpinfo(syr_ind(end),6) syrpumpinfo(syr_ind(end),4)];
    
    
    %% identify when syringe number rolls over and for each 'batch' assign a speed, based on top syringe count:
    
    batch_ind=find(diff(syrchangeinfo(:,5)) < -1);
    batch_ind=[batch_ind; size(syrchangeinfo,1)]; %add on the end entry...
    %     if any(syrchangeinfo(batch_ind,6) ~=3) % only consider cell syringes
    %         keyboard
    %     end
    %
    pumpspeed=nan(size(syrpumpinfo,1),1);
    missing_ind=[];
    for q=1:length(batch_ind)
        if q==1
            ind1=1;
            ind2=syrchangeinfo(batch_ind(q),4);
        else
            ind1=syrchangeinfo(batch_ind(q-1)+1,3);
            ind2=syrchangeinfo(batch_ind(q),4); %assign indexes for easier handling
        end
        
        switch syrchangeinfo(batch_ind(q),5)
            case 50
                pumpspeed(ind1:ind2)=40;
            case 100
                pumpspeed(ind1:ind2)=80;
            case 200
                pumpspeed(ind1:ind2)=160;
            otherwise %record troublesome indexes
                missing_ind=[missing_ind; batch_ind(q) q];
        end
        
    end
    
    %% for any missing batches, use the next closest in time speed, as long as no time delays:
    [a, ~, c]=setxor(missing_ind(:,1),batch_ind);
    for q=1:length(missing_ind)
        %find the smallest, positive number or least negative number as ->
        %favoring back-in-time speeds
        
        test=missing_ind(q,1)-a;
        ii=find(test>0);
        if ~isempty(ii)
            ii=c(ii(end)); %min positive
        else %this is the first block then
            ii=c(1);
        end
        
        %assign indexes for easier handling:
        if missing_ind(q,2)==1, ind1=1; else ind1=syrchangeinfo(batch_ind(missing_ind(q,2)-1)+1,3); end
        ind2=syrchangeinfo(missing_ind(q,1),4);
        
        pumpspeed(ind1:ind2)=pumpspeed(syrchangeinfo(batch_ind(ii),4));
    end
    
    %
    if any(isnan(pumpspeed)) %should be all filled in....
        disp('Pump speeds unaccounted for!')
        keyboard
    end
    
    %% SCREEN SYRINGES FOR ABNORMALLY LONG ACQ TIMES
       
    % Okay, once have a pump speed, examine acquisition times for any long outliers: grouped by pump speed
    
    %find the median, almost max and almost min acq time of each syringe:
    rec=nan(length(syrchangeinfo),4);
    for q=1:length(syrchangeinfo)
        syrind=syrchangeinfo(q,3):syrchangeinfo(q,4);
        syrind=syrind(flag(syrind)==3); %only want the cell records
        if ~isempty(syrind)
            rec(q,1:3)=quantile(acqtime(syrind),[0.05 0.5 0.95]);
            rec(q,4)=unique(pumpspeed(syrind));
        end
    end
    
    for p=[40 80 160]
        %what is the median range within a set speed?
        speed_ind=find(rec(:,4)==p);
        baseline=quantile(rec(speed_ind,3)-rec(speed_ind,2),0.95);
        baseline=baseline+1.5; %tack on an extra 1.5 sec
        %now find syringes that exceed this:
        jj=find(rec(speed_ind,3)-rec(speed_ind,2) > baseline); %find syringes whose max - median is way outside some std dev
        %flag all points that belong to these syringes:
        for j=1:length(jj)
            syrind=syrchangeinfo(speed_ind(jj(j)),3):syrchangeinfo(speed_ind(jj(j)),4);
            flag(syrind)=60;
        end
        
    end
    
    %FIND TROUBLING SYRINGES BETWEEN WASH AND A 'BAD BATCH'
         %%   

    
    %Also change for first records after a syringe refill - acquisition times are longer than the rest of syringe - sheath leaking in?
    %for each syringe, check first and last record against rest of times:
    %%
    for q=1:length(syrchangeinfo)
        syrind=syrchangeinfo(q,3):syrchangeinfo(q,4);
        syrind=syrind(flag(syrind)==3); %only want the cell records
        if ~isempty(syrind)
            if acqtime(syrind(1)) > (max(acqtime(syrind(2:end-1)))+ 2) %if first record is over max by two sec, flag
                flag(syrind(1))=61;
            elseif acqtime(syrind(end)) > (max(acqtime(syrind(2:end-1)))+ 2) %if end record is over max by two sec, flag
                flag(syrind(end))=61;
            end
        end
    end
    
    
else %for rare case of only 1 syringe being processed...
    keyboard
    %OLD CODE HERE:
    P1=find((pump_speed > 70 & pump_speed < 85));
    P2=find((pump_speed > 150 & pump_speed < 165));
    P3=find((pump_speed > 30 & pump_speed < 45));
    
    pumprate=zeros(size(syrpumpinfo,1),1);
    pumprate(P3)=40;
    pumprate(P2)=80;
    pumprate(P1)=160;
end


%% if plotflag is on, plot a final summary for sanity check:
if syrplotflag
    
    figure(16)
    %syringe number and pump speed over time
    subplot(2,1,1,'replace')
    plot(totalstartsec,syrpumpinfo(:,5),'k.-'); %syringe number and time
    hold on
    h1=plot(totalstartsec,pumpspeed,'.','color',[0 0.5 1]); %syringe designated speed and time   
    h2=plot(totalstartsec(missing_ind),pumpspeed(missing_ind),'.','color',[1 0.5 0]); %syringe speed assigned from closest syringe
    legend([h1(1); h2(1)],'Assigned pump speed','Speeds from closest neighbor','location','NorthOutside')
    ylabel('Syringe Number')
    xlabel('Time (sec)')

    %% acquistion time vs. time
    subplot(2,1,2,'replace')
    plot(totalstartsec,acqtime,'.-') %acquistion time vs. time
    hold on
    ii=find(flag==3 | flag==60 | flag ==61);
    h3=plot(totalstartsec(ii),acqtime(ii),'.'); %only cell syringes
    
    jj=find(flag==60); %cell syringes excluded    
    if ~isempty(jj), h4=plot(totalstartsec(jj),acqtime(jj),'o'); end 
    
    mm=find(flag==61); %first or last record of a syringe excluded
    if ~isempty(mm), h5=plot(totalstartsec(mm),acqtime(mm),'o');end
    
    if ~isempty(jj) && ~isempty(mm)
        legend([h3(1); h4(1); h5(1)],'cell records','excl. syr','records excl.','location','NorthOutside')
    elseif ~isempty(jj) && isempty(mm)
        legend([h3(1); h4(1)],'cell records','excl. syr','location','NorthOutside')
    elseif isempty(jj) && ~isempty(mm)
        legend([h3(1); h5(1)],'cell records','records excl.','location','NorthOutside')
    else
        legend(h3(1),'cell syr.','location','NorthOutside')
    end
    ylabel('Acquisition time (sec)')
    xlabel('Time (sec)')
    
    disp('Type dbcont to move to next batch of processing....')
    keyboard
    
end %syrplotflag


%% and now for the magic:

offset = pumpspeed*querytime*(totalvol/maxpos); %offset = steps/sec*(q+d time)*(totalvol/maxpos)
analvol = analvol - offset;

%%
outmatrix = [100*(1:length(totalstartsec))' totalstartsec totalendsec acqtime analvol flag];  %2/24/05 heidi, 100 event records for 12 channels
clear acqtime flag t ind analvol maxpos start stop totalvol

