% add special case only for sequence of files in 2024 with 3 events per record instead of ususal 100
% must be used with fcbreadraw2b
%
%fcbtimeproc = modified from fcbtimeproc to work for 12 channel setup in
%lab, now with 100 events per record (instead of 200 as before) 2/24/05
%3/7/03 fcbtimeproc = modified from cytosubproc2_labalt to handle new file format (with
%header), should work with any combination of port switching during analysis, Heidi
%3/1/03 modified from cytosubproc2 to handle files from SYN expt's in lab
%with two cultures analyzed (switch every half hour); Heidi


%In this script:
% Goal is to calculate time and volume analyzed for each record of each syringe:

% ----Go through and remove syringes that are being refilled, etc. during acquitions
% ----Find syringe pump speed in order to remove correct volume that hasn't been analyzed
% ----Remove syringes where acquisition time of triggers are abnormal - indicating a possible clog or flow problem or sheath intrusion

%flag summary:
%99 - fast syringes
%98 - syr refill at start of acqusition
%97 - syr refill during acquistion
%60 - long acq times within a syr - discard whole syringe
%61 - outlier record (based on acq time) right after refill or just before end of syr, single records discarded
%62 - secondary catch on troublesome syringes near those with flag of 60

% initial info and setup
flag = zeros(size(totalstartsec));

if year2do <= 2004 %just for 2003 & 2004
    querytime=0.0908-0.0275; % time for queries on 100 trigger events - high noise investigation, July 2016
    deadtime=(0.0275/100)*200; %events per record pre-2005 have 200 events in them, adjust deadtime for 200 events
    acqtime = totalendsec - totalstartsec - (querytime+deadtime); %for ts-acq-st  %new query time, July 2016, from high noise bead runs
else
    querytime = 0.0908; %0.0930; %actually this is query+deadtime for 100 tirggers; querytime referring to time it takes to query syringe position
    acqtime = totalendsec - totalstartsec - querytime; %for ts-acq-st  %July 2016, from high noise bead runs
end
if year2do == 2024
    yd_temp = totalstartsec/3600/24;
    ind = yd_temp > 256 & yd_temp < 283.61; 
    querytime=0.0908-0.0275; % time for queries on 100 trigger events - high noise investigation, July 2016
    deadtime=(0.0275/100)*3; %events per record was 3 instead of 100 for a while in 2024
    acqtime(ind) = totalendsec(ind) - totalstartsec(ind) - (querytime+deadtime); %for ts-acq-st  %new query time, July 2016, from high noise bead runs    
    events = 100*ones(size(totalstartsec));
    events(ind) = 3;
end
events = cumsum(events);
%acqtime = totalendsec - totalstartsec - .405;  %for tq-acq-qt?
%acqtime = totalendsec - totalstartsec - .1326; %for vnts-acq-stnv
%acqtime = totalendsec - totalstartsec; %for qt-acq-tq (mr07)
%
start = syrpumpinfo(:,7);
stop = syrpumpinfo(:,8);
totalvol = .25;  %ml vol of syringe
maxpos = 48000; %48000 changed to 48200 for 2014 test, Heidi June 2016
analvol = stop*NaN;

t = find(start - stop >= 0);  %start > stop
analvol(t) = (start(t)-stop(t))/maxpos*totalvol;

%syrpumpinfo legend:
%1=temp,2=humidity,3=start port,4=end port,5=start syr#,6=end syr#,7=start syr pos,8=end syr pos.

% CHECK FOR FLAGS:

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

t = find(syrpumpinfo(:,5) < syrpumpinfo(:,6)); %syringe changes over; refilled during acquisition
analvol(t) = (start(t) + maxpos-stop(t))/maxpos*totalvol;
flag(t) = 97; %add specific flag for these syringes KRHC, 7/13/16
%syringes where finishing one syringe and moving onto another

t = find(syrpumpinfo(:,5) == syrpumpinfo(:,6) & start - stop <= 0); %syringe in middle of refilling at start; June 2015 <=0 (not < 0) for cases in 2014 with stopped syr at top or bottom during v. noisy fast acq
analvol(t) = (maxpos-stop(t))/maxpos*totalvol;
flag(t) = 98; %syringes refilling in middle
%syringes where number has ticked over, but still refilling -> has not reached totalvol yet

%last record before end of set (switch to new valve) when syringe is at end, sometimes these have really long acq times (renegade triggers?)
flag(find(syrpumpinfo(2:end,6) < syrpumpinfo(1:end-1,6) & syrpumpinfo(1:end-1,8) == 10)) = 98;

%Old code, but will keep - the below code catches a large portion of this:
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
% this needs to be substracted (as we didn't measure it). We are also not able to measure cells
% directly after a trigger (dead time), such that volume is being pushed through that
% we cannot measure as well. These numbers really matter when there is a large amount of
% noise in the system, such as years 2014 and 2017.

% From Alexi's measurements (7/12/16) of full noise, we see
%that the average dead time for 100 triggers is .0275s (27.5ms) (without query time)

% Also from Alexi's measurements, at a pump speed of 160 steps/sec,
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
% ---- slope of records/time
% ---- looking at syringe rollover numbers

% Note that pump speed is changed manually - it is on a given speed, until
% the instrument is turned off and reset...

% identify when syringe number advances:
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


% identify when syringe number rolls over and for each 'batch' assign a speed, based on top syringe count:
batch_ind=find(diff(syrchangeinfo(:,5)) < -1);
gap_syr=[];
if ~isempty(j0)
    for qq=1:length(j0)
        gg=find(syrchangeinfo(:,3) >= j0(qq));
        gap_syr=[gap_syr; gg(1)-1];
    end
end%find(diff(syrchangeinfo(:,5)) == 0); %this misses the gaps that do not duplicate syringes!
batch_ind=[batch_ind; gap_syr]; %account for gaps - treat them as batches
batch_ind=unique([batch_ind; size(syrchangeinfo,1)],'sorted'); %add on the end entry...

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
            if ~ismember(batch_ind(q),gap_syr) && batch_ind(q)~=batch_ind(end) %for the rare occassion where a gap interrupts at this syr num!
                pumpspeed(ind1:ind2)=40;
            else
                missing_ind=[missing_ind; batch_ind(q) q];
            end
        case 100
            if ~ismember(batch_ind(q),gap_syr) && batch_ind(q)~=batch_ind(end)
                pumpspeed(ind1:ind2)=80;
            else
                missing_ind=[missing_ind; batch_ind(q) q];
            end
        otherwise
            if syrchangeinfo(batch_ind(q),5) > 100 %for the highest speed, can identify if syringe count goes over 100
                pumpspeed(ind1:ind2)=160;
            else
                %record troublesome indexes
                missing_ind=[missing_ind; batch_ind(q) q];
            end
    end
    
end


%% for any missing batches, use the next closest in time speed, as long as no time delays:
plot_missing_ind=[];
if ~isempty(missing_ind)
    [a, ~, c]=setxor(missing_ind(:,1),batch_ind);
    if ~isempty(a) %meaning that missing_ind and batch_ind are the same - 1 rollove event, likely just a few syringes
        for q=1:size(missing_ind,1)
            %find the smallest, positive number or least negative number as ->
            %favoring back-in-time speeds
            test=missing_ind(q,1)-a; %a is batch_ind that have speeds assigned
            ii=find(test>0);
            if length(ii)==1 %closest batch_ind is first roll over..use start of record
                ii=c(1);
                batch_avgtime=(totalstartsec(syrchangeinfo(batch_ind(ii),3))-totalstartsec(syrchangeinfo(1,3)))./ (syrchangeinfo(batch_ind(ii),5) - syrchangeinfo(1,5));
            elseif isempty(ii) %missing_ind is the first roll over, use next batch
                ii=c(1);
                batch_avgtime=(totalstartsec(syrchangeinfo(batch_ind(ii),3))-totalstartsec(syrchangeinfo(missing_ind(q,1)+3,3)))./ (syrchangeinfo(batch_ind(ii),5) - syrchangeinfo(missing_ind(q,1)+3,5));
            elseif ~isempty(ii)
                ii=c(ii(end)); %min positive
                if ismember(batch_ind(missing_ind(q,2)-1),gap_syr) && missing_ind(q,1)~=batch_ind(end)
                    %gap right before this syr, and not the ending syr, use closest in time:
                    [~, ii]=min(abs(totalstartsec(syrchangeinfo(missing_ind(q,1),3))-totalstartsec(syrchangeinfo(a,3))));
                    ii=c(ii);
                    batch_avgtime=(totalstartsec(syrchangeinfo(batch_ind(ii),3))-totalstartsec(syrchangeinfo(batch_ind(ii-1)+3,3)))./ (syrchangeinfo(batch_ind(ii),5) - syrchangeinfo(batch_ind(ii-1)+3,5));
                else
                    batch_avgtime=(totalstartsec(syrchangeinfo(batch_ind(ii),3))-totalstartsec(syrchangeinfo(batch_ind(ii-1)+3,3)))./ (syrchangeinfo(batch_ind(ii),5) - syrchangeinfo(batch_ind(ii-1)+3,5));
                end
            end
            
            %assign indexes for easier handling: into pumpspeed matrix
            if missing_ind(q,2)==1, ind1=1; else  ind1=syrchangeinfo(batch_ind(missing_ind(q,2)-1)+1,3); end
            ind2=syrchangeinfo(missing_ind(q,1),4);
            
            if syrchangeinfo(missing_ind(q,1),5) <= 3 || (syrpumpinfo(ind2,6)-syrpumpinfo(ind1,5) <= 3) || isnan(batch_avgtime) %happens, if a stopped run very close to rollover number
                %any(missing_ind(q,1)-gap_syr < 5 & missing_ind(q,1)-gap_syr >0)
                %these are records ending at a low syringe number or
                %records ending after a gap and not many syringes have accrued
                %just use closest in time - will not be able to reliably compare avg speeds
                
                pumpspeed(ind1:ind2)=pumpspeed(syrchangeinfo(batch_ind(ii),4));
                
            else %calculate an average missing speed to use:
                %to do an average time, go in a few syringes - avoids restarts of 0's and 1's
                if missing_ind(q,2)==1
                    timeind2=syrchangeinfo(missing_ind(q,1),3);
                    timeind1=syrchangeinfo(1+3,3);
                    syrnum1=syrchangeinfo(1+3,5);
                    syrnum2=syrchangeinfo(missing_ind(q,1),5);
                else
                    timeind2=syrchangeinfo(missing_ind(q,1),3);
                    timeind1=syrchangeinfo(batch_ind(missing_ind(q,2)-1)+3,3);
                    syrnum2=syrchangeinfo(missing_ind(q,1),5);
                    syrnum1=syrchangeinfo(batch_ind(missing_ind(q,2)-1)+3,5);
                end
                
                missing_avgtime=((totalstartsec(timeind2)-totalstartsec(timeind1))./(syrnum2-syrnum1));
                
                %evaluate this average speed:
                %time for rollover divided by number of syringes:
                
                if round(missing_avgtime ./ batch_avgtime) == 2
                    pumpspeed(ind1:ind2)=0.5*pumpspeed(syrchangeinfo(batch_ind(ii),4)); %likely half the speed
                elseif round(missing_avgtime ./ batch_avgtime) == 4
                    pumpspeed(ind1:ind2)=0.25*pumpspeed(syrchangeinfo(batch_ind(ii),4)); %likely quarter the speed
                elseif round(batch_avgtime ./ missing_avgtime) == 2
                    pumpspeed(ind1:ind2)=2*pumpspeed(syrchangeinfo(batch_ind(ii),4)); %likely double the speed
                elseif round(batch_avgtime ./ missing_avgtime) == 4
                    pumpspeed(ind1:ind2)=4*pumpspeed(syrchangeinfo(batch_ind(ii),4)); %likely quadruple the speed
                elseif round(missing_avgtime ./ batch_avgtime) == 1  %for small restarts, where average slope is unreliable
                    pumpspeed(ind1:ind2)=pumpspeed(syrchangeinfo(batch_ind(ii),4));
                else
                    disp('cannot find an appropriate matching speed for this missing index!')
                    keyboard
                end
            end
            %record for plotting later:
            plot_missing_ind=[plot_missing_ind ind1:ind2];
        end
    else %there is only one rollover event...
        missing_avgtime=((totalstartsec(end)-totalstartsec(1))./(syrchangeinfo(end,5)-syrchangeinfo(1,5))); %secs/syringe
        if missing_avgtime > 800 && missing_avgtime < 1400 %probably slow speed
            pumpspeed(1:syrchangeinfo(batch_ind(end),4))=40;
        elseif missing_avgtime > 400 && missing_avgtime <= 800 %probably medium speed
            pumpspeed(1:syrchangeinfo(batch_ind(end),4))=80;
        elseif missing_avgtime > 100 && missing_avgtime <= 400 %probably fast speed
            pumpspeed(1:syrchangeinfo(batch_ind(end),4))=160;
        else
            disp('Too few number of syringes to get a pumpspeed estimate - enter in manually please')
            disp('such as: pumpspeed(1:syrchangeinfo(batch_ind(end),4))= 40, 80 or 160...')
            keyboard
            %pumpspeed(1:syrchangeinfo(batch_ind(end),4))= 40, 80 or 160 ....
        end
    end
end

%Secondary check - are all pumpspeeds accounted for?
if any(isnan(pumpspeed)) %should be all filled in....
    disp('Pump speeds unaccounted for!')
    keyboard
end


%% SCREEN SYRINGES FOR ABNORMALLY LONG ACQ TIMES

% Okay, now, within a pump speed, examine acquisition times for range outliers:

%calculate the median, and almost max and almost min acq time of each syringe:
%The difference between almost max and median seems to do a good job of finding abnormal syringes:
%global look:
rec=nan(length(syrchangeinfo),7);
for q=1:size(syrchangeinfo,1)
    syrind0=(syrchangeinfo(q,3)):(syrchangeinfo(q,4)); %(syrchangeinfo(q,3)+1):(syrchangeinfo(q,4)-1); %do not include end indexes as these are often long aquisitions
    syrind=syrind0(flag(syrind0)==3 | flag(syrind0)==6); %only want the cell records
    if ~isempty(syrind)
        rec(q,1:3)=quantile(acqtime(syrind),[0 0.5 0.96]); %[0 0.5 0.90]
        rec(q,4)=unique(pumpspeed(syrind));
        rec(q,5)=length(syrind);
    end
end

%identify potentially troublesome syringes globally by identifying abnormal acq times within a pumpspeed:
abn_syr=[];
for p=(unique(pumpspeed))'
    speed_ind=find(rec(:,4)==p); %should still be only for cell syringes
    baseline0=mean(rec(speed_ind,3)-rec(speed_ind,2))+5*std(rec(speed_ind,3)-rec(speed_ind,2));
    temp=find(rec(speed_ind,3)-rec(speed_ind,2) > baseline0); %%now find syringes that exceed baseline0
    abn_syr=[abn_syr; speed_ind(temp)]; %abn_syr indexes into syrchangeinfo
end

%% evaluate within a roll over event: calculate baseline and eval from there:
syr_excl=[]; syrflag=0;
for q=1:length(batch_ind)
    
    %first check to make sure this index has a full complement of syringes to use to establish a baseline acqtime:         
    if q==1, temp1=1; else temp=batch_ind(q-1); end
    if  batch_ind(q)-temp < 30
        %batch_ind(q)-batch_ind(q-1) < 30 %nope, need more syringes for a proper base!
        %find closest in time longest batch:
        fullbatch_ind=batch_ind(find(diff(batch_ind)>30)+1);
        if isempty(fullbatch_ind)
            s1=1;
            s2=size(syrchangeinfo,1);
        else
            [~,im]=min(abs(batch_ind(q)-fullbatch_ind));
            s2=fullbatch_ind(im);
            uu=find(batch_ind==fullbatch_ind(im));
            s1=batch_ind(uu-1)+1;
        end
    else
        if q==1, s1=1; else s1=batch_ind(q-1)+1; end
        s2=batch_ind(q);
    end
    
    %exclude any potentially troubling syringes when calculating base line
    syrbase=setdiff(s1:s2,abn_syr);
    temprange=[];
    for s=1:length(syrbase) %get appropriate indices:
        temprange0=(syrchangeinfo(syrbase(s),3):syrchangeinfo(syrbase(s),4))';
        tt=find(flag(temprange0)==3 | flag(temprange0)==6); temprange=[temprange; temprange0(tt)]; %only want cell records!
    end
    
    %   alternative metric:
    %   stdbase=nanstd(acqtime(temprange));
    %   baseline=nanmean(acqtime(temprange))+3*stdbase;
    
    %calculate baseline for comparison:
    stdbase=quantile(acqtime(temprange),[0.15 0.85]);
    stdbase=stdbase(2)-stdbase(1);
    baseline=nanmedian(acqtime(temprange))+3*stdbase;
    if isnan(baseline), baseline=baseline0; end %happens if abn syr in with lots of flagged data
    
    %syringe indexes within this batch:
    if q==1, y1=1; else y1=batch_ind(q-1)+1; end
    y2=batch_ind(q);
    
    for s=y1:y2 %go through each syringe and see if too many records are above the baseline
        
        syrind=syrchangeinfo(s,3):syrchangeinfo(s,4);
        syrind=syrind(find(flag(syrind)==3 | flag(syrind)==6));   %only cell records
        
        if isempty(syrind)
            syrflag=2; %skip!
        elseif year2do==2014  %2014 had noise issues...
             if any(acqtime(syrind)-baseline > 30) %these are given a wide range
                syrflag=1;
            end
        elseif length(syrind) <= 20 %ugh...can't really do stats short syringes
            if any(acqtime(syrind)-baseline > min(10*stdbase,40)) %these are given a wide range
                syrflag=1;
            end
        elseif length(find(acqtime(syrind)-baseline > max(stdbase,3)))/length(syrind) > 0.08 | ismember(s,abn_syr)
            rr=quantile(acqtime(syrind),[0.5 0.95]); % a secondary check to make sure needs to be flagged....
            if rr(2)-rr(1) > 4*stdbase
                syrflag=1;
            elseif any(acqtime(syrind(2:end-1))-baseline > max(4*stdbase,10))  | any(((acqtime(syrind(2:end-1))-baseline)./stdbase) > 20) %just incase
                syrflag=1;
            end
        elseif any(acqtime(syrind(2:end-1))-baseline > max(4*stdbase,10)) | any(((acqtime(syrind(2:end-1))-baseline)./stdbase) > 20) %last catch chance...any obs that is essentially 6 std away - just flag...(excludes slow acq times at beginning of syr - those are caught below)
            syrflag=1;
        end
        
        if syrflag==0 %if syringe wasn't excluded, check for other outliers at beginning and end of syringe:
            if length(syrind) >=10 %don't really try this with syringes with only a few records
                within_syr_metric=mean(acqtime(syrind(3:end-1)))+6*std(acqtime(syrind(3:end-1)));
                if acqtime(syrind(1)) > within_syr_metric
                    flag(syrind(1))=61; %if beginning record is over metric, flag
                    if acqtime(syrind(2)) > within_syr_metric,flag(syrind(2))=61; end  %sometimes, second record are also over this metric
                end
                if acqtime(syrind(end)) > within_syr_metric,flag(syrind(end))=61; end  %if end record is over metric, flag
            end
        elseif syrflag==1 %flag that syringe!
            flag(syrind)=60; %only flag cell records - don't overwrite 97,98's
            syr_excl=[syr_excl; s]; %indexes into syrchangeinfo and rec
        end
        
        syrflag=0; %reset flag
        baseline0=baseline;
    end
end

%% SECONDARY CHECKS:
% a double check to exclude additional syringes not flagged by the above code that are:
% right after a roll over event, and before troublesome syringes
% or inbetween troublesome syringes and not flagged
% or right before an end of a deployment and after troublemsome syringes

for j=1:size(syr_excl,1)
    
    %troublesome syrignes are at start of roll over event -> exclude
    %previous syringes if not included already
    if syrchangeinfo(syr_excl(j),5) == 3 || syrchangeinfo(syr_excl(j),5) == 2
        %             check if previous syringes have been included:
        %keyboard
        count=syrchangeinfo(syr_excl(j),5);
        while count > 1
            if ~ismember(syr_excl(j)-count+1,syr_excl) %if the previous syringe is missing - add!
                %keyboard
                syrind=syrchangeinfo(syr_excl(j)-count+1,3):syrchangeinfo(syr_excl(j)-count+1,4);
                flag(syrind)=62;
            end
            count=count-1;
        end
    end
    
    %syringes often are a problem at the end of a deployment - exclude
    %syringes not already flagged that come after troublesome ones near end
    if any(gap_syr-syr_excl(j) < 3 & gap_syr-syr_excl(j) > 0) %gap is ahead of current index
        gg=find((gap_syr-syr_excl(j) < 3) & (gap_syr - syr_excl(j) >=0)); %identify gap syr
        % keyboard
        to_add=find(ismember(syr_excl(j)+1:gap_syr(gg),syr_excl)==0);
        for i=1:length(to_add)
            % syr_excl=[syr_excl; syr_excl(j,1)+to_add(i) NaN];
            syrind=syrchangeinfo(syr_excl(j)+to_add(i),3):syrchangeinfo(syr_excl(j)+to_add(i),4);
            flag(syrind)=62;
        end
    end
    
    %look ahead: probably sufficient:
    %must use cellind to check indexes of only cell syringes:
    if any(diff(ismember(syr_excl(j)+1:min(syr_excl(j)+4,length(syrchangeinfo)),syr_excl))==1) %there are missing syringes!
        aa=find(ismember(syr_excl(j)+1:min(syr_excl(j)+4,length(syrchangeinfo)),syr_excl));
        for i=1:aa
            if ~ismember(syr_excl(j)+i,syr_excl)
                syrind=syrchangeinfo(syr_excl(j)+i,3):syrchangeinfo(syr_excl(j)+i,4);
                flag(syrind)=62;
                syr_excl=[syr_excl; syr_excl(j)+i];
            end
        end
        
    end
end


%% if plotflag is on, plot a final summary for sanity check:
if syrplotflag
    
    figure(16)
    %syringe number and pump speed over time
    subplot(2,1,1,'replace')
    plot(totalstartsec,syrpumpinfo(:,5),'k.-'); %syringe number and time
    hold on
    h1=plot(totalstartsec,pumpspeed,'.','color',[0 0.5 1]); %syringe designated speed and time
    hleg=legend(h1(1),'Assigned pump speed','location','NorthOutside');
    if ~isempty(plot_missing_ind), h2=plot(totalstartsec(plot_missing_ind),pumpspeed(plot_missing_ind),'.','color',[1 0.5 0]);
        temp=hleg.String; uu=find(strcmp('data1',temp)==1); hleg.String(uu)={'Speeds from closest neighbor'}; end %syringe speed assigned from closest syringe
    
    ylabel('Syringe Number or Syringe Speed')
    xlabel('Time (sec)')
    
    % acquistion time vs. time
    subplot(2,1,2,'replace')
    plot(totalstartsec,acqtime,'.-') %acquistion time vs. time
    hold on
    ii=find(flag==3 | flag==6 | flag==60 | flag==61 | flag==62);
    if ~isempty(ii)
        h3=plot(totalstartsec(ii),acqtime(ii),'.'); %only cell syringes
        hleg=legend(h3(1),'cell records','location','NorthOutside');
    end
    ww=find(flag==60); %cell syringes excluded
    if ~isempty(ww), h4=plot(totalstartsec(ww),acqtime(ww),'o','markersize',4,'linewidth',1.5,'color',[1 0.7 0]);
        temp=hleg.String; uu=find(strcmp('data1',temp)==1); hleg.String(uu)={'syringe excl'}; end
    mm=find(flag==61); %first or last record of a syringe excluded
    if ~isempty(mm), h5=plot(totalstartsec(mm),acqtime(mm),'o','markersize',4,'linewidth',1.5,'color',[0 0 0.7]);
        temp=hleg.String; uu=find(strcmp('data1',temp)==1); hleg.String(uu)={'record excl'}; end
    kk=find(flag==62); %syringe excluded at beginning of a rollover
    if ~isempty(kk), h6=plot(totalstartsec(kk),acqtime(kk),'s','markersize',4,'linewidth',1.5,'color',[0.8 0 0]);
        temp=hleg.String; uu=find(strcmp('data1',temp)==1); hleg.String(uu)={'2^{\circ} catch'}; end
    
    ylabel('Acquisition time (sec)')
    xlabel('Time (sec)')
    
    disp('Check the plots and then type dbcont to move to next batch of processing....')
    disp(['Ending pump speed: ' num2str(pumpspeed(end))])
    keyboard
    
end %syrplotflag


%% and now for the magic:

if year2do <= 2004
    offset = pumpspeed*(deadtime+querytime)*(totalvol/maxpos);
    analvol = analvol - offset;
    outmatrix = [200*(1:length(totalstartsec))' totalstartsec totalendsec acqtime analvol flag];
    
else
    offset = pumpspeed*querytime*(totalvol/maxpos); %offset = steps/sec*(q+d time)*(totalvol/maxpos)
    analvol = analvol - offset;
    outmatrix = [events totalstartsec totalendsec acqtime analvol flag];  %2/24/05 heidi, 100 event records for 12 channels
end

clear acqtime flag t ind analvol maxpos start stop totalvol

