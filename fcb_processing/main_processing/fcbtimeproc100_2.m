%fcbtimeproc = modified from fcbtimeproc to work for 12 channel setup in
%lab, now with 100 events per record (instead of 200 as before) 2/24/05
%3/7/03 fcbtimeproc = modified from cytosubproc2_labalt to handle new file format (with
%header), should work with any combination of port switching during analysis, Heidi
%3/1/03 modified from cytosubproc2 to handle files from SYN expt's in lab
%with two cultures analyzed (switch every half hour); Heidi
%keyboard
%
flag = zeros(size(totalstartsec));
%acqtime = totalendsec - totalstartsec - .405;  %for tq-acq-qt?
%acqtime = totalendsec - totalstartsec - .1326; %for vnts-acq-stnv
querytime = 0.0908; %0.0930;
acqtime = totalendsec - totalstartsec - querytime; %for ts-acq-st  %July 2016, from high noise bead runs

%acqtime = totalendsec - totalstartsec; %for qt-acq-tq (mr07)
flag(:) = 0;  %default all records to Not use
%1=temp,2=humidity,3=start port,4=end port,5=start syr#,6=end syr#,7=start syr pos,8=end syr pos.

%flag(find(syrpumpinfo(:,3) == 6 & syrpumpinfo(:,4) == 6 & syrpumpinfo(:,5) & syrpumpinfo(:,6))) = 1; %these are good records for culture 1
%flag(find(syrpumpinfo(:,3) == 3 & syrpumpinfo(:,4) == 3 & syrpumpinfo(:,5) & syrpumpinfo(:,6))) = 2; %these are good records for culture 2

%good records for cells or beads, skipping cases where syringe refills or
%port (valve) changes, also skipping any fast syringes (i.e., syr# = 0)
%non-zero port, same port at start and end, same syringe # at start and end (i.e. syringe did not refill during acquisition
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

start = syrpumpinfo(:,7);
stop = syrpumpinfo(:,8);
totalvol = .25;  %ml vol of syringe
maxpos = 48000; %48000 changed to 48200 for 2014 test, Heidi June 2016
analvol = stop*NaN;

t = find(start - stop >= 0);  %start > stop
analvol(t) = (start(t)-stop(t))/maxpos*totalvol;

%flags:
%99 - fast syringes
%98 - syr refill at start of acqusition
%97 - syr refill during acquistion
%96 - next two records after syringe refill

%t = find(start - stop < 0);  %start < stop --> syringe refilled
t = find(syrpumpinfo(:,5) < syrpumpinfo(:,6)); %syringe refilled during acquisition
analvol(t) = (start(t) + maxpos-stop(t))/maxpos*totalvol;
flag(t) = 97; %add specific flag for these syringes KRHC, 7/13/16
% sacq=[t+1;t+2]; sacq2=sacq(sacq <= length(syrpumpinfo)); %find 2nd acquistion records after refilling
% flag(sacq2)=96; %flag the next record and record after as likely contains sheath fluid?

t = find(syrpumpinfo(:,5) == syrpumpinfo(:,6) & start - stop <= 0); %syringe in middle of refilling at start; June 2015 <=0 (not < 0) for cases in 2014 with stopped syr at top or bottom during v. noisy fast acq
analvol(t) = (maxpos-stop(t))/maxpos*totalvol;
flag(t) = 98; %syringes refilling in middle
% sacq=[t+1;t+2]; sacq2=sacq(sacq <= length(syrpumpinfo));
% flag(sacq2)=96;

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

%analvol = analvol - 1.692e-4;  %volume offset for query, for mr07 set (3/03)
%analvol = analvol - 5.48e-5;  %older, corresponds to P2 speed for two syringe queries? July 2016


%% FIND VOLUME ANALYZED...turns out, not exactly straightforward...

%Notes:
%There are syringe movements not accounted for in the measurements, such
%that volume is being moved while we are querying for pump position or for
%dead time after a trigger. From Alexi's measurements (7/12/16) of full noise, we see
%that average dead time for 100 triggers is .0275s (27.5ms) (without query time)
%At a pump speed of 160 steps/sec, avg steps taken is 14.89 per 100 records with query time
% At 40 steps/sec, avg steps taken is 3.72. In these cases, since the
% triggers could be considered instantaneous, these steps/time are
% essentially due to deadtime and querytime, which leaves:

%       14.89 steps = (160 steps/s) * (dead time + query time)
%       dead+query = 0.0931 sec
%       3.72 steps = (40 steps/s) * (dead time + query time)
%       dead+query = 0.0930 sec !!!!!

%So, the lost time in query ends up being 0.0655...but, we only really care about the total, so
%multiply pump speed by 0.0930 sec to get steps lost, which gives volume lost!

%3 different speeds of pump:
%P3 - 40 steps/sec -> 20min per syr -> syrnum rolls over at 50
%P2 - 80 steps/sec -> 10min per syr -> syrnum rolls over at 100
%P1 - 160 steps/sec -> 5min per syr -> syrnum rolls over at 200

%%can find pump speed by looking at steps per second, avg syring time, or the most robust metric - slope of records/time:
stepdist=syrpumpinfo(:,7)-syrpumpinfo(:,8);
timediff=totalendsec-totalstartsec;
pump_speed=stepdist./timediff; %rough speed

%     one way to get average syringe time:
%     ds=find(diff(syrpumpinfo(:,5))~=0); %look for syringe changes...
%     avgsyrtime=[totalstartsec(ds(2:end)) (1/60)*(totalstartsec(ds(2:end))-totalstartsec(ds(1:end-1)))]; %time difference inbetween

% a slightly clunkier way to find syringe changes, but perhaps more straight forward:
if length(unique(syrpumpinfo(:,5))) > 1
    count=1;
    ii=1;
    syrnum=syrpumpinfo(1,5);
    %syrchangeinfo=[syr_starttime  syr_endtime  start_index end_index syrnum avgsyrtime syrnum_slope]
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
    
    %keyboard
    
    %find the slope of syringe numbers over time:
    celli=find(syrchangeinfo(:,6)==3); %use only the cell syringes for more reliable slopes
    ii=find(diff(syrchangeinfo(celli,5))<0);
    ro1=celli(ii); %rollover indexes
    ro2=celli(ii+1); %start of new cell syringe
    
    ro=[1; ro1; ro2; size(syrchangeinfo,1)];
    ro=sort(ro);
    
    %     ii=find(diff(syrchangeinfo(:,5))<0);
    %     ii2=find(syrchangeinfo(ii,5)~=1);
    %     ro=ii(ii2); %rollover indexes
    %     ro=[1; ro; size(syrchangeinfo,1)];
    %%
    %check if there are any split syringes?
    % use 3hr gap ~ 10000 sec as cutoff to split the syringe -> case where acq stopped, but resumes syr num?
    ind2add=[];
    for q=1:2:length(ro)-1
        if any(diff(syrchangeinfo(ro(q):ro(q+1),1)) > 10000) %3hr gap ~ 10000 sec - split the syringe -> case where acq stopped, but resumes syr num
            keyboard
            disp(['split syringe?' num2str(q)])
            jj=ro(q):ro(q+1); %easier to handle indexes
            kk=find(diff(syrchangeinfo(jj,1)) > 10000);
            ind2add=[ind2add; jj(kk)'; jj(kk+1)']; %add this split into the rollover indexes
        end
    end
    
    ro=sort([ro; ind2add]);
    
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
    % keyboard
    %%
    %use the syrchangeinfo matrix to populate matrix for each record:
    avgsyrtime=nan(size(syrpumpinfo,1),2); %syringe times and slopes
    for q=1:size(syrchangeinfo,1)
        avgsyrtime(syrchangeinfo(q,3):syrchangeinfo(q,4),1)=(1/60)*(syrchangeinfo(q,2)-syrchangeinfo(q,1));
        avgsyrtime(syrchangeinfo(q,3):syrchangeinfo(q,4),2)=syrchangeinfo(q,8);
    end
    
    P3=find(avgsyrtime(:,2) > 1000 & avgsyrtime(:,2) < 1800); %slope = ~1250 syringes/time ~ rollover
    P2=find(avgsyrtime(:,2) > 500 & avgsyrtime(:,2) < 700); %slope = ~615 syringes/time ~ rollover
    P1=find(avgsyrtime(:,2) > 250 & avgsyrtime(:,2) < 400); %slope = ~314 syringes/time ~ rollover
    
    %% any cells or bead runs unaccounted for?
    ind1=find(flag==3 | flag ==6); %cells and beads
    test0=ismember(ind1,[P1;P2;P3]);
    goodrate=find(test0==1); goodrate=ind1(goodrate); %index still into syrpumpinfo
    test=find(test0==0);
    
    if ~isempty(test)
        jj=ind1(test); %easier to handle indexes
        disp(['bead rate and other rates unaccounted for: ' num2str(length(test)) ' ...using closest speed'])
        
        %     use speed and average syringe time to further help assign a speed:
        %     p2=find((pump_speed(jj) > 70 & pump_speed(jj) < 85) | (avgsyrtime(jj) > 8 & avgsyrtime(jj) < 11));
        %     p1=find((pump_speed(jj) > 150 & pump_speed(jj) < 165) | (avgsyrtime(jj) > 3.5 & avgsyrtime(jj) < 6 ));
        %     p3=find((pump_speed(jj) > 30 & pump_speed(jj) < 45) | (avgsyrtime(jj) > 18 & avgsyrtime(jj) < 22 ));
        
        for i=1:length(jj)
            qq=find(goodrate <= jj(i)); %for looking back in time....
            [~, im]=min(totalstartsec(jj(i))-totalstartsec(goodrate(qq)));  %find closest good bead or cell run in (backwards) time that has a pumprate
            %[~, im]=min(abs(jj(i)-goodrate)); %find closest good bead or cell run that has a pumprate (in index)
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
                    keyboard
                end
            end
        end
    end
    
    pumprate=zeros(size(syrpumpinfo,1),1);
    pumprate(P3)=40;
    pumprate(P2)=80;
    pumprate(P1)=160;
    
    %plots for sanity check:
    if syrplotflag
        subplot(2,1,1,'replace')
        plot(totalstartsec,syrpumpinfo(:,5),'k.-');
        hold on
        if ~isempty(test)
            plot(totalstartsec(ind1(test)),syrpumpinfo(ind1(test),5),'r.');
        end
        ylabel('Syringe Number')
        
        subplot(2,1,2,'replace'), hold on
        plot(totalstartsec,pump_speed,'.-','color',[0.6 0.6 0.6])
        plot(totalstartsec(ind1),pump_speed(ind1),'.','color',[0 0.5 1])
        plot(totalstartsec(ind1),pumprate(ind1),'.','color',[0 0 1])
        if ~isempty(test)
            plot(totalstartsec(ind1(test)),pumprate(ind1(test)),'c.');
        end
        ylabel('Pump rate (steps/sec)')
        ylim([-10 180])
        
        keyboard
        %pause(0.5)
    end
    
else %for rare case of only 1 syringe being processed...
    %keyboard
    P1=find((pump_speed > 70 & pump_speed < 85));
    P2=find((pump_speed > 150 & pump_speed < 165));
    P3=find((pump_speed > 30 & pump_speed < 45));
    
    pumprate=zeros(size(syrpumpinfo,1),1);
    pumprate(P3)=40;
    pumprate(P2)=80;
    pumprate(P1)=160;
end

% So...now we need a metric to flag syringes or records that do not seem to
% follow a constant syringe movement (i.e. a blocked syringe, sheath being
% sucked into a partially clogged syringe, etc.

%This is primarily done with looking at the acquisition time over time:

%Maybe a plan would be to go through, keep a running average of the mean
%acquisition time, then compare, syringe by syringe?

%Hmm...maybe a simple test with median values...and then if more than 0.25
%of the data is outside the expected acquisition time range, discard that
%syringe...

%Ooh- or, better yet - compare to see if roughly follows a normal
%distribution of acquistion times - if not, flag!

%% okay: screen syringes

%go through each syringe, check whether or not follows a normal dist based
%on linear test between ordered data and the normal inverse of data
%positions
%From Andy's notes, plot: y(1) < y(2) < y(3)  (these being the acquistion
%times of a syringe)
%against F^-1(1/(n+1), F^-1(2/n+1), F^-1(3/n+1), where F^-1 is the norm inv
%parameterized from the syringe data

for q=1:length(syrchangeinfo)
    
    %syrchangeinfo has the indexes per syringe at columns 3 and 4:
    inds=syrchangeinfo(q,3):syrchangeinfo(q,4);
    tempacq=acqtime(inds); %dataslices
    tempflag=flag(inds);
    tf=find(tempflag==3); %only evaluate cell records (ask of these, should any records be removed?)
    n=length(tempacq(tf));
    
    if ~isempty(tf) & n~=1 %can't be empty and length can't be 1
          
        finv = (norminv(((1:n)/(n+1)),mean(tempacq(tf)),std(tempacq(tf))))'; %expected value given an underling normal distribution
        [ordered_acq, oind] = sort(tempacq(tf)); %ordered acquisition times              
        ss=sum(sqrt((ordered_acq-finv).^2));
        
        %check deviation from expected:
        
        %first check for some outliers that occur at beginning and end of syringes
        if (oind(end)==1 || oind(end)==n) && (ordered_acq(end)-ordered_acq(end-1) > 1.5*std(tempacq(tf))) %first or last point and greater than 1 std dev away
            
            %remove this point and see if better ss:
            n=n-1;
            ordered_acq2=ordered_acq(1:end-1);
            finv2 = (norminv(((1:n)/(n+1)),mean(ordered_acq2),std(ordered_acq2)))'; %expected value given an underling normal distribution

            ss2=sum(sqrt((ordered_acq2-finv2).^2));
            
            if ss2 < 3*mean(tempacq(tf))
                disp('fixed it! excluding that datapacket from syringe')
                flag(inds(tf(oind(end))))=60; %flag it!
                syringe_acqtime_QC2(finv,ordered_acq,oind,ss,finv2,ordered_acq2,ss2)
                keyboard
            else
                disp('hmmmm...more serious problems with this syringe?')
                syringe_acqtime_QC2(finv,ordered_acq,oind,ss,finv2,ordered_acq2,ss2)
                keyboard
                flag(inds(tf))=61; %flag it!
            end
            
        elseif ss > 4*mean(tempacq(tf))
            disp('hmmmm...more serious problems with this syringe?')
            syringe_acqtime_QC2(finv,ordered_acq,oind,ss)
            keyboard
            flag(inds(tf))=61; %flag it!
        end      
        
        %pause
        
    end
end

%%
%analvol(t) = (start(t)-stop(t))/maxpos*totalvol;
%offset = steps/sec*(q+d time)*(totalvol/maxpos)
offset = pumprate*querytime*(totalvol/maxpos);
analvol = analvol - offset;

%keyboard
%outmatrix = [1:length(totalstartsec) totalstartsec totalendsec acqtime medianinterval flag ];
%outmatrix = [200*(1:length(totalstartsec))' totalstartsec totalendsec acqtime analvol flag];
%%
outmatrix = [100*(1:length(totalstartsec))' totalstartsec totalendsec acqtime analvol flag];  %2/24/05 heidi, 100 event records for 12 channels

clear acqtime flag t ind analvol maxpos start stop totalvol

