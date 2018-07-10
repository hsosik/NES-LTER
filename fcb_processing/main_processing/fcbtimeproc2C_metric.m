%3/7/03 fcbtimeproc = modified from cytosubproc2_labalt to handle new file format (with
%header), should work with any combination of port switching during analysis, Heidi
%3/1/03 modified from cytosubproc2 to handle files from SYN expt's in lab
%with two cultures analyzed (switch every half hour); Heidi

flag = zeros(size(totalstartsec));
%acqtime = totalendsec - totalstartsec - .405;  %for tq-acq-qt?
%acqtime = totalendsec - totalstartsec - .1326; %for vnts-acq-stvn, used until July 2016

querytime=0.0908-0.0275; % time for queries on 100 trigger events - high noise investigation, July 2016
deadtime=(0.0275/100)*200; %events per record pre-2005 have 200 events in them, adjust deadtime for 200 events
acqtime = totalendsec - totalstartsec - (querytime+deadtime); %for ts-acq-st  %new query time, July 2016, from high noise bead runs

%acqtime = totalendsec - totalstartsec; %for qt-acq-tq (mr07)
flag(:) = 0;  %default all records to Not use
%1=temp,2=humidity,3=start port,4=end port,5=start syr#,6=end syr#,7=start syr pos,8=end syr pos.

%flag(find(syrpumpinfo(:,3) == 6 & syrpumpinfo(:,4) == 6 & syrpumpinfo(:,5) & syrpumpinfo(:,6))) = 1; %these are good records for culture 1
%flag(find(syrpumpinfo(:,3) == 3 & syrpumpinfo(:,4) == 3 & syrpumpinfo(:,5) & syrpumpinfo(:,6))) = 2; %these are good records for culture 2

%good records for cells or beads, skipping cases where syringe refills or
%port (valve) changes, also skipping any fast syringes (i.e., syr# = 0)
%keyboard
%non-zero port, same port at start and end, same syringe # at start and end (i.e. syringe did not refill during acquisition
%ind = find((syrpumpinfo(:,3) == syrpumpinfo(:,4)) & syrpumpinfo(:,5) & syrpumpinfo(:,6) & (syrpumpinfo(:,5) == syrpumpinfo(:,6))));
ind = 1 + find((syrpumpinfo(2:end,3) == syrpumpinfo(2:end,4)) & syrpumpinfo(2:end,5) & syrpumpinfo(2:end,6) & (syrpumpinfo(2:end,5) == syrpumpinfo(2:end,6)) & (syrpumpinfo(1:end-1,6) == syrpumpinfo(2:end,5)));
%last part to skip cases where syringe starts to refill during preceeding data transmission (start syringe # ~= previous end syringe #)
%these records often seem to have flow problem or something that makes apparent volume too high
flag(ind) = syrpumpinfo(ind,3);  %set good records to syringe port # (3/03 6=culture 1, 3=culture2, ?=beads)
%follwing finds syringe refills starting mid-record
%ind = find((syrpumpinfo(:,3) == 3 & syrpumpinfo(:,4) == 3) & syrpumpinfo(:,5) & syrpumpinfo(:,6) & (syrpumpinfo(:,5) ~= syrpumpinfo(:,6)));
ind = find(syrpumpinfo(1:end-1,5) & syrpumpinfo(2:end,5) == 0);  %transitions to fast syringes
flag(ind) = 99;

start = syrpumpinfo(:,7);
stop = syrpumpinfo(:,8);
totalvol = .25;  %ml vol of syringe
maxpos = 48000;
analvol = stop*NaN;
t = find(start - stop >= 0);  %start > stop
%tt = start(t)-stop(t);
%tt(tt>48000) = 48000;
%analvol(t) = tt/maxpos*totalvol;
analvol(t) = (start(t)-stop(t))/maxpos*totalvol;
%t = find(start - stop < 0);  %start < stop --> syringe refilled
t = find(syrpumpinfo(:,5) < syrpumpinfo(:,6)); %syringe refilled during acquisition
analvol(t) = (start(t) + maxpos-stop(t))/maxpos*totalvol;
flag(t) = 97; %add specific flag for these syringes KRHC, 7/13/16
% sacq=[t+1;t+2]; sacq2=sacq(sacq <= length(syrpumpinfo)); %find 2nd acquistion records after refilling
% flag(sacq2)=96; %flag the next record and record after as likely contains sheath fluid

t = find(syrpumpinfo(:,5) == syrpumpinfo(:,6) & start - stop < 0); %syringe in middle of refilling at start
analvol(t) = (maxpos-stop(t))/maxpos*totalvol;
flag(t) = 98; %syringes refilling in middle
% sacq=[t+1;t+2]; sacq2=sacq(sacq <= length(syrpumpinfo));
% flag(sacq2)=96;

%last record before end of set (switch to new valve) when syringe is at end, sometimes these have really long acq times (renegade triggers?)
flag(find(syrpumpinfo(2:end,6) < syrpumpinfo(1:end-1,6) & syrpumpinfo(1:end-1,8) == 10)) = 98;


% FIND VOLUME ANALYZED...turns out, not exactly straightforward...

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
%Make a matrix that has all the info on an individual syringe:
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
    
    %Now check for any delayed syringes and adjust matrix to reflect this: delays that simply resume the same syringe will cause problems for slope calculations:
     % use 3hr gap ~ 10000 sec as cutoff to split the syringe -> case where acq stopped, but resumes syr num?
    ds=find(syrchangeinfo(:,2)-syrchangeinfo(:,1) > 10000); %ds for delayed syringe
    if ~isempty(ds)
        for i=1:length(ds)
            disp('delayed syringe?')
            ind_slice=syrchangeinfo(ds(i),3):syrchangeinfo(ds(i),4);
            start_slice=totalstartsec(ind_slice);
            end_slice=totalendsec(ind_slice);
            oo=find(diff(start_slice) > 10000);
            
            syrchangeinfo=[syrchangeinfo; syrchangeinfo(ds(i),:)]; %add extra entry to modify
            syrchangeinfo(ds(i),2)=end_slice(oo); %change the end time
            syrchangeinfo(ds(i),4)=ind_slice(oo); %change the end index
            syrchangeinfo(end,1)=start_slice(oo+1); %change the start time for the new entry
            syrchangeinfo(end,3)=ind_slice(oo+1); %change the index for the new entry
        end
    end
    
    %now, resort so back in chronological order if added any new lines for
    %delayed syringes:
    [~,is]=sort(syrchangeinfo(:,1));
    syrchangeinfo=syrchangeinfo(is,:);
    
    %calculate average syringe time:
    syrchangeinfo=[syrchangeinfo (1/60)*(syrchangeinfo(:,2)-syrchangeinfo(:,1))];
    
    %keyboard
    
    %find the slope of syringe numbers over time:
    %disp([unique(syrpumpinfo(:,3)) unique(syrpumpinfo(:,4))]) %sanity check for year 2003 and 2004
    
    %for 2003 and 2004, there are no syringes labeled as 3, only 6?
    celli=find(syrchangeinfo(:,6)==6); %use only the cell syringes for more reliable slopes
    ii=find(diff(syrchangeinfo(celli,5))<= 0);
    ro1=celli(ii); %rollover indexes
    ro2=celli(ii+1); %start of new cell syringe
    
    if ~isempty(celli)
        ro=[celli(1); ro1; ro2; celli(end)];
    else
        ro=[1; ro1; ro2; size(syrchangeinfo,1)]; 
    end
    
    ro=sort(unique(ro)); %in case file ends on the first cell syringe of a new roll over

    %% calculate the slope for each batch or syringes, in between roll over events:
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
    
    %use the syrchangeinfo matrix to populate a matrix of slopes and syringe times for each record:
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
        
        figure(1)
        subplot(4,1,1,'replace')
        plot(totalstartsec,syrpumpinfo(:,5),'k.-');
        hold on
        if ~isempty(test)
            plot(totalstartsec(ind1(test)),syrpumpinfo(ind1(test),5),'r.');
            legend('all syringes','syringes where slope method failed; using nearest pump rate')
        end

        ylabel('Syringe Number')
        
        subplot(4,1,2,'replace'), hold on
        plot(totalstartsec,pump_speed,'.-','color',[0.6 0.6 0.6])
        plot(totalstartsec(ind1),pump_speed(ind1),'.','color',[0 0.5 1])
        plot(totalstartsec(ind1),pumprate(ind1),'.','color',[0 0 1])
        if ~isempty(test)
            plot(totalstartsec(ind1(test)),pumprate(ind1(test)),'c.');
        end
        legend('all syringes','cell syringes','rate to be used','location','northoutside','orientation','horizontal')
        ylabel('Pump rate (steps/sec)')
        ylim([-10 180])
        
         subplot(4,1,3,'replace'), hold on
         plot(totalstartsec,acqtime,'.-')
         plot(totalstartsec(ind1),acqtime(ind1),'.','color',[0 0.5 1])
         ylabel('Acq. Time (sec)')
        
         subplot(4,1,4,'replace'), hold on
         plot(totalstartsec,syrpumpinfo(:,7),'.-')
         plot(totalstartsec(ind1),syrpumpinfo(ind1,7),'.','color',[0 0.5 1])
         ylabel('Syringe start postion')
         
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

%analvol(t) = (start(t)-stop(t))/maxpos*totalvol;
%offset = steps/sec*(q+d time)*(totalvol/maxpos)
offset = pumprate*(deadtime+querytime)*(totalvol/maxpos);
analvol = analvol - offset;


%analvol = analvol - 1.692e-4;  %volume offset for query, for mr07 set (3/03)
% analvol = analvol - 5.48e-5;  %

%outmatrix = [1:length(totalstartsec) totalstartsec totalendsec acqtime medianinterval flag ];
outmatrix = [200*(1:length(totalstartsec))' totalstartsec totalendsec acqtime analvol flag];

clear acqtime flag t ind analvol maxpos start stop totalvol

