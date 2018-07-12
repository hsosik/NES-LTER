%3/7/03 fcbtimeproc = modified from cytosubproc2_labalt to handle new file format (with
%header), should work with any combination of port switching during analysis, Heidi
%3/1/03 modified from cytosubproc2 to handle files from SYN expt's in lab
%with two cultures analyzed (switch every half hour); Heidi

flag = zeros(size(totalstartsec));
%acqtime = totalendsec - totalstartsec - .405;  %for tq-acq-qt?
acqtime = totalendsec - totalstartsec - .1326; %for vnts-acq-stvn
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
maxpos = 48200;
analvol = stop*NaN;
t = find(start - stop >= 0);  %start > stop
%tt = start(t)-stop(t);
%tt(tt>48000) = 48000;
%analvol(t) = tt/maxpos*totalvol;
analvol(t) = (start(t)-stop(t))/maxpos*totalvol;
%t = find(start - stop < 0);  %start < stop --> syringe refilled
t = find(syrpumpinfo(:,5) < syrpumpinfo(:,6)); %syringe refilled during acquisition
analvol(t) = (start(t) + maxpos-stop(t))/maxpos*totalvol;
t = find(syrpumpinfo(:,5) == syrpumpinfo(:,6) & start - stop < 0); %syringe in middle of refilling at start
analvol(t) = (maxpos-stop(t))/maxpos*totalvol;
flag(t) = 98; %syringes refilling in middle 
%last record before end of set (switch to new valve) when syringe is at end, sometimes these have really long acq times (renegade triggers?)
flag(find(syrpumpinfo(2:end,6) < syrpumpinfo(1:end-1,6) & syrpumpinfo(1:end-1,8) == 10)) = 98; 

%analvol = analvol - 1.692e-4;  %volume offset for query, for mr07 set (3/03)
analvol = analvol - 5.48e-5;  %

%outmatrix = [1:length(totalstartsec) totalstartsec totalendsec acqtime medianinterval flag ];
outmatrix = [200*(1:length(totalstartsec))' totalstartsec totalendsec acqtime analvol flag];

clear acqtime flag t ind analvol maxpos start stop totalvol

