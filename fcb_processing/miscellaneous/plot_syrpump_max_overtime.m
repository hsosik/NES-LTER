% More informative plots on each syringe:

%Okay, need to get a sense of how many times position values go above
%48000...i.e. how systematic this is / regular this is...

%plan:
%identify syringes
%ask, for that syringe, what was the max?
%plot

basepath='\\sosiknas1\Lab_data\MVCO\FCB\';
% basepath='/Volumes/Lab_data/MVCO/FCB/';

max_syrpos=[];
for year2do=2003:2016
    
    switch year2do
        case 2003
            yearlabel='May';
        case 2004
            yearlabel='Apr';
        case 2005
            yearlabel='Apr';
        case 2006
            yearlabel='May';
        case 2007
            yearlabel='Mar';
        otherwise
            yearlabel='Jan';
    end
    
    timepath=[basepath 'MVCO_' yearlabel num2str(year2do) '\data\processed\time\'];
    %     timepath=[basepath 'MVCO_' yearlabel num2str(year2do) '/data/processed/time/'];
    
    filelist=dir([timepath '*.mat']);
    
    for j=1:length(filelist)
        
        filename=filelist(j).name;
        eval(['load ' timepath filename])
        disp(['loaded...' filename])
        
        %for easier handling, declare time info as temptime:
        if year2do <=2005
            eval(['temptime=' filename(1:end-6) ';'])
        else
            tempvar=whos('FCB*');
            eval(['temptime=' tempvar.name ';'])
        end
        
        %Find FCB number:
        if year2do <=2006
            fcbnum=1;
        else
            tempnum=regexp(filename,'FCB(?<fcb>\d{1})_','names');
            tempnum=tempnum.fcb; fcbnum=str2num(tempnum);
        end
        
        syr=find(diff(syrpumpinfo(:,5))~=0);
        
        temprec=nan(length(syr),4);
        for q=1:length(syr)-1
            temp=syrpumpinfo(syr(q)+1:syr(q+1),:); %temporary slice of data corresponding to one syringe
            
            temprec(q,1:2)=[fcbnum temptime(syr(q),2)];
            temprec(q,3:4)=max(temp(:,end-1:end));
        end
        
        max_syrpos=[max_syrpos; temprec];
        
        clear FCB* temptime fcbnum temprec
    end %filelist
    clearvars -except max_syrpos filelist basepath timepath
end %year2do

%% and some plots:

figure, hold on
%plot(max_syrpos(:,2),max_syrpos(:,4),'.')
f1=find(max_syrpos(:,1)==1);
f2=find(max_syrpos(:,1)==2);
plot(max_syrpos(f1,2),max_syrpos(f1,4),'.')
plot(max_syrpos(f2,2),max_syrpos(f2,4),'.')
line(xlim,[48000 48000],'color','r')
line(xlim,[48200 48200],'color','r')
set(gca,'xgrid','on')
ylabel('Max end position reached for each syringe')
xlim([datenum('1-1-03') datenum('12-31-16')])
datetick('x','keeplimits')
ylim([47500 48300])


%% %look at records discarding after every syringe refill:

basepath='\\sosiknas1\Lab_data\MVCO\FCB\';
% basepath='/Volumes/Lab_data/MVCO/FCB/';

max_syrpos=[];
throw_away=[];
for year2do=[2004 2008 2010 2012 2014]
    
    switch year2do
        case 2003
            yearlabel='May';
        case 2004
            yearlabel='Apr';
        case 2005
            yearlabel='Apr';
        case 2006
            yearlabel='May';
        case 2007
            yearlabel='Mar';
        otherwise
            yearlabel='Jan';
    end
    
    timepath=[basepath 'MVCO_' yearlabel num2str(year2do) '\data\processed\time\'];
    %     timepath=[basepath 'MVCO_' yearlabel num2str(year2do) '/data/processed/time/'];
    
    filelist=dir([timepath '*.mat']);
    
    for j=1:length(filelist)
        
        filename=filelist(j).name;
        eval(['load ' timepath filename])
        disp(['loaded...' filename])
        
        %for easier handling, declare time info as temptime:
        if year2do <=2005
            eval(['temptime=' filename(1:end-6) ';'])
        else
            tempvar=whos('FCB*');
            eval(['temptime=' tempvar.name ';'])
        end
        
        %Find FCB number:
        if year2do <=2006
            fcbnum=1;
        else
            tempnum=regexp(filename,'FCB(?<fcb>\d{1})_','names');
            tempnum=tempnum.fcb; fcbnum=str2num(tempnum);
        end
        
        syr=find(diff(syrpumpinfo(:,5))~=0); %syringe changes
        
        temprec=nan(length(syr),4);
        for q=1:length(syr)-1
            temp=syrpumpinfo(syr(q)+1:syr(q+1),:); %temporary slice of data corresponding to one syringe
            temprec(q,1:2)=[fcbnum temptime(syr(q),2)];
            temprec(q,3:4)=max(temp(:,end-1:end));
        end
        
        max_syrpos=[max_syrpos; temprec];
        
        ii=find(temptime(:,end) == 96);
        throw_away=[throw_away; temptime(ii,2) syrpumpinfo(ii,end-1:end)];      
        
        clear FCB* temptime fcbnum temprec
        
    end %filelist
    clearvars -except max_syrpos filelist basepath timepath throw_away
    
end %year2do


%% figures


figure, hold on
plot(max_syrpos(:,2),max_syrpos(:,4),'.')
plot(throw_away(:,1),throw_away(:,2),'.')
%plot(throw_away(:,1),throw_away(:,3),'.')
line(xlim,[48000 48000],'color','r')
line(xlim,[48200 48200],'color','r')
set(gca,'xgrid','on')
ylabel('Max end position reached for each syringe')
xlim([datenum('1-1-03') datenum('12-31-16')])
datetick('x','keeplimits')
ylim([47000 48300])


%% look at all time files and every syringe for a year:

%basepath='\\sosiknas1\Lab_data\MVCO\FCB\';
basepath='/Volumes/Lab_data/MVCO/FCB/';
querytime = 0.0908;
mattime=[];
acqtime=[];
%allsyringes=[];

year2do=2011;
    
    switch year2do
        case 2003
            yearlabel='May';
        case 2004
            yearlabel='Apr';
        case 2005
            yearlabel='Apr';
        case 2006
            yearlabel='May';
        case 2007
            yearlabel='Mar';
        otherwise
            yearlabel='Jan';
    end
    
    %timepath=[basepath 'MVCO_' yearlabel num2str(year2do) '\data\processed\time\'];
    timepath=[basepath 'MVCO_' yearlabel num2str(year2do) '/data/processed/time/'];
    
    filelist=dir([timepath '*.mat']);
    
    for j=1:length(filelist)
        
        filename=filelist(j).name;
        eval(['load ' timepath filename])
        disp(['loaded...' filename])
        
        %for easier handling, declare time info as temptime:
        if year2do <=2005
            eval(['temptime=' filename(1:end-6) ';'])
        else
            tempvar=whos('FCB*');
            eval(['temptime=' tempvar.name ';'])
        end
        
%         mattime=[mattime; temptime(:,2)];
%         acqtime=[acqtime; (24*3600)*(temptime(:,3) - temptime(:,2)) - querytime];
%         
        figure(1)
        acqtime=(24*3600)*(temptime(:,3) - temptime(:,2)) - querytime;
        plot(temptime(:,2),acqtime,'.-')
        datetick('x','mm/dd')
        keyboard
        clear FCB* temptime fcbnum temprec
        
    end %filelist
    


%% older:
% loop over all the years, all the time files and pull out max and min step
% size positions of the FCB syringe pump for just these time chunks:
%
% clear all
% close all
%
% basepath='\\sosiknas1\Lab_data\MVCO\FCB\';
% syr_record=[]; syr_recordB={};
%
% for year2do=2003:2016
%     switch year2do
%         case 2003
%             yearlabel='May';
%         case 2004
%             yearlabel='Apr';
%         case 2005
%             yearlabel='Apr';
%         case 2006
%             yearlabel='May';
%         case 2007
%             yearlabel='Mar';
%         otherwise
%             yearlabel='Jan';
%     end
%
%     timepath=[basepath 'MVCO_' yearlabel num2str(year2do) '\data\processed\time\'];
%
%     filelist=dir([timepath '*.mat']);
%
%     for j=1:length(filelist)
%
%         filename=filelist(j).name;
%         eval(['load ' timepath filename])
%         disp(['loaded...' filename])
%
%         if year2do <=2005
%
%             eval(['temptime=' filename(1:end-6) ';'])
%         else
%             tempvar=whos('FCB*');
%             eval(['temptime=' tempvar.name ';'])
%         end
%
%         find fcbnum:
%         if year2do <=2006
%             fcbnum=1;
%         else
%             tempnum=regexp(filename,'FCB(?<fcb>\d{1})_','names');
%             tempnum=tempnum.fcb; fcbnum=str2num(tempnum);
%         end
%
%         syr_record=[syr_record; fcbnum temptime(1,2) temptime(end,2) min(syrpumpinfo(:,end-1:end)) max(syrpumpinfo(:,end-1:end))];
%         syr_recordB=[syr_recordB; {filename}];
%
%         clear FCB* tempvar temptime fcbnum
%     end
% end
%
% % and for some plots!
%
% syrlabels = {'fcb#';
%     'startime of record';
%     'end time of record';
%     'min start pos';
%     'min end pos';
%     'max start pos';
%     'max end pos'};
%
% %
% plot max end position:
% clf, hold on
% i1=find(syr_record(:,1)==1);
% i2=find(syr_record(:,1)==2);
% plot(syr_record(i1,2),syr_record(i1,end),'k.','markersize',10)
% plot(syr_record(i2,2),syr_record(i2,end),'.','color',[0 0.5 1],'markersize',10)
% plot(syr_record(i1,2),syr_record(12,end),'k.','markersize',10)
% plot(syr_record(i2,2),syr_record(i2,end),'.','color',[0 0.5 1],'markersize',10)
% datetick
% line([datenum('1-1-03') datenum('12-31-16')],[48000 48000],'color','r')
% line([datenum('1-1-03') datenum('12-31-16')],[48200 48200],'color','r')
% ylabel('Syringe Max End Position')
% legend('FCB1','FCB2')

%
% figure
%  hold on
% %i1=find(syr_record(:,1)==1);
% i2=find(syr_record(:,1)==1);
% plot(syr_record(i2,2),syr_record(i2,end),'b.')
% plot(syr_record(i2,2),syr_record(i2,end-1),'r.')
% datetick
% line([datenum('1-1-03') datenum('12-31-16')],[48000 48000],'color','r')
% line([datenum('1-1-03') datenum('12-31-16')],[48200 48200],'color','r')
% ylabel('Syringe Max Position')
% legend('Max End Pos','Max Start Pos')

