clear all
close all

year=2017;
%pathname='/home/kristen/Desktop/MVCO_FCB_time_testfiles/';
%eval(['filelist=dir(''/home/kristen/Desktop/MVCO_FCB_time_testfiles/*' num2str(year) '*.mat'')'])

pathname='/Volumes/Lab_data/MVCO/FCB/MVCO_Jan2017/data/processed/time/';
%pathname='/Volumes/Lab_data/MVCO/FCB/MVCO_Apr2005/data/processed_July2016/time/';
%pathname='/Volumes/Lab_data/MVCO/FCB/MVCO_May2006/data/processed_July2016/time/';

%eval(['filelist=dir(''' pathname '*_*.mat'')'])
eval(['filelist=dir(''' pathname '*' num2str(year) '*.mat'')'])
syrplotflag=1;
%%

clearvars -except year filelist pathname syrplotflag

j=1;

filename=filelist(j).name(1:end-6);
load(strcat(pathname,filelist(j).name))
eval(['temp=' filename ';'])
totalstartsec=24*3600*(temp(:,2)-temp(1,2));
totalendsec=24*3600*(temp(:,3)-temp(1,2));

%% just the current plots:

    figure(16)
    %syringe number and pump speed over time
    subplot(2,1,1,'replace')
    plot(totalstartsec,syrpumpinfo(:,5),'k.-'); %syringe number and time
    hold on
    ylabel('Syringe Number')
    xlabel('Time (sec)')
    
    % acquistion time vs. time
    subplot(2,1,2,'replace')
    plot(totalstartsec,temp(:,4),'.-') %acquistion time vs. time
    hold on
    ii=find(temp(:,6)==3 | temp(:,6)==60 | temp(:,6)==61 | temp(:,6)==62 | temp(:,6)==6);
    h3=plot(totalstartsec(ii),temp(ii,4),'.'); %only cell syringes
    hleg=legend(h3(1),'cell records','location','NorthOutside');
    
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



