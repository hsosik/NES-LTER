%Script that processes solar data form MVCO for input into setup_days for division rate
%model. Data is checked for gaps, and when possible, replaced with data
%from Nantucket buoy offshore.
clear all

solarpath = '\\sosiknas1\Lab_data\MVCO\EnvironmentalData\asit.mininode.CLRohn_rad_2020.csv';
solarsavepath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2020\syn_model\'

if ~exist(solarsavepath)
    mkdir(solarsavepath)
end


X = readtable(solarpath); 

%exclude nan's:
ind = find(~isnan(X.ir_rad_mean));
X = X(ind, :); 

dawnlevel = 5;   %threshold for light

%% Load in good estimates of where dawn/dusk is to detect gaps:

%constructed from raw_data with get_average_dawn.m, 
%load expected UTC values of dawn and dusk over 13 years :)


Solar = X.solar_rad_mean; 
date_met = datenum(X.timestamp)
year2do  = year(X.timestamp(1)); 
Hour = hour(X.timestamp)

load('../setup_and_postprocessing/median_dawn.mat')

%% Find Dawn Hour %%
% 
 unqday = unique(floor(date_met - 5/24)); %approx number of local days
% ; 
% 
%     unqday = unqday(~isnan(unqday));
%     dawnhr=nan(size(unqday));
%     duskhr=nan(size(unqday));
%     for count = 1:length(unqday),
%         ind = find(floor(date_met - 5/24) == unqday(count));
%         ind2 = find(Solar(ind) > dawnlevel);
%         if ~isempty(ind2),
%             %        dawnhr(count) = Hour(ind(ind2(1))) - 1 - 4;  %local time from May onwards through summer, 1 h before dawnlevel?
%             %        duskhr(count) = Hour(ind(ind2(end))) - 4 + 1;  %next hour after end of light
%             dawnhr(count) = Hour(ind(ind2(1))) - 1;  %1 h before dawnlevel?
%             duskhr(count) = Hour(ind(ind2(end))) + 1;  %next hour after end of light
%             if dawnhr(count) == -1, keyboard, end;
%         else
%             dawnhr(count) = NaN;
%             duskhr(count) = NaN;
%         end;
%     end;
% 
%      yrdays=find_yearday(unqday);
%     temp_dawn(yrdays)=dawnhr;
%     temp_dusk(yrdays)=duskhr;
% 
% dawn_median=floor(nanmedian(temp_dawn));
% ii=find(dawn_median <8);
% dawn_median(ii)=dawn_median(ii-1);
% 
% dusk_median=floor(nanmedian(temp_dusk,2));
% %fix for dusk data in middle of summer:
% ii= dusk_median == 1;
% dusk_median(ii)=25; %1 corresponds to first hour of next day :)
% ii=find(dusk_median < 20); 
% dusk_median(ii)=dusk_median(ii-1); %change value to near neighbor
% 

days=(datenum(['1-1-' num2str(year2do)]):datenum(['12-31-' num2str(year2do)]))';
duskhr=nan(length(days),1);
dawnhr=nan(length(days),1);
lighttime_gaps=zeros(length(days),1);

for j=1:length(days)
    
    ind = find(floor(date_met - 5/24) == days(j));
    if ~isempty(ind)
        
        ind2=find(date_met(ind) >= days(j)+(dawn_median(j)-1)/24 & date_met(ind) <= days(j)+(dusk_median(j)+1)/24); %examine light only within the expected dawn-dusk range (sometimes renegard nightime noise!)
        ind3 = find(Solar(ind(ind2)) > dawnlevel); %is there light data?
        
        if ~isempty(ind3) %good, we have light data!
            
            tempdawn=Hour(ind(ind2(ind3(1)))) - 1; %1 h before dawnlevel?
            tempdusk=Hour(ind(ind2(ind3(end)))) + 1; %next hour after end of light
            if tempdusk < 20, tempdusk=tempdusk+24; end %i.e. for hours 1,2 - this is the next day...
            
            %check for gaps around dawn and dusk:           
            if abs(dawn_median(j)-tempdawn) <= 1
                dawnhr(j)=tempdawn; %good, accept this dawn!   
                if ind(ind2(ind3(1)))==ind(1)   %if daylight hour is first hour of day, indicates a gap during the night            
                    date_met = [date_met; days(j)+(tempdawn+1)/24];
                    Solar = [Solar; 0];
                end
            elseif abs(dawn_median(j)-tempdawn) <= 3 %close...small gap, pad with a datapoint for later interpolation!
                 date_met = [date_met; days(j)+(dawn_median(j)+1)/24];
                 Solar = [Solar; 0];
                 dawnhr(j)=dawn_median(j);
                 %the plus 1/24 is important - in setup_days, the
                 %script looks only past the dawn hour, not before...
            else
                lighttime_gaps(j)=1; %gap around dawn
            end
                        
            if abs(dusk_median(j)-tempdusk) <= 1
                duskhr(j)=tempdusk; %good, accept this dusk! 
                %sometimes can have good data close to dusk, but no data
                %further on...check this and then add 0 datapoint if need be:
                if ind(ind2(ind3(end)))==ind(end)
                    date_met = [date_met; days(j)+(tempdusk+0.001)/24]; %add one hour after dusk
                    Solar = [Solar; 0];
                end
            elseif abs(dusk_median(j)-tempdusk) <= 3 %close...small gap, pad with a datapoint for later interpolation!
                 date_met = [date_met; days(j)+(dusk_median(j)+0.001)/24];
                 Solar = [Solar; 0];
                 duskhr(j)=dusk_median(j);
            else
                lighttime_gaps(j)=3; %gap around dusk
            end
            
            %now check for any gaps during light time hours:
            if ~isnan(dawnhr(j)) && ~isnan(duskhr(j))
                if any(find(diff(date_met(ind(ind2))) >= 3/24));
                    lighttime_gaps(j) = 2;
                end;
            end
            
            % a quick double check on days 
            %QC plotting if needed:
%             if ~isnan(lighttime_gaps(j))
%                 hold on
%                 plot(date_met0(ind),Solar0(ind),'ko')
%                 plot(date_met0(ind(ind2)),Solar0(ind(ind2)),'r*')              
%                 keyboard
%             end
        else
           lighttime_gaps(j) = 4;   
        end
    else %if no data...flag this too!
        lighttime_gaps(j) = 4;
    end
end

%in case added any points:
[~,ii]=sort(date_met);
date_met=date_met(ii);
Solar=Solar(ii);

dawn = [days dawnhr];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  check for nighttime noise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

date_met_local=date_met-5/24; %for easier handling

for j=1:length(days)
    
    day=days(j);
    ii=find(date_met_local >= day & date_met_local <= day+1);
    
    %find any light levels during expected nightime hours:
    nn=find(date_met_local(ii) <= day+dawn_median(j)/24-5/24 | date_met_local(ii) >= day+dusk_median(j)/24-5/24); %min normal dawn and max normal dusk

    if ~isempty(nn)
        if any(Solar(ii(nn)) > 20)
            disp(['Found "light" readings during dark period for day: ' num2str(day) ': ' datestr(day) ' correcting...'])
            Solar(ii(nn)) = 0;            
            %             plot(date_met_local(ii(nn)),Solar(ii(nn)),'.-')
            %             xlim([day-3 day+3])
            %             datetick('x','mm dd','keeplimits')
            %             pause
            
        end
    end
    
    
end

%negative noise, simply remove:
jj=find(Solar < 0);
disp(['replacing ' num2str(length(jj)) ' negative values with 0'])
Solar(jj)=0;



%% %%%%%%%%%%%%%%%%%%%%%%   AND PLOT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure
    %colorblock expected nightime:
    for j=1:366
        day=datenum(['1-0-' num2str(year2do)])+j;
        %color in night:
        f1=fill([day-1+dusk_median(j)/24; day+dawn_median(j)/24; day+dawn_median(j)/24; day-1+dusk_median(j)/24],[0 0 1000 1000],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    
    %plot the data!
    h1=plot(date_met,Solar,'-','color',[0.0265 0.6137 0.8135]);
    set(gca, 'layer', 'top')

    % add dawn lines:
    for i=1:length(dawn)
        if isnan(dawn(i,2))
           plot(dawn(i,1)+(dawn_median(i)+10)/24,800,'x','linewidth',2,'color',[0 0 0]) %first matlab default color: [0 0.5 0.8]
        else
           line([dawn(i,1)+dawn(i,2)/24 dawn(i,1)+dawn(i,2)/24],[0 1000],'linewidth',2,'color',[0.8 0.5 0]) %first matlab default color: [0 0.5 0.8]
        end
    end
    
    if exist('buoy_added','var')
        hold on
        h3=plot(buoy_added(:,1),buoy_added(:,2),'o','linewidth',2,'markersize',4,'color',[0 0 0.7]);
       % legend([h1(1); h2(1); h3(1)],'Final data','MVCO data','Original buoy data','location','NorthOutside')
        title([num2str(year2do) ' Solar data: ' num2str(length(find(~isnan(buoy_days)))) ' days added from buoy data'])
    else
       % legend([h1(1); h2(1)],'Final data','MVCO data','location','NorthOutside')
        title(num2str(year2do))
    end
    


xlim([datenum(['1-1-' num2str(year2do)]) datenum(['12-31-' num2str(year2do)])])
ylim([-10 max(Solar)])
datetick('x','keeplimits')



%% %SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


eval(['save ' solarsavepath 'solar' num2str(year2do) '.mat Solar date_met dawn'])















