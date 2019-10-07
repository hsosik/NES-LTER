function [to_use, to_use_SSC, to_use_PE]=exclude_data(cellresultsall,year2do)
%data to exclude for now (and ideally fix later)

%Issues include:
%beads classified as Syn
%merging issues due to baseline wobble
%just bad signatures
%strange acquisiiton times - not sure what exactly these are...
%syringe clogged - bad volume estimates

to_use0=find(~isnan(cellresultsall(:,1)));
cellresultsall=cellresultsall(to_use0,:);

switch year2do
    case 2003
        
        ex_timepts={};     

        ii=find(cellresultsall(:,1)>=datenum('06-08-03 22:00:00') & cellresultsall(:,1)<=datenum('06-08-03 23:00:00'));
        reason=repmat('Syn classified as picoeuks!',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('05-Dec-2003 12:00:00') & cellresultsall(:,1)<=datenum('06-Dec-2003 16:00:00'));
        reason=repmat('smeared bead ssc',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
               
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
    
    case 2004
        %Nothing! :)
        to_use=find(~isnan(cellresultsall(:,1)));
        
    case 2005
        
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('4-20-05 01:00:00') & cellresultsall(:,1)<=datenum('4-20-05 03:00:00'));
        reason=repmat('Beads classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
                
         ii=find(cellresultsall(:,1)>=datenum('10-24-05 22:00:00') & cellresultsall(:,1)<=datenum('10-25-05 00:00:00'));
        reason=repmat('maybe some intrustment start up weirdness',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];

        ii=find(cellresultsall(:,1)>=datenum('10-30-05 23:00:00') & cellresultsall(:,1)<=datenum('10-31-05 00:00:00'));
        reason=repmat('Too long acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];

        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);       

    case  2006
        
        %Nothing! :)
        %Although one weird SSC point for mode on '20-Jun-2006 11:29:55',
        %mean is fine - very, very strange....
        to_use=find(~isnan(cellresultsall(:,1)));
        
    case  2007
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('5-19-07 11:00:00') & cellresultsall(:,1)<=datenum('5-19-07 14:00:00'));
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('12-06-07 20:00:00') & cellresultsall(:,1)<=datenum('12-06-07 21:00:00'));
        reason=repmat('Beads classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];      
        
        ii=find(cellresultsall(:,1)>=datenum('12-13-07 17:00:00') & cellresultsall(:,1)<=datenum('12-13-07 18:00:00'));
        reason=repmat('Missed Syn patch',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
              
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2008
        %Nothing! :)
        to_use=find(~isnan(cellresultsall(:,1)));
        %maybe something suspciuous about restart on 10/26, but not out of
        %the ordinary, cytograms look fine
        
    case 2009
        ex_timepts={};
              
        i1=find(cellresultsall(:,1)>=datenum('June-17-09 22:00:00') & cellresultsall(:,1)<=datenum('June-18-09 06:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('June-18-09 08:00:00') & cellresultsall(:,1)<=datenum('6-19-09 05:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('Jul-1-09 02:00:00') & cellresultsall(:,1)<=datenum('7-1-09 03:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('Jul-7-09 05:00:00') & cellresultsall(:,1)<=datenum('7-8-09 06:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('Jul-11-09 07:00:00') & cellresultsall(:,1)<=datenum('7-8-09 08:00:00'));
        i6=find(cellresultsall(:,1)>=datenum('Jul-12-09 03:00:00') & cellresultsall(:,1)<=datenum('7-12-09 10:00:00'));
        i7=find(cellresultsall(:,1)>=datenum('Jul-19-09 02:00:00') & cellresultsall(:,1)<=datenum('7-19-09 04:00:00'));
        i8=find(cellresultsall(:,1)>=datenum('Jul-21-09 14:00:00') & cellresultsall(:,1)<=datenum('7-21-09 16:00:00'));
        
        ii=[i1;i2;13;i4;i5;i6;i7;i8];
        reason=repmat('Noise classified as Syn :(',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
    
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2010
        
        %Possibly some dying syn around Dec 5-6, where distinct lower PE clusters
        %appear...
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('6-22-10 15:00:00') & cellresultsall(:,1)<=datenum('6-22-10 18:30:00'));
        reason=repmat('A bad syringe left in?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('8-22-10 6:00:00') & cellresultsall(:,1)<=datenum('8-22-10 07:00:00'));
        reason=repmat('Syn classified as picoeuks',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];

        
        i1=find(cellresultsall(:,1)>=datenum('12-17-10 05:00:00') & cellresultsall(:,1)<=datenum('12-17-10 06:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('12-16-10 10:00:00') & cellresultsall(:,1)<=datenum('12-16-10 11:00:00'));
        ii=[i1;i2];
        reason=repmat('Strange acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];

        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2011
        
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('7-12-11 22:00:00') & cellresultsall(:,1)<=datenum('7-14-11 00:00:00'));      
        reason=repmat('Merging issues',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('7-31-11 00:00:00') & cellresultsall(:,1)<=datenum('8-08-11 00:00:00'));
        reason=repmat('Ugh...bad SSC spread coupled with classification isues',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('8-09-11 03:00:00') & cellresultsall(:,1)<=datenum('8-09-11 4:00:00'));      
        i2=find(cellresultsall(:,1)>=datenum('8-10-11 03:00:00') & cellresultsall(:,1)<=datenum('8-10-11 4:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('8-11-11 04:00:00') & cellresultsall(:,1)<=datenum('8-11-11 5:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('8-11-11 14:00:00') & cellresultsall(:,1)<=datenum('8-11-11 17:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('8-12-11 19:00:00') & cellresultsall(:,1)<=datenum('8-12-11 20:00:00'));
        ii=[i1;i2;i3;i4;i5];
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];   
        
        ii=find(cellresultsall(:,1)>=datenum('10-19-11 06:00:00') & cellresultsall(:,1)<=datenum('10-19-11 7:00:00'));
        reason=repmat('Syn patch classified as picoeuks',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('12-15-11 10:00:00') & cellresultsall(:,1)<=datenum('12-15-11 12:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('12-15-11 12:00:00') & cellresultsall(:,1)<=datenum('12-15-11 14:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('12-15-11 15:00:00') & cellresultsall(:,1)<=datenum('12-15-11 18:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('12-15-11 20:00:00') & cellresultsall(:,1)<=datenum('12-15-11 22:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('12-16-11 01:30:00') & cellresultsall(:,1)<=datenum('12-31-11 11:59:00'));
        ii=[i1;i2;i3;i4;i5];
        reason=repmat('Bad SSC data',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2012
        
        ex_timepts={};
        ii=find(cellresultsall(:,1)>=datenum('Jan-09-12 20:00:00') & cellresultsall(:,1)<=datenum('Jan-09-12 22:00:00')); %start until 22:00 (Jan-9-2012)
        reason=repmat('Strang acquisition time pattern?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
               
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2013
        
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('18-Mar-2013 14:00:00') & cellresultsall(:,1)<=datenum('18-Mar-2013 16:00:00'));
        reason=repmat('instrument start weirdness?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
               
        %note- there are more instances where classfier is grabbing noise, these
        %just looked the most severe:
        i1=find(cellresultsall(:,1)>=datenum('3-27-13 19:00:00') & cellresultsall(:,1)<=datenum('3-27-13 22:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('4-3-13 06:00:00') & cellresultsall(:,1)<=datenum('4-3-13 08:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('4-3-13 10:00:00') & cellresultsall(:,1)<=datenum('4-3-13 11:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('4-3-13 20:00:00') & cellresultsall(:,1)<=datenum('4-3-13 22:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('4-4-13 09:00:00') & cellresultsall(:,1)<=datenum('4-4-13 11:00:00'));
        i6=find(cellresultsall(:,1)>=datenum('4-4-13 15:00:00') & cellresultsall(:,1)<=datenum('4-4-13 17:00:00'));
        i7=find(cellresultsall(:,1)>=datenum('4-5-13 09:00:00') & cellresultsall(:,1)<=datenum('4-5-13 11:00:00'));
        i8=find(cellresultsall(:,1)>=datenum('4-5-13 18:00:00') & cellresultsall(:,1)<=datenum('4-5-13 19:00:00'));
        i9=find(cellresultsall(:,1)>=datenum('4-6-13 18:00:00') & cellresultsall(:,1)<=datenum('4-6-13 19:00:00'));
        i10=find(cellresultsall(:,1)>=datenum('4-7-13 19:00:00') & cellresultsall(:,1)<=datenum('4-7-13 22:00:00'));      
        i11=find(cellresultsall(:,1)>=datenum('4-8-13 01:00:00') & cellresultsall(:,1)<=datenum('4-8-13 03:00:00'));
        i12=find(cellresultsall(:,1)>=datenum('4-9-13 09:00:00') & cellresultsall(:,1)<=datenum('4-9-13 10:00:00'));
        i13=find(cellresultsall(:,1)>=datenum('4-9-13 13:00:00') & cellresultsall(:,1)<=datenum('4-9-13 14:00:00'));
        i14=find(cellresultsall(:,1)>=datenum('4-9-13 19:00:00') & cellresultsall(:,1)<=datenum('4-9-13 20:00:00'));        
        i15=find(cellresultsall(:,1)>=datenum('4-12-13 9:00:00') & cellresultsall(:,1)<=datenum('4-12-13 10:00:00'));       
        i16=find(cellresultsall(:,1)>=datenum('4-12-13 18:00:00') & cellresultsall(:,1)<=datenum('4-12-13 20:00:00'));
        i17=find(cellresultsall(:,1)>=datenum('4-13-13 07:00:00') & cellresultsall(:,1)<=datenum('4-13-13 08:00:00'));
        i18=find(cellresultsall(:,1)>=datenum('4-17-13 15:00:00') & cellresultsall(:,1)<=datenum('4-17-13 16:00:00'));
        i19=find(cellresultsall(:,1)>=datenum('4-27-13 16:00:00') & cellresultsall(:,1)<=datenum('4-27-13 17:00:00'));
        i20=find(cellresultsall(:,1)>=datenum('4-28-13 09:00:00') & cellresultsall(:,1)<=datenum('4-28-13 10:00:00'));
        i21=find(cellresultsall(:,1)>=datenum('5-18-13 17:00:00') & cellresultsall(:,1)<=datenum('5-18-13 19:00:00'));
        i22=find(cellresultsall(:,1)>=datenum('5-19-13 11:00:00') & cellresultsall(:,1)<=datenum('5-19-13 13:00:00'));
        i23=find(cellresultsall(:,1)>=datenum('3-29-13 20:00:00') & cellresultsall(:,1)<=datenum('3-27-13 21:00:00'));
        
        ii=[i1;i2;i3;i4;i5;i6;i7;i8;i9;i10;i11;i12;i13;i14;i15;i16;i17;i18;i19;i20;i21;i22;i23];
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
             
        ii=find(cellresultsall(:,1)>=datenum('5-27-13 00:00:00') & cellresultsall(:,1)<=datenum('5-27-13 01:00:00'));
        reason=repmat('Part of Syn patch classified as picoeuks',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('8-02-13 14:00:00') & cellresultsall(:,1)<=datenum('08-02-13 20:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('8-08-13 12:00:00') & cellresultsall(:,1)<=datenum('08-08-13 14:00:00'));
        ii=[i1;i2];
        reason=repmat('Bad syringe left after QC?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];

        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2014
        ex_timepts={};
        
        i1=find(cellresultsall(:,1)>=datenum('1-06-14 17:00:00') & cellresultsall(:,1)<=datenum('1-06-14 20:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('1-10-14 11:00:00') & cellresultsall(:,1)<=datenum('1-10-14 12:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('5-09-14 10:00:00') & cellresultsall(:,1)<=datenum('5-09-14 11:00:00'));
        ii=[i1;i2;i3];
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('1-08-14 00:00:00') & cellresultsall(:,1)<=datenum('1-08-14 20:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('Mar-04-14 19:00:00') & cellresultsall(:,1)<=datenum('Mar-04-14 21:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('May-11-14 03:00:00') & cellresultsall(:,1)<=datenum('May-11-14 04:00:00'));
        ii=[i1;i2;i3];
        reason=repmat('A few bad syringes left in?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('17-Jan-2014 19:00:00') & cellresultsall(:,1)<=datenum('01-Mar-2014 21:00:00')); %from 19:00 - 20:30 - syr pump doesn't reach max
        reason=repmat('This data should probably be added to skip file....',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)]; 
        
        i1=find(cellresultsall(:,1)>=datenum('Apr-03-14 00:00:00') & cellresultsall(:,1)<=datenum('Apr-03-14 01:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('Apr-04-14 23:00:00') & cellresultsall(:,1)<=datenum('Apr-05-14 00:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('Apr-05-14 21:00:00') & cellresultsall(:,1)<=datenum('Apr-05-14 22:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('Apr-06-14 08:00:00') & cellresultsall(:,1)<=datenum('Apr-06-14 09:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('Apr-07-14 05:00:00') & cellresultsall(:,1)<=datenum('Apr-07-14 06:00:00'));
        i6=find(cellresultsall(:,1)>=datenum('Apr-07-14 22:00:00') & cellresultsall(:,1)<=datenum('Apr-07-14 23:00:00'));
        ii=[i1;i2;i3;i4;i5;i6];
        reason=repmat('Bad Syn classification',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('Nov-10-14 19:00:00') & cellresultsall(:,1)<=datenum('Nov-10-14 20:00:00')); %from 19:00 - 20:30 - syr pump doesn't reach max
        reason=repmat('Bad syringe left in?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
               
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
        %PE mode not tracking:        
        i1=find(cellresultsall(:,1)>=datenum('29-Mar-2014 14:00:00') & cellresultsall(:,1)<=datenum('Mar-29-2014 20:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('30-Mar-2014 05:00:00') & cellresultsall(:,1)<=datenum('Mar-30-2014 08:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('31-Mar-2014 07:00:00') & cellresultsall(:,1)<=datenum('Apr-02-2014 12:00:00'));
        ii=[i1;i2;i3];
        reason=repmat('Low mode for PE....too few cells?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
     
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use_PE=find(bt_pts==0);
        
    case 2015
        
        ex_timepts={};
        i1=find(cellresultsall(:,1)>=datenum('3-25-15 15:00:00') & cellresultsall(:,1)<=datenum('3-25-15 17:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('03-Apr-2015 14:00:00') & cellresultsall(:,1)<=datenum('03-Apr-2015 15:00:00'));
        ii=[i1;i2];
        reason=repmat('Beads classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
                   
        ii=find(cellresultsall(:,1)>=datenum('5-17-15 01:00:00') & cellresultsall(:,1)<=datenum('5-17-15 02:00:00'));
        reason=repmat('Bad SSC; bad syringe left in?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('3-30-15 07:00:00') & cellresultsall(:,1)<=datenum('3-30-15 08:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('4-06-15 13:00:00') & cellresultsall(:,1)<=datenum('4-06-15 15:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('6-08-15 07:00:00') & cellresultsall(:,1)<=datenum('6-08-15 08:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('6-11-15 03:00:00') & cellresultsall(:,1)<=datenum('6-11-15 04:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('6-11-15 09:00:00') & cellresultsall(:,1)<=datenum('6-11-15 11:00:00'));
        i6=find(cellresultsall(:,1)>=datenum('6-12-15 11:00:00') & cellresultsall(:,1)<=datenum('6-12-15 13:00:00'));
        i7=find(cellresultsall(:,1)>=datenum('6-12-15 14:00:00') & cellresultsall(:,1)<=datenum('6-12-15 15:00:00'));
        i8=find(cellresultsall(:,1)>=datenum('6-14-15 04:00:00') & cellresultsall(:,1)<=datenum('6-14-15 05:00:00'));
        ii=[i1;i2;i3;i4;i5;i6;i7;i8];
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('30-Jun-2015 01:00:00') & cellresultsall(:,1)<=datenum('30-Jun-2015 06:00:00'));
        reason=repmat('strange ssc drop',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('11-18-15 19:00:00') & cellresultsall(:,1)<=datenum('11-18-15 20:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('11-28-15 19:00:00') & cellresultsall(:,1)<=datenum('11-28-15 20:00:00'));
        ii=[i1;i2];
        reason=repmat('Bad syringes left in?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
        %ugh...abnd estiamtes may be okay, but SSC is just so bad...       
        ii=find(cellresultsall(:,1)>=datenum('20-Jul-2015 00:00:00') & cellresultsall(:,1)<=datenum('11-Sept-2015 00:00:00'));
        reason=repmat('Smeared SSC',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
     
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use_SSC=find(bt_pts==0);
        
        
    case 2016
        
        ex_timepts={};
        
        %Jan 14-15, most bad data removed -some left in good files, remove:
        ii=find(cellresultsall(:,1)>=datenum('1-14-16 22:00:00') & cellresultsall(:,1)<=datenum('1-14-16 23:00:00'));  
        reason=repmat('Just one left over bad syringe',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
                  
        ii=find(cellresultsall(:,1)>=datenum('Apr-5-16 11:00:00') & cellresultsall(:,1)<=datenum('Apr-5-16 12:00:00'));
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('May-9-16 03:00:00') & cellresultsall(:,1)<=datenum('May-9-16 06:00:00'));        
        i2=find(cellresultsall(:,1)>=datenum('5-31-16 02:45:00') & cellresultsall(:,1)<=datenum('5-31-16 04:00:00')); %this one is probably tail end of bad data....
        i3=find(cellresultsall(:,1)>=datenum('6-16-16 12:00:00') & cellresultsall(:,1)<=datenum('6-20-16 00:00:00')); %mess of signatures - THERE MAYBE STILL GOOD FILES IN HERE THOUGH!
        ii=[i1;i2;i3];
        reason=repmat('Syring issues? Strange signatures',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];    

        ii=find(cellresultsall(:,1)>=datenum('Aug-3-2016 19:00:00') & cellresultsall(:,1)<=datenum('Aug-3-2016 23:00:00')); 
        reason=repmat('Strange spread and weird classfications',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];              

        ii=find(cellresultsall(:,1)>=datenum('10-02-16 20:00:00') & cellresultsall(:,1)<=datenum('10-04-16 19:00:00')); 
        reason=repmat('SSC spread - bad signatures',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('10-11-16 05:00:00') & cellresultsall(:,1)<=datenum('10-11-16 15:00:00')); 
        reason=repmat('Bad signatures - spread on both PE and SSC',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];

        ii=find(cellresultsall(:,1)>=datenum('11-08-16 07:00:00') & cellresultsall(:,1)<=datenum('11-08-16 08:00:00')); 
        reason=repmat('Bottom end of Syn cloud classified as picoeuk',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
 
        ii=find(cellresultsall(:,1)>=datenum('12-29-16 20:00:00') & cellresultsall(:,1)<=datenum('12-30-16 10:00:00')); 
        reason=repmat('Small SSC noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
 
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
        
        %ugh...abnd estiamtes should be okay, but SSC is bad...        
        ii=find(cellresultsall(:,1)>=datenum('Sept-22-2016 22:00:00') & cellresultsall(:,1)<=datenum('04-Oct-2016 20:00:00'));
        reason=repmat('Smeared SSC',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
     
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use_SSC=find(bt_pts==0);
        
        
    case 2017
        
        ex_timepts={};
         
        ii=find(cellresultsall(:,1)>=datenum('13-Nov-2017 10:00:00') & cellresultsall(:,1)<=datenum('13-Nov-2017 15:00:00')); 
        reason=repmat('SSC smeared...tending towards 0',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];

        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
                
    case 2018
        %to be determined
        ex_timepts={};
         
         i1=find(cellresultsall(:,1)>=datenum('08-Apr-2018 22:00:00') & cellresultsall(:,1)<=datenum('08-Apr-2018 23:00:00')); 
         i2=find(cellresultsall(:,1)>=datenum('10-Apr-2018 18:00:00') & cellresultsall(:,1)<=datenum('10-Apr-2018 19:00:00'));
         i3=find(cellresultsall(:,1)>=datenum('11-Apr-2018 07:00:00') & cellresultsall(:,1)<=datenum('11-Apr-2018 08:00:00'));
         i4=find(cellresultsall(:,1)>=datenum('11-Apr-2018 11:00:00') & cellresultsall(:,1)<=datenum('11-Apr-2018 14:00:00'));
         i5=find(cellresultsall(:,1)>=datenum('12-Apr-2018 11:00:00') & cellresultsall(:,1)<=datenum('12-Apr-2018 18:00:00'));
         i6=find(cellresultsall(:,1)>=datenum('13-Apr-2018 12:00:00') & cellresultsall(:,1)<=datenum('13-Apr-2018 18:00:00'));
         i7=find(cellresultsall(:,1)>=datenum('14-Apr-2018 01:00:00') & cellresultsall(:,1)<=datenum('14-Apr-2018 05:00:00'));
         i8=find(cellresultsall(:,1)>=datenum('16-Apr-2018 14:00:00') & cellresultsall(:,1)<=datenum('16-Apr-2018 15:00:00'));
         i9=find(cellresultsall(:,1)>=datenum('18-Apr-2018 19:00:00') & cellresultsall(:,1)<=datenum('18-Apr-2018 20:00:00'));
         i10=find(cellresultsall(:,1)>=datenum('19-Apr-2018 00:00:00') & cellresultsall(:,1)<=datenum('19-Apr-2018 01:00:00'));
         i11=find(cellresultsall(:,1)>=datenum('19-Apr-2018 09:00:00') & cellresultsall(:,1)<=datenum('19-Apr-2018 11:00:00'));
         i12=find(cellresultsall(:,1)>=datenum('20-Apr-2018 06:00:00') & cellresultsall(:,1)<=datenum('20-Apr-2018 09:00:00'));
         i13=find(cellresultsall(:,1)>=datenum('20-Apr-2018 19:00:00') & cellresultsall(:,1)<=datenum('21-Apr-2018 04:00:00'));
         i14=find(cellresultsall(:,1)>=datenum('21-Apr-2018 10:00:00') & cellresultsall(:,1)<=datenum('24-Apr-2018 18:00:00'));  %Numbers are just too low for good classification with current algorithm
        ii=[i1;i2;i3;i4;i5;i6;i7;i8;i9;i10;i11;i12;i13;i14];
         reason=repmat('Noise classified as Syn; numbers are just too low for good classfication',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
 
         ii=find(cellresultsall(:,1)>=datenum('Apr-09-2018 12:00:00') & cellresultsall(:,1)<=datenum('Apr-09-2018 13:00:00')); 
        reason=repmat('Bad syringe from restart?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('31-Jul-2018 00:00:00') & cellresultsall(:,1)<=datenum('10-Aug-2018 15:00:00')); 
        reason=repmat('Bad SSC - smeared and tends to 0',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
      
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
        %Abnd estiamtes should be okay, but SSC is just so bad here too...
    
        i1=find(cellresultsall(:,1)>=datenum('July-24-2018 08:00:00') & cellresultsall(:,1)<=datenum('July-24-2018 10:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('July-28-2018 21:00:00') & cellresultsall(:,1)<=datenum('Aug-1-2018 00:00:00'));
        ii=[i1;i2];
        reason=repmat('Noise is bulk of SSC',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
     
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use_SSC=find(bt_pts==0);
end

to_use=to_use0(to_use);

if exist('to_use_SSC','var')
    to_use_SSC=to_use0(to_use_SSC);
else
    to_use_SSC=to_use;
end

if exist('to_use_PE','var')
    to_use_PE=to_use0(to_use_PE);
else
    to_use_PE=to_use;
end

if any(isnan(cellresultsall(to_use,1)))
    keyboard
end
