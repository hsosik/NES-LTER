function [vr_scales] = tsvr_scale(ts_pops,thresholds)

% Created by Silke van Daalen, January 2023, last edit 27-Jan-2023
% Adjusting the timescale-specific variance ratio from returning only two
% timespans, as in Zhao et al, 2020, and Shoemaker et al., 2022
% This function computes a timescale-specific variance ratio, at several scales,
% it requires:
%   - ts_pops should contain a timeseries of all populations, with pop identity
%   across rows and time across columns
%   - thresholds should provide the threshold values by which we demarcate
%   the scales of interest (i.e. [7 30 365] provides a within-week(<7),
%   week-to-month (7-30), month-to-year (30-365), and interannual variance ratio (>365))
%   Currently, the maximum nr. of thresholds is 6

npops=size(ts_pops,2);
ts_tot=sum(ts_pops,2);

[cosp_tot,freq]=cospectrum(ts_tot');
tslength=length(ts_tot);
numer=zeros(tslength-1,1);
numer(:,1)=cosp_tot(1,1,2:end);
freq=freq(:,2:end);
all_spects=zeros(length(ts_pops)-1,npops);
for k=1:npops
    [h_cosp,h_freq]=cospectrum(ts_pops(:,k)');
    all_spects(:,k)=h_cosp(1,1,2:end);
end
denom=sum(all_spects,2);
fsvr_classic=numer./denom; % this is where the rounding errors happen

ts=1./freq;
% code in R:tsvreqclassic
mutot=mean(ts_tot);
cv2_com=numer/(mutot^2);
cv2_comip=denom/(mutot^2);


[cospec,freq]=cospectrum(ts_pops');
lenfreq=length(freq);
freq=freq(:,2:lenfreq);
cospec=cospec(:,:,2:lenfreq);
tsvar=var(ts_pops);
sdg_cospec=zeros(lenfreq-1,1);
for i=1:(lenfreq-1)
    cospec_i=cospec(:,:,i);
    sdg_cospec(i,1)=sum(diag(cospec_i));
end
strength=sdg_cospec/sum(tsvar);


vr_ray=[ts' cv2_com cv2_comip fsvr_classic strength];
output=flip(vr_ray);

nthresh=size(thresholds,2);


% now we get to R:aggts, input = [output, ts_data_long_or_short]
ts_short=output;
ts_short(ts_short(:,1)>=thresholds(:,1),:)=[];
ts_long=output;
ts_long(ts_long(:,1)<thresholds(:,end),:)=[];
% throw out all ts lower than 2, only for short ts
ts_short(ts_short(:,1)<2,:)=[];

ts_medium=zeros(size(output,1),size(output,2),(nthresh-1));


for i=1:(nthresh-1)

    ts_med_temp=output;
    ts_med_temp(ts_med_temp(:,1)<thresholds(:,i),:)=[];
    ts_med_temp(ts_med_temp(:,1)>=thresholds(:,i+1),:)=[];

    ts_medium(1:length(ts_med_temp),:,i)=ts_med_temp;
    
end


% quick fix for now
if nthresh==1

    [~,ind_s]=ismember(ts_short(:,1),output(:,1));
[~,ind_l]=ismember(ts_long(:,1),output(:,1));

T=max(output(:,1));
midind=T/2;

inds_s=round(unique([-ind_s+2*midind; ind_s]));
inds_l=round(unique([-ind_l+2*midind; ind_l]));

ts_short=output(inds_s,:);
vr_short=ts_short(:,4)'*ts_short(:,5)/sum(ts_short(:,5));

ts_long=output(inds_l,:);
vr_long=ts_long(:,4)'*ts_long(:,5)/sum(ts_long(:,5));
vr_scales=[vr_short vr_long];

elseif nthresh==2

    med1=ts_medium(:,:,1);
    med1=med1(any(med1,2),:);

    [~,ind_s]=ismember(ts_short(:,1),output(:,1));
[~,ind_l]=ismember(ts_long(:,1),output(:,1));
[~,ind_m1]=ismember(med1(:,1),output(:,1));

T=max(output(:,1));
midind=T/2;

inds_s=round(unique([-ind_s+2*midind; ind_s]));
inds_l=round(unique([-ind_l+2*midind; ind_l]));
inds_m1=round(unique([-ind_m1+2*midind; ind_m1]));

ts_short=output(inds_s,:);
vr_short=ts_short(:,4)'*ts_short(:,5)/sum(ts_short(:,5));

ts_long=output(inds_l,:);
vr_long=ts_long(:,4)'*ts_long(:,5)/sum(ts_long(:,5));

ts_med1=output(inds_m1,:);
vr_med1=ts_med1(:,4)'*ts_med1(:,5)/sum(ts_med1(:,5));

vr_scales=[vr_short vr_med1 vr_long];

elseif nthresh==3
    med1=ts_medium(:,:,1);
    med1=med1(any(med1,2),:);
        med2=ts_medium(:,:,2);
    med2=med2(any(med2,2),:);

    [~,ind_s]=ismember(ts_short(:,1),output(:,1));
[~,ind_l]=ismember(ts_long(:,1),output(:,1));
[~,ind_m1]=ismember(med1(:,1),output(:,1));
[~,ind_m2]=ismember(med2(:,1),output(:,1));

T=max(output(:,1));
midind=T/2;

inds_s=round(unique([-ind_s+2*midind; ind_s]));
inds_l=round(unique([-ind_l+2*midind; ind_l]));
inds_m1=round(unique([-ind_m1+2*midind; ind_m1]));
inds_m2=round(unique([-ind_m2+2*midind; ind_m2]));

ts_short=output(inds_s,:);
vr_short=ts_short(:,4)'*ts_short(:,5)/sum(ts_short(:,5));

ts_long=output(inds_l,:);
vr_long=ts_long(:,4)'*ts_long(:,5)/sum(ts_long(:,5));

ts_med1=output(inds_m1,:);
vr_med1=ts_med1(:,4)'*ts_med1(:,5)/sum(ts_med1(:,5));

ts_med2=output(inds_m2,:);
vr_med2=ts_med2(:,4)'*ts_med2(:,5)/sum(ts_med2(:,5));

vr_scales=[vr_short vr_med1 vr_med2 vr_long];

elseif nthresh==4

        med1=ts_medium(:,:,1);
    med1=med1(any(med1,2),:);
        med2=ts_medium(:,:,2);
    med2=med2(any(med2,2),:);
        med3=ts_medium(:,:,3);
    med3=med3(any(med3,2),:);

    [~,ind_s]=ismember(ts_short(:,1),output(:,1));
[~,ind_l]=ismember(ts_long(:,1),output(:,1));
[~,ind_m1]=ismember(med1(:,1),output(:,1));
[~,ind_m2]=ismember(med2(:,1),output(:,1));
[~,ind_m3]=ismember(med3(:,1),output(:,1));

T=max(output(:,1));
midind=T/2;

inds_s=round(unique([-ind_s+2*midind; ind_s]));
inds_l=round(unique([-ind_l+2*midind; ind_l]));
inds_m1=round(unique([-ind_m1+2*midind; ind_m1]));
inds_m2=round(unique([-ind_m2+2*midind; ind_m2]));
inds_m3=round(unique([-ind_m3+2*midind; ind_m3]));

ts_short=output(inds_s,:);
vr_short=ts_short(:,4)'*ts_short(:,5)/sum(ts_short(:,5));

ts_long=output(inds_l,:);
vr_long=ts_long(:,4)'*ts_long(:,5)/sum(ts_long(:,5));

ts_med1=output(inds_m1,:);
vr_med1=ts_med1(:,4)'*ts_med1(:,5)/sum(ts_med1(:,5));

ts_med2=output(inds_m2,:);
vr_med2=ts_med2(:,4)'*ts_med2(:,5)/sum(ts_med2(:,5));

ts_med3=output(inds_m3,:);
vr_med3=ts_med3(:,4)'*ts_med3(:,5)/sum(ts_med3(:,5));


vr_scales=[vr_short vr_med1 vr_med2 vr_med3 vr_long];

elseif nthresh==5

        med1=ts_medium(:,:,1);
    med1=med1(any(med1,2),:);
        med2=ts_medium(:,:,2);
    med2=med2(any(med2,2),:);
        med3=ts_medium(:,:,3);
    med3=med3(any(med3,2),:);
        med4=ts_medium(:,:,4);
    med4=med4(any(med4,2),:);

    % the supp info of Zhao et al and the R-script tells me these aliasing frequencies
% below 2 should be included, and so, we shall include them
% but they need to reflect the short/long divide

[~,ind_s]=ismember(ts_short(:,1),output(:,1));
[~,ind_l]=ismember(ts_long(:,1),output(:,1));
[~,ind_m1]=ismember(med1(:,1),output(:,1));
[~,ind_m2]=ismember(med2(:,1),output(:,1));
[~,ind_m3]=ismember(med3(:,1),output(:,1));
[~,ind_m4]=ismember(med4(:,1),output(:,1));

T=max(output(:,1));
midind=T/2;

inds_s=round(unique([-ind_s+2*midind; ind_s]));
inds_l=round(unique([-ind_l+2*midind; ind_l]));
inds_m1=round(unique([-ind_m1+2*midind; ind_m1]));
inds_m2=round(unique([-ind_m2+2*midind; ind_m2]));
inds_m3=round(unique([-ind_m3+2*midind; ind_m3]));
inds_m4=round(unique([-ind_m4+2*midind; ind_m4]));

ts_short=output(inds_s,:);
vr_short=ts_short(:,4)'*ts_short(:,5)/sum(ts_short(:,5));

ts_long=output(inds_l,:);
vr_long=ts_long(:,4)'*ts_long(:,5)/sum(ts_long(:,5));

ts_med1=output(inds_m1,:);
vr_med1=ts_med1(:,4)'*ts_med1(:,5)/sum(ts_med1(:,5));

ts_med2=output(inds_m2,:);
vr_med2=ts_med2(:,4)'*ts_med2(:,5)/sum(ts_med2(:,5));

ts_med3=output(inds_m3,:);
vr_med3=ts_med3(:,4)'*ts_med3(:,5)/sum(ts_med3(:,5));

ts_med4=output(inds_m4,:);
vr_med4=ts_med4(:,4)'*ts_med4(:,5)/sum(ts_med4(:,5));

vr_scales=[vr_short vr_med1 vr_med2 vr_med3 vr_med4 vr_long];

elseif nthresh==6

        med1=ts_medium(:,:,1);
    med1=med1(any(med1,2),:);
        med2=ts_medium(:,:,2);
    med2=med2(any(med2,2),:);
        med3=ts_medium(:,:,3);
    med3=med3(any(med3,2),:);
        med4=ts_medium(:,:,4);
    med4=med4(any(med4,2),:);
        med5=ts_medium(:,:,5);
    med5=med5(any(med5,2),:);

    % the supp info of Zhao et al and the R-script tells me these aliasing frequencies
% below 2 should be included, and so, we shall include them
% but they need to reflect the short/long divide

[~,ind_s]=ismember(ts_short(:,1),output(:,1));
[~,ind_l]=ismember(ts_long(:,1),output(:,1));
[~,ind_m1]=ismember(med1(:,1),output(:,1));
[~,ind_m2]=ismember(med2(:,1),output(:,1));
[~,ind_m3]=ismember(med3(:,1),output(:,1));
[~,ind_m4]=ismember(med4(:,1),output(:,1));
[~,ind_m5]=ismember(med5(:,1),output(:,1));

T=max(output(:,1));
midind=T/2;

inds_s=round(unique([-ind_s+2*midind; ind_s]));
inds_l=round(unique([-ind_l+2*midind; ind_l]));
inds_m1=round(unique([-ind_m1+2*midind; ind_m1]));
inds_m2=round(unique([-ind_m2+2*midind; ind_m2]));
inds_m3=round(unique([-ind_m3+2*midind; ind_m3]));
inds_m4=round(unique([-ind_m4+2*midind; ind_m4]));
inds_m5=round(unique([-ind_m5+2*midind; ind_m5]));

ts_short=output(inds_s,:);
vr_short=ts_short(:,4)'*ts_short(:,5)/sum(ts_short(:,5));

ts_long=output(inds_l,:);
vr_long=ts_long(:,4)'*ts_long(:,5)/sum(ts_long(:,5));

ts_med1=output(inds_m1,:);
vr_med1=ts_med1(:,4)'*ts_med1(:,5)/sum(ts_med1(:,5));

ts_med2=output(inds_m2,:);
vr_med2=ts_med2(:,4)'*ts_med2(:,5)/sum(ts_med2(:,5));

ts_med3=output(inds_m3,:);
vr_med3=ts_med3(:,4)'*ts_med3(:,5)/sum(ts_med3(:,5));

ts_med4=output(inds_m4,:);
vr_med4=ts_med4(:,4)'*ts_med4(:,5)/sum(ts_med4(:,5));

ts_med5=output(inds_m5,:);
vr_med5=ts_med5(:,4)'*ts_med5(:,5)/sum(ts_med5(:,5));

vr_scales=[vr_short vr_med1 vr_med2 vr_med3 vr_med4 vr_med5 vr_long];


else
warning('The number of thresholds is higher than six')
end



return