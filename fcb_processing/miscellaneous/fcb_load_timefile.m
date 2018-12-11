clear all
close all
year=2012;
%pathname='/home/kristen/Desktop/MVCO_FCB_time_testfiles/';
%eval(['filelist=dir(''/home/kristen/Desktop/MVCO_FCB_time_testfiles/*' num2str(year) '*.mat'')'])

pathname='/Volumes/Lab_data/MVCO/FCB/MVCO_Jan2012/data/processed_July2016/time/';
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