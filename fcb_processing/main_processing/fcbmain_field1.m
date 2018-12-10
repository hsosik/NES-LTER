%Main proogram to set up batch processing of FCB data
%this version is for lab data with separate beads, 3/03 Heidi  ??
%this version is for field data with integrated beads, 3/03 Heidi
%fcbmain_field1.m revise from fcbmain_field.m (2011 case) to set up cases for each year

clear all, close all
warning off 

%%USER CHANGE - below here
for year2do = 2005 %[2010 2011 2013 2014] %2003:2004 %2011, 2005
dotime = 1; %0 = NO, 1 = YES
domerge = 0;
doclassify = 0;
doplotgroup = 0;
docells = 1;
dobeads = 0; %ALWAYS MERGE CELLS BEFORE CORRESPONDING BEADS
timeplotflag = 0; %for time: 0 = no plots, 1 = plots
mergeplotflag = 0; %for merge: 0 = no plots, 1 = plots
classplotflag = 0; %for classify: 0 = no plots, 1 = plots
movieflag = 0;
beadmovieflag=0;
syrplotflag=1;
%%USER CHANGE - above here

%defaults
readrawstr = 'fcbreadraw2';
timeprocstr = 'fcbtimeproc100_2';
mergeprocstr = 'fcbmergeproc2';
cellport = 3; 
beadport = 6;
pedist_thre = 9; %10; %8;  %*********CHECK**********
chljunk_coeff = .05; %.1
chljunk_power = 1; %1; %Nov 2015 .9-->1
%synSSCmin = 1;5
%pecutoff_factor = .2;
SSC2PE_cutoff = 50;  %50
cellfiletypelist = ['FCB1_XXXX_0'; 'FCB2_XXXX_0'; 'FCB1_XXXX_1'; 'FCB2_XXXX_1'; 'FCB1_XXXX_2'; 'FCB2_XXXX_2';'FCB1_XXXX_3'; 'FCB2_XXXX_3'];
cellfiletypelist = char(regexprep(cellstr(cellfiletypelist), 'XXXX', num2str(year2do)));

datapath = regexprep('\\sosiknas1\Lab_data\MVCO\FCB\MVCO_JanXXXX\data\', 'XXXX', num2str(year2do));
%datapath = '\\sosiknas1\Lab_data\MVCO\FCB\FCB_tests\docktest25Aug2016\';
%datapath = regexprep('C:\work\MVCO_janXXXX\data\', 'XXXX', num2str(year2do));

switch year2do
    case 2003
        datapath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_May2003\data\';
        cellfiletypelist = ['my1003b'; 'my1003c'; 'my1503c'; 'my1603e'; 'my1703e'; 'my1803a'; 'my1803c'; 'my2803a'; 'jn0503a'; 'jn0903a'; 'jn1603a'; 'jn1703a';...
             'jn2303a'; 'jn2503a'; 'jn2503b'; 'jn2803a'; 'jl0203a'; 'jl1003a'; 'jl2803a'; 'jl2803b'; 'jl2903a'; 'jl3103a'; 'au1403a'; 'au2003a'; 'au2303a';...
             'oc0103a'; 'oc1703a'; 'oc1903a'; 'oc2703a'; 'oc3003a'; 'no0403a'; 'no1103a'; 'no2003a'; 'no2303a'; 'no2903a'; 'de0103a'; 'de0203a'; 'de0503a';...
             'de0503b'; 'de0603a'; 'de0603c'; 'de0703a'];
        %cellfiletypelist = ['au2303a'];
        plotgroupfiletypelist = ['my'; 'jn'; 'jl'; 'au'; 'oc'; 'no'; 'de'];
        readrawstr = 'fcbreadraw1a';
        timeprocstr = 'fcbtimeproc2C_metric';
        mergeprocstr = 'fcbmergeproc1';
        cellport = 6;
        beadport = 1;
        SSC2PE_cutoff = 100;  
    case 2004
        datapath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Apr2004\data\'; 
        cellfiletypelist = ['my2004b'; 'my2504a'; 'my2704a'; 'jn0304b'; 'jn0704a'; 'jn1204b'; 'jn1404a'; 'jn2404a'; 'jn2504a'; 'jn2604b'; 'jn2804a'; 'jl0704b';...
            'jl1904a'; 'jl2304a'; 'jl2304b'; 'au0304a'; 'au1704a'; 'au1904a'; 'au2304a'; 'au2604a'; 'se0304a'; 'se0804b'; 'se1304b'; 'se1404a';'st2204b'; 'se2304a';...
            'se2704a'; 'se2904a'; 'oc0304b'; 'oc0504a'; 'oc0704a'; 'oc0804a'; 'oc1304c'; 'oc1504a'; 'oc1504b'; 'oc2004a'];
        %cellfiletypelist = ['st2204b'];
        plotgroupfiletypelist = ['my'; 'jn'; 'jl'; 'au'; 'se'; 'st'; 'oc'];
        readrawstr = 'fcbreadraw1a';
        timeprocstr = 'fcbtimeproc2B';
        mergeprocstr = 'fcbmergeproc1';
        cellport = 6;
        beadport = 1;
        SSC2PE_cutoff = 100;  
    case 2005
        datapath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Apr2005\data\';
        cellfiletypelist = ['ap1905a'; 'ap2005a'; 'ap2105a'; 'my0405a'; 'my0905c'; 'my1005a'; 'my1205a'; 'my1305a'; 'my2505a'; 'jn0805a'; 'jn1305a'; 'jn1405a';...
            'jn1605a'; 'jn2005a'; 'jn2305a'; 'jl0505a'; 'jl0505b'; 'jl1305a'; 'jl2105a'; 'jl2905a'; 'au0205a'; 'au2305a'; 'au3005a'; 'se0305a'; 'se0305b'; 'se0905a';...
            'se1305a'; 'se1905a'; 'oc2105a'; 'oc2405a'; 'oc2905a'; 'no1505a' ];
        %cellfiletypelist = ['oc2105a'; 'oc2405a'; 'oc2905a'; 'no1505a'];
        plotgroupfiletypelist = ['ap'; 'my'; 'jn'; 'jl'; 'au'; 'se'; 'oc'; 'no'];
        readrawstr = 'fcbreadraw1b';
        SSC2PE_cutoff = 800;
        SSC2PE_cutoff = 200;
    case 2006
        datapath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_May2006\data\';
        SSC2PE_cutoff = 400; 
%        cellfiletypelist = ['FCB1_2006_2']; %['FCB1_2006_1'; 'FCB1_2006_2'; 'FCB1_2006_3'];
    case 2007
        datapath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Mar2007\data\';
        SSC2PE_cutoff = 200; 
        %cellfiletypelist = ['FCB1_2007_1'];
    case 2009
        %cellfiletypelist = ['FCB2_2009_1']; 
        %SSC2PE_cutoff = 75; 
    case 2010
        %cellfiletypelist = ['FCB2_2010_2']; %['FCB1_2013_0'; 'FCB2_2013_0'];
        SSC2PE_cutoff_all = 100*ones(1,8); %for rest
        SSC2PE_cutoff_all(5:6) = [200 10]
    case 2011
        %cellfiletypelist = ['FCB2_2011_2'];
        SSC2PE_cutoff_all = 200*ones(1,8); %for rest
        SSC2PE_cutoff_all(8) = [20]
    case 2012
        SSC2PE_cutoff = 100; 
    case 2013
        %cellfiletypelist = ['FCB2_2013_1'];  
        SSC2PE_cutoff = 100;
    case 2014
    %    cellfiletypelist = ['FCB1_2014_2']; 
        SSC2PE_cutoff = 100;
    case 2015
        %cellfiletypelist = ['FCB2_2015_2'] %; 'FCB2_2015_1'];
        SSC2PE_cutoff = 200;
    case 2016
        SSC2PE_cutoff = 200;
    case 2017
        SSC2PE_cutoff = 200;
    case 2018
        SSC2PE_cutoff = 200;
    end;

    %datapath = regexprep(datapath, '\\\\queenrose\\mvco', 'c:\\work'); %temp for olive
beadfiletypelist = cellfiletypelist;
if ~exist('plotgroupfiletypelist', 'var')
    plotgroupfiletypelist = ['FCB*_' num2str(year2do)];   %will be used with *
end;

temp = [datapath 'processed\']; if ~exist(temp, 'dir'), mkdir(temp), end;
temp = [datapath 'processed\time\']; if ~exist(temp, 'dir'), mkdir(temp), end;
temp = [datapath 'processed\beads\']; if ~exist(temp, 'dir'), mkdir(temp), end;
temp = [datapath 'processed\grouped\']; if ~exist(temp, 'dir'), mkdir(temp), end;
temp = [datapath 'processed\grouped\merged\']; if ~exist(temp, 'dir'), mkdir(temp), end;    
baseprocpath = [datapath 'processed\'];
setsize = 200;  %how many data files to process in one set

if dotime,
    procpath = [baseprocpath 'time\'];
    plotflag = timeplotflag; 
    disp('PROCESSING TIME FILES')
    if docells,
        filetypelist = cellfiletypelist;
        timebatch2
    end;
 end;

if domerge,
    procpath = baseprocpath;
    plotflag = mergeplotflag;
    disp('MERGING')
    if docells,
        filetype = 'cell';  %ALWAYS MERGE CELLS BEFORE CORRESPONDING BEADS
        filetypelist = cellfiletypelist;
        mergebatch2
    end;
end;

if doclassify,
    procpath = baseprocpath;
    timepath = [baseprocpath 'time\'];
    beadpath = [baseprocpath 'beads\']; 
    groupedpath = [baseprocpath 'grouped\']; 
    mergedpath = [baseprocpath 'grouped\merged\'];
    plotflag = classplotflag;
    disp('CLASSIFYING')
    if dobeads,
        datapath = baseprocpath; 
        savepath = [baseprocpath 'beads\']; 
        filetypelist = beadfiletypelist;
        beadbatch6M %temp test for bead movies July 11, 2016
        %beadbatch6_field 
    end;
    if docells,
        filetypelist = cellfiletypelist;
        cellbatch8d_field
        %cellbatch7_field
    end;
end;

if doplotgroup,
    disp('PLOTTING FINAL RESULTS')
    procpath = [baseprocpath 'grouped\'];
    filetypelist = plotgroupfiletypelist;
    plotgroup_field
end;

if movieflag
    addpath \\sosiknas1\Lab_data\MVCO\FCB\MVCO_movies\code\
    groupedpath = [baseprocpath 'grouped\']; 
    mergedpath = [baseprocpath 'grouped\merged\'];
    savepath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_movies\together_movies\';
    mvco_movies_avi_format4(year2do,groupedpath,mergedpath,savepath)
end
end
   