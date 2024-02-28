function [] = test_class_assignments(p)
%test class assignments for a cruise WITHOUT overwriting class data 

% basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20230807_EN706'; 
% assign_class_function = 'assign_class_EN706'; 
% 
% filetype2exclude = {'fcb_bead'; 'FCB_bead'; 'bead'; '(lab test)'; 'test';'Dockwater'; 'Daily'; 'Rinses'; 'discrete'; "Filter config"; "Grazer"; "SFD_AR29_Grazer"; "Cultures"; "08Aug2023"}; %needed for Step2
% OD2setting = 'GL1'; %where was the OD2 filter on this cruise? 'SSC', 'GL1', or 'None' 
% 
% framemaker = 'make_movieframe_density';
% stepsize = 25; %controls resolution of files viewed
% moviechannels = 'late'; 
file2start = '';
while ~isinteger(file2start)
    file2start = input(1,'What file number to view first? (enter for default = 1)');
    if isempty(file2start)
        file2start = 1;
    end
end
stepsize = '';
while ~isinteger(file2start)
    stepsize = input(1,'What step size for file viewing? (enter for default = 25)');
    if isempty(stepsize)
        stepsize = 25;
    end
end

%assign useful paths
fpath = [p.basepath filesep 'FCS' filesep];
outpath = [p.basepath filesep 'bead_calibrated' filesep];
beadfigpath = [outpath filesep 'bead_plots'];
classpath = [outpath 'class' filesep];
load([outpath filesep 'FCSfileinfo.mat'])

if strcmp(p.OD2setting, 'GL1')
    ssc_name = 'GL1';  
else
    ssc_name = 'SSC'; 
end

 %remove filetype2exclude from FCSfileinfo and generate filelist
        %this part has been updated  for Table structure FCS filinfo 
    for iii = 1:length(p.filetype2exclude)    
         t = contains(FCSfileinfo.fcslist, p.filetype2exclude{iii});
         if ~isempty(t)
              FCSfileinfo(t, :) = []; 
         end
    end
    clear t iii  
   
    filelist = FCSfileinfo.fcslist;
     

    %make sure its in chronological order
    matdate = FCSfileinfo.matdate_start;
    QC_flags = FCSfileinfo.QC_flag; 
    [matdate,sort_ind] = sort(matdate);
    filelist = filelist(sort_ind);
    QC_flags = QC_flags(sort_ind); 


    %QC_flags = FCSfileinfo.QC_flowrates(sort_ind,1)<1.1; 

    length(filelist)
    %now go through files of interest, assign classes, and save results
    for count = file2start:stepsize:length(filelist)
        pause(.02)
         if ~rem(count,10)
            disp([num2str(count) ' of ' num2str(length(filelist))])
         end
        filename = [fpath filelist{count}];
        disp(filename)
        [fcsdat,fcshdr] = fca_readfcs(filename);
        [~,fname] = fileparts(filename);
        [class, bounds] = eval([assign_class_function '( fcsdat, fcshdr, 0, fname, FCSfileinfo.QC_flag(count), FCSfileinfo.matdate_start(count) );']); 
        clear fname
   
        
            % assign scattering channels
            ssc_ch = strmatch([ssc_name '-A'], {fcshdr.par.name});
            ssch = strmatch([ssc_name '-H'], {fcshdr.par.name});
            
            % correct for negative SSC-A values
            cf = fitlm(fcsdat(fcsdat(:,ssch)<1000,ssch), fcsdat(fcsdat(:,ssch)<1000, ssc_ch), 'Intercept', false);
            cf = cf.Coefficients.Estimate;
            fcsdat(fcsdat(:,ssc_ch)<0, ssc_ch) = cf*fcsdat(fcsdat(:,ssc_ch)<0, ssch);
            
            % call the function to make the plots and getframe
            eval(['Frame = ', p.framemaker, '(fcsdat, fcshdr, class, moviechannels, QC_flags(count), bounds)']);
     

            set(gcf, 'position', [2.1483e+03 -187 905.3333 590.6667])
    end
end