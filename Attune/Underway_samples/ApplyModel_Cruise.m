
%cd('\\sosiknas1\Lab_data\MVCO\FCB\pico_euk_model\EukWork\Model_1pop\')
% Model now has capacity to fit two subpopulations or just fit 1 population.
%for now, let's fit 2 populations as before. 

%inputspath is link to input.mat files for a given cruise
%for example:
%inputspath = '\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\SlidingWindow\HRS2303'

function ApplyModel_Cruise(inputspath)

cd('\\sosiknas1\Lab_data\MVCO\FCB\pico_euk_model\EukWork\Model_1pop\') %make sure this line connects to wherever the current model scripts are

filelist = dir([inputspath filesep '*input.mat']);


%start while loop 
i = length(filelist); 
while i >= 1%length(filelist)
     savename = [filelist(i).name(1:end-9) 'output.mat']; 

     if exist([filelist(i).folder filesep savename])
         i = i-1; %switched it to go backwards, also changed i-1 to i+1 in informedguess linse 41 and 42 
     else
    %load days data 
    eval(['load ' filelist(i).folder filesep filelist(i).name])

    if exist('eukvolbins')
        volbins = eukvolbins; 
    else
        volbins = synvolbins; 
    end
        volbins(volbins == 0) = []; 

    if size(N_dist, 2) == 25
        
        %Interpolate light data
        time=0:(1/6):25;
        nnind = find(~isnan(Edata(:,2)));
        Edata=Edata(nnind,:);
        [unqE, eind]=unique(Edata(:,1));
        Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
        Einterp(isnan(Einterp)) = 0;
       
        day = floor(datenum(daystarttime)); 
        
        if i < length(filelist) && exist([filelist(i).folder filesep filelist(i+1).name(1:end-9) 'output.mat']) %get previous result if available
          load([filelist(i).folder filesep filelist(i+1).name(1:end-9) 'output.mat'], 'modelresults')
         informedguess = modelresults(2:15); 
         clear modelresults
            [modelresults, modelfits, allstarts, simPROPS, simCONC] =  OneDayFastStart2(day, Einterp, volbins, N_dist, 24, ts, informedguess);     
        else
            [modelresults, modelfits, allstarts, simPROPS, simCONC] =  OneDayFastStart2(day, Einterp, volbins, N_dist, 24, ts);     
        end
                  
        %save results
        save([filelist(i).folder filesep savename], 'N_dist', 'Vhists', 'cellsperml', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp', 'ts')
               
    else
        disp([filelist(i).name(1:9) 'not enough hours'])
    end %if we have 25 hours
    
   
   clearvars('-except', 'Writerobj1', 'i', 'filelist', 'filepath', 'savepath', 'ts'); 
   
i = i - 1; 
     end

end

cd('\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\Scripts')

end