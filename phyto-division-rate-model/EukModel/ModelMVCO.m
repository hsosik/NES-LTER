

%modified to make it a function rather than a script. 

%input looks something like: modelpath =
%'\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2020\euk_model\dawnstart_inputs'; 

%will look for eukvolbins or synvolbins within input file to decide if euk
%or syn model should be fit.


function ModelMVCO(modelpath) 


filelist = dir([modelpath '*data.mat']);
savepath = regexprep(modelpath, 'inputs', 'outputs'); 

if ~exist(savepath)
    mkdir(savepath); 
end

%start while loop 
i = 1; 
while i <= length(filelist)
     
    %load days data 
    eval(['load ' modelpath filesep filelist(i).name])
    savename = [filelist(i).name(1:9) 'output.mat']; 
    
    if exist('eukvolbins')
        volbins = eukvolbins; 
        ts = 0; 
    else
        volbins = synvolbins; 
        ts = 6; 
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
       
        eval('dayval = day'); %stupid variable name causing trouble again

        %apply model
        if i > 1 && exist([filelist(i).folder filesep filelist(i-1).name(1:end-9) 'output.mat']) %get previous result if available
          load([filelist(i).folder filesep filelist(i-1).name(1:end-9) 'output.mat'], 'modelresults')
         informedguess = modelresults(2:15); 
         clear modelresults
            [modelresults, modelfits, allstarts, simPROPS, simCONC] =  OneDayFastStart2(dayval, Einterp, volbins, N_dist, 24, ts, informedguess);     
        else
            [modelresults, modelfits, allstarts, simPROPS, simCONC] =  OneDayFastStart2(dayval, Einterp, volbins, N_dist, 24, ts);     
        end
        

        %save results
        save([savepath filesep savename], 'N_dist', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp', 'ts')
               
    else
        disp([filelist(i).name(1:9) 'not enough hours'])
    end %if we have 25 hours
    
   
   clearvars('-except', 'Writerobj1', 'i', 'filelist', 'modelpath', 'savepath'); 
   
i = i + 1; 

end

end
