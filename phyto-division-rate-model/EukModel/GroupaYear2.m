%% Summarize Results for a Year


%outputpath = directory where the model outputs are 
% e.g. outputpath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2020\syn_model\dawnstart_outputs'

%obs_conc = Observed Concentration over day
% pred_conc = predicted / simulated concnetration for each hour of simulation
%obs_hr_mu = Observed net grwoth rate for each hour to hour
%pred_hr_mu = Predicted hourly growth rate for each hour to hour in perday 
%hr_loss = difference between net growth and predicted growth 

%modelresults are the same as usual 
% vector = [day xmin fmin mu mu1 mu2 exitflag length(modelfits)]
%where xmin is optimal theta value. Therefore Theta = modelresults(2:15)
%fmin is negative log likelihood given these optimal parameters
% mu mu1 mu2 are results of growth_rate function for these paramters
% exitflag is the exitflag for that optimal run, and length(modelfits)
% gives us a sense of how long it took the function to find the answer

function GroupaYear2(outputpath)

inputpath = regexprep(outputpath, 'outputs', 'inputs'); 

%%make a savename and Variable for grouped results
savename = [outputpath filesep 'AllModelOutputs.mat']; 
AllResults = []; 

filelist = dir([outputpath filesep '*output.mat']);

%set up for making video
Writerobj1 = VideoWriter([outputpath filesep 'ModelOutputs.avi']);
open(Writerobj1); 

obs_conc = []; 
obs_hr_mu = [];
pred_hr_mu = obs_hr_mu; 
hr_loss = obs_hr_mu;
pred_conc = obs_conc; 

for i = 1:length(filelist);
    eval(['load ' outputpath filesep filelist(i).name])
    
    %for 2006 only
    %allstarts = allstarts{:}' ; 
    %modelfits = modelfits{:}'; 
    
    
    AllResults = [modelresults; AllResults]; %Save model results to grouped results
    

    %load the input variables too, they are usefull 
    day = modelresults(1); 
    eval(['load ' inputpath filesep 'day' num2str(day) 'data.mat'])

    %Interpolate light data
        time=0:(1/6):25;
        nnind = find(~isnan(Edata(:,2)));
        Edata=Edata(nnind,:);
        [unqE, eind]=unique(Edata(:,1));
        Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
        Einterp(isnan(Einterp)) = 0;
        
    %make the figure
    Modelfit_Frame
    F1=getframe(gcf);%and add this frame to video
    
    writeVideo(Writerobj1, F1);
    clf
    
    if ~exist('CONC')
        CONC = Vhists.*cellsperml; 
    end

    sumconc = sum(CONC(:, 1:24)); 
    sumsim = sum(simCONC); 
    obs_conc(i,:) = sumconc; 
    obs_hr_mu(i,:) = 24*log(sumconc(2:end)./sumconc(1:(end-1))); 
    pred_hr_mu(i,:) = 24*log(sumsim(2:end)./sumsim(1:(end-1))); 
    pred_conc(i,:) = sumsim; 
    hr_loss(i,:) = pred_hr_mu(i,:) - obs_hr_mu(i,:); 
        

end

close(Writerobj1) 

%Flip it around cuz its upside down 
AllResults = AllResults(end:-1:1,:);

save(savename, 'AllResults', 'obs_conc', 'obs_hr_mu', 'pred_hr_mu', 'pred_conc', 'hr_loss'); 


