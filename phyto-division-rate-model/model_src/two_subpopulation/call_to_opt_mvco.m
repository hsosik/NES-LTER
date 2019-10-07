%main script to process field data with 14-parameter model assuming a dirichlet multinomial distribution
%for likelihood calculation:

%script runs batches of random startpoints with multistart as
%this seems to be faster than running fmincon one solver run at time to

% clear all
% close all

%some front matter:
hr1=7; hr2=25; %time window for the model

restitles={'day';'gmax1';'b1';'E*1';'dmax1';'gmax2';'b2';'E*2';'dmax2';'proportion';'m1';'m2';'sigma1';'sigma2';'s';'-logL';'mu';'mu1';'mu2';'ending proportion 1';'ending proportion 2'; 'exitflag';'number solver runs'};
notes='E* bounds are from 0 to max(Einterp), run with scaled s (by a factor of 100)';

ms=MultiStart('Display','off','TolX',1e-5,'UseParallel','always','StartPointsToRun','bounds');
opts=optimset('Display','off','TolX',1e-8,'Algorithm','interior-point','UseParallel','always','MaxIter', 3000,'MaxFunEvals',10000);
icsTol=0.2;
tolvec=[0.01 0.01 100 0.005 0.01 0.01 100 0.005 0.01 0.5 0.5 0.5 0.5 10];

%both savepath and pathname should be specified before running this code!
%This is done in divmodle_setup!

%Which ddata's are we running?

switch datatype
    
    case 'field'
        %Check first to see if data already exists:
        if exist(fullfile(savepath, ['mvco_14par_dmn_' num2str(year2do) '.mat']),'file') == 2 %meaning already a file, probably was interrupted...
            
            disp('found a file with this name; load and continue?')
            keyboard
            
            load(fullfile(savepath, ['mvco_14par_dmn_' num2str(year2do) '.mat']))
            jj=find(modelresults(:,1)==0);
            
            if unique(diff(jj))==1 %meaning that there are zeros in date column, but they are sequential
                start_file_num=jj(1); %begin from
            else
                disp('Check to make sure models runs will be done for correct days in this file')
                keyboard
            end
            
        else %fresh start!
            
            disp(num2str(year2do))
            filelist = dir([pathname 'day*data.mat']); %find the input files
            modelresults=zeros(length(filelist),23);
            allmodelruns=cell(length(filelist),2);
            start_file_num=1;
            
        end
        
    case 'redo'
        
        disp(num2str(year2do))
        filelist = dir([pathname 'day*data.mat']); %find the input files
        temp=char({filelist(:).name}'); temp=temp(:,4:9); temp=str2num(temp);
        qq=find(ismember(temp,days2redo));
        filelist=filelist(qq);
        
        modelresults=zeros(length(filelist),23);
        allmodelruns=cell(length(filelist),2);
        start_file_num=1;
        
    case 'benchtop'
        %indexes of 'good days' to run:
        load(fullfile(codepath,'validation/FCB3_benchtop/model_setup_and_postprocessing/indexes_new.mat'))
        
        %find both port 2 and port 3 files, and check that correct files are found:
        filelist2 = dir([pathname 'day*_2_*.mat']); %find the input files
        temp_daylist2=str2num(char(cellfun(@(x) x(4:9), {filelist2(indtemp2).name}','UniformOutput',false)));
        if ~isequal(culture_date(c2m2(indtemp2)),temp_daylist2), disp('Uh oh - not finding correct benchtop files!'); keyboard, end
        
        filelist3 = dir([pathname 'day*_3_*.mat']); %find the input files
        temp_daylist3=str2num(char(cellfun(@(x) x(4:9), {filelist3(indtemp3).name}','UniformOutput',false)));
        if ~isequal(culture_date(c2m3(indtemp3)),temp_daylist3), disp('Uh oh - not finding correct benchtop files!'); keyboard, end
        
        filelist=[filelist2(indtemp2); filelist3(indtemp3)];
        
        modelresults=zeros(length(filelist),23);
        allmodelruns=cell(length(filelist),2);
        start_file_num=1;
end

%%
for filenum=start_file_num:length(filelist)
    
    filename=filelist(filenum).name;
    day=str2num(filename(4:9));
    
    disp(['optimizing day: ' num2str(day) ' file#: ' num2str(filenum) ' out of ' num2str(length(filelist))])
    eval(['load ' pathname filename])
    
    %special fix for shorter days, but still have enough hours to fit the model
    if size(N_dist,2) < 25
        m=size(N_dist,2);
        N_dist=[nan(57,25-m) N_dist];
        Vhists=[nan(57,25-m) Vhists];
    end
    
    %temporary - should fix in getSolar for benchtop data:
    if exist('valtype','var')
        ind= Edata(:,2) < 6; %deal with extra noise (dark current values)
        noise=mean(Edata(ind,2));
        Edata(:,2)=Edata(:,2)-noise;
        ind = Edata(:,2) < 0;
        Edata(ind,2) = 0;
        ind=find(Edata(:,2) < 4);
        threshold=mean(Edata(ind,2)) + std(Edata(ind,2));
        ii=find(Edata(ind,2) > threshold); %find values that are still too noisy
        Edata(ind(ii),2)=0;
        indh=find(Edata(:,1) < 1 | Edata(:,1) > 16.5); %for spurious light values at end of day...
        ind2=find(Edata(indh,2) > threshold);
        Edata(indh(ind2),2)=0;
    end
    
    %Fix and Interpolate Light Data:
    time=0:(1/6):25;
    nnind = find(~isnan(Edata(:,2)));
    Edata=Edata(nnind,:);
    if Edata(1,1) >= 1, Edata=[0 0; Edata]; end %rare case for gaps around dawn that are not caught by getSolar...i.e.Jun-10-2003
    [unqE eind]=unique(Edata(:,1));
    Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
    Einterp(find(isnan(Einterp))) = 0;
    
    lb=-[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 5 5 1 1 1e-4]; %parameter bounds
    ub=[1 15 max(Einterp) 1 1 15 max(Einterp) 1 0.5 50 50 15 15 1e4];
    
    a1=-1*eye(14); %set parameter bounds to be interpretted by fmincon
    a2=eye(14);
    A=zeros(28,14);
    A(1:2:27,:)=a1;
    A(2:2:28,:)=a2;
    
    B=zeros(27,1);
    B(1:2:27)=lb;
    B(2:2:28)=ub;
    
    %starting conditions:
    x0=[0.2*rand 6*rand max(Einterp)*rand 0.1*rand 0.2*rand 6*rand max(Einterp)*rand 0.1*rand 0.5*rand 30*rand+20 30*rand+20 10*rand+2 10*rand+2 1e4*rand];
    %random start points:
    tpoints = CustomStartPointSet([0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 0.5*rand(40,1) 30*rand(40,1)+20 30*rand(40,1)+20 10*rand(40,1)+2 10*rand(40,1)+2 1e4*rand(40,1)]);
    
    problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) negloglike_calc(Einterp,N_dist,theta,volbins,hr1,hr2),'Aineq',A,'bineq',B,'options',opts);
    [xmin,fmin,exitflag,~,soln] = run(ms,problem,tpoints);
    
    %open up the soln structure:
    temp=zeros(40,17);
    start_points=zeros(40,14);
    temp=zeros(40,17);
    c=1;
    
    for j=1:length(soln)
        %check to see if all start points led to an individual solution or
        %not (MultiStart will only return unique solutions)
        g=cell2mat(soln(j).X0);
        if length(g)==14 %only one start_point led to that solution
            start_points(c,:)=g;
            temp(c,1:14)=soln(j).X;
            temp(c,15)=soln(j).Fval;
            temp(c,16)=growth_rate(Einterp,volbins,N_dist,temp(c,1:13),hr1,hr2);
            temp(c,17)=soln(j).Exitflag;
            c=c+1;
        else
            num=length(g)/14;
            start_points(c:c+num-1,:)=squeeze(reshape(g',1,14,num))';
            temp(c:c+num-1,1:14)=repmat(soln(j).X,num,1);
            temp(c:c+num-1,15)=repmat(soln(j).Fval,num,1);
            temp(c:c+num-1,16)=repmat(growth_rate(Einterp,volbins,N_dist,temp(c,1:13),hr1,hr2),num,1);
            temp(c:c+num-1,17)=repmat(soln(j).Exitflag,num,1);
            c=c+num;
        end
    end
    %just in case have rows left as zeros
    qq=find(temp(:,1)~=0);
    temp=temp(qq,:);
    
    largepopn=zeros(size(temp,1),7);
    smallpopn=zeros(size(temp,1),7);
    for h=1:size(temp,1)
        if temp(h,10) > temp(h,11)
            largepopn(h,:) = temp(h,[1:4 9 10 12]);
            smallpopn(h,:) = [temp(h,5:8) 1-temp(h,9) temp(h,[11 13])];
        else %11 > 10
            largepopn(h,:) = [temp(h,5:8) 1-temp(h,9) temp(h,[11 13])];
            smallpopn(h,:) = temp(h,[1:4 9 10 12]);
        end
    end
    
    modelfits=[smallpopn(:,1:4) largepopn(:,1:4) smallpopn(:,5) smallpopn(:,6) largepopn(:,6) smallpopn(:,7) largepopn(:,7) temp(:,14:end)];
    start_points=start_points(qq,:);
    allstarts=start_points;
    
    %let's now ask, in the first batch run, did the solver "converge"?
    [sortlogL ii]=sort(modelfits(:,15));
    if abs(sortlogL(min(5,size(sortlogL,1)))-sortlogL(1)) < icsTol
        flag1 = 0;
    else
        disp(num2str(sortlogL(1:min(5,size(sortlogL,1)))))
        flag1 = 1;
    end;
    
    partol=max(modelfits(ii(1:min(5,size(sortlogL,1))),1:14))-min(modelfits(ii(1:min(5,size(sortlogL,1))),1:14));
    if sum(abs(partol) < tolvec)==14 || sum((abs(partol./modelfits(ii(1),1:14)) < 0.05))==14 %either the modelfits are within an absolute tolerance or within a relative tolerance
        flag2 = 0;
    else
        flag2 = 1;
    end
    
    disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2) ' for day ' datestr(day) ': ' num2str(filenum) ' out of ' num2str(length(filelist))])
    
    %Do more model runs until both stopping criteria are met:
    %In case of model redo's, do at least 200 model runs:  
    
    k=1; %batch number
    if strcmp(datatype,'redo'), w=5; else w=0; end %separate counter to ensure 200 model runs are done
    
    while (((flag1 || flag2) || size(modelfits,1) <= 20) && k <= 5) || w>0
        
        disp(['k: ' num2str(k)])
        disp(['w: ' num2str(w)])
        k=k+1;
        w=w-1;
        
        x0=[0.2*rand 6*rand max(Einterp)*rand 0.1*rand 0.2*rand 6*rand max(Einterp)*rand 0.1*rand 0.5*rand 30*rand+20 30*rand+20 10*rand+2 10*rand+2 1e4*rand]; %random start point
        
        tpoints = CustomStartPointSet([0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 0.5*rand(40,1) 30*rand(40,1)+20 30*rand(40,1)+20 10*rand(40,1)+2 10*rand(40,1)+2 1e4*rand(40,1)]);
        
        problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) negloglike_calc(Einterp,N_dist,theta,volbins,hr1,hr2),'Aineq',A,'bineq',B,'options',opts);
        [xmin,fmin,exitflag,~,soln] = run(ms,problem,tpoints);
        
        %open up the soln sturcutre:
        temp=zeros(40,17);
        start_points=zeros(40,14);
        temp=zeros(40,17);
        c=1;
        
        for j=1:length(soln)
            %check to see if all start points led to an individual solution or
            %not (MultiSTart will only return unique solutions)
            g=cell2mat(soln(j).X0);
            if length(g)==14 %only one start_point led to that solution
                start_points(c,:)=g;
                temp(c,1:14)=soln(j).X;
                temp(c,15)=soln(j).Fval;
                temp(c,16)=growth_rate(Einterp,volbins,N_dist,temp(c,1:13),hr1,hr2);
                temp(c,17)=soln(j).Exitflag;
                c=c+1;
            else
                num=length(g)/14;
                start_points(c:c+num-1,:)=squeeze(reshape(g',1,14,num))';
                temp(c:c+num-1,1:14)=repmat(soln(j).X,num,1);
                temp(c:c+num-1,15)=repmat(soln(j).Fval,num,1);
                temp(c:c+num-1,16)=repmat(growth_rate(Einterp,volbins,N_dist,temp(c,1:13),hr1,hr2),num,1);
                temp(c:c+num-1,17)=repmat(soln(j).Exitflag,num,1);
                c=c+num;
            end
        end
        %just in case have rows left as zeros
        qq=find(temp(:,1)~=0);
        temp=temp(qq,:);
        start_points=start_points(qq,:);
        
        largepopn=zeros(size(temp,1),7); %large population has the larger mean starting volume bin
        smallpopn=zeros(size(temp,1),7);
        for h=1:size(temp,1)
            if temp(h,10) > temp(h,11)
                largepopn(h,:) = temp(h,[1:4 9 10 12]);
                smallpopn(h,:) = [temp(h,5:8) 1-temp(h,9) temp(h,[11 13])];
            else %11 > 10
                largepopn(h,:) = [temp(h,5:8) 1-temp(h,9) temp(h,[11 13])];
                smallpopn(h,:) = temp(h,[1:4 9 10 12]);
            end
        end
        
        modelfits=[modelfits; smallpopn(:,1:4) largepopn(:,1:4) smallpopn(:,5) smallpopn(:,6) largepopn(:,6) smallpopn(:,7) largepopn(:,7) temp(:,14:end)];
        allstarts=[allstarts; start_points];
        
        %okay, now see after this batch run, did the solver "converge"?
        [sortlogL ii]=sort(modelfits(:,15));
        
        if abs(sortlogL(5)-sortlogL(1)) < icsTol
            flag1 = 0;
        else
            disp(num2str(sortlogL(1:5))) %should be 5, but occassionally get less than 5 solver runs returned...
            flag1 = 1;
        end;
        
        partol=max(modelfits(ii(1:5),1:14))-min(modelfits(ii(1:5),1:14));
        if sum(abs(partol) < tolvec)==14 || sum((abs(partol./modelfits(ii(1),1:14)) < 0.05))==14 %either the modelfits are within an absolute tolerance or within a relative tolerance
            flag2 = 0;
        else
            flag2 = 1;
        end
        disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2) ' for day ' datestr(day) ': ' num2str(filenum) ' out of ' num2str(length(filelist))])
        
    end  %while loop
    
    [s jj]=sort(modelfits(:,15));
    xmin=modelfits(jj(1),1:14);
    fmin=modelfits(jj(1),15);
    exitflag=modelfits(jj(1),17);
    
    [mu mu1 mu2 p1 p2]=growth_rate(Einterp,volbins,N_dist,xmin(1:13),hr1,hr2);
    
    modelresults(filenum,:)=[day xmin fmin mu mu1 mu2 p1 p2 exitflag length(modelfits)];
    allmodelruns{filenum,1}=modelfits;
    allmodelruns{filenum,2}=allstarts;
    
    switch datatype
        case 'field'
            eval(['save ' savepath 'mvco_14par_dmn_' num2str(year2do) ' modelresults allmodelruns'])
        case 'redo'
            eval(['save ' savepath 'mvco_14par_dmn_' num2str(year2do) '_redos modelresults allmodelruns'])
        case 'benchtop'
            eval(['save ' savepath 'benchtop_14par_dmn_modelfits modelresults allmodelruns'])
    end
    
end

