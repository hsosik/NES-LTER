%Goal: Find best fit of model to 1 day of data. 

%Inputs: 
%day = matlab date number for the day of interest 
%Solar is a two collumn light input, with time and irradiance, in LOCAL time.
% must contain but not necessarily be limited to the day of interst. 
%volbins = designation of size classes of cells in population
%dayCONC = matrix with the number of cells in each size class (rows) as
%specificed by volbins, for each hour (collumns) of day of interest

%Outputs: 
%modelresults is 1 x 21 vector = [day xmin fmin mu mu1 mu2 exitflag length(modelfits)]
%where xmin is optimal theta value. Therefore Theta = modelresults(2:15)
%fmin is negative log likelihood given these optimal parameters
% mu mu1 mu2 are results of growth_rate function for these paramters
% exitflag is the exitflag for that optimal run, and length(modelfits)
% gives us a sense of how long it took the function to find the answer

%modelfits stores [ theta negloglike mu exitflag ] for each result of fmincon
%allstarts stores initial parameter value for every run
%daySimPROPS is matrix with the expected proportion of cells in each size class
%(rows) as specified by volbins, for each hour (collumns) according to the
%simulation with best fit parameters 

%Einterp is the interpolated light data for the particular day of interest 

%to allow for day starting at different times of day, and still limiting
%division for syn after dawn. ts can be a 1x2 value with min and max
% so division is limited if hr >= ts(1) & <= ts(2) 
%should be hour numbers relative to "start of day" whether or not that's dawn. 


% we're going to stop if division rate converges to within 3 decimal places
% .001 
%changed ictol.9 



function [modelresults, modelfits, allstarts, simPROPS, simCONC, Einterp] = OneDayFastStart2(day, Einterp, volbins, COUNTS, hr2, ts, informedguess)

    
hr1 = 1; %choose start hour for model fitting
 
%% Fit Model to Data

ms=MultiStart('Display','off','TolX',1e-5,'UseParallel', false ,'StartPointsToRun','bounds-ineqs');
opts=optimset('Display','off','TolX',1e-8,'Algorithm','interior-point','UseParallel','always','MaxIter', 3000,'MaxFunEvals',10000);
icsTol=0.9;
tolvec=[0.01 0.01 100 0.005 0.01 0.01 100 0.005 0.01 0.5 0.5 0.5 0.5 10];

%Define parameter bounds 
lb=-[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 5 5 .5 .5 1e-4]; %negative lower bounds
ub=[1 15 1000 1 1 15 1000 1 0.5 50 50 15 15 1e4]; %upper bounds 

    a1=-1*eye(14); %set parameter bounds to be interpretted by fmincon
    a2=eye(14);
    A=zeros(28,14);
    A(1:2:27,:)=a1;
    A(2:2:28,:)=a2;

    B=zeros(27,1);
    B(1:2:27)=lb;
    B(2:2:28)=ub;
    
%Create random initial parameter values within bounds 
x0 = zeros(1,14);
num_starts = 40; %choose the number of starting values we want, number of trials of fmincon
x2 = zeros(num_starts, 14); %matrix to create tpoints
for i = 1:14
    x0(i) = -lb(i) + (ub(i) + lb(i))*rand(1);  %remember lb has negative lower bounds. so this is lower bound + rand * difference 
    x2(:,i) = -lb(i) + (ub(i) + lb(i))*rand(num_starts,1);
end
clear i
x2(:,9) = 1/(2*num_starts)*(1:num_starts); %disperse Theta(9) evenly between 0 and .5, rather than randomly 

if exist('informedguess')
   informedguess(9) = min(informedguess(9), 1-informedguess(9)); 
    x2(1,:) = informedguess; 
end
tpoints = CustomStartPointSet(x2); 


problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) negloglike_calc(Einterp,COUNTS,theta,volbins,hr1,hr2, ts),'Aineq',A,'bineq',B,'options',opts);
[~,~,~,~,soln] = run(ms,problem,tpoints);
    


    %open up the soln structure:
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
            temp(c,16)=growth_rate(Einterp,volbins,COUNTS,temp(c,1:13),hr1,hr2, ts);
            temp(c,17)=soln(j).Exitflag;
            c=c+1;
        else
            num=length(g)/14;
            start_points(c:c+num-1,:)=squeeze(reshape(g',1,14,num))';
            temp(c:c+num-1,1:14)=repmat(soln(j).X,num,1);
            temp(c:c+num-1,15)=repmat(soln(j).Fval,num,1);
            temp(c:c+num-1,16)=repmat(growth_rate(Einterp,volbins,COUNTS,temp(c,1:13),hr1,hr2, ts),num,1);
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
    [sortlogL, ii]=sort(modelfits(:,15));
    if abs(sortlogL(min(5,size(sortlogL,1)))-sortlogL(1)) < icsTol  %Did likelihood converge within tolerance
        flag1 = 0; %Yes it did converge
    else
        disp(num2str(sortlogL(1:min(5,size(sortlogL,1)))))
        flag1 = 1; %no it didn't 
    end

       %did parameters converge within a tolerance 
    partol=max(modelfits(ii(1:min(5,size(sortlogL,1))),1:14))-min(modelfits(ii(1:min(5,size(sortlogL,1))),1:14));
    if sum(abs(partol) < tolvec)==14 || sum((abs(partol./modelfits(ii(1),1:14)) < 0.05))==14 %either the modelfits are within an absolute tolerance or within a relative tolerance
        flag2 = 0; %yes they did converge 
    else
        flag2 = 1; %no they didn't converge 
    end

    disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])

    
    %check if division rates converge
        [mui] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(1),1:14),hr1,hr2, ts);
        [mui2] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(2),1:14),hr1,hr2, ts);
        [mui3] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(3),1:14),hr1,hr2, ts);
        [mui4] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(4),1:14),hr1,hr2, ts);
        [mui5] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(5),1:14),hr1,hr2, ts);
        
        
        flag3 = (range([mui mui2 mui3 mui4 mui5]) < .001 ) ;
        
        
    k=1; %batch number
    % need negloglike within .9 and divrate within .001 and at least 20 runs or
    % max out k. 
    while ((flag1 && ~flag3) || size(modelfits,1) <= 20) && k <= 5 %run MultiStart 4 more times D: 

        disp(['k: ' num2str(k)])
        k=k+1;

        %Create random initial parameter values within bounds 
        x0 = zeros(1,14);
        num_starts = 40; %choose the number of starting values we want, number of trials of fmincon
        x2 = zeros(num_starts, 14); %matrix to create tpoints
        for i = 1:14
             x0(i) = -lb(i) + (ub(i) + lb(i))*rand(1);  %remember lb has negative lower bounds. so this is lower bound + rand * difference 
             x2(:,i) = -lb(i) + (ub(i) + lb(i))*rand(num_starts,1);
        end
        clear i
        x2(:,9) = 1/(2*num_starts)*(1:num_starts); %disperse Theta(9) evenly between 0 and .5, rather than randomly 
        tpoints = CustomStartPointSet(x2); 

        problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) negloglike_calc(Einterp,COUNTS,theta,volbins,hr1,hr2, ts),'Aineq',A,'bineq',B,'options',opts);
        [~,~,~,~,soln] = run(ms,problem,tpoints);

        %open up the soln structure:
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
                temp(c,16)=growth_rate(Einterp,volbins,COUNTS,temp(c,1:13),hr1,hr2, ts);
                temp(c,17)=soln(j).Exitflag;
                c=c+1;
            else
                num=length(g)/14;
                start_points(c:c+num-1,:)=squeeze(reshape(g',1,14,num))';
                temp(c:c+num-1,1:14)=repmat(soln(j).X,num,1);
                temp(c:c+num-1,15)=repmat(soln(j).Fval,num,1);
                temp(c:c+num-1,16)=repmat(growth_rate(Einterp,volbins,COUNTS,temp(c,1:13),hr1,hr2, ts),num,1);
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

        %modelfits stores the proportion of population with smaller mean for theta(9), 
        modelfits=[modelfits; smallpopn(:,1:4) largepopn(:,1:4) smallpopn(:,5) smallpopn(:,6) largepopn(:,6) smallpopn(:,7) largepopn(:,7) temp(:,14:end)];
        allstarts=[allstarts; start_points];

        %okay, now see after this batch run, did the solver "converge"?
        [sortlogL, ii]=sort(modelfits(:,15));

        if abs(sortlogL(5)-sortlogL(1)) < icsTol
            flag1 = 0;
        else
            disp(num2str(sortlogL(1:5))) %should be 5, but occassionally get less than 5 solver runs returned...
            flag1 = 1;
        end

        partol=max(modelfits(ii(1:5),1:14))-min(modelfits(ii(1:5),1:14));
        if sum(abs(partol) < tolvec)==14 || sum((abs(partol./modelfits(ii(1),1:14)) < 0.05))==14 %either the modelfits are within an absolute tolerance or within a relative tolerance
            flag2 = 0;
        else
            flag2 = 1;
        end
        disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])

        
        
                
        %check if division rates converge
        [mui] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(1),1:14),hr1,hr2, ts);
        [mui2] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(2),1:14),hr1,hr2, ts);
        [mui3] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(3),1:14),hr1,hr2, ts);
        [mui4] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(4),1:14),hr1,hr2, ts);
        [mui5] = growth_rate(Einterp,volbins,COUNTS,modelfits(ii(5),1:14),hr1,hr2, ts);
        
        
        flag3 = (range([mui mui2 mui3 mui4 mui5]) < .001 ) ;
        
        
    end  %while loop

    
    %Assign Outputs
    [~, jj]=sort(modelfits(:,15));
    xmin=modelfits(jj(1),1:14);
    fmin=modelfits(jj(1),15);
    exitflag=modelfits(jj(1),17);
    
    %calculate growth rate according to ideal theta
    [mu, mu1, mu2] = growth_rate(Einterp,volbins,COUNTS,xmin(1:13),hr1,hr2, ts);

    modelresults=[day xmin fmin mu mu1 mu2 exitflag length(modelfits) ts]; %this is how the optimization results are stored in model output files 

    
% %Step 5: Simulate according to results of optimization 
% 
[simCONC, simPROPS] = simdata(Einterp,COUNTS,xmin,volbins,hr1,hr2, ts);
simPROPS = [NaN*ones(length(volbins), hr1-1) simPROPS]; %fill in null hours before model kicked in 
simCONC = [NaN*ones(length(volbins), hr1-1) simCONC]; %fill in null hours before model kicked in 


end 




















