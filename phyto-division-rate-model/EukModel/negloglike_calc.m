function [prob]=loglike_DMN_14params_phours_plateau(Einterp,N_dist,theta,volbins,hr1,hr2,ts)
%July 27 , 2020 adjusted so that if Theta is sized for one_pop, will assign
%appropriate parameter values. 

%Calculates negative log likelihood of a set of parameters given a day of
%observations (counts of cells in each size class for each hour). Model version is a two subpopulation model structure
% as described in Hunter-Cevera et al. 2014. Subpopulations are
% distinguished based on starting mean volume. The subpopulation with the
% smaller volume is referred to subpopn 1

%The inputs are as follows:

%Einterp - interpolated light data for every 10 min of the day (W/m2)
%N_dist = number of counts of cells in each size class as specified by volbins
%volbins - cell size classes (micrometers cubed)
%hr1 and hr2 refer to the starting and ending hour of the portion of day
%you want to fit. In the paper, we've used hr1=7 hours after dawn and
%hr2=25 (run till end of day).
%theta - set of parameters, described below:

%gmax1=theta(1); %max fraction of cells growing into next size class, subpopn 1
%b1=theta(2);  %shape parameter division function, subpopn 1
%E_star1=theta(3); %shape parameter of growth function (point where function switches from linear to constant), subpopn 1
%dmax1=theta(4); %max fraction of cells able to divide in a given size class, subpopn 1
%gmax2=theta(5); %max fraction of cells growing into next size class, subpopn 2
%b2=theta(6); %shape parameter division function, subpopn 2
%E_star2=theta(7); %shape parameter of growth function (point where function switches from linear to constant), subpopn 2
%dmax2=theta(8); %max fraction of cells able to divide in a given size class, subpopn 2
%f=theta(9); %proportion parameter, specifies starting fraction of subpopn 1
%m1=theta(10); %mean volume for starting cell size distribution, subpopn 1
%m2=theta(11); %mean volume for starting cell size distribution, subpopn 2
%sigma1=theta(12); %variance parameter for starting cell size distributions for popn 1
%sigma2=theta(13); %variance parameter for starting cell size distributions for popn 2
%s=theta(14); %overdispersion parameter for the Dirichlet-multinomial distribution

if length(theta) >= 13
    
gmax1=theta(1);
b1=theta(2);
E_star1=theta(3);
dmax1=theta(4);
gmax2=theta(5);
b2=theta(6);
E_star2=theta(7);
dmax2=theta(8);
f=theta(9); %fraction of starting distribution
m1=theta(10);
m2=theta(11);
sigma1=theta(12);
sigma2=theta(13);
s=100*theta(14); %scaled to help solver shrink distance between small numbers and very large numbers encountered in theta vector

elseif length(theta) == 7
    
gmax1=theta(1);
b1=theta(2);
E_star1=theta(3);
dmax1=theta(4);
gmax2=NaN;
b2=NaN;
E_star2=NaN;
dmax2=NaN;
f=1; %fraction of starting distribution
m1=theta(5);
m2=NaN;
sigma1=theta(6);
sigma2=NaN;
s=100*theta(7);
end


q=hr2-hr1;
numbins = length(volbins);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create all B matrices for the hours:

B_day1=zeros(numbins,numbins,q);
B_day2=zeros(numbins,numbins,q);

for t=(hr1-1):(hr2-2)
     B1=matrix_const(t,Einterp,volbins,b1,dmax1,E_star1,gmax1,ts);
     B_day1(:,:,t-hr1+2)=B1;
     B2=matrix_const(t,Einterp,volbins,b2,dmax2,E_star2,gmax2,ts);
     B_day2(:,:,t-hr1+2)=B2;
end


%% Project forward each subcomponent: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1=normpdf(1:numbins,m1,sigma1);
y2=normpdf(1:numbins,m2,sigma2);


%this just determines the initial size of the two simulated populations,
%which really doesn't matter since we normalize later
if ~isnan(sum(N_dist(:,hr1)))
    Nt1=f*sum(N_dist(:,hr1))*(y1./sum(y1));
    Nt2=(1-f)*sum(N_dist(:,hr1))*(y2./sum(y2));
else 
   a = find(~isnan(sum(N_dist)), 1); 
   Nt1=f*sum(N_dist(:,a))*(y1./sum(y1));
   Nt2=(1-f)*sum(N_dist(:,a))*(y2./sum(y2));
end

Nt1=Nt1';
Nt2=Nt2';

    
if ~isnan(Nt2)
    simdist(:,1)=(Nt1+Nt2)./sum(Nt1+Nt2);  %Only if starting hour isn't all zeros!
else
    simdist(:,1)= Nt1./sum(Nt1);  %Only if starting hour isn't all zeros!
end
        
for t=1:q

    Nt1(:,t+1)=B_day1(:,:,t)*Nt1(:,t);           %project forward with the numbers
    Nt2(:,t+1)=B_day2(:,:,t)*Nt2(:,t);
    
    
    if ~isnan(Nt2)
    simdist(:,t+1)=(Nt1(:,t+1)+Nt2(:,t+1))./sum(Nt1(:,t+1)+Nt2(:,t+1)); %normalize to get distribution for likelihood
    else
    simdist(:,t+1)=(Nt1(:,t+1))./sum(Nt1(:,t+1)); 
    end
    
    
    if any(isnan(simdist(:,t+1))) %just in case
        disp(['DMN 14param...simdist has a nan? theta:' num2str(theta)])
    end

end

%% Now calculate the log likelihood using the Dirichlet Multinomial distribution: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specifiy the expected distribution:
alpha=s*simdist;
TotN=sum(N_dist); 

logL=zeros(q+1,1);
for t=1:q+1
    if ~isnan(TotN(:,t+hr1-1))  
        C = gammaln(s) - gammaln(TotN(:,t+hr1-1)+s); %constant out in front
        logL(t)=C+sum(gammaln(N_dist(:,t+hr1-1)+alpha(:,t)) - gammaln(alpha(:,t)));
    end
end

prob=-sum(logL); %negative for the minimization routine fmincon

end 

