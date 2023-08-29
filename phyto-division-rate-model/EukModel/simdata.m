%this function simulates a day of cell observations based on projections from
%matrix model. 

%Outputs: 
%simPROPS is a matrix with the expected proportion of cells in each size class
%(rows) as specified by volbins, for each hour (collumns) according to the
%simulation. 
%Vt1 and Vt2 are basically versions of simPROPS for each subpopulation 
%sampleCOUNTS is a matrix with the number of cells in each size class (rows)
% for each hour (collumns) as sampled from distributions according to simPROPS
% sample size is exactly that of the observed data for each hour. 

%Inputs: 
%Einterp is a 1 x 151 vector with radiation data in W/m^2 every 10 minutes 
%COUNTS is a matrix with the number of cells in each size class (rows) as
%specificed by volbins, for each hour (collumns).
% volbins is the designation of size classes that the observations are
% broken into 
% hr 1 is when we want our model to start
% hr 2 is when we want our model to end 

% theta is 1 x 14 vector of parameters as defined below 
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


function [simCONC,simPROPS,Nt1,Nt2]=simdata(Einterp,N_dist,theta,volbins,hr1,hr2, ts)

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
s=100*theta(14);

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


numbins = length(volbins) ; 
q=hr2-hr1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create all B matrices for the hours:

B_day1=zeros(numbins,numbins,q);
B_day2=zeros(numbins,numbins,q);

for t=(hr1-1):(hr2-2)
     B1=matrix_const(t,Einterp,volbins,b1,dmax1,E_star1,gmax1, ts);
     B_day1(:,:,t-hr1+2)=B1;
     B2=matrix_const(t,Einterp,volbins,b2,dmax2,E_star2,gmax2, ts);
     B_day2(:,:,t-hr1+2)=B2;
end


%% Project forward each subcomponent: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%relative distribution of each subpopulation assuming normal distributions 
y1=normpdf(1:numbins,m1,sigma1); 
y2=normpdf(1:numbins,m2,sigma2);

%counts of cells in each size bin, each hour for each subpopulation 
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
Vt1=Nt1./sum(Nt1+Nt2);
Vt2=Nt2./sum(Nt1+Nt2);

if ~isnan(Nt2)
simPROPS(:,1)=(Nt1+Nt2)./sum(Nt1+Nt2);  %Only if starting hour has no zeros!
else 
    simPROPS(:,1)=(Nt1)./sum(Nt1);  %Only if starting hour has no zeros!
end

for t=1:q

    Nt1(:,t+1)=B_day1(:,:,t)*Nt1(:,t);           %project forward with the numbers
    Nt2(:,t+1)=B_day2(:,:,t)*Nt2(:,t);
    
    if ~isnan(Nt2)
    simPROPS(:,t+1)=(Nt1(:,t+1)+Nt2(:,t+1))./sum(Nt1(:,t+1)+Nt2(:,t+1)); %normalize to get distribution for likelihood
    else
     simPROPS(:,t+1)=(Nt1(:,t+1))./sum(Nt1(:,t+1)); %normalize to get distribution for likelihood
    end
            
    if any(isnan(simPROPS(:,t+1))) %just in case
        disp(['DMN 13param...simdist has a nan? theta:' num2str(theta)])
%         keyboard
    end

end

if ~isnan(Nt2)
simCONC = Nt1 + Nt2; 
else
simCONC = Nt1; 
end
