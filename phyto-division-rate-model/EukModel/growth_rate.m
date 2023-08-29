%This function calculates the division rate of the population for each hour
%according to the model's projection, 
%given the parameters, one day's light data, and an initial count of 
%cells in the population at the begining of a day. 

%Inputs: 
%Einterp is a 1 x 151 vector with radiation data in W/m^2 every 10 minutes 

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

% volbins is the designation of size classes that we want our population to
% be broken into 

% hr 1 is when we want our model to start
% hr 2 is when we want our model to end 

%N_Dist is a matrix with the number of cells in each size class (rows) as
%specificed by volbins, for each hour (collumns). Though in fact, this
%function only uses the first hour of COUNTS to get an initial population
%size. 

%Outputs: 
% mu_mat is the matrix of hourly division rates
% first row total population, 2nd and 3rd row subpopulations 
% pp is vector of subpopulation proportions of the total population at
% each hour 

function [mu mu1 mu2]=growth_rate_phours_13params_plateau(Einterp,volbins,N_dist,theta,hr1,hr2,ts)

%label individual parameters:
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


%useful variables 
q=hr2-hr1; 
numbins = length(volbins);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%create all B matrices for the hours:
B_day1=zeros(numbins,numbins,q);
B_day2=zeros(numbins,numbins,q); 
for t=(hr1-1):(hr2-2)
     [B1] =matrix_const(t,Einterp,volbins,b1,dmax1,E_star1,gmax1,ts);
     B_day1(:,:,t-hr1+2)=B1;
     [B2] =matrix_const(t,Einterp,volbins,b2,dmax2,E_star2,gmax2,ts);
     B_day2(:,:,t-hr1+2)=B2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create distributions according to parameters for two subpops 
y1=normpdf(1:numbins,m1,sigma1);
y2=normpdf(1:numbins,m2,sigma2);

%distribute initial # of cells into two pops according to f 
%start to build two matrices for the day, of counts of cells in each size bin
%(for each hour/collumn) one for each subpopulation 
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

%simPROPS is the matrix for the proportions of cells in each size bin (for each 
%hour) for the whole population
if ~isnan(Nt2)
 simPROPS(:,1)=(Nt1+Nt2)./sum(Nt1+Nt2);  %Only if starting hour has no zeros!
else
 simPROPS(:,1)=(Nt1)./sum(Nt1);  %Only if starting hour has no zeros!
end

for t=1:q %go through all the hours of interest

    Nt1(:,t+1)=B_day1(:,:,t)*Nt1(:,t); %project forward with the numbers
    Nt2(:,t+1)=B_day2(:,:,t)*Nt2(:,t);
    
    if ~isnan(Nt2)
    simPROPS(:,t+1)=(Nt1(:,t+1)+Nt2(:,t+1))./sum(Nt1(:,t+1)+Nt2(:,t+1)); %normalize to get distribution for likelihood
    elseif isnan(Nt2)
    simPROPS(:,t+1)=(Nt1(:,t+1))./sum(Nt1(:,t+1)); %normalize to get distribution for likelihood
    end
    
    if any(isnan(simPROPS(:,t+1))) %just in case
        disp('hmmm...simdist has a nan?')
         keyboard
    end
    
    %calculate hourly growth rate, equation is log(N(24)/N(0))*24 since unit of interest is 1 day
    %mu_mat = [mu_mat [24*log(sum(Nt1(:,t+1) + Nt2(:,t+1)) / sum(Nt1(:,t) + Nt2(:,t))); 24*log(sum(Nt1(:,t+1)) / sum(Nt1(:,t))); 24*log(sum(Nt2(:,t+1)) / sum(Nt2(:,t)))]]; 
        
end

if ~isnan(Nt2)
mu = (log(sum((Nt1(:,end)+Nt2(:,end)))./sum(Nt1(:,1)+Nt2(:,1)))); %division rate for whole day 
mu1=log(sum((Nt1(:,end)))./sum(Nt1(:,1)));
mu2=log(sum((Nt2(:,end)))./sum(Nt2(:,1)));
else
mu =log(sum(Nt1(:,end))./sum(Nt1(:,1))); %division rate for whole day
mu1=log(sum(Nt1(:,end))./sum(Nt1(:,1)));
mu2=NaN;
end

end


