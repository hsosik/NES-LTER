function [dirsample, simdist,Vt1,Vt2]=simdata_dirichlet_sample_plt(Einterp,N_dist,theta,volbins,hr1,hr2, ts)

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



q=hr2-hr1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create all B matrices for the hours:
numbins = length(volbins);

B_day1=zeros(numbins,numbins,q);
B_day2=zeros(numbins,numbins,q);

for t=(hr1-1):(hr2-2)
     B1=matrix_const(t,Einterp,volbins,b1,dmax1,E_star1,gmax1, ts);
     B_day1(:,:,t-hr1+2)=B1;
     B2=matrix_const(t,Einterp,volbins,b2,dmax2,E_star2,gmax2, ts);
     B_day2(:,:,t-hr1+2)=B2;
end


%% Project forward each subcomponent: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1=normpdf(1:numbins,m1,sigma1); %relative distribution of each subpopulation assuming normal distributions 
y2=normpdf(1:numbins,m2,sigma2);

%distribute initial # of cells into two pops according to f 
%uses N_dist to determine initial concentration of cells, but we normalize
%later
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
simdist(:,1)=(Nt1+Nt2)./sum(Nt1+Nt2);  %Only if starting hour has no zeros!
else 
    simdist(:,1)=(Nt1)./sum(Nt1);  %Only if starting hour has no zeros!
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
        disp(['DMN 13param...simdist has a nan? theta:' num2str(theta)])
%         keyboard
    end

end

%simdist is then used to generate sample from the Dirichlet distribution:

dirsample=zeros(numbins,q);
for i=1:q+1
dirsample(:,i)=gamrnd(s*simdist(:,i),1);
dirsample(:,i)=dirsample(:,i)./sum(dirsample(:,i));
dirsample(:,i)=mnrnd(round(sum(N_dist(:,hr1-1+i))),dirsample(:,i));
end
end
