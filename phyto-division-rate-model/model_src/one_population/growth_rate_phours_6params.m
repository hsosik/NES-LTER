function [mu]=growth_rate_phours_6params(Einterp,volbins,N_dist,theta,hr1,hr2)
%really, only the first hour of N_dist or Vhists is needed

gmax=theta(1);
b=theta(2);
E_star=theta(3);
dmax=theta(4);
m1=theta(5);
sigma=theta(6);

q=hr2-hr1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
 %create all B matrices for the hours:

B_day=zeros(57,57,q);

for t=(hr1-1):(hr2-2)
     B=matrix_const(t,Einterp,volbins,b,dmax,E_star,gmax);
     B_day(:,:,t-hr1+2)=B;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1=normpdf(1:57,m1,sigma);
Nt=sum(N_dist(:,hr1))*(y1./sum(y1));
Nt=Nt';

simdist=Nt./sum(Nt);  %Only if starting hour has no zeros!

for t=1:q

    Nt(:,t+1)=B_day(:,:,t)*Nt(:,t);           %project forward with the numbers
    simdist(:,t+1)=Nt(:,t+1)./sum(Nt(:,t+1)); %normalize to get distribution for likelihood

    if any(isnan(simdist(:,t+1))) %just in case
        disp('hmmm...simdist has a nan?')
         keyboard
    end

end


mu=(log(sum((Nt(:,end)))./sum(Nt(:,1))));
