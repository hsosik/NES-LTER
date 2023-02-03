function phi_x = loreau_vr(x)
%loreau_vr A community variance ratio.
%   This function takes as its only input a matrix x with columns (one for
%   each species) representing some population level measure (like
%   population size).  It returns a variance ratio, phi_x, defined by 
%   defined in eq. 9 of Loreau and de Mazancourt (Am Nat, 2008, 
%   DOI: 10.1086/589746).
%
%   Last edited: 2023-01-31 2:26pm by MGN

    x_tau = sum(x,2);   % aggregate quantity 
                        %  (i.e. total abundance if x is population sizes).
    
    phi_x = var(x_tau)/(sum(std(x))^2);
    
end