function [ts_phi_x,freq] = ts_loreau_vr(x)
%ts_loreau_vr A timescale specific community variance ratio.
%   This function takes as its only input a matrix x with columns (one for
%   each species) representing some population level measure (like
%   population size).  It returns a variance ratio, phi_x, defined by 
%   defined in eq. 9 of Loreau and de Mazancourt (Am Nat, 2008, 
%   DOI: 10.1086/589746).
%
%   Last edited: 2023-01-31 2:26pm by MGN

    npops=size(x,2);
    x_tau = sum(x,2);   % aggregate quantity 
                        %  (i.e. total abundance if x is population sizes).

    [cosp_tot,freq]=cospectrum(x_tau'); % transform total timeseries into cospectrum
    tslength=length(x_tau);
    numer=zeros(tslength-1,1);
    numer(:,1)=cosp_tot(1,1,2:end); % enumerator
    freq=freq(:,2:end);

    denom = zeros(size(numer));
    for i = 1:npops
        [cosp,freq] = cospectrum(x(:,i)');
        denom(:,1) = denom(:,1) + sqrt(squeeze(cosp(1,1,2:end)));
    end
    denom = denom.^2;

    ts_phi_x = numer./denom;
    freq=freq(:,2:end);
end