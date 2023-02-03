function [phi_sigma,freq] = tsvr_all_sigma(ts_pops)

% Created by Silke van Daalen, January 2023
% Based on Eqs in Fig 2 and appendix from Zhao et al., 2020
% This function computes the variance ratio at all possible timescales 
% (as determined by the fft frequencies), it requires:
%   - ts_pops should contain a timeseries of all populations, with population
% identity across columns and time across rows 

npops=size(ts_pops,2);
ts_tot=sum(ts_pops,2); % sum subpopulations into aggregate

% individual spectral densities for total and subpops
[cosp_tot,freq]=cospectrum(ts_tot'); % transform total timeseries into cospectrum
tslength=length(ts_tot);
numer=zeros(tslength-1,1);
numer(:,1)=cosp_tot(1,1,2:end); % enumerator
freq=freq(:,2:end);
all_spects=zeros(length(ts_pops)-1,npops);
for k=1:npops
    [h_cosp,h_freq]=cospectrum(ts_pops(:,k)');
    all_spects(:,k)=h_cosp(1,1,2:end);
end
denom=sum(all_spects,2); % denominator
phi_sigma=numer./denom; % I think this is phi(sigma)

return