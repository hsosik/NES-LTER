function [phi_sigma,freq] = tsvr_all_sigma(ts_pops)

% This function computes a short and long timescale variance ratio, it
% requires:
%   - ts_pops should contain a timeseries of all populations, with pop identity
%   across rows and time across columns
%   - threshold determines the divide between the short and long timescale 

npops=size(ts_pops,1);
ts_tot=sum(ts_pops); % sum subpopulations into aggregate

% individual spectral densitiies for total and subpops
[cosp_tot,freq]=cospectrum(ts_tot); % transform total timeseries into cospectrum
tslength=length(ts_tot);
numer=zeros(tslength-1,1);
numer(:,1)=cosp_tot(1,1,2:end);
freq=freq(:,2:end);
all_spects=zeros(length(ts_pops)-1,npops);
for k=1:npops
    [h_cosp,h_freq]=cospectrum(ts_pops(k,:));
    all_spects(:,k)=h_cosp(1,1,2:end);
end
denom=sum(all_spects,2);
phi_sigma=numer./denom; % I think this is phi(sigma)

return