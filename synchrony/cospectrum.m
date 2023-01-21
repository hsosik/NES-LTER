function [cosp,freqs]=cospectrum(X_ts)
% X_ts is a timeseries of n populations on the rows, and time series data
% across columns

tslength=length(X_ts);
freqs=[0:1/tslength:(1-1/tslength)];
npops=size(X_ts,1);
allffts=zeros(tslength,npops);
for i=1:npops
    allffts(:,i)=fft(X_ts(i,:)');
end
cosp=zeros(npops,npops,tslength);
for i=1:npops
    for j=1:npops
        cosp(i,j,:)=real(conj(allffts(:,i)).*allffts(:,j))/tslength/(tslength-1);
    end
end

return