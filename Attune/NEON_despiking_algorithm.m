%NEON's Despiking Algorithm
%function [spike_res]=spike_test(window,step_size,data,timeseries)
w = 10; %window size, if possible should always have more than nine points
step_size = 1; %step size
data = ones(100,1); %data set
data(51) = 100;
n = length(data);
q = 1; %MAD threshold
k = 1.4826;%constant scale factor
spike = zeros(n,1);

if n > 9
    b = n/(n-0.8);
elseif n == 4
      b=  1.363;
elseif n == 5
    b = 1.206;
elseif n == 6
    b = 1.200;
elseif n == 7
    b = 1.140;
elseif n == 8
    b = 1.129;
elseif n == 9
    b = 1.107;
elseif n < 4
    display("n must be greater than 4")
end
        
% Median Absolute Deviation

for i=1:s:n-10
    window = data(i:i+w-1,:);
    MAD = median(abs(data - median(window)));
    MAD_adjusted = b*q*k*MAD;
    if  ((data(i) <=(median(window)-MAD_adjusted))||(data(i) >=( median(window)+MAD_adjusted)))==1
        display('Spike Detected')
        display(i)
        spike(i) == 1;
    end
end
%end
