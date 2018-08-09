function syringe_acqtime_QC2(finv,ordered_acq,oind,ss,varargin)

n=length(ordered_acq);
subplot(1,2,1,'replace')
plot(finv,ordered_acq,'o'), hold on

if oind(end)==1 || oind(end)==n
    plot(finv(end),ordered_acq(end),'ro','markerface','r')
end

line([ordered_acq(1) ordered_acq(end)],[ordered_acq(1) ordered_acq(end)]) %expectation
title(['sum of squares:' num2str(ss)])

subplot(1,2,2,'replace')

if ~isempty(varargin)
    finv2=varargin{1};
    ordered_acq2=varargin{2};
    ss2=varargin{3};
    plot(finv2,ordered_acq2,'o')
    line([ordered_acq2(1) ordered_acq2(end)],[ordered_acq2(1) ordered_acq2(end)])    
    title(['sum of squares:' num2str(ss2) ' | std: ' num2str(std(ordered_acq))])
    
end
end %syringe acqtime plotting function


